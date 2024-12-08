/*
 * Copyright (C) 2024 Yifeng Li (tomli@tomli.me)
 * Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "engine_tiling.h"
#include "tools/denormal.h"
#include "tools/tiling/tbb/paraexec.h"
//#include "tools/tiling/tbb/first_touch.h"

//! \brief construct an Engine_Tiling instance
//! it's the responsibility of the caller to free the returned pointer
Engine_Tiling* Engine_Tiling::New(const Operator_Tiling* op)
{
	std::cout << "Create FDTD engine ";
	std::cout << "(compressed SSE + multi-threading + spatial/temporal tiling)";
	std::cout << std::endl;

	Engine_Tiling* e = new Engine_Tiling(op);
	e->Init();
	return e;
}

Engine_Tiling::Engine_Tiling(const Operator_Tiling* op) : Engine(op)
{
	m_type = TILING_V1;
	Op = op;

	// speed up the calculation of denormal floating point values (flush-to-zero)
	Denormal::Disable();
}

void Engine_Tiling::Init()
{
	Engine::Init();

	// This engine uses its own SIMD arrays to represent E&M fields, so
	// free the arrays from the base class.
	delete volt_ptr;
	volt_ptr = NULL;
	delete curr_ptr;
	curr_ptr = NULL;

	numVectors = ceil((double) numLines[2] / (double) veclen);

	// allocate our own storage
	// See engine_tiling.h for comments
	m_volt_v = new ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>(
		"volt_v", {numLines[0], numLines[1], numVectors}
	);

	m_curr_v = new ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>(
		"curr_v", {numLines[0], numLines[1], numVectors}
	);

	// TODO: modify ArrayLib to include a "no memset()" option
	// TODO: explicit first touch for NUMA

	for (const Tiling::Plan3D& plan : Op->m_perNodePlan)
		InitializeTiling(plan);
}

// Before simulation starts, we need to pass the pre-calculated tiling
// plan to each engine extension, so that each extension has prior
// knowledge about which tiles require extension, and which tiles do not.
void Engine_Tiling::InitializeTiling(const Tiling::Plan3D plan)
{
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		std::cerr << "Using ";
		std::cerr << m_Eng_exts.at(n)->GetExtensionName();
		std::cerr << std::endl;

		m_Eng_exts.at(n)->InitializeTiling(plan);
	}
}

bool Engine_Tiling::IterateTS(unsigned int iterTS)
{
	size_t numBatches = (iterTS * 2) / Op->m_tileHalfTs;
	size_t numRem     = (iterTS * 2) % Op->m_tileHalfTs;

	if (numRem > 0)
		numBatches += 1;

	for (size_t batchId = 0; batchId < numBatches; batchId++)
	{
		for (size_t stage = 0; stage < Op->m_perNodePlan[0].size(); stage++)
		{
			for (size_t nodeId = 0; nodeId < ParaExec::numNodes(); nodeId++)
			{
				const Tiling::TileList3D& tileList = Op->m_perNodePlan[nodeId][stage];

				ParaExec::paraFor(nodeId, tileList.size(), [&](size_t tileId)
				{
					const Tiling::Tile3D& tile = tileList[tileId];

					for (const Tiling::Subtile3D& subtile : tile)
					{
						size_t totalHt;
						if (numRem > 0 && batchId == numBatches - 1)
							totalHt = numRem;
						else
							totalHt = subtile.size();

						for (size_t ht = 0; ht < totalHt; ht++)
						{
							if (ht % 2 == 0)
							{
								const Tiling::Range3D<size_t>& vRange = subtile[ht];

								DoPreVoltageUpdates(numTS + ht / 2, vRange);
								updateVoltageRange(vRange);
								DoPostVoltageUpdates(numTS + ht / 2, vRange);
								Apply2Voltages(numTS + ht / 2, vRange);
							}
							else
							{
								const Tiling::Range3D<size_t>& cRange = subtile[ht];

								DoPreCurrentUpdates(numTS + ht / 2, cRange);
								updateCurrentRange(cRange);
								DoPostCurrentUpdates(numTS + ht / 2, cRange);
								Apply2Current(numTS + ht / 2, cRange);
							}
						}
					}
				});
			}
			ParaExec::barrier();
		}

		if (numRem > 0 && batchId == numBatches - 1)
			numTS += numRem / 2;
		else
			numTS += Op->m_tileHalfTs / 2;
	}

	return true;
}

void Engine_Tiling::updateVoltageRange(const Tiling::Range3D<> range)
{
	// The K dimension is vectorized but the loop range (first_k, last_k)
	// is often at the middle of some vectors. Thus, we split the inner
	// loop into three ranges: head, body, tail. The body is exactly aligned
	// to a vector, but the misaligned head and tail need special treatments.
	std::pair<uint32_t, uint32_t> headRange, bodyRange, tailRange;
	constexpr const uint32_t invalid = UINT32_MAX;

	if (range.last[2] - range.first[2] + 1 < veclen * 2 &&
	    range.first[2] % veclen != 0 && (range.last[2] + 1) % veclen != 0)
	{
		headRange = {range.first[2], range.last[2]};
		bodyRange = {invalid, invalid};
		tailRange = {invalid, invalid};
	}
	else {
		if (range.first[2] % veclen == 0)
		{
			headRange = {invalid, invalid};
			bodyRange.first = range.first[2];
		}
		else
		{
			headRange.first = range.first[2];
			bodyRange.first = range.first[2] +
			                  (veclen - (range.first[2] % veclen));
			headRange.second = bodyRange.first - 1;
		}

		if ((range.last[2] + 1) % veclen == 0)
		{
			tailRange = {invalid, invalid};
			bodyRange.second = range.last[2];
		}
		else
		{
			bodyRange.second = range.last[2] - (range.last[2] % veclen) - 1;
			tailRange.first = bodyRange.second + 1;
			tailRange.second = range.last[2];
		}
	}

	for (uint32_t i = range.first[0]; i <= range.last[0]; i++)
	{
		for (uint32_t j = range.first[1]; j <= range.last[1]; j++)
		{
			// prologue: head
			if (headRange.first != invalid)
			{
				updateVoltageScalarKernel(
					i, j,
					headRange.first, headRange.second
				);
			}

			// body
			if (bodyRange.first == 0)
			{
					updateVoltageVectorKernel<true>(
						i, j,
						0, veclen - 1
					);
					if (bodyRange.second - bodyRange.first > veclen - 1)
					{
						updateVoltageVectorKernel<false>(
							i, j,
							veclen, bodyRange.second
						);
					}
			}
			else if (bodyRange.first != invalid)
			{
				updateVoltageVectorKernel<false>(
					i, j,
					bodyRange.first, bodyRange.second
				);
			}

			// epilogue: tail
			if (tailRange.first != invalid)
			{
				updateVoltageScalarKernel(
					i, j,
					tailRange.first, tailRange.second
				);
			}
		}
	}
}

inline void
Engine_Tiling::updateVoltageScalarKernel(
	uint32_t i,
	uint32_t j,
	uint32_t first_k,
	uint32_t last_k
)
{
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& volt_v = *m_volt_v;
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& curr_v = *m_curr_v;
	const Operator_Tiling& op = *Op;

	uint32_t pi = i > 0 ? i - 1 : 0;
	uint32_t pj = j > 0 ? j - 1 : 0;

	for (uint32_t k = first_k; k <= last_k; k++)
	{
		// 3 FP32 loads
		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;

		float volt0_ci_cj_ck = volt_v(0, i, j, vk)[vk_offset];
		float volt1_ci_cj_ck = volt_v(1, i, j, vk)[vk_offset];
		float volt2_ci_cj_ck = volt_v(2, i, j, vk)[vk_offset];

		// 3 FP32 loads
		float vv0_ci_cj_ck = op.vv_v(0, i, j, vk)[vk_offset];
		float vv1_ci_cj_ck = op.vv_v(1, i, j, vk)[vk_offset];
		float vv2_ci_cj_ck = op.vv_v(2, i, j, vk)[vk_offset];

		// 3 FP32 loads
		float vi0_ci_cj_ck = op.vi_v(0, i, j, vk)[vk_offset];
		float vi1_ci_cj_ck = op.vi_v(1, i, j, vk)[vk_offset];
		float vi2_ci_cj_ck = op.vi_v(2, i, j, vk)[vk_offset];

		// 9 FP32 loads
		assert(k > 0);
		uint32_t pk = k - 1;
		uint32_t vpk = pk / veclen;
		uint32_t vpk_offset = pk % veclen;

		float curr0_ci_cj_ck = curr_v(0, i,  j,  vk )[vk_offset];
		float curr1_ci_cj_ck = curr_v(1, i,  j,  vk )[vk_offset];
		float curr2_ci_cj_ck = curr_v(2, i,  j,  vk )[vk_offset];
		float curr0_ci_cj_pk = curr_v(0, i,  j,  vpk)[vpk_offset];
		float curr1_ci_cj_pk = curr_v(1, i,  j,  vpk)[vpk_offset];
		float curr0_ci_pj_ck = curr_v(0, i,  pj, vk )[vk_offset];
		float curr2_ci_pj_ck = curr_v(2, i,  pj, vk )[vk_offset];
		float curr1_pi_cj_ck = curr_v(1, pi, j,  vk )[vk_offset];
		float curr2_pi_cj_ck = curr_v(2, pi, j,  vk )[vk_offset];

		// 6 FLOPs, for x polarization
		volt0_ci_cj_ck *= vv0_ci_cj_ck;
		volt0_ci_cj_ck +=
			vi0_ci_cj_ck * (
				curr2_ci_cj_ck -
				curr2_ci_pj_ck -
				curr1_ci_cj_ck +
				curr1_ci_cj_pk
			);

		// 6 FLOPs, for y polarization
		volt1_ci_cj_ck *= vv1_ci_cj_ck;
		volt1_ci_cj_ck +=
			vi1_ci_cj_ck * (
				curr0_ci_cj_ck -
				curr0_ci_cj_pk -
				curr2_ci_cj_ck +
				curr2_pi_cj_ck
			);

		// 6 FLOPs, for z polarization
		volt2_ci_cj_ck *= vv2_ci_cj_ck;
		volt2_ci_cj_ck +=
			vi2_ci_cj_ck * (
				curr1_ci_cj_ck -
				curr1_pi_cj_ck -
				curr0_ci_cj_ck +
				curr0_ci_pj_ck
			);

		volt_v(0, i, j, vk)[vk_offset] = volt0_ci_cj_ck;
		volt_v(1, i, j, vk)[vk_offset] = volt1_ci_cj_ck;
		volt_v(2, i, j, vk)[vk_offset] = volt2_ci_cj_ck;
	}
}

template <bool boundary>
inline void
Engine_Tiling::updateVoltageVectorKernel(
	uint32_t i,
	uint32_t j,
	uint32_t first_k,
	uint32_t last_k
)
{
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& volt_v = *m_volt_v;
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& curr_v = *m_curr_v;
	const Operator_Tiling& op = *Op;

	static_assert(veclen == 4, "Only float4 vectors are implemented for now");

	uint32_t pi = i > 0 ? i - 1 : 0;
	uint32_t pj = j > 0 ? j - 1 : 0;

	uint32_t first_vk = first_k / veclen;
	uint32_t last_vk = last_k / veclen;

	for (uint32_t vk = first_vk; vk <= last_vk; vk += 1)
	{
		Simd<float, veclen> volt0_ci_cj_ck = volt_v(0, i, j, vk);
		Simd<float, veclen> volt1_ci_cj_ck = volt_v(1, i, j, vk);
		Simd<float, veclen> volt2_ci_cj_ck = volt_v(2, i, j, vk);

		// 12 (3 x 4) FP32 loads
		Simd<float, veclen> curr0_ci_cj_ck = curr_v(0, i, j, vk);
		Simd<float, veclen> curr1_ci_cj_ck = curr_v(1, i, j, vk);
		Simd<float, veclen> curr2_ci_cj_ck = curr_v(2, i, j, vk);

		// 16 (4 x 4) FP32 loads
		Simd<float, veclen> curr0_ci_pj_ck = curr_v(0, i,  pj, vk);
		Simd<float, veclen> curr2_ci_pj_ck = curr_v(2, i,  pj, vk);
		Simd<float, veclen> curr1_pi_cj_ck = curr_v(1, pi, j,  vk);
		Simd<float, veclen> curr2_pi_cj_ck = curr_v(2, pi, j,  vk);

		// 2 misaligned FP32 loads
		Simd<float, veclen> curr0_ci_cj_pk;
		curr0_ci_cj_pk[1] = curr0_ci_cj_ck[0];
		curr0_ci_cj_pk[2] = curr0_ci_cj_ck[1];
		curr0_ci_cj_pk[3] = curr0_ci_cj_ck[2];

		Simd<float, veclen> curr1_ci_cj_pk;
		curr1_ci_cj_pk[1] = curr1_ci_cj_ck[0];
		curr1_ci_cj_pk[2] = curr1_ci_cj_ck[1];
		curr1_ci_cj_pk[3] = curr1_ci_cj_ck[2];

		// "boundary" is a templated Boolean constant, no branching occurs
		if (boundary)
		{
			// outside boundary
			curr0_ci_cj_pk[0] = curr0_ci_cj_ck[0];
			curr1_ci_cj_pk[0] = curr1_ci_cj_ck[0];
		}
		else
		{
			curr0_ci_cj_pk[0] = curr_v(0, i, j, vk - 1)[3];
			curr1_ci_cj_pk[0] = curr_v(1, i, j, vk - 1)[3];
		}

		// 24 (6 x 4) FP32 loads
		Simd<float, veclen> vv0_ci_cj_ck = op.vv_v(0, i, j, vk);
		Simd<float, veclen> vv1_ci_cj_ck = op.vv_v(1, i, j, vk);
		Simd<float, veclen> vv2_ci_cj_ck = op.vv_v(2, i, j, vk);
		Simd<float, veclen> vi0_ci_cj_ck = op.vi_v(0, i, j, vk);
		Simd<float, veclen> vi1_ci_cj_ck = op.vi_v(1, i, j, vk);
		Simd<float, veclen> vi2_ci_cj_ck = op.vi_v(2, i, j, vk);

		// x-polarization
		volt0_ci_cj_ck *= vv0_ci_cj_ck;
		volt0_ci_cj_ck +=
			vi0_ci_cj_ck * (
				curr2_ci_cj_ck -
				curr2_ci_pj_ck -
				curr1_ci_cj_ck +
				curr1_ci_cj_pk
			);

		// y-polarization
		volt1_ci_cj_ck *= vv1_ci_cj_ck;
		volt1_ci_cj_ck +=
			vi1_ci_cj_ck * (
				curr0_ci_cj_ck -
				curr0_ci_cj_pk -
				curr2_ci_cj_ck +
				curr2_pi_cj_ck
			);

		// z-polarization
		volt2_ci_cj_ck *= vv2_ci_cj_ck;
		volt2_ci_cj_ck +=
			vi2_ci_cj_ck * (
				curr1_ci_cj_ck -
				curr1_pi_cj_ck -
				curr0_ci_cj_ck +
				curr0_ci_pj_ck
			);

		// 12 (3 x 4) FP32 stores
		volt_v(0, i, j, vk) = volt0_ci_cj_ck;
		volt_v(1, i, j, vk) = volt1_ci_cj_ck;
		volt_v(2, i, j, vk) = volt2_ci_cj_ck;
	}
}

void Engine_Tiling::updateCurrentRange(const Tiling::Range3D<> range)
{
	// The K dimension is vectorized but the loop range (first_k, last_k)
	// is often at the middle of some vectors. Thus, we split the inner
	// loop into three ranges: head, body, tail. The body is exactly aligned
	// to a vector, but the misaligned head and tail need special treatments.
	std::pair<uint32_t, uint32_t> headRange, bodyRange, tailRange;
	uint32_t invalid = UINT32_MAX;

	if (range.last[2] - range.first[2] + 1 < veclen * 2 &&
	    range.first[2] % veclen != 0 && (range.last[2] + 1) % veclen != 0)
	{
		headRange = {range.first[2], range.last[2]};
		bodyRange = {invalid, invalid};
		tailRange = {invalid, invalid};
	}
	else {
		if (range.first[2] % veclen == 0)
		{
			headRange = {invalid, invalid};
			bodyRange.first = range.first[2];
		}
		else
		{
			headRange.first = range.first[2];
			bodyRange.first = range.first[2] +
			                  (veclen - (range.first[2] % veclen));
			headRange.second = bodyRange.first - 1;
		}

		if ((range.last[2] + 1) % veclen == 0)
		{
			tailRange = {invalid, invalid};
			bodyRange.second = range.last[2];
		}
		else
		{
			bodyRange.second = range.last[2] - (range.last[2] % veclen) - 1;
			tailRange.first = bodyRange.second + 1;
			tailRange.second = range.last[2];
		}
	}

	for (uint32_t i = range.first[0]; i <= range.last[0]; i++)
	{
		for (uint32_t j = range.first[1]; j <= range.last[1]; j++)
		{
			// prologue: head
			if (headRange.first != invalid)
			{
				updateCurrentScalarKernel(
					i, j,
					headRange.first, headRange.second
				);
			}

			// body
			if (bodyRange.second == numLines[2] - 1)
			{
				if (bodyRange.second - bodyRange.first > veclen - 1)
				{
					updateCurrentVectorKernel<false>(
						i, j,
						bodyRange.first, bodyRange.second - veclen
					);
					updateCurrentVectorKernel<true>(
						i, j,
						bodyRange.second - veclen + 1, bodyRange.second
					);
				}
				else {
					updateCurrentVectorKernel<true>(
						i, j,
						bodyRange.first, bodyRange.second
					);
				}
			}
			else if (bodyRange.first != invalid)
			{
				updateCurrentVectorKernel<false>(
					i, j,
					bodyRange.first, bodyRange.second
				);
			}

			// epilogue: tail
			if (tailRange.first != invalid)
			{
				updateCurrentScalarKernel(
					i, j,
					tailRange.first, tailRange.second
				);
			}
		}
	}
}

inline void
Engine_Tiling::updateCurrentScalarKernel(
	uint32_t i,
	uint32_t j,
	uint32_t first_k,
	uint32_t last_k
)
{
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& curr_v = *m_curr_v;
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& volt_v = *m_volt_v;
	const Operator_Tiling& op = *Op;

	uint32_t pi = i > 0 ? i - 1 : 0;
	uint32_t pj = j > 0 ? j - 1 : 0;

	for (uint32_t k = first_k; k <= last_k; k++)
	{
		// 3 FP32 loads
		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;

		float curr0_ci_cj_ck = curr_v(0, i, j, vk)[vk_offset];
		float curr1_ci_cj_ck = curr_v(1, i, j, vk)[vk_offset];
		float curr2_ci_cj_ck = curr_v(2, i, j, vk)[vk_offset];

		// 3 FP32 loads
		float ii0_ci_cj_ck = op.ii_v(0, i, j, vk)[vk_offset];
		float ii1_ci_cj_ck = op.ii_v(1, i, j, vk)[vk_offset];
		float ii2_ci_cj_ck = op.ii_v(2, i, j, vk)[vk_offset];

		// 3 FP32 loads
		float iv0_ci_cj_ck = op.iv_v(0, i, j, vk)[vk_offset];
		float iv1_ci_cj_ck = op.iv_v(1, i, j, vk)[vk_offset];
		float iv2_ci_cj_ck = op.iv_v(2, i, j, vk)[vk_offset];

		// 9 FP32 loads
		uint32_t nk = k + 1;
		uint32_t vnk = nk / veclen;
		uint32_t vnk_offset = nk % veclen;

		float volt0_ci_cj_ck = volt_v(0, i,   j,   vk )[vk_offset];
		float volt1_ci_cj_ck = volt_v(1, i,   j,   vk )[vk_offset];
		float volt0_ci_cj_nk = volt_v(0, i,   j,   vnk)[vnk_offset];
		float volt1_ci_cj_nk = volt_v(1, i,   j,   vnk)[vnk_offset];
		float volt2_ci_cj_ck = volt_v(2, i,   j,   vk )[vk_offset];
		float volt0_ci_nj_ck = volt_v(0, i,   j+1, vk )[vk_offset];
		float volt2_ci_nj_ck = volt_v(2, i,   j+1, vk )[vk_offset];
		float volt1_ni_cj_ck = volt_v(1, i+1, j,   vk )[vk_offset];
		float volt2_ni_cj_ck = volt_v(2, i+1, j,   vk )[vk_offset];

		// 6 FLOPs, for x polarization
		curr0_ci_cj_ck *= ii0_ci_cj_ck;
		curr0_ci_cj_ck +=
			iv0_ci_cj_ck * (
				volt2_ci_cj_ck -
				volt2_ci_nj_ck -
				volt1_ci_cj_ck +
				volt1_ci_cj_nk
			);

		// 6 FLOPs, for y polarization
		curr1_ci_cj_ck *= ii1_ci_cj_ck;
		curr1_ci_cj_ck +=
			iv1_ci_cj_ck * (
				volt0_ci_cj_ck -
				volt0_ci_cj_nk -
				volt2_ci_cj_ck +
				volt2_ni_cj_ck
			);

		// 6 FLOPs, for z polarization
		curr2_ci_cj_ck *= ii2_ci_cj_ck;
		curr2_ci_cj_ck +=
			iv2_ci_cj_ck * (
				volt1_ci_cj_ck -
				volt1_ni_cj_ck -
				volt0_ci_cj_ck +
				volt0_ci_nj_ck
			);

		// 3 FP32 stores
		curr_v(0, i, j, vk)[vk_offset] = curr0_ci_cj_ck;
		curr_v(1, i, j, vk)[vk_offset] = curr1_ci_cj_ck;
		curr_v(2, i, j, vk)[vk_offset] = curr2_ci_cj_ck;
	}
}

template <bool boundary>
inline void
Engine_Tiling::updateCurrentVectorKernel(
	uint32_t i,
	uint32_t j,
	uint32_t first_k,
	uint32_t last_k
)
{
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& curr_v = *m_curr_v;
	ArrayLib::ArrayNIJK<Simd<float, veclen>, size_t>& volt_v = *m_volt_v;
	const Operator_Tiling& op = *Op;

	static_assert(veclen == 4, "Only float4 vectors are implemented for now");

	uint32_t pi = i > 0 ? i - 1 : 0;
	uint32_t pj = j > 0 ? j - 1 : 0;

	uint32_t first_vk = first_k / veclen;
	uint32_t last_vk = last_k / veclen;

	for (uint32_t vk = first_vk; vk <= last_vk; vk += 1)
	{
		// 12 (3 x 4) FP32 loads
		Simd<float, veclen> curr0_ci_cj_ck = curr_v(0, i,     j,     vk);
		Simd<float, veclen> curr1_ci_cj_ck = curr_v(1, i,     j,     vk);
		Simd<float, veclen> curr2_ci_cj_ck = curr_v(2, i,     j,     vk);

		// 12 (3 x 4) FP32 loads
		Simd<float, veclen> volt0_ci_cj_ck = volt_v(0, i,     j,     vk);
		Simd<float, veclen> volt1_ci_cj_ck = volt_v(1, i,     j,     vk);
		Simd<float, veclen> volt2_ci_cj_ck = volt_v(2, i,     j,     vk);

		// 16 (4 x 4) FP32 loads
		Simd<float, veclen> volt0_ci_nj_ck = volt_v(0, i,     j + 1, vk);
		Simd<float, veclen> volt2_ci_nj_ck = volt_v(2, i,     j + 1, vk);
		Simd<float, veclen> volt1_ni_cj_ck = volt_v(1, i + 1, j,     vk);
		Simd<float, veclen> volt2_ni_cj_ck = volt_v(2, i + 1, j,     vk);

		// 2 misaligned FP32 loads
		Simd<float, veclen> volt0_ci_cj_nk;
		volt0_ci_cj_nk[0] = volt0_ci_cj_ck[1];
		volt0_ci_cj_nk[1] = volt0_ci_cj_ck[2];
		volt0_ci_cj_nk[2] = volt0_ci_cj_ck[3];

		Simd<float, veclen> volt1_ci_cj_nk;
		volt1_ci_cj_nk[0] = volt1_ci_cj_ck[1];
		volt1_ci_cj_nk[1] = volt1_ci_cj_ck[2];
		volt1_ci_cj_nk[2] = volt1_ci_cj_ck[3];

		// "boundary" is a templated Boolean constant, no branching occurs
		if (boundary)
		{
			// outside boundary
			volt0_ci_cj_nk[3] = 0;
			volt1_ci_cj_nk[3] = 0;
		}
		else
		{
			volt0_ci_cj_nk[3] = volt_v(0, i, j, vk + 1)[0];
			volt1_ci_cj_nk[3] = volt_v(1, i, j, vk + 1)[0];
		}

		// 24 (6 x 4) FP32 loads
		Simd<float, veclen> ii0_ci_cj_ck = op.ii_v(0, i, j, vk);
		Simd<float, veclen> ii1_ci_cj_ck = op.ii_v(1, i, j, vk);
		Simd<float, veclen> ii2_ci_cj_ck = op.ii_v(2, i, j, vk);
		Simd<float, veclen> iv0_ci_cj_ck = op.iv_v(0, i, j, vk);
		Simd<float, veclen> iv1_ci_cj_ck = op.iv_v(1, i, j, vk);
		Simd<float, veclen> iv2_ci_cj_ck = op.iv_v(2, i, j, vk);

		// x-polarization
		curr0_ci_cj_ck *= ii0_ci_cj_ck;
		curr0_ci_cj_ck +=
			iv0_ci_cj_ck * (
				volt2_ci_cj_ck -
				volt2_ci_nj_ck -
				volt1_ci_cj_ck +
				volt1_ci_cj_nk
			);

		// y-polarization
		curr1_ci_cj_ck *= ii1_ci_cj_ck;
		curr1_ci_cj_ck +=
			iv1_ci_cj_ck * (
				volt0_ci_cj_ck -
				volt0_ci_cj_nk -
				volt2_ci_cj_ck +
				volt2_ni_cj_ck
			);

		// z-polarization
		curr2_ci_cj_ck *= ii2_ci_cj_ck;
		curr2_ci_cj_ck +=
			iv2_ci_cj_ck * (
				volt1_ci_cj_ck -
				volt1_ni_cj_ck -
				volt0_ci_cj_ck +
				volt0_ci_nj_ck
			);

		// 12 (3 x 4) FP32 stores
		curr_v(0, i, j, vk) = curr0_ci_cj_ck;
		curr_v(1, i, j, vk) = curr1_ci_cj_ck;
		curr_v(2, i, j, vk) = curr2_ci_cj_ck;
	}
}

void Engine_Tiling::DoPreVoltageUpdates(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
		m_Eng_exts.at(n)->DoPreVoltageUpdates(timestep, range);
}

void Engine_Tiling::DoPostVoltageUpdates(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->DoPostVoltageUpdates(timestep, range);
}

void Engine_Tiling::Apply2Voltages(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->Apply2Voltages(timestep, range);
}

void Engine_Tiling::DoPreCurrentUpdates(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
		m_Eng_exts.at(n)->DoPreCurrentUpdates(timestep, range);
}

void Engine_Tiling::DoPostCurrentUpdates(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->DoPostCurrentUpdates(timestep, range);
}

void Engine_Tiling::Apply2Current(int timestep, const Tiling::Range3D<> range)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->Apply2Current(timestep, range);
}

void Engine_Tiling::Reset()
{
	Engine::Reset();

	delete m_volt_v;
	delete m_curr_v;
	m_volt_v = NULL;
	m_curr_v = NULL;
}

Engine_Tiling::~Engine_Tiling()
{
	Reset();
}
