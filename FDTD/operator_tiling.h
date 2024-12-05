/*
 * Copyright (C) 2024 Yifeng Li (tomli@tomli.me)
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
#ifndef OPERATOR_TILING_H
#define OPERATOR_TILING_H

#include <cstdint>
#include <boost/program_options.hpp>
#include "operator.h"
#include "engine_tiling.h"
#include "tools/global.h"
#include "tools/arraylib/array_nijk.h"
#include "tools/arraylib/simd.h"
#include "tools/tiling/tiling.h"

using ArrayLib::Simd;
constexpr static const size_t veclen = 4;

class Operator_Tiling : public Operator
{
	friend class Engine_Interface_FDTD;

public:
	static boost::program_options::options_description optionDesc();

	//! Create a new operator
	static Operator_Tiling* New(uint32_t numThreads=0);
	//static Operator_Tiling* New();
	virtual ~Operator_Tiling();

	virtual Engine* CreateEngine();

	// Scalar getters and setters APIs, used by the Operator base class
	// for simulation setups, so they must always be implemented.
	virtual inline
	float GetVV(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const
	{
		return GetOpr(OP_TYPE::VV, n, i, j, k);
	}

	virtual inline
	float GetVI(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const
	{
		return GetOpr(OP_TYPE::VI, n, i, j, k);
	}

	virtual inline
	float GetII(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const
	{
		return GetOpr(OP_TYPE::II, n, i, j, k);
	}

	virtual inline
	float GetIV(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const
	{
		return GetOpr(OP_TYPE::IV, n, i, j, k);
	}

	virtual inline
	void SetVV(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val)
	{
		return SetOpr(OP_TYPE::VV, n, i, j, k, val);
	}

	virtual inline
	void SetVI(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val)
	{
		return SetOpr(OP_TYPE::VI, n, i, j, k, val);
	}

	virtual inline
	void SetII(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val)
	{
		return SetOpr(OP_TYPE::II, n, i, j, k, val);
	}

	virtual inline
	void SetIV(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val)
	{
		return SetOpr(OP_TYPE::IV, n, i, j, k, val);
	}

	// Vector getters and setters, used by the Tiling engine internally
	inline
	const Simd<float, veclen>& vv_v(uint32_t n, uint32_t i, uint32_t j, uint32_t vk) const
	{
		return GetOprVec(OP_TYPE::VV, n, i, j, vk);
	}

	inline
	const Simd<float, veclen>& vi_v(uint32_t n, uint32_t i, uint32_t j, uint32_t vk) const
	{
		return GetOprVec(OP_TYPE::VI, n, i, j, vk);
	}

	inline
	const Simd<float, veclen>& iv_v(uint32_t n, uint32_t i, uint32_t j, uint32_t vk) const
	{
		return GetOprVec(OP_TYPE::IV, n, i, j, vk);
	}

	inline
	const Simd<float, veclen>& ii_v(uint32_t n, uint32_t i, uint32_t j, uint32_t vk) const
	{
		return GetOprVec(OP_TYPE::II, n, i, j, vk);
	}

	std::vector<Tiling::Plan3D> m_perNodePlan;
	unsigned int m_tileHalfTs;

protected:
	//! use New() for creating a new Operator
	Operator_Tiling();

	virtual void Init();
	void Delete();
	virtual void Reset();
	virtual void setNumThreads(unsigned int numThreads);
	virtual void InitOperator();
	virtual int CalcECOperator(DebugFlags debugFlags=None);

	void CompressOperators(
		std::vector<Tiling::Plan3D> perNodePlan,
		uint32_t oprType
	);

	Tiling::Plan3D makePlan(size_t tileHalfTs);
	virtual void applyOptions();

	std::array<int, 3>  m_tileSize;
	std::array<char, 3> m_tileType;
	bool m_numa_enable = false;
	unsigned int m_numa_node = 0;
	unsigned int m_numThreads;

	enum OP_TYPE{
		VV = 0,
		VI = 1,
		II = 2,
		IV = 3,
		LAST = 3,
	};

private:
	// underlying implementation of the scalar getter
	inline
	float GetOpr(uint32_t type, uint32_t n, uint32_t i, uint32_t j, uint32_t k) const
	{

		uint32_t vecIdx    = k / veclen;
		uint32_t vecOffset = k % veclen;

		if (m_compOprArr[type] != NULL) {
			ArrayLib::ArrayNIJK<uint32_t, uint32_t, 4>& index = *m_index;

			uint32_t idx = index(OP_TYPE::VV, i, j, vecIdx);
			return m_compOprArr[type][idx * 3 + n][vecOffset];
		}
		else {
			ArrayLib::ArrayNIJK<Simd<float, veclen>>& oprArr = *m_rawOprArr[type];
			return oprArr(n, i, j, vecIdx)[vecOffset];
		}
	}

	// underlying implementation of the scalar setter
	inline
	void SetOpr(uint32_t type, uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val) const
	{
		if (m_compOprArr[type] != NULL) {
			throw std::runtime_error(
				"operators are already compressed and read-only "
			);
		}
		else {
			uint32_t vecIdx    = k / veclen;
			uint32_t vecOffset = k % veclen;

			ArrayLib::ArrayNIJK<Simd<float, veclen>>& oprArr = *m_rawOprArr[type];
			oprArr(n, i, j, vecIdx)[vecOffset] = val;
		}
	}

	// underlying implementation of the vector getter/setter
	inline
	const Simd<float, veclen>& GetOprVec(uint32_t type, uint32_t n, uint32_t i, uint32_t j, uint32_t vk) const
	{
		ArrayLib::ArrayNIJK<uint32_t, uint32_t, 4>& index = *m_index;

		uint32_t idx = index(type, i, j, vk);
		return m_compOprArr[type][idx * 3 + n];
	}

	// raw uncompressed operators, used during setup
	std::array<
		// 4D array of operator values.
		// To get a value, use the standard FDTD coordinates.
		ArrayLib::ArrayNIJK<Simd<float, veclen>>*,
		4
	> m_rawOprArr; // vv, vi, iv, ii

	// deduplicated operators stored in a compressed format,
	// used during actual simulation
	std::array<
		// 1D array of operator values
		// To get a value, use the index lookup table.
		Simd<float, veclen>*,
		4
	> m_compOprArr; // vv, vi, iv, ii

	// index lookup table
	// Input is FDTD coordinates (i, j, vk, OP_TYPE)
	// Output: an unique index
	ArrayLib::ArrayNIJK<
		uint32_t,
		uint32_t,
		4
	>* m_index; // vv, vi, iv, ii

	unsigned int numVectors;
};

#endif // OPERATOR_TILING_H
