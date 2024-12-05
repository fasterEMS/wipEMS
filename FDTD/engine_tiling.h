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

#ifndef ENGINE_TILING_H
#define ENGINE_TILING_H

#include "engine.h"
#include "extensions/engine_extension.h"
#include "operator_tiling.h"
#include "tools/array_ops.h"
#include "tools/arraylib/simd.h"
#include "tools/tiling/tiling.h"

using ArrayLib::Simd;

static const uint32_t veclen = 4;

class Operator_Tiling;

class Engine_Tiling : public Engine
{
public:
	static Engine_Tiling* New(const Operator_Tiling* op);
	virtual ~Engine_Tiling();

	virtual void Init();
	virtual void Reset();

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

	inline virtual float
	GetVolt(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const override
	{
		ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v = *m_volt_v;

		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;
		return volt_v(n, i, j, k / veclen)[vk_offset];
	}

	inline virtual float
	GetVolt(uint32_t n, const uint32_t pos[3]) const override
	{
		return GetVolt(n, pos[0], pos[1], pos[2]);
	}

	inline virtual float
	GetCurr(uint32_t n, uint32_t i, uint32_t j, uint32_t k) const override
	{
		ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v = *m_curr_v;

		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;
		return curr_v(n, i, j, vk)[vk_offset];
	}

	inline virtual float
	GetCurr(uint32_t n, const uint32_t pos[3]) const override
	{
		return GetCurr(n, pos[0], pos[1], pos[2]);
	}

	inline virtual void
	SetVolt(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val) override
	{
		ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v = *m_volt_v;

		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;
		volt_v(n, i, j, vk)[vk_offset] = val;
	}

	inline virtual void
	SetVolt(uint32_t n, const uint32_t pos[3], float val) override
	{
		SetVolt(n, pos[0], pos[1], pos[2], val);
	}

	inline virtual void
	SetCurr(uint32_t n, uint32_t i, uint32_t j, uint32_t k, float val) override
	{
		ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v = *m_curr_v;

		uint32_t vk = k / veclen;
		uint32_t vk_offset = k % veclen;
		curr_v(n, i, j, vk)[vk_offset] = val;
	}

	inline virtual void
	SetCurr(uint32_t n, const uint32_t pos[3], float val) override
	{
		SetCurr(n, pos[0], pos[1], pos[2], val);
	}

	virtual bool IterateTS(unsigned int iterTS);

protected:
	void InitializeTiling(const Tiling::Plan3D plan);

	void DoPreVoltageUpdates(int timestep, const Tiling::Range3D<> range);
	void DoPostVoltageUpdates(int timestep, const Tiling::Range3D<> range);
	void Apply2Voltages(int timestep, const Tiling::Range3D<> range);
	void DoPreCurrentUpdates(int timestep, const Tiling::Range3D<> range);
	void DoPostCurrentUpdates(int timestep, const Tiling::Range3D<> range);
	void Apply2Current(int timestep, const Tiling::Range3D<> range);

	Engine_Tiling(const Operator_Tiling* op);
	const Operator_Tiling* Op;

	// Benchmark showed that size_t indices is 1.5x faster than
	// uint32_t indices on 64-bit machines, likely due to reasons
	// related to compiler codegen.
	ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>* m_volt_v;
	ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>* m_curr_v;

private:
	inline void
	updateVoltageScalarKernel(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const Operator_Tiling& op,
		uint32_t i,
		uint32_t j,
		uint32_t first_k,
		uint32_t last_k
	);

	template <bool boundary>
	inline void
	updateVoltageVectorKernel(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const Operator_Tiling& op,
		uint32_t i,
		uint32_t j,
		uint32_t first_k,
		uint32_t last_k
	);

	void updateVoltageRange(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const Operator_Tiling& op,
		const Tiling::Range3D<> range
	);

	inline void
	updateCurrentScalarKernel(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const Operator_Tiling& op,
		uint32_t i,
		uint32_t j,
		uint32_t first_k,
		uint32_t last_k
	);

	template <bool boundary>
	inline void
	updateCurrentVectorKernel(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const Operator_Tiling& op,
		uint32_t i,
		uint32_t j,
		uint32_t first_k,
		uint32_t last_k
	);

	void updateCurrentRange(
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& curr_v,
		const ArrayLib::ArrayNIJK<Simd<float, 4>, size_t>& volt_v,
		const Operator_Tiling& op,
		const Tiling::Range3D<> range
	);
};

#endif // ENGINE_TILING_H
