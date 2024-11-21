#ifndef ARRAYLIB_ARRAY_NIJK_H
#define ARRAYLIB_ARRAY_NIJK_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>

#include "impl/array_base.h"

namespace ArrayLib
{
	template <typename T, typename IndexType=uint32_t, size_t extentN=3>
	class ArrayNIJK;
};

template <typename T, typename IndexType, size_t extentN>
class ArrayLib::ArrayNIJK :
	public ArrayBase<ArrayNIJK<T, IndexType>, T, 4, IndexType>
{
private:
	std::array<IndexType, 2> m_stride;

public:
	using Base = ArrayBase<ArrayNIJK<T, IndexType>, T, 4, IndexType>;
	using Base::operator();

	// 2-phase initialization: user declares a dummy array object first,
	// and calls Init() later. Ugly but needed because openEMS's FDTD
	// engine itself uses the 2-phase initialization pattern, we can't
	// determine simulation domain size in the constructor.
	ArrayNIJK() {}

	void Init(std::string name, std::array<IndexType, 3> extent)
	{
		if (this->m_ptr != NULL)
			Base::AllocatorType::free(this->m_ptr, this->m_size);

		this->m_name = name;
		this->m_size = extent[0] * extent[1] * extent[2] * extentN;
		this->m_bytes = sizeof(T) * this->m_size;
		this->m_ptr = Base::AllocatorType::alloc(this->m_size);

		this->m_extent = {extentN, extent[0], extent[1], extent[2]};
		this->m_stride[0] = extent[1] * extent[2] * extentN;
		this->m_stride[1] = extent[2] * extentN;
	}

	void Init(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	// standard RAII initialization
	ArrayNIJK(std::string name, std::array<IndexType, 3> extent)
	{
		Init(name, extent);
	}

	ArrayNIJK(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	// Access syntax is (n, i, j, k) but its memory layout is (i, j, k, n),
	// in which n is the stride-1 dimension. In openEMS, performance is
	// 20% higher.
	IndexType linearIndex(std::array<IndexType, 4> tupleIndex) const
	{
		return tupleIndex[1] * m_stride[0] + \
		       tupleIndex[2] * m_stride[1] +
			   tupleIndex[3] * extentN +
			   tupleIndex[0];
	}

	T& operator() (IndexType n, IndexType i, IndexType j, IndexType k) const
	{
		return this->m_ptr[linearIndex({n, i, j, k})];
	}

	~ArrayNIJK()
	{
		Base::AllocatorType::free(this->m_ptr, this->m_size);
	}
};

#endif // ARRAYLIB_ARRAY_NIJK_H
