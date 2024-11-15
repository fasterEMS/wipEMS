#ifndef ARRAYLIB_SUBSCRIPT_H
#define ARRAYLIB_SUBSCRIPT_H

#include <cstddef>
#include <array>
#include <stdexcept>

// Subscript is a wrapper class that enables access an ArrayType& via
// arr[i][j][k] syntax, with one index per operator[]. This is implemented
// via template metaprogramming to either recursively return another
// Subscript, or return an actual reference T& when recursion reaches
// the base case.
//
// Performance of nested operator[] should be the same as as flat operator()
// (the last time I've check, with identical assembly code generation even
// on MSVC).
namespace ArrayLib
{
	template <typename ArrayType, size_t maxDim, size_t dim>
	class Subscript;
};

template <typename ArrayType, size_t maxDim, size_t dim>
class ArrayLib::Subscript
{
public:
	Subscript(ArrayType& array) : m_array(array) {}

	Subscript(
		const ArrayType& array,
		std::array<typename ArrayType::IndexType, maxDim> sub
	) :
		m_array(array),
		m_sub(sub)
	{
	}

	template <typename ArrayType_=ArrayType> // dummy param for enable_if
	typename std::enable_if<
		/* enable when */ dim + 1 < maxDim,
		/* return type */ Subscript<ArrayType_, maxDim, dim + 1>
	>::type
	operator[] (size_t nextTupleIdx)
	{
		m_sub[dim] = nextTupleIdx;
		return Subscript<ArrayType, maxDim, dim + 1>(m_array, m_sub);
	}

	template <typename ArrayType_=ArrayType> // dummy param for enable_if
	typename std::enable_if<
		/* enable when */ dim + 1 == maxDim,
		/* return type */ typename ArrayType_::ValueType&
	>::type
	operator[] (size_t finalTupleIdx)
	{
		m_sub[dim] = finalTupleIdx;
		return m_array(m_sub);
	}

	// catch common misuses
	void operator=(typename ArrayType::ValueType val)
	{
		(void) val;

		// equivalent to static_assert(false) with friendly error messages
		static_assert(
			maxDim == dim, // evaluates to false
			"Too few [] operators are given to a multi-dimensional array. "
			"Indexing a pointer without dereference can cause the same error."
			" (e.g. ptr[x][y][z] should be (*ptr)[x][y][z])"
		);
	}

	operator typename ArrayType::ValueType() const
	{
		// equivalent to static_assert(false) with friendly error messages
		static_assert(
			maxDim == dim, // evaluates to false
			"Too few [] operators are given to a multi-dimensional array. "
			"Indexing a pointer without dereference can cause the same error "
			" (e.g. ptr[x][y][z] should be (*ptr)[x][y][z])"
		);
		throw std::logic_error("bad template usage");
	}

private:
	const ArrayType& m_array;
	std::array<typename ArrayType::IndexType, maxDim> m_sub;
};

#endif // ARRAYLIB_SUBSCRIPT_H
