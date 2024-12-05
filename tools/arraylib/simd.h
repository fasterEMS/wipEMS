/*
 * Copyright (C) 2024 Yifeng Li (tomli@tomli.me)
 * Copyright (C) 2010, 2019 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

// simd.h: new SIMD wrapper class
//
// Class simd.h is a new wrapper for manipulating SIMD/SSE vector
// data. In comparison to the original f4vector class, it has two
// advantages:
// 
// 1. It's templated so that one can create a vector of arbitrary data
// type and length, such as:
// 
//     Simd<int, 4> v1, v2;
//     Sime<float, 8> v3;
// 
// 2. Individual elements can be accessed by the natural operator[]
// syntax, such as:
// 
//     v1[2]
//     v2[7]
// 
// 3. The underlying union is private, operations are applied directly
// on a vector via overloaded operators (without the member variable
// ".v"), such as:
// 
//     v1 += v2;
//
// 4. Implement operator== using memcmp, so manual vector comparison
// is not necessary.
// 
// 5. It relies on vector extension in GCC/clang. Manual hardcoded
// fallbacks for each data size and type are needed for MSVC. Only
// float4 is implemnted for now. Since clang can be used as the MSVC
// complier, it's probably a good idea to make the switch?

#ifndef ARRAYLIB_SIMD_H
#define ARRAYLIB_SIMD_H

#include <cstddef>
#include <array>

#ifndef __GNUC__
#include <emmintrin.h>
#endif

namespace ArrayLib
{
	template <typename Tsca_, size_t len_>
	class Simd;

	using float4 = Simd<float, 4>;
};

template <typename Tsca_, size_t len_>
class ArrayLib::Simd
{
public:
	constexpr static size_t len()	{ return len_;                  }
	constexpr static size_t bytes() { return sizeof(Tsca_) * len_;  }

	typedef Tsca_ Tsca;

#ifdef __GNUC__
	// Define Tvec to be a GCC vector derived from Tsca, GCC
	// automatically generates SIMD code on multiple architectures
	typedef Tsca Tvec __attribute__ ((vector_size(bytes())));
#else
	// Fallback for MSVC
	static_assert(
		std::is_same<Tsca, float>::value && len_ == 4,
		"Only float4 vector is implemented on MSVC."
	);

	using Tvec =
		typename std::enable_if<
			std::is_same<Tsca, float>::value && len_ == 4,
			__m128
		>::type;
#endif

	Simd()
	{
	}

	Simd(std::array<Tsca, len_> val)
	{
		for (size_t i = 0; i < len_; i++) {
			elem[i] = val[i];
		}
	}

#ifdef __GNUC__
	// Enable only on GCC. On MSVC/x86, initialization is ambiguous
	// since both std:array and __m128 can be brace-initialized
	Simd(Tvec v)
	{
		all = v;
	}
#endif

	Tsca&		operator[] (size_t idx)         { return elem[idx]; }
	const Tsca&	operator[] (size_t idx) const   { return elem[idx]; }

	bool operator== (const Simd<Tsca, len_>& that) const
	{
		return memcmp(this, &that, bytes()) == 0;
	}

	void operator+= (const Simd<Tsca, len_>& that)
	{
#ifdef __GNUC__
		this->all += that.all;
#else
		// one "if" is always true on MSVC, no actual braching occurs
		if (std::is_same<Tvec, __m128>::value)
			this->all = _mm_add_ps(this->all, that.all);
		//if (std::is_same<Tvec, __m256>::value)
		//      TODO: support more types in the future
#endif
	}

	void operator-= (const Simd<Tsca, len_>& that)
	{
#ifdef __GNUC__
		this->all -= that.all;
#else
		// one "if" is always true on MSVC, no actual braching occurs
		if (std::is_same<Tvec, __m128>::value)
			this->all = _mm_sub_ps(this->all, that.all);
#endif
	}

	void operator*= (const Simd<Tsca, len_>& that)
	{
#ifdef __GNUC__
		this->all *= that.all;
#else
		// one "if" is always true on MSVC, no actual braching occurs
		if (std::is_same<Tvec, __m128>::value)
			this->all = _mm_mul_ps(this->all, that.all);
#endif
	}

	void operator/= (const Simd<Tsca, len_>& that)
	{
#ifdef __GNUC__
		this->all /= that.all;
#else
		// one "if" is always true on MSVC, no actual braching occurs
		if (std::is_same<Tvec, __m128>::value)
			this->all = _mm_div_ps(this->all, that.all);
#endif
	}

private:
	// Type punning via union is undefined by the C++ standard, but both
	// GCC's vector extension and MSVC's SSE/AVX intrinsic code uses it to
	// access elements in a vector.
	union
	{
		Tsca elem[len_];
		Tvec all;
	};
};

template <typename Tsca, size_t len>
ArrayLib::Simd<Tsca, len> operator+ (
	const ArrayLib::Simd<Tsca, len>& a,
	const ArrayLib::Simd<Tsca, len>& b
)
{
	ArrayLib::Simd<Tsca, len> tmp = a;
	tmp += b;
	return tmp;
}

template <typename Tsca, size_t len>
ArrayLib::Simd<Tsca, len> operator- (
	const ArrayLib::Simd<Tsca, len>& a,
	const ArrayLib::Simd<Tsca, len>& b
)
{
	ArrayLib::Simd<Tsca, len> tmp = a;
	tmp -= b;
	return tmp;
}

template <typename Tsca, size_t len>
ArrayLib::Simd<Tsca, len> operator* (
	const ArrayLib::Simd<Tsca, len>& a,
	const ArrayLib::Simd<Tsca, len>& b
)
{
	ArrayLib::Simd<Tsca, len> tmp = a;
	tmp *= b;
	return tmp;
}

template <typename Tsca, size_t len>
ArrayLib::Simd<Tsca, len> operator/ (
	const ArrayLib::Simd<Tsca, len>& a,
	const ArrayLib::Simd<Tsca, len>& b
)
{
	ArrayLib::Simd<Tsca, len> tmp = a;
	tmp /= b;
	return tmp;
}

#endif // ARRAYLIB_SIMD_H
