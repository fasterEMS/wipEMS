#ifndef ARRAYLIB_ALLOCATOR_H
#define ARRAYLIB_ALLOCATOR_H

#include <cstddef>

namespace ArrayLib
{
	template <typename T>
	class SimpleAllocator;

	template <typename T>
	class AlignedAllocator;
};

template <typename T>
class ArrayLib::SimpleAllocator
{
public:
	static T* alloc(size_t numelem)
	{
		T* buf = malloc(numelem * sizeof(T));
		if (buf == NULL)
		{
			std::cerr << "Failed to allocate aligned memory" << std::endl;
			throw std::bad_alloc();
		}
		memset(buf, 0, numelem * sizeof(T));
		return new (buf) T[numelem];
	}

	static void free(T* ptr)
	{
		if (ptr)
		{
			ptr->~T();
			free(ptr);
		}
	}
};

#ifdef WIN32
#include <malloc.h>
#endif

template <typename T>
class ArrayLib::AlignedAllocator
{
public:
	static T* alloc(size_t numelem)
	{
		// Align to 1 MiB, a common hugepage size on many systems.
		// It allows us to enable hugepage in the future. It's also
		// a (large) multiple of all SIMD data types.
		constexpr size_t alignment = 1 * 1024 * 1024;

		T* buf;
#ifdef WIN32
		buf = _mm_malloc(numelem * sizeof(T), alignment);
		if (buf == NULL)
		{
			std::cerr << "Failed to allocate aligned memory" << std::endl;
			throw std::bad_alloc();
		}
#else
		int retval = posix_memalign((void**) &buf, alignment, numelem * sizeof(T));
		if (retval != 0)
		{
			std::cerr << "Failed to allocate aligned memory" << std::endl;
			throw std::bad_alloc();
		}
#endif
		memset(buf, 0, numelem * sizeof(T));
		return new (buf) T[numelem];

	}

	static void free(T* ptr)
	{
		if (ptr)
		{
			ptr->~T();
#ifdef WIN32
			_mm_free(ptr);
#else
			free(ptr);
#endif
		}

	}
};

#endif // ARRAYLIB_ALLOCATOR_H
