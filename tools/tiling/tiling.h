/*
 * BSD Zero Clause License
 *
 * Copyright (C) 2024 Yifeng Li
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

// tiling - FDTD temporial tiling algorithm implementation.
//
// FDTD simulation is a form of stencil computation, its order
// of calculations can be changed in such a way to partition
// the simulation domain into multiple small regions, in which
// multiple timesteps can be calculated at a time. It increases
// simulation speed significantly due to better cache locality.
// This technique is known as tiling or blocking. More
// specifically, temporal tiling.
//
// Correct tiling is tricky to implement and test, so the plan,
// thus, tiling.h is a library for pre-calculate the necessary
// coordinates of all regions according to the parallelogram and
// trapezoid tiling algorithms.
//
// This is used by the Tiling engine to generate a "plan" before
// starting the simulation.
//
// For more information, see:
//
// * Fukaya, T., & Iwashita, T. (2018). Time-space tiling with
//   tile-level parallelism for the 3D FDTD method. Proceedings of
//   the International Conference on High Performance Computing
//   in Asia-Pacific Region - HPC Asia 2018.
//
//   doi: 10.1145/3149457.3149478
//
// * The original implementation of this library, with notes and
//   diagnostic tools: https://github.com/biergaizi/project-diamond

#ifndef TILING_H
#define TILING_H

#include <cstdint>
#include <cstring>
#include <array>
#include <vector>

namespace Tiling {
	using size_t = std::size_t;

	template <typename T>
	struct WrappedVector
	{
	public:
		using iterator = typename std::vector<T>::iterator;
		using const_iterator = typename std::vector<T>::const_iterator;
		iterator       begin()        { return m_vector.begin();    }
		iterator       end()          { return m_vector.end();      }
		const_iterator begin()  const { return m_vector.begin();    }
		const_iterator end()    const { return m_vector.end();      }
		const_iterator cbegin() const { return m_vector.cbegin();   }
		const_iterator cend()   const { return m_vector.cend();     }
		size_t         size()   const { return m_vector.size();     }

		void     reserve(size_t numelem)       { m_vector.reserve(numelem); }
		void     push_back(T elem)             { m_vector.push_back(elem);  }
		T&       operator[] (size_t idx)       { return m_vector[idx];      }
		const T& operator[] (size_t idx) const { return m_vector[idx];      }

	protected:
		std::vector<T> m_vector;
	};

	template <typename T>
	struct Range1D
	{
		T first, last;
	};

	struct Tile1D : WrappedVector<Range1D<size_t>>
	{
	public:
		Tile1D(size_t id) : m_id(id) {}
		size_t id() const { return m_id; }

	private:
		size_t m_id;
	};

	using TileList1D = std::vector<Tile1D>;

	struct Plan1D : WrappedVector<TileList1D>
	{
	public:
		Plan1D(size_t elems)
		{
			m_vector.resize(elems);
		}

		size_t numTiles() const
		{
			size_t retval = 0;
			for (const TileList1D& tileList : m_vector) {
				retval += tileList.size();
			}

			return retval;
		}
	};

	Plan1D
	computeParallelogramTiles(
		size_t totalWidth, size_t tileWidth,
		size_t halfTimesteps
	);

	Plan1D
	computeTrapezoidTiles(
		size_t totalWidth, size_t tileWidth,
		size_t halfTimesteps
	);

	template <typename T=size_t>
	struct Range3D
	{
		std::array<T, 3> first;
		std::array<T, 3> last;

		bool operator== (const Range3D& other) const
		{
			return this->first == other.first && this->last == other.last;
		}
	};

	struct Subtile3D : WrappedVector<Range3D<size_t>>
	{
	public:
		Subtile3D() {}
		Subtile3D(size_t id) : m_id(id) {}
		size_t id() const { return m_id; }

		void push_back(Range3D<size_t> range)
		{
			m_vector.push_back(range);
			for (size_t i = 0; i < 3; i++) {
				first[i] = std::min(range.first[i], first[i]);
				 last[i] = std::max(range.last[i],   last[i]);
			}
		}

		std::array<size_t, 3> first = {SIZE_MAX, SIZE_MAX, SIZE_MAX};
		std::array<size_t, 3> last  = {0, 0, 0};

	private:
		size_t m_id;
	};

	struct Tile3D : WrappedVector<Subtile3D>
	{
	public:
		Tile3D(std::array<size_t, 3> id) : m_id(id) {}
		std::array<size_t, 3> id() const { return m_id; }

	private:
		std::array<size_t, 3> m_id;
	};

	using TileList3D = std::vector<Tile3D>;

	struct Plan3D : WrappedVector<TileList3D>
	{
		Plan3D() {}

		Plan3D(std::array<size_t, 3> numTiles, size_t numStages)
		{
			m_numTiles = numTiles;
			m_vector.resize(numStages);
		}

		std::array<size_t, 3> numTiles() { return m_numTiles; }

	private:
		std::array<size_t, 3> m_numTiles;
	};

	Plan3D
	combineTilesTTT(const Plan1D& i, const Plan1D& j, const Plan1D& k);

	Plan3D
	combineTilesTTP(const Plan1D& i, const Plan1D& j, const Plan1D& k);

	Plan3D
	toLocalCoords(Plan3D plan);

	std::vector<Plan3D>
	divideNumaWorkload(Plan3D& plan, size_t numSocket);

	void visualizeTiles(
		const Plan1D& plan,
		size_t totalWidth, size_t tileWidth,
		size_t halfTimesteps
	);

	template <typename T=size_t>
	class Range3DHasher
	{
	public:
		size_t operator() (const Range3D<T>& range) const
		{
			// Copy the raw bytes within a Range3D<T> into a
			// std::u32string buffer, which is hashable by
			// std::hash.
			std::string bytestream;
			bytestream.resize(sizeof(range));
			std::memcpy(&bytestream[0], &range, sizeof(range));
			return std::hash<std::string>{}(bytestream);;
		}
	};
}

// hack: make it header-only for quick tests without build systems
#ifdef TILING_HEADER_ONLY
#include "tiling.cpp"
#endif

#endif // TILING_H
