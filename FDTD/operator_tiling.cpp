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
#include "engine_tiling.h"
#include "operator_tiling.h"
#include "tools/tiling/tbb/paraexec.h"
//#include "tools/tiling/tbb/first_touch.h"

Operator_Tiling* Operator_Tiling::New(unsigned int numThreads)
//Operator_Tiling* Operator_Tiling::New()
{
	std::cout << "Create FDTD operator ";
	std::cout << "(compressed SSE + multi-threading + spatial/temporal tiling)";
	std::cout << std::endl;

	Operator_Tiling* op = new Operator_Tiling();
	op->setNumThreads(numThreads);
	op->Init();
	return op;
}

Engine* Operator_Tiling::CreateEngine()
{
	m_Engine = Engine_Tiling::New(this);
	return m_Engine;
}

Operator_Tiling::Operator_Tiling() : Operator()
{
	for (uint32_t oprType = 0; oprType <= OP_TYPE::LAST; oprType++)
	{
		m_rawOprArr[oprType] = NULL;
		m_compOprArr[oprType] = NULL;
	}
	m_index = NULL;
}

void Operator_Tiling::setNumThreads(unsigned int numThreads)
{
	m_numThreads = numThreads;
}

void Operator_Tiling::Init()
{
	Delete();
	Operator::Init();

	applyOptions();

	for (uint32_t oprType = 0; oprType <= OP_TYPE::LAST; oprType++)
	{
		m_rawOprArr[oprType] = NULL;
		m_compOprArr[oprType] = NULL;
	}
	m_index = NULL;

	if (!m_numa_enable)
	    ParaExec::init(m_numThreads, false, false);
	else if (m_numa_node == 0)
	    ParaExec::init(m_numThreads, true, true);
	else if (m_numa_node == 1)
	    ParaExec::init(m_numThreads, false, true);

	ParaExec::printSystemInfo();
}

namespace po = boost::program_options;

po::options_description
Operator_Tiling::optionDesc()
{
	po::options_description optdesc(
		"Tiling engine options (needs --engine=tiling)"
	);
	optdesc.add_options()
		("tile-size",   po::value<std::string>()->default_value("17t,17t,200p"),
		                "Work unit size for Tiling engine \n"
		                "Syntax: it,jt,kt or it,jt,kp, suffix 't' "
		                "denotes trazepoid tiling, suffix 'p' denotes "
		                "parallelogram tiling, e.g. 20t,20t,200t or "
		                "10t,10t,200p. Optimal tile size is crucial "
		                "for performance but difficult to determine. "
		                "Smaller size increases CPU cache hits, "
		                "bigger size increases parallelism. See "
		                "manual for tuning guidelines.")
		("tile-height", po::value<unsigned int>()->default_value(14),
		                "Half timestep batch size for Tiling engine "
		                "Must be even and smaller than specified tile "
		                "size's smallest dimensions. See manual.")
		("numa",        "Enable NUMA or L3 cache awareness. See manual")
		("numa-node",   po::value<unsigned int>()->default_value(0),
		                "Enable NUMA or L3 cache awareness. See manual");

	return optdesc;
}

void Operator_Tiling::applyOptions()
{
	if (g_settings.hasOption("tile-size"))
	{
		std::array<std::string, 3> tileArgString;
		char* tileArgC = strdup(
			g_settings.getOption("tile-size").as<std::string>().c_str()
		);
		if (!tileArgC)
			throw std::bad_alloc();

		tileArgString[0] = strtok(tileArgC, ",");
		tileArgString[1] = strtok(NULL, ",");
		tileArgString[2] = strtok(NULL, ",");

		for (size_t dim = 0; dim < 3; dim++)
		{
			std::string& arg = tileArgString[dim];

			if (arg[arg.size() - 1] != 't' && arg[arg.size() - 1] != 'p')
			{
				throw std::invalid_argument(
					std::string("tile suffix must be 't' or 'p', got {}") +
					            arg[arg.size() - 1]
				);
			}

			m_tileType[dim] = arg[arg.size() - 1];
			arg[arg.size() - 1] = '\0';
			m_tileSize[dim] = atoi(arg.c_str());
		}

		if (m_tileType[0] != 't' || m_tileType[1] != 't')
		{
			throw std::invalid_argument(
				"dimension i and j only support trapezoid tiling (suffix t)"
			);
		}

		free(tileArgC);
	}

	if (g_settings.hasOption("tile-height"))
		m_tileHalfTs = g_settings.getOption("tile-height").as<unsigned int>();

	if (g_settings.hasOption("numa"))
		m_numa_enable = true;
	m_numa_node = g_settings.getOption("numa-node").as<unsigned int>();

	if (!m_numa_enable && m_numa_node)
		throw std::runtime_error("To use --numa-node, --numa must be enabled.");
}

void Operator_Tiling::InitOperator()
{
	Delete();

	numVectors = ceil((double) numLines[2] / (double) veclen);

	m_rawOprArr[OP_TYPE::VV] = new ArrayLib::ArrayNIJK<Simd<float, 4>>(
		"vv", {numLines[0], numLines[1], numVectors}
	);
	m_rawOprArr[OP_TYPE::VI] = new ArrayLib::ArrayNIJK<Simd<float, 4>>(
		"vi", {numLines[0], numLines[1], numVectors}
	);
	m_rawOprArr[OP_TYPE::IV] = new ArrayLib::ArrayNIJK<Simd<float, 4>>(
		"iv", {numLines[0], numLines[1], numVectors}
	);
	m_rawOprArr[OP_TYPE::II] = new ArrayLib::ArrayNIJK<Simd<float, 4>>(
		"ii", {numLines[0], numLines[1], numVectors}
	);
	m_index = new ArrayLib::ArrayNIJK<uint32_t, uint32_t, 4>(
		"index", {numLines[0], numLines[1], numVectors}
	);

	Tiling::Plan3D plan = makePlan(m_tileHalfTs);
	m_perNodePlan = divideNumaWorkload(plan, ParaExec::numNodes());
}

Tiling::Plan3D Operator_Tiling::makePlan(size_t tileHalfTs)
{
	std::array<uint32_t, 3> gridSize = {
		numLines[0], numLines[1], numLines[2]
	};

	printf("tile\t\t" "%04u x %04u x %04u\n",
		   m_tileSize[0], m_tileSize[1], m_tileSize[2]);
	printf("timestep\t\t" "%02u\n", m_tileHalfTs);

	Tiling::Plan1D i = Tiling::computeTrapezoidTiles(
		gridSize[0], m_tileSize[0], m_tileHalfTs
	);
	Tiling::Plan1D j = Tiling::computeTrapezoidTiles(
		gridSize[1], m_tileSize[1], m_tileHalfTs
	);

	if (m_tileType[2] == 'p')
	{
		Tiling::Plan1D k = Tiling::computeParallelogramTiles(
			gridSize[2], m_tileSize[2], m_tileHalfTs
		);
		Tiling::Plan3D plan = Tiling::combineTilesTTP(i, j, k);
		return plan;
	}
	else if (m_tileType[2] == 't')
	{
		Tiling::Plan1D k = Tiling::computeTrapezoidTiles(
			gridSize[2], m_tileSize[2], m_tileHalfTs
		);
		Tiling::Plan3D plan = Tiling::combineTilesTTT(i, j, k);
		return plan;
	}
	else
	{
		throw std::invalid_argument(
			std::string("tile suffix must be 't' or 'p', got ") + m_tileType[2]
		);
	}
}

int Operator_Tiling::CalcECOperator(DebugFlags debugFlags)
{
	int err = Operator::CalcECOperator(debugFlags);

	if (g_settings.GetVerboseLevel() > 0)
		cout << "Compressing the FDTD operator... this may take a while..." << endl;

	for (uint32_t oprType = 0; oprType <= OP_TYPE::LAST; oprType++)
		CompressOperators(m_perNodePlan, oprType);

	// uncompressed raw operators are never used once simulation starts
	for (uint32_t oprType = 0; oprType <= OP_TYPE::LAST; oprType++)
	{
		delete m_rawOprArr[oprType];
		m_rawOprArr[oprType] = NULL;
	}

	return err;
}

// Input: A VectorGroup, which is 3 vectors (X, Y, Z polarization) from
// the same region in the simulation domain.
//
// Output: Hash code used by std::unordered_map, for detecting duplicated
// VectorGroups.
template <typename Tvec>
struct VectorGroupHasher
{
	size_t operator() (std::array<Tvec, 3> vg) const
	{
		// Copy the raw bytes within a SIMD vector into a std::u32string
		// buffer, which is hashable by std::hash.
		std::u32string bytestream;
		bytestream.resize(Tvec::len() * 3);

		size_t idx = 0;
		for (size_t n = 0; n < 3; n++) {
			Tvec v = vg[n];

			for (size_t i = 0; i < Tvec::len(); i++) {
				// Assume float can be reinterpreted as uint32_t, which
				// should work on all CPUs with IEEE 754.
				static_assert(sizeof(typename Tvec::Tsca) == sizeof(uint32_t));

				// Can also be (uint32_t*) &v[i], but its behavior
				// is undefined. To do safe type punning, memcpy() is the
				// only way, its result is implementation-defined. In
				// practice, compiler will usually optimize it out, so
				// it's as efficient as a pointer cast.
				memcpy(&bytestream[idx], &v[i], sizeof(typename Tvec::Tsca));
				idx++;
			}
		}

		return std::hash<std::u32string>{}(bytestream);
	}
};

void Operator_Tiling::CompressOperators(
	std::vector<Tiling::Plan3D> perNodePlan,
	uint32_t oprType
)
{
	// Pass 1: For each VectorGroup in the operator array, check if
	// it's a duplicate that we've already seen before. If so, we
	// obtain its index. If not, an index is assigned to this unique
	// vector.
	//
	// Duplication check uses an unordered_map. The key is the hashcode
	// computed from raw FP32 bytes, the value is its unique index.
	std::unordered_map<
		std::array<Simd<float, 4>, 3>, uint32_t,
		VectorGroupHasher<Simd<float, 4>>
	> uniqIdxMap;

	ArrayLib::ArrayNIJK<Simd<float, 4>>& rawOprArr = *m_rawOprArr[oprType];
	ArrayLib::ArrayNIJK<uint32_t, uint32_t, 4>& index = *m_index;

	for (size_t i = 0; i < numLines[0]; i++) {
		for (size_t j = 0; j < numLines[1]; j++) {
			for (size_t vk = 0; vk < numVectors; vk++) {
				std::array<Simd<float, 4>, 3> vecGroup = {
					rawOprArr(0, i, j, vk),
					rawOprArr(1, i, j, vk),
					rawOprArr(2, i, j, vk)
				};

				// not seen before, insert it into the hash table
				if (uniqIdxMap.find(vecGroup) == uniqIdxMap.end()) {
					// The size of the hash table is the number of unique
					// VectorGroups, which increases monotonically so it's
					// also used as the unique index of the current Vector
					// Group.
					uint32_t uniqIdx = uniqIdxMap.size();
					uniqIdxMap[vecGroup] = uniqIdx;
				}

				// already seen before or just inserted, find the unique
				// index of this vector group, and insert it into the
				// coordinate-to-index lookup table.
				index(oprType, i, j, vk) = uniqIdxMap[vecGroup];
			}
		}
	}

	m_compOprArr[oprType] = ArrayLib::AlignedAllocator<Simd<float, 4>>::alloc(
		uniqIdxMap.size() * 3
	);

	fprintf(stderr, "%s: %u -> %lu operators\n",
			m_rawOprArr[oprType]->name().c_str(),
			m_rawOprArr[oprType]->size(),
			uniqIdxMap.size()
	);

	// Pass 2: For each VectorGroup in the operator array, copy its
	// values into the Compressed Unique Operator Array. The index
	// of the VectorGroup determines its memory location in the destination
	// array (obtained from the previously-built unordered_map).
	//
	// On multi-node systems, the first memory write determines the owner
	// of memory, and cross-node memory access is expensive. Thus, the
	// copying process should closely follow the logic of the parallel
	// simulation code to implement first-touch. Implicitly, this maximizes
	// the chance that the same memory is copied and used on the same node.
	for (size_t stage = 0; stage < perNodePlan[0].size(); stage++) {
		for (size_t nodeId = 0; nodeId < ParaExec::numNodes(); nodeId++) {
			const Tiling::TileList3D& tileList = perNodePlan[nodeId][stage];

			ParaExec::paraFor(nodeId, tileList.size(), [&](size_t tileId) {
				const Tiling::Tile3D& tile = tileList[tileId];

				for (const Tiling::Subtile3D& subtile : tile) {
					const std::array<size_t, 3>& first = subtile.first;
					const std::array<size_t, 3>& last = subtile.last;

					for (size_t i = first[0]; i < last[0]; i++) {
						for (size_t j = first[1]; j < last[1]; j++) {
							// The K axis is always compressed entirely
							// regardless of the tiling plan.
							for (size_t vk = 0; vk < rawOprArr.extent(3); vk++) {
								std::array<Simd<float, 4>, 3> vg = {
									rawOprArr(0, i, j, vk),
									rawOprArr(1, i, j, vk),
									rawOprArr(2, i, j, vk)
								};

								uint32_t idx = uniqIdxMap.at(vg);

								m_compOprArr[oprType][idx * 3 + 0] = vg[0];
								m_compOprArr[oprType][idx * 3 + 1] = vg[1];
								m_compOprArr[oprType][idx * 3 + 2] = vg[2];
							}
						}
					}
				}
			});
		}
		ParaExec::barrier();
	}
}

Operator_Tiling::~Operator_Tiling()
{
	Delete();
}

void Operator_Tiling::Delete()
{
	for (uint32_t oprType = 0; oprType <= OP_TYPE::LAST; oprType++)
	{
		delete m_rawOprArr[oprType];
		delete m_compOprArr[oprType];

		m_rawOprArr[oprType] = NULL;
		m_compOprArr[oprType] = NULL;
	}
	delete m_index;
	m_index = NULL;
}

void Operator_Tiling::Reset()
{
	Delete();
	Operator::Reset();
}
