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
#ifndef TILING_MACHTOPO_H
#define TILING_MACHTOPO_H

#include <vector>
#include <hwloc.h>

// Collect and hold information about machine's CPU topology, including the
// number of nodes and the corresponding CPU core IDs under that node.
// Cross-platform via the hwloc library.
//
// The definition of a node is "a group of CPU cores with a dedicated L3
// cache." This includes traditional NUMA nodes on different CPU sockets,
// and also single-socket CPUs with multiple dies.
//
// Currently only one instance is created and held internally by parexec.hpp.
class MachTopo
{
public:
	inline MachTopo(size_t numCpus=0, bool useMultipleNodes=true, bool nodeAffinity=true);

	inline size_t numNodes();
	inline hwloc_topology_t topology();
	inline std::vector<size_t> cpuList(size_t nodeIdx);
	inline std::vector<hwloc_cpuset_t> cpuSet(size_t nodeIdx);

private:
	size_t m_numNodes;
	hwloc_topology_t m_topology;
	std::vector<std::vector<size_t>> m_cpuList;
	std::vector<std::vector<hwloc_cpuset_t>> m_cpuSet;
};

MachTopo::MachTopo(size_t numCpus, bool useMultipleNodes, bool nodeAffinity)
{
	hwloc_topology_init(&m_topology);
	hwloc_topology_load(m_topology);

	size_t numNodes;

	if (useMultipleNodes && nodeAffinity)
	{
		numNodes = hwloc_get_nbobjs_by_type(m_topology, HWLOC_OBJ_L3CACHE);
		m_numNodes = numNodes;
	}
	else if (useMultipleNodes && !nodeAffinity)
	{
		numNodes = hwloc_get_nbobjs_by_type(m_topology, HWLOC_OBJ_L3CACHE);
		m_numNodes = 1;
	}
	else
	{
		numNodes = 1;
		m_numNodes = 1;
	}

	m_cpuList.resize(numNodes);
	m_cpuSet.resize(numNodes);

	for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++)
	{
		hwloc_obj_t node = hwloc_get_obj_by_type(
			m_topology, HWLOC_OBJ_L3CACHE, nodeIdx
		);

		size_t systemNumCpus = hwloc_get_nbobjs_inside_cpuset_by_type(
			m_topology, node->cpuset, HWLOC_OBJ_PU
		);

		if (numCpus != 0 && systemNumCpus > numCpus)
			systemNumCpus = numCpus;

		for (size_t cpuIdx = 0; cpuIdx < systemNumCpus; cpuIdx++)
		{
			hwloc_obj_t cpu = hwloc_get_obj_inside_cpuset_by_type(
				m_topology, node->cpuset, HWLOC_OBJ_PU, cpuIdx
			);

			if (nodeAffinity)
			{
				m_cpuList[nodeIdx].push_back(cpu->os_index);
				m_cpuSet[nodeIdx].push_back(cpu->cpuset);
			}
			else
			{
				m_cpuList[0].push_back(cpu->os_index);
				m_cpuSet[0].push_back(cpu->cpuset);
			}
		}
	}
}

size_t MachTopo::numNodes()
{
	return m_numNodes;
}

hwloc_topology_t MachTopo::topology()
{
	return m_topology;
}

std::vector<size_t>
MachTopo::cpuList(size_t nodeIdx)
{
	return m_cpuList[nodeIdx];
}

std::vector<hwloc_cpuset_t>
MachTopo::cpuSet(size_t nodeIdx)
{
	return m_cpuSet[nodeIdx];
}

#endif // TILING_MACHTOPO_H
