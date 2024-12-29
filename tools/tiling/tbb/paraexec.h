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
#ifndef TILING_PARAEXEC_H
#define TILING_PARAEXEC_H

#include <functional>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <tbb/task_scheduler_observer.h>

#include "machtopo.h"

// This header-only library uses Intel Threading Bloking Block (Intel TBB)
// to implement a parallel_for function (similar to OpenMP's #pragma omp
// parallel for or SYCL's parallel_for). However, instead of using a global
// thread pool, it creates one thread pool per CPU node using TBB's task_arena
// and task_group feature for NUMA awareness.
namespace ParaExec
{
	inline MachTopo* machine;
	inline std::vector<tbb::task_arena> arenaList;
	inline std::vector<tbb::task_group*> taskGroupList;
	inline std::vector<double> runtimeList;

	class nodePinner : public tbb::task_scheduler_observer
	{
	public:
		inline nodePinner(tbb::task_arena &a, size_t nodeId);
		inline void on_scheduler_entry(bool worker) override;
		inline void on_scheduler_exit(bool worker) override;

	private:
		size_t m_nodeId;
	};
	inline std::vector<nodePinner*> pinningObserversList;

	inline void init(size_t numCpus=0, bool useMultipleNodes=true, bool nodeAwareness=true);
	inline void printSystemInfo();

	inline size_t numNodes();

	inline void runOnNode(
		size_t nodeId,
		const std::function<void (void)>& func
	);

	inline void paraFor(
		size_t nodeId, size_t size, size_t granularity,
		const std::function<void (size_t)>& func
	);

	inline void paraFor(
		size_t nodeId, size_t size,
		const std::function<void (size_t)>& func
	);

	inline void barrier();
}

void ParaExec::init(size_t numCpus, bool useMultipleNodes, bool nodeAwareness)
{
	machine = new MachTopo(numCpus, useMultipleNodes, nodeAwareness);

	arenaList.resize(machine->numNodes());
	runtimeList.resize(machine->numNodes());

	for (size_t nodeId = 0; nodeId < ParaExec::numNodes(); nodeId++)
	{
		arenaList[nodeId].initialize(
			machine->cpuList(nodeId).size()
		);

		// cat't use non-pointer because it's not copy-constructable.
		taskGroupList.push_back(new tbb::task_group());

		if (nodeAwareness)
		{
			pinningObserversList.push_back(
				new nodePinner(arenaList[nodeId], nodeId)
			);
		}
	}
}

void ParaExec::printSystemInfo()
{
	fprintf(stderr, "total L3 cache node: %zu\n", numNodes());
	for (size_t i = 0; i < numNodes(); i++)
	{
		const std::vector<size_t>& cpuList = ParaExec::machine->cpuList(i);

		fprintf(stderr, "node %zu (%zu CPUs): ", i, cpuList.size());
		for (size_t cpuId : cpuList)
			fprintf(stderr, "%zu ", cpuId);
		fprintf(stderr, "\n");
	}
}

size_t ParaExec::numNodes()
{
	return ParaExec::machine->numNodes();
}

void ParaExec::runOnNode(
	size_t nodeId,
	const std::function<void (void)>& func
)
{
	arenaList[nodeId].execute([=]() {
		taskGroupList[nodeId]->run([=]()
		{
			auto t1 = std::chrono::high_resolution_clock().now();
			func();
			auto t2 = std::chrono::high_resolution_clock().now();

			std::chrono::duration<double, std::micro> elapsed = t2 - t1;
			runtimeList[nodeId] = elapsed.count();
		});
	});
}

void ParaExec::paraFor(
	size_t nodeId, size_t size, size_t granularity,
	const std::function<void (size_t)>& func
)
{
	assert(arenaList.size() == taskGroupList.size());

	runOnNode(nodeId, [=]()
	{
		tbb::parallel_for(
			tbb::blocked_range<size_t>{0, size, granularity},
			[=](const tbb::blocked_range<size_t> &r) {
				for (size_t i = r.begin(); i < r.end(); i++) {
					func(i);
				}
			}
		);
	});
}

void ParaExec::paraFor(
	size_t nodeId, size_t size,
	const std::function<void (size_t)>& func
)
{
	paraFor(nodeId, size, 1, func);
}

void ParaExec::barrier()
{
	for (size_t nodeId = 0; nodeId < numNodes(); nodeId++)
	{
		arenaList[nodeId].execute([=]() {
			taskGroupList[nodeId]->wait();
		});
	}
}

// nodePinner: pin all TBB tasks within a task arena onto a specific CPU node.
//
// Unlike the traditional definition of NUMA, our own definition of a CPU node
// is "groups of CPU sharing the same L3 cache," so this includes traditional
// NUMA nodes on different CPU sockets, and also single-socket CPUs with
// multiple dies.
ParaExec::nodePinner::nodePinner(tbb::task_arena &a, size_t nodeId)
:
	tbb::task_scheduler_observer(a),
	m_nodeId(nodeId)
{
	observe(true);
}

void ParaExec::nodePinner::on_scheduler_entry(bool worker)
{
	size_t tid = tbb::this_task_arena::current_thread_index();

	hwloc_set_cpubind(
		ParaExec::machine->topology(),
		ParaExec::machine->cpuSet(m_nodeId)[tid],
		HWLOC_CPUBIND_THREAD
	);
}

void ParaExec::nodePinner::on_scheduler_exit(bool worker)
{
	hwloc_cpuset_t fullSet = hwloc_bitmap_alloc();
	if (fullSet == NULL)
	{
		throw std::runtime_error(
			"ParaExec::nodePinner: hwloc_bitmap_alloc() failed!"
		);
	}

	hwloc_bitmap_fill(fullSet);
	hwloc_set_cpubind(
		ParaExec::machine->topology(),
		fullSet,
		HWLOC_CPUBIND_THREAD
	);
	hwloc_bitmap_free(fullSet);
}

#endif // TILING_PARAEXEC_H
