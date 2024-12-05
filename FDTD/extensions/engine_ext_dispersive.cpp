/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "engine_ext_dispersive.h"
#include "operator_ext_dispersive.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Dispersive::Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp) : Engine_Extension(op_ext_disp)
{
	m_TilingSupported = true;

	m_Op_Ext_Disp = op_ext_disp;
	int order = m_Op_Ext_Disp->m_Order;
	curr_ADE = new FDTD_FLOAT**[order];
	volt_ADE = new FDTD_FLOAT**[order];
	for (int o=0;o<order;++o)
	{
		curr_ADE[o] = new FDTD_FLOAT*[3];
		volt_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			if (m_Op_Ext_Disp->m_curr_ADE_On[o]==true)
			{
				curr_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					curr_ADE[o][n][i]=0.0;
			}
			else
				curr_ADE[o][n] = NULL;
			if (m_Op_Ext_Disp->m_volt_ADE_On[o]==true)
			{
				volt_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					volt_ADE[o][n][i]=0.0;
			}
			else
				volt_ADE[o][n] = NULL;
		}
	}
}

Engine_Ext_Dispersive::~Engine_Ext_Dispersive()
{
	if (curr_ADE==NULL && volt_ADE==NULL)
		return;

	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		for (int n=0; n<3; ++n)
		{
			delete[] curr_ADE[o][n];
			delete[] volt_ADE[o][n];
		}
		delete[] curr_ADE[o];
		delete[] volt_ADE[o];
	}
	delete[] curr_ADE;
	curr_ADE=NULL;

	delete[] volt_ADE;
	volt_ADE=NULL;

	m_cellmap.clear();
}

void Engine_Ext_Dispersive::InitializeTiling(const Tiling::Plan3D& plan)
{
	for (const Tiling::TileList3D& stage : plan)
	{
		for (const Tiling::Tile3D& tile : stage)
		{
			for (const Tiling::Subtile3D& subtile : tile)
			{
				for (const Tiling::Range3D<>& range : subtile)
				{
					InitializeTilingImpl(range);
				}
			}
		}
	}

	if (g_settings.getOption("verbose").as<unsigned int>() > 0)
	{
		for (int o = 0; o < m_Op_Ext_Disp->m_Order; o++)
		{
			std::cout << "Engine_Ext_Dispersive::InitializeTiling: ";
			std::cout << "order " << o << " has "
			          << m_cellmap[o].size() << " cells\n";
		}
	}
}

// Before simulation starts, Engine_Tiling tells us how the space is
// partitioned into ranges of tiles via a data structure called "Plan".
// Here we check all ranges in the plan if this range overlaps with
// the range controlled by us, we append the affected ranges into a
// hash table at m_cellmap[o], otherwise the table is empty. It allows
// us to immediately know where to update when Engine_Tiling passses
// a tile to us.
void Engine_Ext_Dispersive::InitializeTilingImpl(Tiling::Range3D<> range)
{
	// m_cellmap itself is a vector with a length determined by the
	// material physical model's order. The hash tables are located
	// at m_cellmap[o].
	m_cellmap.resize(m_Order);

	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (!m_Op_Ext_Disp->m_volt_ADE_On[o] &&
		    !m_Op_Ext_Disp->m_curr_ADE_On[o])
		{
			continue;
		}

		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
		{
			if (InsideTile(range, pos, i))
				m_cellmap[o][range].push_back(
					{i, pos[0][i], pos[1][i], pos[2][i]}
				);
		}
	}
}

// Whether the target cell is inside the tile that is currently
// being processed.
bool Engine_Ext_Dispersive::InsideTile(
	Tiling::Range3D<> range,
	unsigned int** target,
	unsigned int i
)
{
	for (unsigned int n = 0; n < 3; n++)
	{
		if (target[n][i] < range.first[n] || target[n][i] > range.last[n])
			return false;
	}

	return true;
}

template <typename EngineType>
void Engine_Ext_Dispersive::Apply2VoltagesImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> volt_ADE
)
{
	eng->EngineType::SetVolt(0, pos[0], pos[1], pos[2],
		eng->EngineType::GetVolt(0, pos[0], pos[1], pos[2]) - volt_ADE[0]
	);
	eng->EngineType::SetVolt(1, pos[0], pos[1], pos[2],
		eng->EngineType::GetVolt(1, pos[0], pos[1], pos[2]) - volt_ADE[1]
	);
	eng->EngineType::SetVolt(2, pos[0], pos[1], pos[2],
		eng->EngineType::GetVolt(2, pos[0], pos[1], pos[2]) - volt_ADE[2]
	);
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_Dispersive::Apply2Voltages()
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
		{
			ENG_DISPATCH_ARGS(
				Apply2VoltagesImpl,
				{pos[0][i], pos[1][i], pos[2][i]},
				{volt_ADE[o][0][i], volt_ADE[o][1][i], volt_ADE[o][2][i]}
			);
		}
	}
}

// Only used by Tiling engine.
void Engine_Ext_Dispersive::Apply2Voltages(int timestep, Tiling::Range3D<> range)
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false)
			continue;
		if (m_cellmap[o].count(range) == 0)
			continue;

		for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
		{
			unsigned int idx = val[0];
			std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

			ENG_DISPATCH_ARGS(
				Apply2VoltagesImpl,
				pos,
				{volt_ADE[o][0][idx], volt_ADE[o][1][idx], volt_ADE[o][2][idx]}
			);
		}
	}
}

template <typename EngineType>
void Engine_Ext_Dispersive::Apply2CurrentImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> curr_ADE
)
{
	eng->EngineType::SetCurr(0, pos[0], pos[1], pos[2],
		eng->EngineType::GetCurr(0, pos[0], pos[1], pos[2]) - curr_ADE[0]
	);
	eng->EngineType::SetCurr(1, pos[0], pos[1], pos[2],
		eng->EngineType::GetCurr(1, pos[0], pos[1], pos[2]) - curr_ADE[1]
	);
	eng->EngineType::SetCurr(2, pos[0], pos[1], pos[2],
		eng->EngineType::GetCurr(2, pos[0], pos[1], pos[2]) - curr_ADE[2]
	);
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_Dispersive::Apply2Current()
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
		{
			ENG_DISPATCH_ARGS(
				Apply2CurrentImpl,
				{pos[0][i], pos[1][i], pos[2][i]},
				{curr_ADE[o][0][i], curr_ADE[o][1][i], curr_ADE[o][2][i]}
			);
		}
	}
}

// Only used by Tiling engine.
void Engine_Ext_Dispersive::Apply2Current(int timestep, Tiling::Range3D<> range)
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false)
			continue;
		if (m_cellmap[o].count(range) == 0)
			continue;

		for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
		{
			unsigned int idx = val[0];
			std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

			ENG_DISPATCH_ARGS(
				Apply2CurrentImpl,
				pos,
				{curr_ADE[o][0][idx], curr_ADE[o][1][idx], curr_ADE[o][2][idx]}
			);
		}
	}
}
