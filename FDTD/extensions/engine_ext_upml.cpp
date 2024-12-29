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

#include "engine_ext_upml.h"
#include "operator_ext_upml.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"

Engine_Ext_UPML::Engine_Ext_UPML(Operator_Ext_UPML* op_ext) : Engine_Extension(op_ext)
{
	m_Op_UPML = op_ext;

	//this ABC extension should be executed first!
	m_Priority = ENG_EXT_PRIO_UPML;
	m_TilingSupported = true;

	volt_flux.Init("volt_flux", m_Op_UPML->m_numLines);
	curr_flux.Init("curr_flux", m_Op_UPML->m_numLines);

	SetNumberOfThreads(1);
}

Engine_Ext_UPML::~Engine_Ext_UPML()
{
}

void Engine_Ext_UPML::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_Op_UPML->m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

template <typename EngineType>
void Engine_Ext_UPML::DoPreVoltageUpdatesImpl(
	EngineType* eng,
	Tiling::Range3D<unsigned int> range
)
{
	if (m_Eng == NULL)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (loc_pos[0] = range.first[0]; loc_pos[0] <= range.last[0]; loc_pos[0]++)
	{
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1] = range.first[1]; loc_pos[1] <= range.last[1]; loc_pos[1]++)
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2] = range.first[2]; loc_pos[2] <= range.last[2]; loc_pos[2]++)
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = m_Op_UPML->vv[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetVolt(0,pos)
						 - m_Op_UPML->vvfo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetVolt(0,pos, volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

				f_help = m_Op_UPML->vv[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetVolt(1,pos)
						 - m_Op_UPML->vvfo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetVolt(1,pos, volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

				f_help = m_Op_UPML->vv[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetVolt(2,pos)
						 - m_Op_UPML->vvfo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetVolt(2,pos, volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
			}
		}
	}
}

void Engine_Ext_UPML::DoPreVoltageUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range3D<unsigned int> range = {
		.first = {m_start.at(threadID), 0, 0},
		.last = {
			m_start.at(threadID) + m_numX.at(threadID) - 1,
			m_Op_UPML->m_numLines[1] - 1,
			m_Op_UPML->m_numLines[2] - 1
		}
	};

	ENG_DISPATCH_ARGS(DoPreVoltageUpdatesImpl, range);
}

void Engine_Ext_UPML::DoPreVoltageUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range3D<unsigned int> range;

	bool overlap = ToLocalCoords(tileRange, range);
	if (overlap)
	{
		DoPreVoltageUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, range
		);
	}
}

template <typename EngineType>
void Engine_Ext_UPML::DoPostVoltageUpdatesImpl(
	EngineType* eng,
	Tiling::Range3D<unsigned int> range
)
{
	if (m_Eng == NULL)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (loc_pos[0] = range.first[0]; loc_pos[0] <= range.last[0]; loc_pos[0]++)
	{
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1] = range.first[1]; loc_pos[1] <= range.last[1]; loc_pos[1]++)
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2] = range.first[2]; loc_pos[2] <= range.last[2]; loc_pos[2]++)
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetVolt(0,pos);
				eng->EngineType::SetVolt(0,pos, f_help + m_Op_UPML->vvfn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

				f_help = volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetVolt(1,pos);
				eng->EngineType::SetVolt(1,pos, f_help + m_Op_UPML->vvfn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

				f_help = volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetVolt(2,pos);
				eng->EngineType::SetVolt(2,pos, f_help + m_Op_UPML->vvfn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
			}
		}
	}
}

void Engine_Ext_UPML::DoPostVoltageUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range3D<unsigned int> range = {
		.first = {m_start.at(threadID), 0, 0},
		.last = {
			m_start.at(threadID) + m_numX.at(threadID) - 1,
			m_Op_UPML->m_numLines[1] - 1,
			m_Op_UPML->m_numLines[2] - 1
		}
	};

	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, range);
}

void Engine_Ext_UPML::DoPostVoltageUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range3D<unsigned int> range;

	bool overlap = ToLocalCoords(tileRange, range);
	if (overlap)
	{
		DoPostVoltageUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, range
		);
	}
}

template <typename EngineType>
void Engine_Ext_UPML::DoPreCurrentUpdatesImpl(
	EngineType* eng,
	Tiling::Range3D<unsigned int> range
)
{
	if (m_Eng == NULL)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (loc_pos[0] = range.first[0]; loc_pos[0] <= range.last[0]; loc_pos[0]++)
	{
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1] = range.first[1]; loc_pos[1] <= range.last[1]; loc_pos[1]++)
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2] = range.first[2]; loc_pos[2] <= range.last[2]; loc_pos[2]++)
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = m_Op_UPML->ii[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetCurr(0,pos)
						 - m_Op_UPML->iifo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetCurr(0,pos, curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

				f_help = m_Op_UPML->ii[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetCurr(1,pos)
						 - m_Op_UPML->iifo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetCurr(1,pos, curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

				f_help = m_Op_UPML->ii[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng->EngineType::GetCurr(2,pos)
						 - m_Op_UPML->iifo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				eng->EngineType::SetCurr(2,pos, curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
				curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

			}
		}
	}
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_UPML::DoPreCurrentUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range3D<unsigned int> range = {
		.first = {m_start.at(threadID), 0, 0},
		.last = {
			m_start.at(threadID) + m_numX.at(threadID) - 1,
			m_Op_UPML->m_numLines[1] - 1,
			m_Op_UPML->m_numLines[2] - 1
		}
	};

	ENG_DISPATCH_ARGS(DoPreCurrentUpdatesImpl, range);
}

// Only used by Tiling engine.
void Engine_Ext_UPML::DoPreCurrentUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range3D<unsigned int> range;

	bool overlap = ToLocalCoords(tileRange, range);
	if (overlap)
	{
		DoPreCurrentUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, range
		);
	}
}

template <typename EngineType>
void Engine_Ext_UPML::DoPostCurrentUpdatesImpl(
	EngineType* eng,
	Tiling::Range3D<unsigned int> range
)
{
	if (m_Eng == NULL)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (loc_pos[0] = range.first[0]; loc_pos[0] <= range.last[0]; loc_pos[0]++)
	{
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1] = range.first[1]; loc_pos[1] <= range.last[1]; loc_pos[1]++)
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2] = range.first[2]; loc_pos[2] <= range.last[2]; loc_pos[2]++)
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetCurr(0,pos);
				eng->EngineType::SetCurr(0,pos, f_help + m_Op_UPML->iifn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

				f_help = curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetCurr(1,pos);
				eng->EngineType::SetCurr(1,pos, f_help + m_Op_UPML->iifn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

				f_help = curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
				curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng->EngineType::GetCurr(2,pos);
				eng->EngineType::SetCurr(2,pos, f_help + m_Op_UPML->iifn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
			}
		}
	}
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_UPML::DoPostCurrentUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range3D<unsigned int> range = {
		.first = {m_start.at(threadID), 0, 0},
		.last = {
			m_start.at(threadID) + m_numX.at(threadID) - 1,
			m_Op_UPML->m_numLines[1] - 1,
			m_Op_UPML->m_numLines[2] - 1
		}
	};

	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, range);
}

// Only used by Tiling engine.
void Engine_Ext_UPML::DoPostCurrentUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range3D<unsigned int> range;

	bool overlap = ToLocalCoords(tileRange, range);

	if (overlap)
	{
		DoPostCurrentUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, range
		);
	}
}

// When the tiling engine is used the global 3D space is divided into small
// tiles for processing. We need to know whether the current tile range
// overlaps with UPML. If there's no overlap, return immediately. If there's
// overlap, return the local UPML coordinates of the overlapped region.
bool
Engine_Ext_UPML::ToLocalCoords(
	const Tiling::Range3D<>& tileRange,
	Tiling::Range3D<unsigned int>& overlapRange
)
{
	Tiling::Range3D<unsigned int> pmlRange = {
		.first = {
			m_Op_UPML->m_StartPos[0],
			m_Op_UPML->m_StartPos[1],
			m_Op_UPML->m_StartPos[2]
		},
		.last = {
			m_Op_UPML->m_StartPos[0] + m_Op_UPML->m_numLines[0] - 1,
			m_Op_UPML->m_StartPos[1] + m_Op_UPML->m_numLines[1] - 1,
			m_Op_UPML->m_StartPos[2] + m_Op_UPML->m_numLines[2] - 1
		}
	};


	for (int n = 0; n < 3; n++)
	{
		if (tileRange.first[n] <= pmlRange.last[n] &&
		    tileRange.last[n] >= pmlRange.first[n])
		{
			// find the overlapped region between the current tile and UPML
			overlapRange.first[n] = std::max(
				(uint32_t) tileRange.first[n], pmlRange.first[n]
			);
			overlapRange.last[n] = std::min(
				(uint32_t) tileRange.last[n], pmlRange.last[n]
			);

			// transform to local coordinate used by PML arrays
			overlapRange.first[n] -= pmlRange.first[n];
			overlapRange.last[n] -= pmlRange.first[n];
		}
		else
		{
			// Tile and UPML do not overlap.
			return false;
		}
	}

	return true;
}
