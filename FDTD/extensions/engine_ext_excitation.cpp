/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_excitation.h"
#include "operator_ext_excitation.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Excitation::Engine_Ext_Excitation(Operator_Ext_Excitation* op_ext) : Engine_Extension(op_ext)
{
	m_Op_Exc = op_ext;
	m_Priority = ENG_EXT_PRIO_EXCITATION;
	m_TilingSupported = true;
}

Engine_Ext_Excitation::~Engine_Ext_Excitation()
{

}

// Whether the excited cell is inside the tile that is currently
// being processed.
//
// TODO: It checks every excited cell and is inefficient. It's
// better to pre-determine all cells involved before simulation
// in InitializeTiling(). But it works for now since most simulations
// only have tiny ports.
bool Engine_Ext_Excitation::InsideTile(
	Tiling::Range3D<> range,
	unsigned int exc[3]
)
{
	for (unsigned int n = 0; n < 3; n++)
	{
		if (exc[n] < range.first[n] || exc[n] > range.last[n])
			return false;
	}

	return true;
}

template <typename EngineType, bool tiling>
void Engine_Ext_Excitation::Apply2VoltagesImpl(
	EngineType *eng,
	int numTS,
	Tiling::Range3D<> range
)
{
	//soft voltage excitation here (E-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_volt =  m_Op_Exc->m_Exc->GetVoltageSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
	{
		pos[0] = m_Op_Exc->Volt_index[0][n];
		pos[1] = m_Op_Exc->Volt_index[1][n];
		pos[2] = m_Op_Exc->Volt_index[2][n];

		if constexpr (tiling)
		{
			if (!InsideTile(range, pos))
				continue;
		}

		exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Volt_dir[n];
		eng->EngineType::SetVolt(ny,pos, eng->EngineType::GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
	}
}

void Engine_Ext_Excitation::Apply2Voltages()
{
	ENG_DISPATCH_ARGS(
		Apply2VoltagesImpl,
		m_Eng->GetNumberOfTimesteps(),
		Tiling::Range3D<>()  // dummy value
	);
}

void Engine_Ext_Excitation::Apply2Voltages(int timestep, Tiling::Range3D<> range)
{
	Apply2VoltagesImpl<Engine_Tiling, true>(
		(Engine_Tiling*) m_Eng,
		timestep,
		range
	);
}

template <typename EngineType, bool tiling>
void Engine_Ext_Excitation::Apply2CurrentImpl(
	EngineType *eng,
	int numTS,
	Tiling::Range3D<> range
)
{
	//soft current excitation here (H-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_curr =  m_Op_Exc->m_Exc->GetCurrentSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
	{
		pos[0] = m_Op_Exc->Curr_index[0][n];
		pos[1] = m_Op_Exc->Curr_index[1][n];
		pos[2] = m_Op_Exc->Curr_index[2][n];

		if constexpr (tiling)
		{
			if (!InsideTile(range, pos))
				continue;
		}

		exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Curr_dir[n];
		eng->EngineType::SetCurr(ny,pos, eng->EngineType::GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
	}
}

void Engine_Ext_Excitation::Apply2Current()
{
	ENG_DISPATCH_ARGS(
		Apply2CurrentImpl,
		m_Eng->GetNumberOfTimesteps(),
		Tiling::Range3D<>()  // dummy value
	);
}

void Engine_Ext_Excitation::Apply2Current(int timestep, Tiling::Range3D<> range)
{
	Apply2CurrentImpl<Engine_Tiling, true>(
		(Engine_Tiling*) m_Eng,
		timestep,
		range
	);
}
