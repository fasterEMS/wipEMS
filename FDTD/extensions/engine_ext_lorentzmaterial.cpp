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

#include "engine_ext_lorentzmaterial.h"
#include "operator_ext_lorentzmaterial.h"
#include "FDTD/engine_sse.h"

Engine_Ext_LorentzMaterial::Engine_Ext_LorentzMaterial(Operator_Ext_LorentzMaterial* op_ext_lorentz) : Engine_Ext_Dispersive(op_ext_lorentz)
{
	m_TilingSupported = true;

	m_Op_Ext_Lor = op_ext_lorentz;
	m_Order = m_Op_Ext_Lor->GetDispersionOrder();
	int order = m_Op_Ext_Lor->m_Order;

	curr_Lor_ADE = new FDTD_FLOAT**[order];
	volt_Lor_ADE = new FDTD_FLOAT**[order];
	for (int o=0;o<order;++o)
	{
		curr_Lor_ADE[o] = new FDTD_FLOAT*[3];
		volt_Lor_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			if (m_Op_Ext_Lor->m_curr_Lor_ADE_On[o]==true)
			{
				curr_Lor_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Lor->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count[o]; ++i)
					curr_Lor_ADE[o][n][i]=0.0;
			}
			else
				curr_Lor_ADE[o][n] = NULL;

			if (m_Op_Ext_Lor->m_volt_Lor_ADE_On[o]==true)
			{
				volt_Lor_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Lor->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count[o]; ++i)
					volt_Lor_ADE[o][n][i]=0.0;
			}
			else
				volt_Lor_ADE[o][n] = NULL;
		}
	}
}

Engine_Ext_LorentzMaterial::~Engine_Ext_LorentzMaterial()
{
	if (curr_Lor_ADE==NULL && volt_Lor_ADE==NULL)
		return;

	for (int o=0;o<m_Op_Ext_Lor->m_Order;++o)
	{
		for (int n=0; n<3; ++n)
		{
			delete[] curr_Lor_ADE[o][n];
			delete[] volt_Lor_ADE[o][n];
		}
		delete[] curr_Lor_ADE[o];
		delete[] volt_Lor_ADE[o];
	}
	delete[] curr_Lor_ADE;
	curr_Lor_ADE=NULL;

	delete[] volt_Lor_ADE;
	volt_Lor_ADE=NULL;

	m_cellmap.clear();
}

void Engine_Ext_LorentzMaterial::InitializeTilingImpl(Tiling::Range3D<> range)
{
	for (int o=0;o<m_Order;++o)
	{
		if (!m_Op_Ext_Lor->m_volt_ADE_On[o] &&
		    !m_Op_Ext_Lor->m_curr_ADE_On[o])
		{
			continue;
		}

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
		{
			if (InsideTile(range, pos, i))
				m_cellmap[o][range].push_back(
					{i, pos[0][i], pos[1][i], pos[2][i]}
				);
		}
	}
}

void Engine_Ext_LorentzMaterial::InitializeTiling(const Tiling::Plan3D& plan)
{
	Engine_Ext_Dispersive::InitializeTiling(plan);

	m_cellmap.resize(m_Order);

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
		for (int o=0;o<m_Order;++o)
		{
			std::cout << "Engine_Ext_LorentzMaterial::InitializeTiling: ";
			std::cout << "order " << o << " has "
			          << m_cellmap[o].size() << " cells\n";
		}
	}
}

template <typename EngineType>
void Engine_Ext_LorentzMaterial::DoPreVoltageUpdatesImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> v_int_ADE,
	std::array<float, 3> v_ext_ADE,
	std::array<float, 3> v_Lor_ADE,
	std::array<float*, 3> volt_Lor_ADE,
	std::array<float*, 3> volt_ADE
)
{
	*volt_Lor_ADE[0] += v_Lor_ADE[0] * (*volt_ADE[0]);
	*volt_ADE[0] *= v_int_ADE[0];
	*volt_ADE[0] += v_ext_ADE[0] * (eng->EngineType::GetVolt(0,pos[0],pos[1],pos[2]) - (*volt_Lor_ADE[0]));

	*volt_Lor_ADE[1] += v_Lor_ADE[1] * (*volt_ADE[1]);
	*volt_ADE[1] *= v_int_ADE[1];
	*volt_ADE[1] += v_ext_ADE[1] * (eng->EngineType::GetVolt(1,pos[0],pos[1],pos[2]) - (*volt_Lor_ADE[1]));

	*volt_Lor_ADE[2] += v_Lor_ADE[2] * (*volt_ADE[2]);
	*volt_ADE[2] *= v_int_ADE[2];
	*volt_ADE[2] += v_ext_ADE[2] * (eng->EngineType::GetVolt(2,pos[0],pos[1],pos[2]) - (*volt_Lor_ADE[2]));
}

template <typename EngineType>
void Engine_Ext_LorentzMaterial::DoPreVoltageUpdatesImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> v_int_ADE,
	std::array<float, 3> v_ext_ADE,
	std::array<float*, 3> volt_ADE
)
{
	*volt_ADE[0] *= v_int_ADE[0];
	*volt_ADE[0] += v_ext_ADE[0] * eng->EngineType::GetVolt(0,pos[0],pos[1],pos[2]);

	*volt_ADE[1] *= v_int_ADE[1];
	*volt_ADE[1] += v_ext_ADE[1] * eng->EngineType::GetVolt(1,pos[0],pos[1],pos[2]);

	*volt_ADE[2] *= v_int_ADE[2];
	*volt_ADE[2] += v_ext_ADE[2] * eng->EngineType::GetVolt(2,pos[0],pos[1],pos[2]);
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_LorentzMaterial::DoPreVoltageUpdates()
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_volt_ADE_On[o]==false)
			continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_volt_Lor_ADE_On[o])
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
			{
				std::array<float, 3> v_int_ADE, v_ext_ADE, v_Lor_ADE;
				std::array<float*, 3> Lor_ADE, ADE;

				for (int n = 0; n < 3; n++)
				{
					v_int_ADE[n] = m_Op_Ext_Lor->v_int_ADE[o][n][i];
					v_ext_ADE[n] = m_Op_Ext_Lor->v_ext_ADE[o][n][i];
					v_Lor_ADE[n] = m_Op_Ext_Lor->v_Lor_ADE[o][n][i];
					Lor_ADE[n] = &volt_Lor_ADE[o][n][i];
					ADE[n] = &volt_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					{pos[0][i], pos[1][i], pos[2][i]},
					v_int_ADE, v_ext_ADE, v_Lor_ADE, Lor_ADE, ADE
				);
			}
		}
		else
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
			{
				std::array<float, 3> v_int_ADE, v_ext_ADE;
				std::array<float*, 3> ADE;

				for (int n = 0; n < 3; n++)
				{
					v_int_ADE[n] = m_Op_Ext_Lor->v_int_ADE[o][n][i];
					v_ext_ADE[n] = m_Op_Ext_Lor->v_ext_ADE[o][n][i];
					ADE[n] = &volt_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					{pos[0][i], pos[1][i], pos[2][i]},
					v_int_ADE, v_ext_ADE, ADE
				);
			}
		}
	}
}

// Only used by Tiling engine.
void Engine_Ext_LorentzMaterial::DoPreVoltageUpdates(int timestep, Tiling::Range3D<> range)
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_volt_ADE_On[o]==false)
			continue;
		if (m_cellmap[o].count(range) == 0)
			continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_volt_Lor_ADE_On[o])
		{
			for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
			{
				unsigned int i = val[0];
				std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

				std::array<float, 3> v_int_ADE, v_ext_ADE, v_Lor_ADE;
				std::array<float*, 3> Lor_ADE, ADE;

				for (int n = 0; n < 3; n++)
				{
					v_int_ADE[n] = m_Op_Ext_Lor->v_int_ADE[o][n][i];
					v_ext_ADE[n] = m_Op_Ext_Lor->v_ext_ADE[o][n][i];
					v_Lor_ADE[n] = m_Op_Ext_Lor->v_Lor_ADE[o][n][i];
					Lor_ADE[n] = &volt_Lor_ADE[o][n][i];
					ADE[n] = &volt_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					pos, v_int_ADE, v_ext_ADE, v_Lor_ADE, Lor_ADE, ADE
				);
			}
		}
		else
		{
			for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
			{
				unsigned int i = val[0];
				std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

				std::array<float, 3> v_int_ADE, v_ext_ADE;
				std::array<float*, 3> ADE;

				for (int n = 0; n < 3; n++)
				{
					v_int_ADE[n] = m_Op_Ext_Lor->v_int_ADE[o][n][i];
					v_ext_ADE[n] = m_Op_Ext_Lor->v_ext_ADE[o][n][i];
					ADE[n] = &volt_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					pos, v_int_ADE, v_ext_ADE, ADE
				);
			}
		}
	}
}

template <typename EngineType>
void Engine_Ext_LorentzMaterial::DoPreCurrentUpdatesImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> i_Lor_ADE,
	std::array<float, 3> i_int_ADE,
	std::array<float, 3> i_ext_ADE,
	std::array<float*, 3> curr_Lor_ADE,
	std::array<float*, 3> curr_ADE
)
{
	*curr_Lor_ADE[0] += i_Lor_ADE[0] * (*curr_ADE[0]);
	*curr_ADE[0] *= i_int_ADE[0];
	*curr_ADE[0] += i_ext_ADE[0] * (eng->EngineType::GetCurr(0,pos[0],pos[1],pos[2]) - (*curr_Lor_ADE[0]));

	*curr_Lor_ADE[1] += i_Lor_ADE[1] * (*curr_ADE[1]);
	*curr_ADE[1] *= i_int_ADE[1];
	*curr_ADE[1] += i_ext_ADE[1] * (eng->EngineType::GetCurr(1,pos[0],pos[1],pos[2]) - (*curr_Lor_ADE[1]));

	*curr_Lor_ADE[2] += i_Lor_ADE[2] * (*curr_ADE[2]);
	*curr_ADE[2] *= i_int_ADE[2];
	*curr_ADE[2] += i_ext_ADE[2] * (eng->EngineType::GetCurr(2,pos[0],pos[1],pos[2]) - (*curr_Lor_ADE[2]));
}

template <typename EngineType>
void Engine_Ext_LorentzMaterial::DoPreCurrentUpdatesImpl(
	EngineType* eng,
	std::array<unsigned int, 3> pos,
	std::array<float, 3> i_int_ADE,
	std::array<float, 3> i_ext_ADE,
	std::array<float*, 3> curr_ADE
)
{
	*curr_ADE[0] *= i_int_ADE[0];
	*curr_ADE[0] += i_ext_ADE[0] * eng->EngineType::GetCurr(0,pos[0],pos[1],pos[2]);

	*curr_ADE[1] *= i_int_ADE[1];
	*curr_ADE[1] += i_ext_ADE[1] * eng->EngineType::GetCurr(1,pos[0],pos[1],pos[2]);

	*curr_ADE[2] *= i_int_ADE[2];
	*curr_ADE[2] += i_ext_ADE[2] * eng->EngineType::GetCurr(2,pos[0],pos[1],pos[2]);
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_LorentzMaterial::DoPreCurrentUpdates()
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_curr_ADE_On[o]==false)
			continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_curr_Lor_ADE_On[o])
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
			{
				std::array<float, 3> i_int_ADE, i_ext_ADE, i_Lor_ADE;
				std::array<float*, 3> Lor_ADE, ADE;

				for (int n = 0; n < 3; n++)
				{
					i_int_ADE[n] = m_Op_Ext_Lor->i_int_ADE[o][n][i];
					i_ext_ADE[n] = m_Op_Ext_Lor->i_ext_ADE[o][n][i];
					i_Lor_ADE[n] = m_Op_Ext_Lor->i_Lor_ADE[o][n][i];
					Lor_ADE[n] = &curr_Lor_ADE[o][n][i];
					ADE[n] = &curr_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					{pos[0][i], pos[1][i], pos[2][i]},
					i_int_ADE, i_ext_ADE, i_Lor_ADE, Lor_ADE, ADE
				);
			}
		}
		else
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
			{
				std::array<float, 3> i_int_ADE, i_ext_ADE;
				std::array<float*, 3> ADE;

				for (int n = 0; n < 3; n++)
				{
					i_int_ADE[n] = m_Op_Ext_Lor->i_int_ADE[o][n][i];
					i_ext_ADE[n] = m_Op_Ext_Lor->i_ext_ADE[o][n][i];
					ADE[n] = &curr_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreVoltageUpdatesImpl,
					{pos[0][i], pos[1][i], pos[2][i]},
					i_int_ADE, i_ext_ADE, ADE
				);
			}
		}
	}
}

// Only used by Tiling engine.
void Engine_Ext_LorentzMaterial::DoPreCurrentUpdates(int timestep, Tiling::Range3D<> range)
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_curr_ADE_On[o]==false)
			continue;
		if (m_cellmap[o].count(range) == 0)
			continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_curr_Lor_ADE_On[o])
		{
			for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
			{
				unsigned int i = val[0];
				std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

				std::array<float, 3> i_int_ADE, i_ext_ADE, i_Lor_ADE;
				std::array<float*, 3> Lor_ADE, ADE;

				for (int n = 0; n < 3; n++)
				{
					i_int_ADE[n] = m_Op_Ext_Lor->i_int_ADE[o][n][i];
					i_ext_ADE[n] = m_Op_Ext_Lor->i_ext_ADE[o][n][i];
					i_Lor_ADE[n] = m_Op_Ext_Lor->i_Lor_ADE[o][n][i];
					Lor_ADE[n] = &curr_Lor_ADE[o][n][i];
					ADE[n] = &curr_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreCurrentUpdatesImpl,
					pos, i_int_ADE, i_ext_ADE, i_Lor_ADE, Lor_ADE, ADE
				);
			}
		}
		else
		{
			for (const std::array<unsigned int, 4>& val : m_cellmap[o].at(range))
			{
				unsigned int i = val[0];
				std::array<unsigned int, 3> pos = {val[1], val[2], val[3]};

				std::array<float, 3> i_int_ADE, i_ext_ADE;
				std::array<float*, 3> ADE;

				for (int n = 0; n < 3; n++)
				{
					i_int_ADE[n] = m_Op_Ext_Lor->i_int_ADE[o][n][i];
					i_ext_ADE[n] = m_Op_Ext_Lor->i_ext_ADE[o][n][i];
					ADE[n] = &curr_ADE[o][n][i];
				}

				ENG_DISPATCH_ARGS(
					DoPreCurrentUpdatesImpl,
					pos, i_int_ADE, i_ext_ADE, ADE
				);
			}
		}
	}
}
