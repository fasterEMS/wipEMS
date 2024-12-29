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

#include <boost/format.hpp>
#include "engine_ext_mur_abc.h"
#include "operator_ext_mur_abc.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "operator_ext_excitation.h"

Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC(Operator_Ext_Mur_ABC* op_ext) :
	Engine_Extension(op_ext),
	m_Mur_Coeff_dir2 (op_ext->m_Mur_Coeff_dir2),
	m_Mur_Coeff_dir3(op_ext->m_Mur_Coeff_dir3),
	m_volt_dir2 ("volt_dir2",  op_ext->m_numLines),
	m_volt_dir3("volt_dir3", op_ext->m_numLines)
{
	m_Op_mur = op_ext;
	m_numLines[0] = m_Op_mur->m_numLines[0];
	m_numLines[1] = m_Op_mur->m_numLines[1];
	m_dir1 = m_Op_mur->m_dir1;
	m_dir2 = m_Op_mur->m_dir2;
	m_dir3 = m_Op_mur->m_dir3;
	m_LineNr = m_Op_mur->m_LineNr;
	m_LineNr_Shift = m_Op_mur->m_LineNr_Shift;

	//find if some excitation is on this mur-abc and find the max length of this excite, so that the abc can start after the excitation is done...
	int maxDelay=-1;
	Operator_Ext_Excitation* Exc_ext = m_Op_mur->m_Op->GetExcitationExtension();
	for (unsigned int n=0; n<Exc_ext->GetVoltCount(); ++n)
	{
		if ( ((Exc_ext->Volt_dir[n]==m_dir2) || (Exc_ext->Volt_dir[n]==m_dir3)) && (Exc_ext->Volt_index[m_dir1][n]==m_LineNr) )
		{
			if ((int)Exc_ext->Volt_delay[n]>maxDelay)
				maxDelay = (int)Exc_ext->Volt_delay[n];
		}
	}
	m_start_TS = 0;
	if (maxDelay>=0)
	{
		m_start_TS = maxDelay + m_Op_mur->m_Op->GetExcitationSignal()->GetLength() + 10; //give it some extra timesteps, for the excitation to travel at least one cell away
		cerr << "Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC: Warning: Excitation inside the Mur-ABC #" <<  m_dir1 << "-" << (int)(m_LineNr>0) << " found!!!!  Mur-ABC will be switched on after excitation is done at " << m_start_TS << " timesteps!!! " << endl;
	}

	SetNumberOfThreads(1);

	m_TilingSupported = true;
}

Engine_Ext_Mur_ABC::~Engine_Ext_Mur_ABC()
{
}

void Engine_Ext_Mur_ABC::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

void Engine_Ext_Mur_ABC::InitializeTiling(const Tiling::Plan3D& plan)
{
	for (const Tiling::TileList3D& stage : plan)
	{
		for (const Tiling::Tile3D& tile : stage)
		{
			for (const Tiling::Subtile3D& subtile : tile)
			{
				for (size_t i = 0; i < subtile.size(); i++)
				{
					// Mur ABC only uses electric field, so the magnetic
					// field's range (at odd half timesteps) must NOT be
					// checked. The magnetic mesh is one cell smaller than
					// the electric mesh, so it appears to be out of range
					// to us, causing a spurious range check failure in
					// InitializeTilingImpl().
					if (i % 2 == 1)
						continue;

					const Tiling::Range3D<>& range = subtile[i];
					InitializeTilingImpl(range);
				}
			}
		}
	}
}

// Before simulation starts, Engine_Tiling tells us how the space is
// partitioned into ranges of tiles via a data structure called "Plan".
// Here we check all ranges in the plan if this range overlaps with
// the range controlled by us, we mark the range as "true" in the hash
// table m_needRun, otherwise it's marked False. It allows us to quickly
// skip irrelevant tiles without range comparisons.
void Engine_Ext_Mur_ABC::InitializeTilingImpl(Tiling::Range3D<> range)
{
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_dir1] = m_LineNr;
	pos_shift[m_dir1] = m_LineNr_Shift;

	for (unsigned int i = 0; i < m_numLines[0]; i++)
	{
		pos[m_dir2] = i;
		pos_shift[m_dir3] = i;

		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos[m_dir3] = j;
			pos_shift[m_dir3] = j;

			if (InsideTile(range, pos) && InsideTile(range, pos_shift))
				// Mur ABC cells (controlled by us) is in this tile
				m_needRun[range] = true;
			else if (!InsideTile(range, pos) && !InsideTile(range, pos_shift))
				// Mur ABC cells is NOT in this tile
				m_needRun[range] = false;
			else
			{
				// Mur ABC cells is only partially in this tile
				boost::format fmt = boost::format(
					"Unsupported tiling partitioning for Mur ABC detected! "
					"Cells (%u, %u, %u) and (%u, %u, %u) must be both inside "
					"tile (%lu-%lu, %lu-%lu, %lu-%lu), not one inside, one "
					"outside."
				) % pos[0] % pos[1] % pos[2]
				  % pos_shift[0] % pos_shift[1] % pos_shift[2]
				  % range.first[0] % range.last[0]
				  % range.first[1] % range.last[1]
				  % range.first[2] % range.last[2];

				throw std::runtime_error(boost::str(fmt));
			}
		}
	}
}

bool Engine_Ext_Mur_ABC::InsideTile(
	Tiling::Range3D<> range,
	unsigned int target[3]
)
{
	for (unsigned int n = 0; n < 3; n++)
	{
		if (target[n] < range.first[n] || target[n] > range.last[n])
			return false;
	}

	return true;
}

template <typename EngineType, bool tiling>
void Engine_Ext_Mur_ABC::DoPreVoltageUpdatesImpl(
	EngineType* eng,
	Tiling::Range2D<> abcRange,
	Tiling::Range3D<> tileRange
)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_dir1] = m_LineNr;
	pos_shift[m_dir1] = m_LineNr_Shift;

	for (unsigned int i = abcRange.first[0]; i <= abcRange.last[0]; i++)
	{
		pos[m_dir2] = i;
		pos_shift[m_dir3] = i;

		for (unsigned int j = abcRange.first[1]; j <= abcRange.last[1]; j++)
		{
			pos[m_dir3] = j;
			pos_shift[m_dir3] = j;

			if (tiling && !InsideTile(tileRange, pos) && !InsideTile(tileRange, pos_shift))
				continue;

			m_volt_dir2[i][j] = eng->EngineType::GetVolt(m_dir2, pos_shift) -
						m_Op_mur->m_Mur_Coeff_dir2[i][j] *
						eng->EngineType::GetVolt(m_dir2, pos);

			m_volt_dir3[i][j] = eng->EngineType::GetVolt(m_dir3, pos_shift) -
						m_Op_mur->m_Mur_Coeff_dir3[i][j] *
						eng->EngineType::GetVolt(m_dir3, pos);
		}
	}
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_Mur_ABC::DoPreVoltageUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range2D<> abcRange;
	abcRange.first = {m_start.at(threadID), 0};
	abcRange.last =
	{
		m_start.at(threadID) + m_numX.at(threadID) - 1,
		m_numLines[1] - 1,
	};

	ENG_DISPATCH_ARGS(
		DoPreVoltageUpdatesImpl,
		abcRange,
		Tiling::Range3D<>()  // dummy value
	);
}

// Only used by Tiling engine.
void Engine_Ext_Mur_ABC::DoPreVoltageUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range2D<> abcRange;
	abcRange.first = {0, 0};
	abcRange.last = {m_numLines[0] - 1, m_numLines[1] - 1};

	//if (m_tileMap.find(tileRange))

	if (m_needRun[tileRange])
	{
		DoPreVoltageUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, abcRange, tileRange
		);
	}
}

template <typename EngineType, bool tiling>
void Engine_Ext_Mur_ABC::DoPostVoltageUpdatesImpl(
	EngineType* eng,
	Tiling::Range2D<> abcRange,  // dir2, dir3
	Tiling::Range3D<> tileRange
)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_dir1] = m_LineNr;
	pos_shift[m_dir1] = m_LineNr_Shift;

	for (unsigned int i = abcRange.first[0]; i <= abcRange.last[0]; i++)
	{
		pos[m_dir2] = i;
		pos_shift[m_dir3] = i;

		for (unsigned int j = abcRange.first[1]; j <= abcRange.last[1]; j++)
		{
			pos[m_dir3] = j;
			pos_shift[m_dir3] = j;

			if (tiling && !InsideTile(tileRange, pos) && !InsideTile(tileRange, pos_shift))
				continue;

			m_volt_dir2[i][j] +=
				m_Op_mur->m_Mur_Coeff_dir2[i][j] *
				eng->EngineType::GetVolt(m_dir2, pos_shift);

			m_volt_dir3[i][j] +=
				m_Op_mur->m_Mur_Coeff_dir3[i][j] *
				eng->EngineType::GetVolt(m_dir3, pos_shift);
		}
	}
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_Mur_ABC::DoPostVoltageUpdates(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range2D<> abcRange;
	abcRange.first = {m_start.at(threadID), 0};
	abcRange.last =
	{
		m_start.at(threadID) + m_numX.at(threadID) - 1,
		m_numLines[1] - 1,
	};

	ENG_DISPATCH_ARGS(
		DoPostVoltageUpdatesImpl,
		abcRange,
		Tiling::Range3D<>()  // dummy value
	);
}

// Only used by Tiling engine.
void Engine_Ext_Mur_ABC::DoPostVoltageUpdates(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range2D<> abcRange;
	abcRange.first = {0, 0};
	abcRange.last = {m_numLines[0] - 1, m_numLines[1] - 1};

	if (m_needRun[tileRange])
	{
		DoPostVoltageUpdatesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, abcRange, tileRange
		);
	}
}

template <typename EngineType, bool tiling>
void Engine_Ext_Mur_ABC::Apply2VoltagesImpl(
	EngineType* eng,
	Tiling::Range2D<> abcRange,  // dir2, dir3
	Tiling::Range3D<> tileRange
)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	pos[m_dir1] = m_LineNr;

	for (unsigned int i = abcRange.first[0]; i <= abcRange.last[0]; i++)
	{
		pos[m_dir2] = i;

		for (unsigned int j = abcRange.first[1]; j <= abcRange.last[1]; j++)
		{
			pos[m_dir3] = j;

			if (tiling && !InsideTile(tileRange, pos))
				continue;

			eng->EngineType::SetVolt(m_dir2, pos, m_volt_dir2[i][j]);
			eng->EngineType::SetVolt(m_dir3, pos, m_volt_dir3[i][j]);
		}
	}
}

// Used for Basic, SSE, and SSE_Compressed engines.
void Engine_Ext_Mur_ABC::Apply2Voltages(int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	Tiling::Range2D<> abcRange;
	abcRange.first = {m_start.at(threadID), 0};
	abcRange.last =
	{
		m_start.at(threadID) + m_numX.at(threadID) - 1,
		m_numLines[1] - 1,
	};

	ENG_DISPATCH_ARGS(
		Apply2VoltagesImpl,
		abcRange,
		Tiling::Range3D<>()  // dummy value
	);
}

// Only used by Tiling engine.
void Engine_Ext_Mur_ABC::Apply2Voltages(int timestep, Tiling::Range3D<> tileRange)
{
	Tiling::Range2D<> abcRange;
	abcRange.first = {0, 0};
	abcRange.last = {m_numLines[0] - 1, m_numLines[1] - 1};

	if (m_needRun[tileRange])
	{
		Apply2VoltagesImpl<Engine_Tiling>(
			(Engine_Tiling*) m_Eng, abcRange, tileRange
		);
	}
}
