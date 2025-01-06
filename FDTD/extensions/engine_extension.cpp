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

#include "engine_extension.h"
#include "operator_extension.h"

#include "FDTD/engine.h"

Engine_Extension::Engine_Extension(Operator_Extension* op_ext)
{
	m_Op_ext = op_ext;
	m_Eng = NULL;
	m_Priority = ENG_EXT_PRIO_DEFAULT;
	m_NrThreads = 1;
	m_TilingSupported = false;
}

Engine_Extension::~Engine_Extension()
{
}

void Engine_Extension::SetNumberOfThreads(int nrThread)
{
	if (nrThread<1)
		return;
	m_NrThreads=nrThread;
}

string Engine_Extension::GetExtensionName() const
{
	if (m_Op_ext)
		return m_Op_ext->GetExtensionName();
	else
		return "Unknown Extension";
}

void Engine_Extension::TilingUnsupportedError(void)
{
	throw std::runtime_error(
		GetExtensionName() + " does not support the tiling engine."
	);
}

void Engine_Extension::InitializeTiling(const Tiling::Plan3D& plan)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::DoPreVoltageUpdates(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		DoPreVoltageUpdates();
}

void Engine_Extension::DoPreVoltageUpdates(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::DoPostVoltageUpdates(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		DoPostVoltageUpdates();
}

void Engine_Extension::DoPostVoltageUpdates(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::Apply2Voltages(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		Apply2Voltages();
}

void Engine_Extension::Apply2Voltages(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::DoPreCurrentUpdates(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		DoPreCurrentUpdates();
}

void Engine_Extension::DoPreCurrentUpdates(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::DoPostCurrentUpdates(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		DoPostCurrentUpdates();
}

void Engine_Extension::DoPostCurrentUpdates(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

void Engine_Extension::Apply2Current(int threadID)
{
	//if this method gets called the derived extension obviously doesn't support multithrading, calling non-MT method...
	if (threadID==0)
		Apply2Current();
}

void Engine_Extension::Apply2Current(int timestep, Tiling::Range3D<> range)
{
	// If this method gets called, the derived extension either
	// doesn't need to do anything in this method or the extension
	// doesn't support spatial and temporial tiling. In the former
	// case, it's a no-op. In the latter case, the program terminates.
	if (!m_TilingSupported)
		TilingUnsupportedError();
}

bool Engine_Extension::operator< (const Engine_Extension& other)
{
	return (GetPriority()<other.GetPriority());
}
