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

#ifndef ENGINE_EXT_EXCITATION_H
#define ENGINE_EXT_EXCITATION_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_Excitation;

class Engine_Ext_Excitation : public Engine_Extension
{
public:
	Engine_Ext_Excitation(Operator_Ext_Excitation* op_ext);
	virtual ~Engine_Ext_Excitation();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

	virtual void Apply2Voltages(
		int timestep,
		Tiling::Range3D<> range
	);

	virtual void Apply2Current(
		int timestep,
		Tiling::Range3D<> range
	);

protected:
	template <typename EngineType, bool tiling=false>
	void Apply2VoltagesImpl(
		EngineType *eng,
		int numTS,
		Tiling::Range3D<> range
	);

	template <typename EngineType, bool tiling=false>
	void Apply2CurrentImpl(
		EngineType *eng,
		int numTS,
		Tiling::Range3D<> range
	);

	bool InsideTile(
		Tiling::Range3D<> range,
		unsigned int exc[3]
	);

	Operator_Ext_Excitation* m_Op_Exc;
};

#endif // ENGINE_EXT_EXCITATION_H
