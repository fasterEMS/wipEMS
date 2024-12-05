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

#ifndef ENGINE_EXT_MUR_ABC_H
#define ENGINE_EXT_MUR_ABC_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"
#include "tools/arraylib/array_ij.h"
#include "tools/tiling/tiling.h"

class Operator_Ext_Mur_ABC;

class Engine_Ext_Mur_ABC : public Engine_Extension
{
public:
	Engine_Ext_Mur_ABC(Operator_Ext_Mur_ABC* op_ext);
	virtual ~Engine_Ext_Mur_ABC();

	virtual void InitializeTiling(const Tiling::Plan3D& plan);

	virtual void SetNumberOfThreads(int nrThread);

	virtual void DoPreVoltageUpdates() {Engine_Ext_Mur_ABC::DoPreVoltageUpdates(0);}
	virtual void DoPreVoltageUpdates(int threadID);
	virtual void DoPostVoltageUpdates() {Engine_Ext_Mur_ABC::DoPostVoltageUpdates(0);}
	virtual void DoPostVoltageUpdates(int threadID);
	virtual void Apply2Voltages() {Engine_Ext_Mur_ABC::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);

	virtual void DoPreVoltageUpdates(int timestep, Tiling::Range3D<> range);
	virtual void DoPostVoltageUpdates(int timestep, Tiling::Range3D<> range);
	virtual void Apply2Voltages(int timestep, Tiling::Range3D<> range);

protected:
	bool InsideTile(Tiling::Range3D<> range, unsigned int target[3]);

	template <typename EngineType, bool tiling=false>
	void DoPreVoltageUpdatesImpl(
		EngineType* eng,
		Tiling::Range3D<> abcRange,
		Tiling::Range3D<> tileRange
	);

	template <typename EngineType, bool tiling=false>
	void DoPostVoltageUpdatesImpl(
		EngineType* eng,
		Tiling::Range3D<> abcRange,
		Tiling::Range3D<> tileRange
	);

	template <typename EngineType, bool tiling=false>
	void Apply2VoltagesImpl(
		EngineType* eng,
		Tiling::Range3D<> abcRange,
		Tiling::Range3D<> tileRange
	);

	Operator_Ext_Mur_ABC* m_Op_mur;

	inline bool IsActive() {if (m_Eng->GetNumberOfTimesteps()<m_start_TS) return false; return true;}
	unsigned int m_start_TS;

	int m_ny;
	int m_nyP,m_nyPP;
	unsigned int m_LineNr;
	int m_LineNr_Shift;
	unsigned int m_numLines[2];

	vector<unsigned int> m_start;
	vector<unsigned int> m_numX;

	ArrayLib::ArrayIJ<FDTD_FLOAT>& m_Mur_Coeff_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>& m_Mur_Coeff_nyPP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> m_volt_nyP; //n+1 direction
	ArrayLib::ArrayIJ<FDTD_FLOAT> m_volt_nyPP; //n+2 direction

private:
	void InitializeTilingImpl(Tiling::Range3D<> range);

	std::unordered_map<
		Tiling::Range3D<>, bool, Tiling::Range3DHasher<>
	> m_needRun;
};

#endif // ENGINE_EXT_MUR_ABC_H
