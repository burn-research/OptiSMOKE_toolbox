/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_ShockTubeReactor_InitialConditions_H
#define	OpenSMOKE_ShockTubeReactor_InitialConditions_H

#include "math/OpenSMOKEVector.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"

namespace OpenSMOKE
{
	class GasStream
	{
	public:

		GasStream() {};
		void Set(const double T, const double P, const double v, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
		{
			T_ = T;
			P_ = P;
			v_ = v;
			omega_ = omega;
		}

		double T() const { return T_; }
		double P() const { return P_; }
		double v() const { return v_; }
		const OpenSMOKE::OpenSMOKEVectorDouble& omega() const { return omega_; }

	private:
		double T_;
		double P_;
		double v_;
		OpenSMOKE::OpenSMOKEVectorDouble omega_;
	};

	//!  A class for simulating perfectly stirred reactors with constant pressure in adiabatic conditions
	/*!
		 The purpose of this class is to simulate a perfectly stirred reactor with constant pressure, in adiabatic conditions
		 The conservation equations of species and energy are solved in terms of mass fractions and temperature, respectively.
	*/

	class ShockTubeReactor_InitialConditions
	{
	public:
		ShockTubeReactor_InitialConditions(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap);
	
		void SetName(const std::string name);
		
		void Set_IncidentShock_CaseA(const double _UShock, const double _T1, const double _p1, OpenSMOKE::OpenSMOKEVectorDouble& _omega1);
		void Set_IncidentShock_CaseB(const double _UShock, const double _T2, const double _p2, OpenSMOKE::OpenSMOKEVectorDouble& _omega2);
		void Set_ReflectedShock_CaseA(const double _T5, const double _p5, OpenSMOKE::OpenSMOKEVectorDouble& _omega5);
		void Set_ReflectedShock_CaseB(GasStream &shockStream, const double _UShock, const double _Urs, const double _T1, const double _p1, OpenSMOKE::OpenSMOKEVectorDouble& _omega1);
		void Set_ReflectedShock_CaseC(GasStream &shockStream, const double _UShock, const double _Urs, const double _T2, const double _p2, OpenSMOKE::OpenSMOKEVectorDouble& _omega2);

		void Solve_IncidentShock_CaseA(GasStream &shockStream);
		void Solve_IncidentShock_CaseB(GasStream &shockStream);
		void Solve_ReflectedShock_CaseA(GasStream &shockStream);
		void Solve_ReflectedShock_CaseB(GasStream &shockStream);
		void Solve_ReflectedShock_CaseC(GasStream &shockStream);

		void	VideoSummary();
		double rho1() const { return rho1_; }
		double	LengthBoundaryLayerCorrection(const double d, const double mu300 = 1.85e-5);


	private:

		double u1;
		double p1;
		double T1;
		double rho1_;
		double h1;
		double pm1;
		double Cp1;
		double gamma1;
		double M1;
		double group1;
	//	double mu1;
		OpenSMOKE::OpenSMOKEVectorDouble omega1;
		OpenSMOKE::OpenSMOKEVectorDouble x1;

		double u2;
		double p2;
		double T2;
		double rho2;
		double h2;
		double pm2;
		double Cp2;
		double gamma2;
		double M2;
		double group2;
	//	double mu2;
		OpenSMOKE::OpenSMOKEVectorDouble omega2;
		OpenSMOKE::OpenSMOKEVectorDouble x2;

		double Urs;
		double u5;
		double p5;
		double T5;
		double rho5;
		double h5;
		double pm5;
		double Cp5;
		double gamma5;
		double M5;
	//	double mu5;
		OpenSMOKE::OpenSMOKEVectorDouble omega5;
		OpenSMOKE::OpenSMOKEVectorDouble x5;


	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
		double UShock;
		double Beta;
		bool iIncidentShock;
		double ALFAMAX;

	private:

		std::string name_object;
		void ErrorMessage(const std::string message);
		void WarningMessage(const std::string message);

	private:

		double AlfaFirstGuessIncidentShock_CaseA();
		double AlfaFirstGuessIncidentShock_CaseB();
		double AlfaFirstGuessReflectedShock_CaseB();
		double AlfaMaxIncidentShock_CaseA();
		double AlfaMaxIncidentShock_CaseB();
		double AlfaMaxReflectedShock_CaseB();

	public:

		double IncidentShock_CaseA(double alfa);
		double IncidentShock_CaseB(double alfa);
		double ReflectedShock_CaseB(double alfa);
	};
}

#include "ShockTubeReactor_InitialConditions.hpp"

#endif	/* OpenSMOKE_ShockTubeReactor_InitialConditions_H */

