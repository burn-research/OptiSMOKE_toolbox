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

#include "math/native-nlf-solvers/NonLinearFunctionSolver_Robust.h"

namespace OpenSMOKE
{
	ShockTubeReactor_InitialConditions *pt_ics;

	void ShockTubeReactor_InitialConditions::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  OpenSMOKE_ShockTube_InitialConditions"	<< std::endl;
		std::cout << "Object: " << name_object							<< std::endl;
		std::cout << "Error:  " << message								<< std::endl;
		std::cout << "Press a key to continue... "						<< std::endl;
		getchar();
		exit(-1);
	}

	void ShockTubeReactor_InitialConditions::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:   OpenSMOKE_ShockTube_InitialConditions"	<< std::endl;
		std::cout << "Object:  " << name_object							<< std::endl;
		std::cout << "Warning: " << message								<< std::endl;
		std::cout << "Press a key to continue... "						<< std::endl;
		getchar();
	}

	ShockTubeReactor_InitialConditions::ShockTubeReactor_InitialConditions
	(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap) :
	thermodynamicsMap_(thermodynamicsMap)
	{
		name_object	= "Default name";	// Object Name
		ALFAMAX		= 30.;				// Maximum temperature ratio
		pt_ics		= this;


		h1 = u1 = p1 = T1 = rho1_ = pm1 = Cp1 = gamma1 = M1 = group1 = 0.;
		h2 = u2 = p2 = T2 = rho2 = pm2 = Cp2 = gamma2 = M2 = group2 = 0.;
		h5 = u5 = p5 = T5 = rho5 = pm5 = Cp5 = gamma5 = M5 = 0.;
	}

	void ShockTubeReactor_InitialConditions::SetName(const std::string name)
	{
		name_object	= name;				// Object Name
	}

	double ShockTubeReactor_InitialConditions::AlfaFirstGuessIncidentShock_CaseA()
	{
		double alfa	 =	( gamma1*(M1*M1)-0.50*(gamma1-1.) ) * ( 0.50*(gamma1-1.)*(M1*M1)+1. ) / 
						( boost::math::pow<2>(0.50*(gamma1+1.)*M1) );
		return alfa;
	}

	double ShockTubeReactor_InitialConditions::AlfaFirstGuessIncidentShock_CaseB()
	{
		const double gamma	= 1.40;
		const double M		= 1.;
		const double alfa	=	( gamma*(M*M)-0.50*(gamma-1.) ) * ( 0.50*(gamma-1.)*(M*M)+1. ) / 
								( boost::math::pow<2>(0.50*(gamma+1.)*M) );
		return alfa;
	}

	double ShockTubeReactor_InitialConditions::AlfaFirstGuessReflectedShock_CaseB()
	{
		const double eta		= rho2/rho1_;
		const double csi		= eta*(M1*M1*(eta-1.)*gamma1+eta)/(M1*M1*(eta-1.)*(gamma1-1.)+eta);
		const double T5overT1	= 1.+(M1*M1)*(gamma1-1.)*(csi-1.)*(eta-1.)/(eta*(csi-eta));
		const double alfa		= T5overT1*T1/T2;
		return alfa;
	}

	double ShockTubeReactor_InitialConditions::AlfaMaxIncidentShock_CaseA()
	{
		return (group1+1.)*(group1+1.)/4./group1;
	}

	double ShockTubeReactor_InitialConditions::AlfaMaxIncidentShock_CaseB()
	{
		const double csi = rho2*u1*u1/p2;
		
		if (csi >= 4.)	return ALFAMAX;
		else			return (sqrt(csi)+2.)/sqrt(csi)/(4.-csi);
	}

	double ShockTubeReactor_InitialConditions::AlfaMaxReflectedShock_CaseB()
	{
		if (Urs == 0.)	
		{
			return ALFAMAX;
		}
		else
		{
			const double u2p	= Urs + u1 - u2;
			const double coeff	= rho2*u2p*u2p/p2;
			return (1.+coeff)*(1.+coeff)/4./coeff; 
		}
	}

	void ShockTubeReactor_InitialConditions::Set_IncidentShock_CaseA(const double _UShock, const double _T1, const double _p1, OpenSMOKE::OpenSMOKEVectorDouble& _omega1)
	{
		iIncidentShock = true;
	
		UShock	= _UShock;
		T1		= _T1;
		p1		= _p1;
		omega1	= _omega1;
		u1		= UShock;
	
		ChangeDimensions(omega1.Size(), &x1, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x1.GetHandle(), pm1, omega1.GetHandle());
		rho1_	= p1*pm1/PhysicalConstants::R_J_kmol/T1;

		thermodynamicsMap_.SetTemperature(T1);
		thermodynamicsMap_.SetPressure(p1);

		Cp1 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x1.GetHandle());
		Cp1 /= pm1;

		h1 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x1.GetHandle());
		h1 /= pm1;
		
		gamma1	= Cp1*pm1/(Cp1*pm1-PhysicalConstants::R_J_kmol);
		M1		= u1*std::sqrt(rho1_/gamma1/p1);

		group1  = rho1_*u1*u1/p1;

		omega2	= omega1;
		pm2		= pm1;
		x2		= x1;
	}

	void ShockTubeReactor_InitialConditions::Set_IncidentShock_CaseB(const double _UShock, const double _T2, const double _p2, OpenSMOKE::OpenSMOKEVectorDouble& _omega2)
	{
		iIncidentShock = true;
	
		UShock	= _UShock;
		T2		= _T2;
		p2		= _p2;
		omega2	= _omega2;
		u1		= UShock;
	
		ChangeDimensions(omega2.Size(), &x2, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x2.GetHandle(), pm2, omega2.GetHandle());
		rho2	= p2*pm2/PhysicalConstants::R_J_kmol/T2;

		thermodynamicsMap_.SetTemperature(T2);
		thermodynamicsMap_.SetPressure(p2);

		Cp2 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x2.GetHandle());
		Cp2 /= pm2;

		h2 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x2.GetHandle());
		h2 /= pm2;
		
		gamma2	= Cp2*pm2/(Cp2*pm2-PhysicalConstants::R_J_kmol);

		omega1	= omega2;
		x1		= x2;
		pm1		= pm2;
	}

	void ShockTubeReactor_InitialConditions::Set_ReflectedShock_CaseA(const double _T5, const double _p5, OpenSMOKE::OpenSMOKEVectorDouble& _omega5)
	{
		iIncidentShock = false;

		Urs		= 0.;
		u5		= 0.;
		M5		= 0.;

		T5		= _T5;
		p5		= _p5;
		omega5	= _omega5;

		ChangeDimensions(omega5.Size(), &x5, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x5.GetHandle(), pm5, omega5.GetHandle());
		rho5	= p5*pm5/PhysicalConstants::R_J_kmol/T5;

		thermodynamicsMap_.SetTemperature(T5);
		thermodynamicsMap_.SetPressure(p5);

		Cp5 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x5.GetHandle());
		Cp5 /= pm5;

		h5 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x5.GetHandle());
		h5 /= pm5;

		gamma5	= Cp5*pm5/(Cp5*pm5-PhysicalConstants::R_J_kmol);
	}

	void ShockTubeReactor_InitialConditions::Set_ReflectedShock_CaseB(GasStream& shockStream, const double _UShock, const double _Urs, const double _T1, const double _p1, OpenSMOKE::OpenSMOKEVectorDouble& _omega1)
	{
		Set_IncidentShock_CaseA(_UShock, _T1, _p1, _omega1);
		Solve_IncidentShock_CaseA(shockStream);

		iIncidentShock = false;

		Urs		= _Urs;
		omega5	= omega1;
		x5		= x1;
		pm5		= pm1;
	}

	void ShockTubeReactor_InitialConditions::Set_ReflectedShock_CaseC(GasStream &shockStream, const double _UShock, const double _Urs, const double _T2, const double _p2, OpenSMOKE::OpenSMOKEVectorDouble& _omega2)
	{
		Set_IncidentShock_CaseB(_UShock, _T2, _p2, _omega2);
		Solve_IncidentShock_CaseB(shockStream);

		iIncidentShock = false;

		Urs		= _Urs;
		omega5	= omega1;
		x5		= x1;
		pm5		= pm1;
	}

	double NLS_IncidentShock_CaseA(double alfa)
	{
		return pt_ics->IncidentShock_CaseA(alfa);
	}

	double NLS_IncidentShock_CaseB(double alfa)
	{
		return pt_ics->IncidentShock_CaseB(alfa);
	}

	double NLS_ReflectedShock_CaseB(double alfa)
	{
		return pt_ics->ReflectedShock_CaseB(alfa);
	}

	double ShockTubeReactor_InitialConditions::IncidentShock_CaseA(double alfa)
	{
		Beta = 0.50*( (1.+group1)+sqrt(boost::math::pow<2>(1.+group1) - 4.*alfa*group1) );

		T2	= T1*alfa;
		p2	= p1*Beta;

		thermodynamicsMap_.SetTemperature(T2);
		thermodynamicsMap_.SetPressure(p2);

		ChangeDimensions(omega2.Size(), &x2, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x2.GetHandle(), pm2, omega2.GetHandle());
		h2 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x2.GetHandle());
		h2 /= pm2;

		double f = h1 + u1*u1/2.*(1.-alfa*alfa/Beta/Beta) - h2;

		return f;
	}

	double ShockTubeReactor_InitialConditions::IncidentShock_CaseB(double alfa)
	{
		double csi = rho2*u1*u1/p2;

		Beta = 0.50 * ( (1.+csi*alfa) + sqrt( boost::math::pow<2>(1.+csi*alfa) -4.*csi*alfa*alfa) );

		T1	= T2/alfa;
		p1	= p2/Beta;

		thermodynamicsMap_.SetTemperature(T1);
		thermodynamicsMap_.SetPressure(p1);

		ChangeDimensions(omega1.Size(), &x1, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x1.GetHandle(), pm1, omega1.GetHandle());
		h1 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x1.GetHandle());
		h1 /= pm1;
	
		double f = h1 + u1*u1/2.*(1.-alfa*alfa/Beta/Beta) - h2;

		return f;
	}

	void ShockTubeReactor_InitialConditions::Solve_IncidentShock_CaseA(GasStream &shockStream)
	{
		#if OPENSMOKE_USE_BZZMATH == 0
		OpenSMOKE::NonLinearFunctionSolver_Robust m(AlfaFirstGuessIncidentShock_CaseA(), NLS_IncidentShock_CaseA, 1., AlfaMaxIncidentShock_CaseA());
		#else
		BzzFunctionRootRobust m(AlfaFirstGuessIncidentShock_CaseA(), NLS_IncidentShock_CaseA, 1., AlfaMaxIncidentShock_CaseA());
		#endif

		m.OnlyOneRoot();
		m();
	
		double alfa		= m.GetTSolution();
		double residual	= m.GetYSolution();

		T2		= alfa*T1;
		p2		= Beta*p1;
		rho2	= p2*pm2/PhysicalConstants::R_J_kmol/T2;
		u2		= rho1_*u1/rho2;

		ChangeDimensions(omega2.Size(), &x2, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x2.GetHandle(), pm2, omega2.GetHandle());
		thermodynamicsMap_.SetTemperature(T2);
		thermodynamicsMap_.SetPressure(p2);
		Cp2 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x2.GetHandle());
		Cp2 /= pm2;
		h2 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x2.GetHandle());
		h2 /= pm2;

		gamma2	= Cp2*pm2/(Cp2*pm2-PhysicalConstants::R_J_kmol);
		M2		= u2*std::sqrt(rho2/gamma2/p2);
	
		group2  = rho2*u2*u2/p2;

	//	mix->SpeciesViscosityFromFitting(T1); 
	//	mu1		= mix->MixViscosity_FromMolarFractions(x1);

	//	mix->SpeciesViscosityFromFitting(T2); 
	//	mu2		= mix->MixViscosity_FromMolarFractions(x2);

		shockStream.Set(T2, p2, u2, omega2);
	}

	void ShockTubeReactor_InitialConditions::Solve_IncidentShock_CaseB(GasStream &shockStream)
	{
		#if OPENSMOKE_USE_BZZMATH == 0
		OpenSMOKE::NonLinearFunctionSolver_Robust m(AlfaFirstGuessIncidentShock_CaseB(), NLS_IncidentShock_CaseB, 1., AlfaMaxIncidentShock_CaseB());
		#else
		BzzFunctionRootRobust m(AlfaFirstGuessIncidentShock_CaseB(), NLS_IncidentShock_CaseB, 1., AlfaMaxIncidentShock_CaseB());
		#endif

		m.OnlyOneRoot();
		m();
	
		double alfa		= m.GetTSolution();
		double residual	= m.GetYSolution();

		T1		= T2/alfa;
		p1		= p2/Beta;
		rho1_	= p1*pm1/PhysicalConstants::R_J_kmol/T1;

		ChangeDimensions(omega1.Size(), &x1, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x1.GetHandle(), pm1, omega1.GetHandle());
		thermodynamicsMap_.SetTemperature(T1);
		thermodynamicsMap_.SetPressure(p1);
		Cp1 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x1.GetHandle());
		Cp1 /= pm1;
		h1 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x1.GetHandle());
		h1 /= pm1;

		gamma1	= Cp1*pm1/(Cp1*pm1-PhysicalConstants::R_J_kmol);

		u2		= rho1_*u1/rho2;
		M1		= u1*std::sqrt(rho1_/gamma1/p1);
		M2		= u2*std::sqrt(rho2/gamma2/p2);
	
		group1  = rho1_*u1*u1/p1;
		group2  = rho2*u2*u2/p2;

	//	mix->SpeciesViscosityFromFitting(T1); 
	//	mu1		= mix->MixViscosity_FromMolarFractions(x1);

	//	mix->SpeciesViscosityFromFitting(T2); 
	//	mu2		= mix->MixViscosity_FromMolarFractions(x2);

		shockStream.Set(T2, p2, u2, omega2);
	}

	double ShockTubeReactor_InitialConditions::ReflectedShock_CaseB(double alfa)
	{
		double f;
		double u2p;

		double eta = rho2*rho1_;
		double csi = eta*(M1*M1*(eta-1.)*gamma1+eta)/(M1*M1*(eta-1.)*(gamma1-1.)+eta);

		if (Urs == 0.)
		{
			double coeff	= 1.+rho2*(u1-u2)*(u1-u2)/p2 + alfa;
			Beta			= 0.50*(coeff+sqrt(coeff*coeff-4.*alfa));
			f				= h2+0.50*(u1-u2)*(u1-u2)*(1.+alfa/Beta)/(1.-alfa/Beta) - h5;
		}
		else
		{
			u2p		= Urs + u1 - u2;
			double coeff	= rho2*u2p*u2p/p2;
			Beta			= 0.50*( (1.+coeff) + sqrt(boost::math::pow<2>(1.+coeff) -4.*alfa*coeff));
		}

		T5	= alfa*T2;
		p5	= Beta*p2;

		thermodynamicsMap_.SetTemperature(T5);
		thermodynamicsMap_.SetPressure(p5);

		ChangeDimensions(omega5.Size(), &x5, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x5.GetHandle(), pm5, omega5.GetHandle());
		h5 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x5.GetHandle());
		h5 /= pm5;
	
		if (Urs == 0.)
			f				= h2+0.50*(u1-u2)*(u1-u2)*(1.+alfa/Beta)/(1.-alfa/Beta) - h5;
		else
			f				= h2+0.50*u2p*u2p*(1.-alfa*alfa/Beta/Beta)-h5;

		return f;
	}

	void ShockTubeReactor_InitialConditions::Solve_ReflectedShock_CaseA(GasStream &shockStream)
	{
		shockStream.Set(T5, p5, u5, omega5);
	}

	void ShockTubeReactor_InitialConditions::Solve_ReflectedShock_CaseB(GasStream &shockStream)
	{
		#if OPENSMOKE_USE_BZZMATH == 0
		OpenSMOKE::NonLinearFunctionSolver_Robust m(AlfaFirstGuessReflectedShock_CaseB(), NLS_ReflectedShock_CaseB, 1., AlfaMaxReflectedShock_CaseB());
		#else
		BzzFunctionRootRobust m(AlfaFirstGuessReflectedShock_CaseB(), NLS_ReflectedShock_CaseB, 1., AlfaMaxReflectedShock_CaseB());
		#endif

		m.OnlyOneRoot();
		m();
	
		double alfa		= m.GetTSolution();
		double residual	= m.GetYSolution();

		T5		= alfa*T2;
		p5		= Beta*p2;
		rho5	= p5*pm5/PhysicalConstants::R_J_kmol/T5;

		ChangeDimensions(omega5.Size(), &x5, true);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x5.GetHandle(), pm5, omega5.GetHandle());
		thermodynamicsMap_.SetTemperature(T5);
		thermodynamicsMap_.SetPressure(p5);
		Cp5 = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x5.GetHandle());
		Cp5 /= pm5;
		h5 = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x5.GetHandle());
		h5 /= pm5;

		gamma5	= Cp5*pm5/(Cp5*pm5-PhysicalConstants::R_J_kmol);

		if (Urs == 0.)	Urs = alfa/Beta*(u1-u2)/(1.-alfa/Beta);

		u5		= rho2/rho5*(Urs+u1-u2);
		M5		= u5*std::sqrt(rho5/gamma5/p5);

	//	mix->SpeciesViscosityFromFitting(T5); 
	//	mu5		= mix->MixViscosity_FromMolarFractions(x5);

		shockStream.Set(T5, p5, u5, omega5);
	}

	void ShockTubeReactor_InitialConditions::Solve_ReflectedShock_CaseC(GasStream &shockStream)
	{
		Solve_ReflectedShock_CaseB(shockStream);
	}

	void ShockTubeReactor_InitialConditions::VideoSummary()
	{
		std::cout.setf(std::ios::scientific);

		if (iIncidentShock == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------"	<< std::endl;
			std::cout << std::setw(14) << std::left << "Shock" << std::setw(20) << std::left << "Before" << "After"	<< std::endl;
			std::cout << "-----------------------------------------------------------------------------"	<< std::endl;
			std::cout << std::setw(14) << std::left << "T [K]"			<< std::setw(20) << std::left << T1				<< T2			<< std::endl;
			std::cout << std::setw(14) << std::left << "P [atm]"		<< std::setw(20) << std::left << p1/101325.		<< p2/101325.	<< std::endl;
			std::cout << std::setw(14) << std::left << "rho [kg/m3]"	<< std::setw(20) << std::left << rho1_			<< rho2			<< std::endl;
			std::cout << std::setw(14) << std::left << "u [m/s]"		<< std::setw(20) << std::left << u1				<< u2			<< std::endl;
			std::cout << std::setw(14) << std::left << "M [-]"			<< std::setw(20) << std::left << M1				<< M2			<< std::endl;
			std::cout << std::setw(14) << std::left << "h [J/kg]"		<< std::setw(20) << std::left << h1				<< h2			<< std::endl;
			std::cout << std::setw(14) << std::left << "Cp [J/kg/K]"	<< std::setw(20) << std::left << Cp1			<< Cp2			<< std::endl;
			std::cout << std::setw(14) << std::left << "gamma [-]"		<< std::setw(20) << std::left << gamma1			<< gamma2		<< std::endl;
		}
		else
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------"				<< std::endl;
			std::cout << std::setw(14) << std::left << "Shock" << std::setw(20) << std::left << "Before" << std::setw(20) << std::left << "After" << "Reflected"	<< std::endl;
			std::cout << "-----------------------------------------------------------------------------"				<< std::endl;
			std::cout << std::setw(14) << std::left << "T [K]"			<< std::setw(20) << std::left << T1				<< std::setw(20) << std::left << T2			<< T5			<< std::endl;
			std::cout << std::setw(14) << std::left << "P [atm]"		<< std::setw(20) << std::left << p1/101325.		<< std::setw(20) << std::left << p2/101325.	<< p5/101325.	<< std::endl;
			std::cout << std::setw(14) << std::left << "rho [kg/m3]"	<< std::setw(20) << std::left << rho1_			<< std::setw(20) << std::left << rho2		<< rho5			<< std::endl;
			std::cout << std::setw(14) << std::left << "u [m/s]"		<< std::setw(20) << std::left << u1				<< std::setw(20) << std::left << u2			<< u5			<< std::endl;
			std::cout << std::setw(14) << std::left << "M [-]"			<< std::setw(20) << std::left << M1				<< std::setw(20) << std::left << M2			<< M5			<< std::endl;
			std::cout << std::setw(14) << std::left << "h [J/kg]"		<< std::setw(20) << std::left << h1				<< std::setw(20) << std::left << h2			<< h5			<< std::endl;
			std::cout << std::setw(14) << std::left << "Cp [J/kg/K]"	<< std::setw(20) << std::left << Cp1			<< std::setw(20) << std::left << Cp2		<< Cp5			<< std::endl;
			std::cout << std::setw(14) << std::left << "gamma [-]"		<< std::setw(20) << std::left << gamma1			<< std::setw(20) << std::left << gamma2		<< gamma5		<< std::endl;
		}
	}

	double viscosity(const double mu0, const double beta, const double T)
	{
		return mu0 * std::pow(T/300.,beta);
	}

	double ShockTubeReactor_InitialConditions::LengthBoundaryLayerCorrection(const double d, const double mu300)
	{
		const double mu1 = viscosity(mu300, 0.6756, T1);
		const double mu2 = viscosity(mu300, 0.6756, T2);
		
		const double rhoWall = rho1_*p2/p1;
		const double muWall  = mu1;

		const double C = std::pow( (rho2/rhoWall)*(mu2/muWall), 0.37);
		const double W = rho2/rho1_;
		
		double Z = (gamma1+1.)/(gamma1-1.);
		if (W>Z) Z = W;
	
		const double B = 1.59*C*(1.+(1.796+0.802*W)/(Z*W-1.));
		const double lm = boost::math::pow<2>(d*rho2/4./B/rhoWall)/(W-1.)*(u2/(muWall/rhoWall));

		std::cout << std::endl;
		std::cout << " * Gas viscosity before the shock: " << mu1 << " [kg/m/s]" << std::endl;
		std::cout << " * Boundary layer parameter beta:  " << B << " [-]" << std::endl;
		std::cout << " * Limiting separation between shock and contact surface:  " << lm << " m" << std::endl;
		std::cout << std::endl;

		return lm;
	}
}
