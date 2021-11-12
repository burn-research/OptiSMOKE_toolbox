/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alessandro Stagni <alessandro.stagni@polimi.it>               |
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
|   Copyright(C) 2016, 2015 Alessandro Stagni & Alberto Cuoci             |
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

/**
 \class EoSModel_mix
 \brief Abstract base class to define eos properties for mixtures

 */

#ifndef OPENSMOKEPP_EOSMODEL_MIX_H
#define OPENSMOKEPP_EOSMODEL_MIX_H

class EosModel_mix 
{
	public:

		EosModel_mix(	const std::vector<double>& Tc, const std::vector<double>& Pc, 
						const std::vector<double>& omega, const std::vector<double>& MW);
  
		inline double a() const { return a_; }
		inline double b() const { return b_; }
		inline double A() const { return A_; }
		inline double B() const { return B_; }
		inline const std::vector<double>& ZR() const { return ZR_; }
  
		double Zmin() const;
		double Zmax() const;
 
  
	protected:
  
		double a_mix(const std::vector<double>& x);
		double b_mix(const std::vector<double>& x);
  

	private:
		
		friend class pengrobinson_mix;
    
		const std::vector<double>& Tc_;
		const std::vector<double>& Pc_;
		const std::vector<double>& omega_;
		const std::vector<double>& MW_;
		unsigned int NS_;

		double R; // gas constant [J/(mol K)]

		std::vector<double> ZR_;
		std::vector<double> Tr; //Reduced temperature
		std::vector<double> Pr; //Reduced pressure

		std::vector<double> m;
		std::vector<double> alfa;
  
		std::vector<double> a_species_;
		std::vector<double> b_species_;
		std::vector<double> A_species_;
		std::vector<double> B_species_;

		double a_;
		double b_;
		double A_;
		double B_;
  
		double a1;
		double a2;
		double a3;
  
		virtual void Solve(const double T, const double P, const std::vector<double>& x) {};
};

#include "EosModel_mix.hpp"

#endif	// OPENSMOKEPP_EOSMODEL_MIX_H
