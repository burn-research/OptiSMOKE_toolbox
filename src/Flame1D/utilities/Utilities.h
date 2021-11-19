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

#ifndef OpenSMOKE_PremixedLaminarFlame1D_Utilities_H
#define OpenSMOKE_PremixedLaminarFlame1D_Utilities_H

	/**
	*@brief Reads a solution from a backup file
	*@param path_file path to the backup file
	*@param x axial coordinate [m]
	*@param T temperatures [K]
	*@param P pressure profile [Pa]
	*@param m mass flow rate profile [kg/m2/s]
	*@param omega species mass fractions profiles
	*@param names_species names of species
	*/
	void ReadFromBackupFile(const boost::filesystem::path path_file, std::vector<double>& x, std::vector<double>& T, std::vector<double>& P, 
							std::vector<double>& m, std::vector< std::vector<double> >& omega, std::vector<std::string>& names_species);

	/**
	*@brief Reads a solution from a backup file
	*@param path_file path to the backup file
	*@param thermodynamicsMap thermodynamic map from which names of species can be extracted
	*@param x axial coordinate [m]
	*@param T temperatures [K]
	*@param P pressure profile [Pa]
	*@param m mass flow rate profile [kg/m2/s]
	*@param omega species mass fractions profiles
	*/
	void ReadFromBackupFile(const boost::filesystem::path path_file, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
							std::vector<double>& x, std::vector<double>& T, std::vector<double>& P, std::vector<double>& m, std::vector< std::vector<double> >& omega);

#include "Utilities.hpp"

#endif