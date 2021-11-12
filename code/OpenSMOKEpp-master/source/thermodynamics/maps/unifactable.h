/*-----------------------------------------------------------------------*\
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
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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

/*! \class UnifacTable
    \brief Class to Manage Unifac table coefficients 

    Class to Manage Unifac table coefficients
*/

#ifndef UNIFACTABLE_H
#define UNIFACTABLE_H

#include "Eigen/Dense"

namespace fs = boost::filesystem;
namespace lc = libconfig;

class UnifacTable
{
public:
  UnifacTable();
  
  //Access functions
  inline const vector<string>& groupname() const                    {return groupname_;};
  inline const vector<int>& maingroup() const                       {return maingroup_;};
  inline const vector<int>& subgroup() const                        {return subgroup_;};
  inline const vector<double>& R() const                            {return R_;};
  inline const vector<double>& Q() const                            {return Q_;};
  inline const Eigen::MatrixXd& a() const                           {return a_;};
  
  void SetupTable(const fs::path& map_path);
  void SetupInteractionCoefficients(const fs::path& map_path);
private:
  vector<string> groupname_;
  vector<int> maingroup_;
  vector<int> subgroup_;
  vector<double> R_;
  vector<double> Q_;
  
  Eigen::MatrixXd a_;
  
  
  
  void TableAllocation(int ngroups);
};

#include "unifactable.hpp"



#endif
