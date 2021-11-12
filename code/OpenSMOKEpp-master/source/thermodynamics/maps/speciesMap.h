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

/*! \class speciesMap
    \brief Map of species from which properties of pure compounds can be created, and mixtures can be set

    Map of species from which properties of pure compounds can be created, and mixtures can be set
*/

#ifndef SPECIESMAP_H_
#define SPECIESMAP_H_

#include <iostream>
#include <string>
#include <sstream>
#include "boost/filesystem.hpp"
#include "libconfig.h++"
#include "math/OpenSMOKEUtilities.h"

//Property models
#include "thermodynamics/pure/eos/eos.h"
#include "thermodynamics/pure/densityL/densityL.h"
#include "thermodynamics/pure/vaporpressure/vaporpressure.h"
#include "thermodynamics/pure/hov/heatOfVaporization.h"
#include "thermodynamics/pure/viscosity/viscosity.h"
#include "thermodynamics/pure/conductivity/conductivity.h"
#include "thermodynamics/pure/heatcapacity/heatcapacity.h"
#include "thermodynamics/pure/surfacetension/surfacetension.h"
#include "thermodynamics/pure/fugacity/fugacity.h"

#include "thermodynamics/maps/unifactable.h"

namespace fs = boost::filesystem;
namespace lc = libconfig;

class speciesMap
{
public:
  speciesMap(const fs::path& config_folder);
  ~speciesMap();


  double rhoL(string species, double T, double P) const;
  double d_rhoL_over_dT(string species, double T, double P) const;
  double pVap(string species, double T, double P) const;
  double dpVapdT(string species, double T, double P) const;
  double deltaHv(string species, double T) const;
  double etaL(string species, double T) const;
  double lambdaL(string species, double T) const;
  double cpL(string species, double T) const;
  double cpG(string species, double T) const;
  double etaG(string species, double T) const;
  double lambdaG(string species, double T) const;
  double sigma(string species, double T) const;

  //Access functions to independent properties
  inline int NS() const                                               {return NS_;};
  inline const std::vector<std::string>& name()  const                      {return name_;};
  inline const std::vector<double>& MW() const                              {return MW_;};
  inline const std::vector<double>& Tc() const                              {return Tc_;};
  inline const std::vector<double>& Pc()  const                             {return Pc_;};
  inline const std::vector<double>& Rhoc() const                            {return Rhoc_;};
  inline const std::vector<double>& Tnbp()    const                         {return Tnbp_;};
  inline const std::vector<double>& dipolMoment()         const             {return dipolMoment_;};
  inline const std::vector<double>& omega()    const                        {return omega_;};
  
  //Access functions to models
  inline const std::vector<DensityModel*>& dM() const                       {return dM_;};
  inline const std::vector<VaporPressureModel*>& vpM() const                {return vpM_;};
  
  //Access functions to UNIFAC model
  inline const std::vector<std::vector<double> >& unifacgroup_R()  const            {return unifacgroup_R_;};
  inline const std::vector<std::vector<double> >& unifacgroup_Q()  const            {return unifacgroup_Q_;};
  inline const std::vector<std::vector<int> >&    unifacgroup_subgroup()  const     {return unifacgroup_subgroup_;};
  inline const std::vector<std::vector<int> >&    unifacgroup_maingroup()  const    {return unifacgroup_maingroup_;};
  inline const std::vector<std::vector<int> >&    unifacgroup_number()  const       {return unifacgroup_number_;};
  const Eigen::MatrixXd& interaction_table() const                          {return unifactable_->a();};
  
  double specific_volume_stp(string species)    const;
  inline const std::vector<double>& specific_volume_stp()    const          {return specific_volume_stp_;};

private:

  void ReadSpeciesNames(const fs::path& map_path);
  void MemoryAllocation();
  void ReadSpeciesProperties(const fs::path& map_path);
  void SetupProperties();

  int NS_;
  const fs::path& config_folder_;
  

  //Constant properties
  std::vector<std::string> name_;
  std::vector<double> MW_;
  std::vector<double> Tc_;
  std::vector<double> Pc_;
  std::vector<double> Rhoc_;
  std::vector<double> Tnbp_;
  std::vector<double> dipolMoment_;
  std::vector<double> omega_;
  std::vector<double> specific_volume_stp_;

  //Dependent properties: correlation indices
  std::vector<int> pVap_index_;
  std::vector<int> dHv_index_;
  std::vector<int> etaL_index_;
  std::vector<int> lambdaL_index_;
  std::vector<int> cpG_index_;
  std::vector<int> cpL_index_;
  std::vector<int> lambdaG_index_;
  std::vector<int> etaG_index_;
  std::vector<int> rhoL_index_;
  std::vector<int> sigma_index_;

  //Model objects
  std::vector<DensityModel*> dM_;
  std::vector<VaporPressureModel*> vpM_;
  std::vector<HeatOfVaporizationModel*> hovM_;
  std::vector<ViscosityModel*> etaLM_;
  std::vector<ConductivityModel*> lambdaLM_;
  std::vector<CapacityModel*> cpGM_;
  std::vector<CapacityModel*> cpLM_;
  std::vector<ViscosityModel*> etaGM_;
  std::vector<ConductivityModel*> lambdaGM_;
  std::vector<SurfaceTensionModel*> sigmaM_;
  UnifacTable* unifactable_;
  
  //Unifac model
  std::vector<std::vector<double> > unifacgroup_R_;
  std::vector<std::vector<double> > unifacgroup_Q_;
  std::vector<std::vector<int> >    unifacgroup_subgroup_;
  std::vector<std::vector<int> >    unifacgroup_maingroup_;
  std::vector<std::vector<int> >    unifacgroup_number_;


  //Coefficients
  std::vector<std::vector<double> > rhoLc_;
  std::vector<std::vector<double> > pVapc_;
  std::vector<std::vector<double> > hovc_;
  std::vector<std::vector<double> > etaLc_;
  std::vector<std::vector<double> > lambdaLc_;
  std::vector<std::vector<double> > cpGc_;
  std::vector<std::vector<double> > cpLc_;
  std::vector<std::vector<double> > etaGc_;
  std::vector<std::vector<double> > lambdaGc_;
  std::vector<std::vector<double> > sigmac_;
};

#include "speciesMap.hpp"

#endif
