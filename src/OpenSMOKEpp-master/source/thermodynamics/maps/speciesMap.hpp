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

speciesMap::speciesMap(const fs::path& config_folder) :
config_folder_(config_folder)
{
  fs::path input_path = config_folder_ / "speciesSettings.cfg";

  ReadSpeciesNames(input_path);
  MemoryAllocation();
  ReadSpeciesProperties(input_path);
  SetupProperties();
}

speciesMap::~speciesMap()
{
  for (int i = 0; i < NS_; i++)
    {
      delete dM_[i];
      delete vpM_[i];
      delete hovM_[i];
      delete etaLM_[i];
      delete lambdaLM_[i];
      delete cpGM_[i];
      delete cpLM_[i];
      delete etaGM_[i];
      delete lambdaGM_[i];
      delete sigmaM_[i];
    }

  delete unifactable_;
}

void speciesMap::ReadSpeciesNames(const fs::path& map_path)
{
  lc::Config doc;

  try
    {
      doc.readFile(map_path.string().c_str());
      const lc::Setting &species = doc.lookup("species");
      NS_ = species.getLength();
      name_.resize(NS_);
      for (int i = 0; i < NS_; i++)
        {
          species[i].lookupValue("speciesName", name_[i]);
        }
    }
  catch (const lc::SettingNotFoundException &nfex)
    {
      std::cerr << "Setting error at " << nfex.getPath() << std::endl;
    }
  catch (const lc::ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    }
  catch (const lc::FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file:" << fioex.what() << std::endl;
    }
}

void speciesMap::MemoryAllocation()
{
  //Independent properties
  MW_.resize(NS_);
  Tc_.resize(NS_);
  Pc_.resize(NS_);
  Rhoc_.resize(NS_);
  Tnbp_.resize(NS_);
  dipolMoment_.resize(NS_);
  omega_.resize(NS_);
  specific_volume_stp_.resize(NS_);

  //Dependent properties: correlation indices
  pVap_index_.resize(NS_);
  dHv_index_.resize(NS_);
  etaL_index_.resize(NS_);
  lambdaL_index_.resize(NS_);
  cpG_index_.resize(NS_);
  cpL_index_.resize(NS_);
  lambdaG_index_.resize(NS_);
  etaG_index_.resize(NS_);
  rhoL_index_.resize(NS_);
  sigma_index_.resize(NS_);


  //Model objects
  dM_.resize(NS_);
  vpM_.resize(NS_);
  hovM_.resize(NS_);
  etaLM_.resize(NS_);
  lambdaLM_.resize(NS_);
  cpGM_.resize(NS_);
  cpLM_.resize(NS_);
  etaGM_.resize(NS_);
  lambdaGM_.resize(NS_);
  sigmaM_.resize(NS_);

  //Coefficients
  rhoLc_.resize(NS_);
  pVapc_.resize(NS_);
  hovc_.resize(NS_);
  etaLc_.resize(NS_);
  lambdaLc_.resize(NS_);
  cpGc_.resize(NS_);
  cpLc_.resize(NS_);
  etaGc_.resize(NS_);
  lambdaGc_.resize(NS_);
  sigmac_.resize(NS_);

  //Unifac model
  unifacgroup_R_.resize(NS_);
  unifacgroup_Q_.resize(NS_);
  unifacgroup_subgroup_.resize(NS_);
  unifacgroup_maingroup_.resize(NS_);
  unifacgroup_number_.resize(NS_);

  for (int i = 0; i < NS_; i++)
    {
      rhoLc_[i].resize(4);
      pVapc_[i].resize(7);
      hovc_[i].resize(5);
      etaLc_[i].resize(5);
      lambdaLc_[i].resize(5);
      cpGc_[i].resize(5);
      cpLc_[i].resize(5);
      etaGc_[i].resize(5);
      lambdaGc_[i].resize(5);
      sigmac_[i].resize(2);
    }
}

void speciesMap::ReadSpeciesProperties(const fs::path& map_path)
{
  lc::Config doc;
  doc.readFile(map_path.string().c_str());

  try
    {
      const lc::Setting &molarMassCfg = doc.lookup("molarMass");
      const lc::Setting &TcritCfg = doc.lookup("Tcrit");
      const lc::Setting &pCritCfg = doc.lookup("pCrit");
      const lc::Setting &rhoCritCfg = doc.lookup("rhoCrit");
      const lc::Setting &TnbpCfg = doc.lookup("Tnbp");
      const lc::Setting &acentricCfg = doc.lookup("acentric");
      const lc::Setting &dipolMomentCfg = doc.lookup("dipolMoment");

      for (int i = 0; i < NS_; i++)
        {
          molarMassCfg.lookupValue(name_[i], MW_[i]);
          TcritCfg.lookupValue(name_[i], Tc_[i]);
          pCritCfg.lookupValue(name_[i], Pc_[i]);
          rhoCritCfg.lookupValue(name_[i], Rhoc_[i]);
          TnbpCfg.lookupValue(name_[i], Tnbp_[i]);
          dipolMomentCfg.lookupValue(name_[i], dipolMoment_[i]);
          acentricCfg.lookupValue(name_[i], omega_[i]);
        }

      lc::Setting &species = doc.lookup("species");
      for (int i = 0; i < NS_; i++)
        {
          species[i].lookupValue("pVapEq", pVap_index_[i]);
          species[i].lookupValue("deltaHvEq", dHv_index_[i]);
          species[i].lookupValue("etaLEq", etaL_index_[i]);
          species[i].lookupValue("lambdaLEq", lambdaL_index_[i]);
          species[i].lookupValue("cpGEq", cpG_index_[i]);
          species[i].lookupValue("cpLEq", cpL_index_[i]);
          species[i].lookupValue("lambdaGEq", lambdaG_index_[i]);
          species[i].lookupValue("etaGEq", etaG_index_[i]);
          species[i].lookupValue("rhoLEq", rhoL_index_[i]);
          species[i].lookupValue("sigmaEq", sigma_index_[i]);
        }

      {
        // Liquid density
        lc::Setting &rhoLCfg = doc.lookup("rhoL");
        for (int i = 0; i < NS_; i++)
          if (rhoL_index_[i] > 0)
            for (int j = 0; j < 4; j++)
              rhoLc_[i][j] = rhoLCfg[name_[i].c_str()][j];
      }


      {
        //Vapor pressure
        lc::Setting &vpMCfg = doc.lookup("pVap");
        for (int i = 0; i < NS_; i++)
          if (pVap_index_[i] > 0)
            for (int j = 0; j < 7; j++)
              pVapc_[i][j] = vpMCfg[name_[i].c_str()][j];
      }

      {
        //Heat of vaporization
        lc::Setting &hovMCfg = doc.lookup("deltaHv");
        for (int i = 0; i < NS_; i++)
          if (dHv_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              hovc_[i][j] = hovMCfg[name_[i].c_str()][j];
      }

      {
        //Liquid dynamic viscosity
        lc::Setting &etaLMCfg = doc.lookup("etaL");
        for (int i = 0; i < NS_; i++)
          if (etaL_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              etaLc_[i][j] = etaLMCfg[name_[i].c_str()][j];
      }

      {
        //Liquid conductivity
        lc::Setting &lambdaLMCfg = doc.lookup("lambdaL");
        for (int i = 0; i < NS_; i++)
          if (lambdaL_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              lambdaLc_[i][j] = lambdaLMCfg[name_[i].c_str()][j];
      }

      {
        //Gas heat capacity
        lc::Setting &cpGMCfg = doc.lookup("cpG");
        for (int i = 0; i < NS_; i++)
          if (cpG_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              cpGc_[i][j] = cpGMCfg[name_[i].c_str()][j];
      }

      {
        //Liquid heat capacity
        lc::Setting &cpLMCfg = doc.lookup("cpL");
        for (int i = 0; i < NS_; i++)
          if (cpL_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              cpLc_[i][j] = cpLMCfg[name_[i].c_str()][j];
      }

      {
        //Gas dynamic viscosity
        lc::Setting &etaGMCfg = doc.lookup("etaG");
        for (int i = 0; i < NS_; i++)
          if (etaG_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              etaGc_[i][j] = etaGMCfg[name_[i].c_str()][j];
      }

      {
        //Gas conductivity
        lc::Setting &lambdaGMCfg = doc.lookup("lambdaG");
        for (int i = 0; i < NS_; i++)
          if (lambdaG_index_[i] > 0)
            for (int j = 0; j < 5; j++)
              lambdaGc_[i][j] = lambdaGMCfg[name_[i].c_str()][j];
      }

      {
        //Surface tension
        lc::Setting &sigmaMCfg = doc.lookup("sigma");
        for (int i = 0; i < NS_; i++)
          for (int j = 0; j < 2; j++)
            sigmac_[i][j] = sigmaMCfg[name_[i].c_str()][j];
      }

      {
        //Unifac activity model
        unifactable_ = new UnifacTable();
        const fs::path table_path = config_folder_ / "unifacTableCoeffs.cfg";
        unifactable_->SetupTable(table_path);
        const fs::path interaction_coefficients_path = config_folder_ / "interactionCoefficients.cfg";
        unifactable_->SetupInteractionCoefficients(interaction_coefficients_path);
        lc::Setting &unifacCfg = doc.lookup("UNIFAC");
        for (int i = 0; i < NS_; i++)
          {
            lc::Setting &unifacDict = unifacCfg[name_[i].c_str()];
            unsigned int n_groups = unifacDict.getLength();
            unifacgroup_R_[i].resize(n_groups);
            unifacgroup_Q_[i].resize(n_groups);
            unifacgroup_subgroup_[i].resize(n_groups);
            unifacgroup_maingroup_[i].resize(n_groups);
            unifacgroup_number_[i].resize(n_groups);
            for (unsigned int j = 0; j < n_groups; j++)
              {
                string nameStruct;
                unsigned int numberStruct;
                unifacDict[j].lookupValue("name", nameStruct);
                unifacDict[j].lookupValue("number", numberStruct);
                unsigned int group_index = OpenSMOKE::Index(nameStruct, unifactable_->groupname());
                unifacgroup_R_[i][j] = unifactable_->R()[group_index];
                unifacgroup_Q_[i][j] = unifactable_->Q()[group_index];
                unifacgroup_subgroup_[i][j] = unifactable_->subgroup()[group_index];
                unifacgroup_maingroup_[i][j] = unifactable_->maingroup()[group_index];
                unifacgroup_number_[i][j] = numberStruct;
              }
          }
      }


    }
  catch (const lc::SettingNotFoundException &nfex)
    {
      std::cout << "Setting error at " << nfex.getPath() << std::endl;
      exit(-1);
    }
  catch (const lc::ParseException &pex)
    {
      std::cout << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
      exit(-1);
    }
  catch (const lc::FileIOException &fioex)
    {
      std::cout << "I/O error while reading file:" << fioex.what() << std::endl;
      exit(-1);
    }

}

void speciesMap::SetupProperties()
{
  for (int i = 0; i < NS_; i++)
    {
      switch (rhoL_index_[i])
        {
          case 0:
            dM_[i] = NULL;
            break;
          case 1:
            dM_[i] = new RhoVDI(&rhoLc_[i][0]);
            break;
          case 2:
            dM_[i] = new RhoYaws(&rhoLc_[i][0]);
            break;
          case 3:
            dM_[i] = new RhoLPengRobinson(Tc_[i], Pc_[i], omega_[i], MW_[i]);
            break;
          default:
            cout << "Equation type for rhoL of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (pVap_index_[i])
        {
          case 0:
            vpM_[i] = NULL;
            break;
          case 1:
            vpM_[i] = new PVapVDI(&pVapc_[i][0], Tc_[i], Pc_[i]);
            break;
          case 2:
            vpM_[i] = new PVapYaws(&pVapc_[i][0]);
            break;
          case 3:
            vpM_[i] = new PVapPrausnitz(&pVapc_[i][0]);
            break;
          default:
            cout << "Equation type for pVap of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (dHv_index_[i])
        {
          case 0:
            hovM_[i] = NULL;
            break;
          case 1:
            hovM_[i] = new DeltaHvVDI(&hovc_[i][0], Tc_[i]);
            break;
          case 2:
            hovM_[i] = new DeltaHvYaws(&hovc_[i][0], Tc_[i], MW_[i]);
            break;
          default:
            cout << "Equation type for deltaHv of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (etaL_index_[i])
        {
          case 0:
            etaLM_[i] = NULL;
            break;
          case 1:
            etaLM_[i] = new EtaLVDI(&etaLc_[i][0]);
            break;
          case 2:
            etaLM_[i] = new EtaLYaws(&etaLc_[i][0]);
            break;
          case 3:
            etaLM_[i] = new EtaLPrausnitz(&etaLc_[i][0], Tc_[i]);
            break;
          default:
            cout << "Equation type for etaL of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (lambdaL_index_[i])
        {
          case 0:
            lambdaLM_[i] = NULL;
            break;
          case 1:
            lambdaLM_[i] = new LambdaLVDI(&lambdaLc_[i][0], Tc_[i]);
            break;
          case 2:
            lambdaLM_[i] = new LambdaLYaws(&lambdaLc_[i][0], Tc_[i]);
            break;
          default:
            cout << "Equation type for lambdaL of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (cpG_index_[i])
        {
          case 0:
            cpGM_[i] = NULL;
            break;
          case 1:
            cpGM_[i] = new CpGVDI(&cpGc_[i][0]);
            break;
          case 2:
            cpGM_[i] = new CpGYaws(&cpGc_[i][0], MW_[i]);
            break;
          default:
            cout << "Equation type for cpG of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (cpL_index_[i])
        {
          case 0:
            cpLM_[i] = NULL;
            break;
          case 1:
            cpLM_[i] = new CpLVDI(&cpLc_[i][0]);
            break;
          case 2:
            cpLM_[i] = new CpLYaws(&cpLc_[i][0], MW_[i]);
            break;
          default:
            cout << "Equation type for cpL of " << name_[i] << " unkonwn!" << endl;
            exit(-1);
        }
      switch (etaG_index_[i])
        {
          case 0:
            etaGM_[i] = NULL;
            break;
          case 1:
            etaGM_[i] = new EtaGVDI(&etaGc_[i][0]);
            break;
          case 2:
            etaGM_[i] = new EtaGYaws(&etaGc_[i][0]);
            break;
          default:
            cout << "Equation type for etaG of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }
      switch (lambdaG_index_[i])
        {
          case 0:
            lambdaGM_[i] = NULL;
            break;
          case 1:
            lambdaGM_[i] = new LambdaGVDI(&lambdaGc_[i][0]);
            break;
          case 2:
            lambdaGM_[i] = new LambdaGYaws(&lambdaGc_[i][0]);
            break;
          case 3:
            lambdaGM_[i] = new LambdaSutherland(cpGM_[i], etaGM_[i], MW_[i]);
            break;
          default:
            cout << "Equation type for lambdaG of " << name_[i] << " unknown!" << endl;
            exit(-1);
        }

      sigmaM_[i] = new SurfaceTensionModel(&sigmac_[i][0], Tc_[i]);

      //Calculating specific volume at STP
      if (rhoL_index_[i] > 0)
        specific_volume_stp_[i] = 1. / rhoL(name_[i], 273., 101325.);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Liquid  density in [kg/m^3]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::rhoL(string species, double T, double P) const
{
  if (rhoL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Liquid density not available for " << species << "." << endl;
      exit(-1);
    }

  return dM_[OpenSMOKE::Index(species, name_)]->rho(T, P);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Derivative of Liquid Gas density in [kg/m^3]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::d_rhoL_over_dT(string species, double T, double P) const
{
  if (rhoL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Derivative of density not available for " << species << "." << endl;
      exit(-1);
    }

  return dM_[OpenSMOKE::Index(species, name_)]->d_rhoL_over_dT(T, P);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Vapor pressure in [Pa]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::pVap(string species, double T, double P) const
{
  if (pVap_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Vapor pressure not available for " << species << "." << endl;
      exit(-1);
    }

  int species_index = OpenSMOKE::Index(species, name_);
  if (T > 0.999 * Tc_[species_index]) return P;
  else return vpM_[species_index]->pVap(T);
}

double speciesMap::dpVapdT(string species, double T, double P) const
{
  if (pVap_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Derivative of vapor pressure not available for " << species << "." << endl;
      exit(-1);
    }

  int species_index = OpenSMOKE::Index(species, name_);
  return vpM_[species_index]->dpVapdT(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Heat of vaporization in in [J/kg]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::deltaHv(string species, double T) const
{
  if (dHv_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Derivative of vapor pressure not available for " << species << "." << endl;
      exit(-1);
    }

  return hovM_[OpenSMOKE::Index(species, name_)]->deltaHv(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Liquid dynamic viscosity in [W/mK]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::etaL(string species, double T) const
{
  if (etaL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Liquid viscosity not available for " << species << "." << endl;
      exit(-1);
    }

  return etaLM_[OpenSMOKE::Index(species, name_)]->eta(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Liquid thermal conductivity in [W/mK]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::lambdaL(string species, double T) const
{
  if (lambdaL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Liquid viscosity not available for " << species << "." << endl;
      exit(-1);
    }

  return lambdaLM_[OpenSMOKE::Index(species, name_)]->lambda(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Liquid heat capacity at constant pressure in [J/kg K]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::cpL(string species, double T) const
{
  if (cpL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Liquid specific heat not available for " << species << "." << endl;
      exit(-1);
    }

  return cpLM_[OpenSMOKE::Index(species, name_)]->cp(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Gas heat capacity at constant pressure in [J/kg K]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::cpG(string species, double T) const
{
  if (cpG_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Gas specific heat not available for " << species << "." << endl;
      exit(-1);
    }

  return cpGM_[OpenSMOKE::Index(species, name_)]->cp(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Gas dynamic viscosity in [Pa s]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::etaG(string species, double T) const
{
  if (etaG_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Gas viscosity not available for " << species << "." << endl;
      exit(-1);
    }

  return etaGM_[OpenSMOKE::Index(species, name_)]->eta(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Gas thermal conductivity in [W/m K]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::lambdaG(string species, double T) const
{
  if (lambdaG_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Gas conductivity not available for " << species << "." << endl;
      exit(-1);
    }

  return lambdaGM_[OpenSMOKE::Index(species, name_)]->lambda(T);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Surface tension in [N/m]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double speciesMap::sigma(string species, double T) const
{
  if (lambdaG_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Gas conductivity not available for " << species << "." << endl;
      exit(-1);
    }

  return sigmaM_[OpenSMOKE::Index(species, name_)]->sigma(T);
}

double speciesMap::specific_volume_stp(string species) const
{
  if (rhoL_index_[OpenSMOKE::Index(species, name_)] == 0)
    {
      cout << "Fatal Error: Specific volume at STP not available for  " << species << "." << endl;
      exit(-1);
    }

  return specific_volume_stp_[OpenSMOKE::Index(species, name_)];
}
