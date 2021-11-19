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

mixtureL::mixtureL(const std::vector<std::string>& species, speciesMap& map) :
species_(species),
map_(map)
{
	NS_ = species_.size();

	// Mixing rules
	density_rule_ = 1;		// Density via mixture PR
	capacity_rule_ = 1;		// Cp via ideal mixing
	viscosity_rule_ = 1;	// Viscosity via Grunberg and Nisan (1949) rule - no interaction parameters
	conductivity_rule_ = 1; // Conductivity via Vredeveld method
	dinf_rule_ = 1;			// Siddiqi-Lucas correlation
	activity_rule_ = 1;		// Ideal mixture
	fugacity_rule_ = 2;		// 1: Ideal gas. 2: Indirect method. 3: Direct Method

	diffusion_correction_ = 1.;

	MemoryAllocation();
	IndependentProperties();
	SetupProperties();
}

mixtureL::mixtureL(const std::vector<std::string>& species, speciesMap& map, OpenSMOKE::OpenSMOKE_Dictionary& dictionary) :
species_(species),
map_(map)
{
	NS_ = species_.size();

	// Mixing rules
	density_rule_ = 1;		// Density via mixture PR
	capacity_rule_ = 1;		// Cp via ideal mixing
	viscosity_rule_ = 1;	// Viscosity via Grunberg and Nisan (1949) rule - no interaction parameters
	conductivity_rule_ = 1; // Conductivity via Vredeveld method
	dinf_rule_ = 1;			// Siddiqi-Lucas correlation
	activity_rule_ = 1;		// UNIFAC model
	fugacity_rule_ = 2;		// 1: Ideal gas. 2: Indirect method. 3: Direct Method

	diffusion_correction_ = 1.;
	conductivity_correction_ = 1.;

	SetupFromDictionary(dictionary);

	MemoryAllocation();
	IndependentProperties();
	SetupProperties();
}

mixtureL::~mixtureL()
{
	delete rhoLM_;
	delete cpLM_;
	delete etaLM_;
	delete lambdaLM_;
	delete DinfM_;

	delete smM_;
	delete fugacityM_;
}

void mixtureL::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
{
	Grammar_LiquidMixture species_grammar;
	dictionary.SetGrammar(species_grammar);

	//Non compulsory options
	if (dictionary.CheckOption("@Density"))
	{
		std::string value;
		dictionary.ReadString("@Density", value);

		if (value == "Peng-Robinson")
			density_rule_ = 1;
		else if (value == "AdditiveVolumes")
			density_rule_ = 2;
		else
			OpenSMOKE::FatalErrorMessage("@Density options are: Peng-Robinson | AdditiveVolumes");
	}

	if (dictionary.CheckOption("@HeatCapacity"))
	{
		std::string value;
		dictionary.ReadString("@HeatCapacity", value);
		
		if (value == "Ideal")
			capacity_rule_ = 1;
		else
			OpenSMOKE::FatalErrorMessage("@HeatCapacity options are: Ideal");
	}

	if (dictionary.CheckOption("@Viscosity"))
	{
		std::string value;
		dictionary.ReadString("@Viscosity", value);
		
		if (value == "GrunbergNissanIdeal")
			viscosity_rule_ = 1;
		else
			OpenSMOKE::FatalErrorMessage("@Viscosity options are: GrunbergNissanIdeal");
	}

	if (dictionary.CheckOption("@Conductivity"))
	{
		std::string value;
		dictionary.ReadString("@Conductivity", value);
		
		if (value == "Vredeveld")
			conductivity_rule_ = 1;
		else
			OpenSMOKE::FatalErrorMessage("@Conductivity options are: Vredeveld");
	}

	if (dictionary.CheckOption("@Dinf"))
	{
		std::string value;
		dictionary.ReadString("@Dinf", value);
		
		if (value == "Siddiqi-Lucas")
			dinf_rule_ = 1;
		else
			OpenSMOKE::FatalErrorMessage("@Dinf options are: Siddiqi-Lucas");
	}

	if (dictionary.CheckOption("@Activity"))
	{
		std::string value;
		dictionary.ReadString("@Activity", value);

		if (value == "Ideal")
			activity_rule_ = 1;
		else if (value == "Unifac")
			activity_rule_ = 2;
		else
			OpenSMOKE::FatalErrorMessage("@Activity options are: Ideal || Unifac");
	}

	if (dictionary.CheckOption("@Fugacity"))
	{
		std::string value;
		dictionary.ReadString("@Fugacity", value);

		if (value == "Raoult")
			fugacity_rule_ = 1;
		else if (value == "PengRobinsonIndirect")
			fugacity_rule_ = 2;
		else if (value == "PengRobinsonDirect")
			fugacity_rule_ = 3;
		else
			OpenSMOKE::FatalErrorMessage("@Fugacity options are: Raoult || PengRobinsonIndirect || PengRobinsonDirect");
	}

	if (dictionary.CheckOption("@DiffusionCorrection"))
		dictionary.ReadDouble("@DiffusionCorrection", diffusion_correction_);

	if (dictionary.CheckOption("@ConductivityCorrection"))
		dictionary.ReadDouble("@ConductivityCorrection", conductivity_correction_);
}

void mixtureL::MemoryAllocation()
{
	//Physical properties
	mapindex_.resize(NS_);
	Tc_.resize(NS_);
	Pc_.resize(NS_);
	MW_.resize(NS_);
	Rhoc_.resize(NS_);
	Tnbp_.resize(NS_);
	dipolMoment_.resize(NS_);
	omega_.resize(NS_);
	specific_volume_stp_.resize(NS_);
	dM_.resize(NS_);
	vpM_.resize(NS_);

	//UNIFAC properties
	unifacgroup_R_.resize(NS_);
	unifacgroup_Q_.resize(NS_);
	unifacgroup_subgroup_.resize(NS_);
	unifacgroup_maingroup_.resize(NS_);
	unifacgroup_number_.resize(NS_);
}

void mixtureL::IndependentProperties()
{
	for (unsigned int i = 0; i < NS_; i++)
	{
		mapindex_[i] = OpenSMOKE::Index(species_[i], map_.name());

		//Physical properties
		Tc_[i] = map_.Tc()[mapindex_[i]];
		Pc_[i] = map_.Pc()[mapindex_[i]];
		MW_[i] = map_.MW()[mapindex_[i]];
		Rhoc_[i] = map_.Rhoc()[mapindex_[i]];
		Tnbp_[i] = map_.Tnbp()[mapindex_[i]];
		dipolMoment_[i] = map_.dipolMoment()[mapindex_[i]];
		omega_[i] = map_.omega()[mapindex_[i]];
		specific_volume_stp_[i] = map_.specific_volume_stp()[mapindex_[i]];

		dM_[i] = map_.dM()[mapindex_[i]];
		vpM_[i] = map_.vpM()[mapindex_[i]];

		//Unifac properties
		unifacgroup_R_[i] = map_.unifacgroup_R()[mapindex_[i]];
		unifacgroup_Q_[i] = map_.unifacgroup_Q()[mapindex_[i]];
		unifacgroup_subgroup_[i] = map_.unifacgroup_subgroup()[mapindex_[i]];
		unifacgroup_maingroup_[i] = map_.unifacgroup_maingroup()[mapindex_[i]];
		unifacgroup_number_[i] = map_.unifacgroup_number()[mapindex_[i]];
	}
}

void mixtureL::SetupProperties()
{
	{
		switch (density_rule_)
		{
			case 1:
				rhoLM_ = new rhoLPengRobinson_mix(NS_, Tc_, Pc_, omega_, MW_);
				break;
			case 2:
				rhoLM_ = new rhoLAdditiveVolumes_mix(NS_, dM_, MW_);	
				break;		
			default:
				cout << "Mixture rhoL model unknown!" << endl;
				exit(-1);
		}
			
		switch (capacity_rule_)
		{
			case 1:
				cpLM_ = new cpLIdealmix(NS_, MW_);
				break;
			
			default:
				cout << "Mixture cp model unknown!" << endl;
				exit(-1);
		}
			
		switch (viscosity_rule_)
		{
			case 1:
				etaLM_ = new etaL_GNI_mix(NS_);
				break;
			
			default:
				cout << "Mixture viscosity model unknown!" << endl;
				exit(-1);
		}
			
		switch (conductivity_rule_)
		{
			case 1:
				lambdaLM_ = new lambdaLV_mix(NS_, MW_);
				lambdaLM_->SetConductivityCorrection(conductivity_correction_);
				break;

			default:
				cout << "Mixture conductivity model unknown!" << endl;
				exit(-1);
		}
		
		switch (dinf_rule_)
		{
			case 1:
				DinfM_ = new DinfSL(NS_, MW_, Tnbp_, dM_);
				DinfM_->SetDiffusionCorrection(diffusion_correction_);
				break;
			
			default:
				cout << "Infinite diffusivity model unknown!" << endl;
				exit(-1);
		}
		
		switch (activity_rule_)
		{
			case 1:
				activityM_ = new IdealActivityModel(NS_);
				break;
			
			case 2:
				activityM_ = new UnifacModel(NS_, unifacgroup_R_, unifacgroup_Q_,
				unifacgroup_subgroup_, unifacgroup_maingroup_, unifacgroup_number_,
				map_.interaction_table());
				break;
			
			default:
				cout << "Activity model type unknown!" << endl;
				exit(-1);
		}
			
		switch (fugacity_rule_)
		{
			case 1:
				fugacityM_ = new FugacityLIdealGas_mix(NS_, vpM_);
				break;
			
			case 2:
				fugacityM_ = new FugacityLPengRobinsonIndirect_mix(NS_, Tc_, Pc_, omega_, MW_, dM_, vpM_, *activityM_);
				break;
			
			case 3:
				fugacityM_ = new FugacityLPengRobinsonDirect_mix(NS_, Tc_, Pc_, omega_, MW_);
				break;
			
			default:
				cout << "Fugacity model type unknown!" << endl;
				exit(-1);
		}
	}

	smM_ = new stefanmaxwell(*DinfM_, *activityM_);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Mixture Liquid  density in [kg/m^3]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double mixtureL::rhoL_mix(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);
	return rhoLM_->rho(T, P, x);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Mixture derivative of Liquid  density
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double mixtureL::d_rhoL_over_dT_mix(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);
	return rhoLM_->d_rhoL_over_dT(T, P, x);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Mixture Specific heat in [W/kgK]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double mixtureL::cpL_mix(const double T, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	std::vector<double> cp_species(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		cp_species[i] = map_.cpL(species_[i], T);

	return cpLM_->cp(cp_species, x);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Mixture Liquid viscosity in [Pa s]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double mixtureL::etaL_mix(const double T, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	std::vector<double> etaL_species(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		etaL_species[i] = map_.etaL(species_[i], T);

	return etaLM_->etaL(etaL_species, x);
}

double mixtureL::lambdaL_mix(const double T, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	std::vector<double> lambdaL_species(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		lambdaL_species[i] = map_.lambdaL(species_[i], T);

	return lambdaLM_->lambdaL(lambdaL_species, x);
}

const Eigen::MatrixXd mixtureL::Dinf(const double T, const double P) const
{
  std::vector<double> etaL_species(NS_);
  for (unsigned int i = 0; i < NS_; i++)
    etaL_species[i] = map_.etaL(species_[i], T);

  DinfM_->UpdateInfiniteDiffusivities(etaL_species, T, P);

  return DinfM_->Dinf();
}

const Eigen::MatrixXd mixtureL::D_sm(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	Dinf(T, P);
	smM_->UpdateBinaryDiffusivities(x);
  
	return smM_->D_sm();
}

const Eigen::MatrixXd mixtureL::B(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	D_sm(T, P, x);
	smM_->UpdateB(T, P, x);
	
	return smM_->B();
}

const Eigen::MatrixXd mixtureL::thermodynamicmatrix_sm(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	smM_->UpdateThermodynamicMatrix(T, P, x);

	return smM_->tm();
}

const std::vector<double> mixtureL::gamma(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	activityM_->UpdateGamma(T, P, x);
  
	return activityM_->gamma();
}

void mixtureL::SetStefanMaxwellLastSpecies(const std::string species)
{
	smM_->SetLastIndex(OpenSMOKE::Index(species, species_));
}

const Eigen::VectorXd mixtureL::molar_diffusion_fluxes(	const double T, const double P,
														const std::vector<double>& x, const double ctot, 
														const std::vector<double>& dx) const
{
	Dinf(T, P);

	return smM_->molar_diffusion_fluxes(T, P, x, ctot, dx);
}

const std::vector<double> mixtureL::fugacity(const double T, const double P, const std::vector<double>& x_original) const
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	return fugacityM_->fugacity(T, P, x);
}
