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

mixtureG::mixtureG(	const std::vector<string>& species, speciesMap& map) :
species_(species),
map_(map)
{
	NS_ = species_.size();
	CheckBasicSpecies();

	fugacity_rule_ = 1;		// 1. Ideal gas. 2: Peng-Robinson equation of state

	MemoryAllocation();
	IndependentProperties();
	SetupProperties();
}

mixtureG::mixtureG(const std::vector<string>& species, speciesMap& map, OpenSMOKE::OpenSMOKE_Dictionary& dictionary) :
species_(species),
map_(map)
{
	NS_ = species_.size();
	CheckBasicSpecies();

	fugacity_rule_ = 1;		// 1. Ideal gas. 2: Peng-Robinson equation of state

	SetupFromDictionary(dictionary);

	MemoryAllocation();
	IndependentProperties();
	SetupProperties();
}

void mixtureG::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
{
	Grammar_GasMixture species_grammar;
	dictionary.SetGrammar(species_grammar);

	//Non compulsory options 
	if (dictionary.CheckOption("@Fugacity"))
	{
		string value;
		dictionary.ReadString("@Fugacity", value);
		if (value == "Raoult")
		fugacity_rule_ = 1;
		else if (value == "PengRobinson")
		fugacity_rule_ = 2;
		else
		OpenSMOKE::FatalErrorMessage("@Fugacity options are: Raoult || PengRobinson");
	}
}

mixtureG::~mixtureG()
{
	delete fugacityM_;
}

void mixtureG::CheckBasicSpecies() 
{
}

void mixtureG::MemoryAllocation()
{
	mapindex_.resize(NS_);
	Tc_.resize(NS_);
	Pc_.resize(NS_);
	MW_.resize(NS_);
	Rhoc_.resize(NS_);
	Tnbp_.resize(NS_);
	dipolMoment_.resize(NS_);
	omega_.resize(NS_);
	specific_volume_stp_.resize(NS_);
}

void mixtureG::IndependentProperties()
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
	}
}

void mixtureG::SetupProperties()
{
	switch (fugacity_rule_)
	{
		case 1:
			fugacityM_ = new FugacityGIdealGas_mix(NS_);
			break;

		case 2:
			fugacityM_ = new FugacityGPengRobinson_mix(NS_, Tc_, Pc_, omega_, MW_);
			break;

		default:
			cout << "Mixture rhoL model unknown!" << endl;
			exit(-1);
	}
}

const std::vector<double> mixtureG::fugacity(const double T, const double P, const std::vector<double>& x) const
{
	return fugacityM_->fugacity(T, P, x);
}
