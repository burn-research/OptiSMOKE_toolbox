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

#include "Grammar_SurfaceOnTheFlyROPA.h"

namespace OpenSMOKE
{
	void ReorderROPACoefficients(const std::vector<unsigned int>& positive_indices, const std::vector<unsigned int>& negative_indices,
		const std::vector<double>& positive_coefficients, const std::vector<double>& negative_coefficients,
		std::vector<int>& reordered_positive_indices, std::vector<double>& reordered_positive_coefficients,
		std::vector<int>& reordered_negative_indices, std::vector<double>& reordered_negative_coefficients);

	SurfaceOnTheFlyROPA::SurfaceOnTheFlyROPA(	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap,
												OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsMap) :
		
		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	
	{
		merge_reaction_rates_ = false;
		is_active_ = false;
		threshold_contributions_ = 0.03;
		compact_mode_ = true;
		threshold_formation_rate_ = 1.e-32;
	}

	void SurfaceOnTheFlyROPA::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output)
	{
		is_active_ = true;

		Grammar_SurfaceOnTheFlyROPA grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@MergeForwardAndBackwardReactions") == true)
			dictionary.ReadBool("@MergeForwardAndBackwardReactions", merge_reaction_rates_);

		if (dictionary.CheckOption("@CompactOutput") == true)
			dictionary.ReadBool("@CompactOutput", compact_mode_);

		if (dictionary.CheckOption("@Threshold") == true)
			dictionary.ReadDouble("@Threshold", threshold_contributions_);

		if (dictionary.CheckOption("@Species") == true)
		{
			dictionary.ReadOption("@Species", list_species_);
			if (list_species_[0] == "ALL" || list_species_[0] == "all")
			{
				list_species_.resize(thermodynamicsMap_.NumberOfSpecies());
				list_species_ = thermodynamicsMap_.NamesOfSpecies();
			}
		}

		Setup(path_kinetics_output);
	}
		
	void SurfaceOnTheFlyROPA::Setup(const boost::filesystem::path& path_kinetics_output)
	{
		boost::filesystem::path file_name = path_kinetics_output / "surface_reaction_names.xml";
		OpenSMOKE::ImportReactionNames(file_name, kineticsMap_.NumberOfReactions(), reaction_names_);

		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &R_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &P_, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &D_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rf_, true);
		OpenSMOKE::ChangeDimensions(kineticsMap_.NumberOfReactions(), &rb_, true);
	}

	void SurfaceOnTheFlyROPA::WriteHead(std::ofstream& fOut, const std::string reactor_type)
	{
		OpenSMOKE::OpenSMOKE_logo(fOut, "Rate of Production Analysis");

		fOut << "Frequency factors:      [kmol,m,s]" << std::endl;
		fOut << "Activation energies:    [cal/mol]" << std::endl;
		fOut << "Formation rates:        [kmol,m,s]" << std::endl;
		fOut << "Reaction rates:         [kmol,m,s]" << std::endl;
		fOut << "Compact mode:           " << compact_mode_ << std::endl;
		fOut << "Merge forward/backward: " << merge_reaction_rates_ << std::endl;
		fOut << "Threshold:              " << threshold_contributions_ * 100 << "%" << std::endl;
		fOut << std::endl;

		fOut << "Reactor: " << reactor_type << std::endl;
		fOut << "Date:    " << OpenSMOKE::GetCurrentDate() << std::endl;
		fOut << "Time:    " << OpenSMOKE::GetCurrentTime() << std::endl;
		fOut << std::endl;

		fOut << "Number of species (total):     " << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		fOut << "Number of reactions (surface): " << kineticsMap_.NumberOfReactions() << std::endl;
		fOut << std::endl;
	}

	void SurfaceOnTheFlyROPA::Analyze(	std::ofstream& fOut, const int iteration, const double t, const double T, const double P, 
										const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEVectorDouble& omega, 
										const OpenSMOKE::OpenSMOKEVectorDouble& Z, const OpenSMOKE::OpenSMOKEVectorDouble& a, const OpenSMOKE::OpenSMOKEVectorDouble& gamma)
	{
		{
			// Concentration [kmol/m3]
			const double cTot = c.SumElements();

			// Molecular weight [kg/kmol]
			double MW = 0.;
			for (unsigned int i = 1; i <= thermodynamicsMap_.number_of_gas_species(); i++)
				MW += c[i] * thermodynamicsMap_.MW(i-1);
			MW /= cTot;

			// Density [kg/m3]
			const double rho = P*MW / PhysicalConstants::R_J_kmol / T;

			// Calculates thermodynamic properties
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P);

			// Calculates kinetics
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P);

			// Reaction rates
			kineticsMap_.KineticConstants();
			kineticsMap_.ReactionRates(c.GetHandle(), Z.GetHandle(), a.GetHandle(), gamma.GetHandle());
			kineticsMap_.ProductionAndDestructionRates(P_.GetHandle(), D_.GetHandle());

			// Rate of Production Analysis
			OpenSMOKE::ROPA_Data ropa;
			if (merge_reaction_rates_ == true)
			{
				kineticsMap_.FormationRates(R_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa);
			}
			else
			{
				kineticsMap_.GetForwardReactionRates(rf_.GetHandle());
				kineticsMap_.GetBackwardReactionRates(rb_.GetHandle());
				kineticsMap_.RateOfProductionAnalysis(ropa, rf_.GetHandle(), rb_.GetHandle());
			}
			
			fOut << "***********************************************************************************************************" << std::endl;
			fOut << std::setw(14) << std::left << "Time[s]";
			fOut << std::setw(14) << std::left << "T[K]";
			fOut << std::setw(14) << std::left << "P[atm]";
			fOut << std::setw(14) << std::left << "rho[kg/m3]";
			fOut << std::setw(14) << std::left << "MW[kg/kmol]";
			fOut << std::setw(14) << std::left << "Gamma[kmol/m2]";
			fOut << std::endl;
			fOut << std::setw(14) << std::left << std::scientific << std::setprecision(2) << t;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << T;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << P / 101325.;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(4) << rho;
			fOut << std::setw(14) << std::left << std::fixed << std::setprecision(2) << MW;
			fOut << std::setw(14) << std::left << std::scientific << std::setprecision(4) << gamma[1];
			fOut << std::endl;
			fOut << "***********************************************************************************************************" << std::endl;
			fOut << std::endl;

			fOut << std::setw(30) << std::left << "GasSpecies";
			fOut << std::setw(16) << std::left << "Conc.[kmol/m3]";
			fOut << std::setw(16) << std::left << "Mole fract.";
			fOut << std::setw(16) << std::left << "Mass fract.";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 1; i <= thermodynamicsMap_.number_of_gas_species(); i++)
			{
				if (compact_mode_ == false)
				{
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[i - 1];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i] / cTot;
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << omega[i];
					fOut << std::endl;
				}
				else
				{
					if (omega[i] >= 1.e-16)
					{
						fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[i - 1];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << c[i] / cTot;
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << omega[i];
						fOut << std::endl;
					}
				}
			}
			fOut << std::endl;
			fOut << std::endl;

			fOut << std::setw(30) << std::left << "SurfaceSpecies";
			fOut << std::setw(16) << std::left << "Fraction[-]";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 1; i <= thermodynamicsMap_.number_of_site_species(); i++)
			{
				unsigned int j = thermodynamicsMap_.number_of_gas_species() + i;
				if (compact_mode_ == false)
				{
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j - 1];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << Z[i];
					fOut << std::endl;
				}
				else
				{
					if (Z[i] >= 1.e-16)
					{
						fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j - 1];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << Z[i];
						fOut << std::endl;
					}
				}
			}
			fOut << std::endl;
			fOut << std::endl;

			fOut << std::setw(30) << std::left << "BulkSpecies";
			fOut << std::setw(16) << std::left << "Activity[-]";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 1; i <= thermodynamicsMap_.number_of_bulk_species(); i++)
			{
				unsigned int j = thermodynamicsMap_.number_of_gas_species() + thermodynamicsMap_.number_of_site_species() + i;
				if (compact_mode_ == false)
				{
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j - 1];
					fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << a[i];
					fOut << std::endl;
				}
				else
				{
					if (a[i] >= 1.e-16)
					{
						fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j - 1];
						fOut << std::setw(16) << std::left << std::scientific << std::setprecision(5) << a[i];
						fOut << std::endl;
					}
				}
			}
			fOut << std::endl;
			fOut << std::endl;

			// Loop over all the species
			for (unsigned int i = 0; i < list_species_.size(); i++)
			{
				const unsigned int index_of_species = thermodynamicsMap_.IndexOfSpecies(list_species_[i]) - 1;

				// Reorder species
				std::vector<int> reordered_positive_indices;
				std::vector<double> reorder_positive_coefficients;
				std::vector<int> reordered_negative_indices;
				std::vector<double> reorder_negative_coefficients;
				ReorderROPACoefficients(ropa.production_reaction_indices[index_of_species],
					ropa.destruction_reaction_indices[index_of_species],
					ropa.production_coefficients[index_of_species], ropa.destruction_coefficients[index_of_species],
					reordered_positive_indices, reorder_positive_coefficients,
					reordered_negative_indices, reorder_negative_coefficients);

				if (P_[index_of_species + 1] >= threshold_formation_rate_ || D_[index_of_species + 1] >= threshold_formation_rate_)
				{
					// Write 
					fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[index_of_species];
					fOut << "Formation: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << P_[index_of_species + 1];
					fOut << "Consumption: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << -D_[index_of_species + 1];
					fOut << "Net: " << std::setw(12) << std::left << std::setprecision(2) << std::scientific << P_[index_of_species + 1] - D_[index_of_species + 1];
					fOut << std::endl;
					fOut << "-----------------------------------------------------------------------------------------------------------" << std::endl;

					// Production contribution
					const double sum_positive_coefficients = std::accumulate(reorder_positive_coefficients.begin(), reorder_positive_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_positive_indices.size(); i++)
					{
						const double percentage = reorder_positive_coefficients[i] / sum_positive_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_positive_indices[i];
							fOut << std::setw(6) << std::left << j;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << reorder_positive_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << percentage*100. << "%";
//							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
//							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
//							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j)* PhysicalConstants::R_cal_mol;
							fOut << std::left << "   " << reaction_names_[j] << std::endl;
						}
					}

					// Consumption contributions
					const double sum_negative_coefficients = std::accumulate(reorder_negative_coefficients.begin(), reorder_negative_coefficients.end(), 0.);
					for (unsigned int i = 0; i < reordered_negative_indices.size(); i++)
					{
						const double percentage = reorder_negative_coefficients[i] / sum_negative_coefficients;
						if (percentage > threshold_contributions_)
						{
							unsigned int j = reordered_negative_indices[i];
							fOut << std::setw(6) << std::left << j;
							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << -reorder_negative_coefficients[i];
							fOut << std::setw(10) << std::right << std::setprecision(1) << std::fixed << -percentage*100. << "%";
//							fOut << std::setw(10) << std::right << std::setprecision(2) << std::scientific << kineticsMap_.A(j);
//							fOut << std::setw(8) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.Beta(j);
//							fOut << std::setw(12) << std::right << std::setprecision(2) << std::fixed << kineticsMap_.E_over_R(j)* PhysicalConstants::R_cal_mol;
							fOut << std::left << "   " << reaction_names_[j] << std::endl;
						}
					}

					fOut << std::endl;
				}
			}
		}
	}
	
}

