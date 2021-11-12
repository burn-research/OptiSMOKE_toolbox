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
|	License                                                           |
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


void ReadFromBackupFile(const boost::filesystem::path path_file, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
	std::vector<double>& x, std::vector<double>& T, std::vector<double>& P, std::vector<double>& m, std::vector< std::vector<double> >& omega)
{
	std::vector < std::vector<double> > omega_backup;
	std::vector<std::string> names_species_backup;

	ReadFromBackupFile(path_file, x, T, P, m, omega_backup, names_species_backup);

	std::vector<int> index_previous(thermodynamicsMap.NumberOfSpecies());
	for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
		index_previous[i] = -1;

	unsigned int count = 0;
	for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
	{
		bool iFound = false;
		for (unsigned int j = 0; j < names_species_backup.size(); j++)
		if (thermodynamicsMap.NamesOfSpecies()[i] == names_species_backup[j])
		{
			index_previous[i] = j;
			count++;
			break;
		}
	}

	std::cout << " * Number of species in the previous kinetic mechanism: " << names_species_backup.size() << std::endl;
	std::cout << " * Number of species in the current kinetic mechanism:  " << thermodynamicsMap.NumberOfSpecies() << std::endl;
	std::cout << " * Number of species which are available from backup:   " << count << std::endl;

	omega.resize(thermodynamicsMap.NumberOfSpecies());

	unsigned np = x.size();
	for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
	{
		omega[i].resize(np);
		if (index_previous[i] != -1)
		{
			for (unsigned int j = 0; j < np; j++)
				omega[i][j] = omega_backup[index_previous[i]][j];
		}
		else
		{
			for (unsigned int j = 0; j < np; j++)
				omega[i][j] = 0.;
		}
	}

	const double epsilon = 1.e-4;
	unsigned int iInert = 0;

	// Looking for the inert species
	{
		double sum_max = 0.;
		for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
		{
			double sum = 0.;
			for (unsigned int j = 0; j < np - 1; j++)
				sum += 0.50*(omega[i][j] + omega[i][j + 1])*(x[j + 1] - x[j]);
			if (sum > sum_max)
			{
				iInert = i;
				sum_max = sum;
			}
		}
	}
	for (unsigned int j = 0; j < np; j++)
	{
		double sum = 0.;
		for (unsigned int i = 0; i < thermodynamicsMap.NumberOfSpecies(); i++)
			sum += omega[i][j];

		if (sum > 1. + epsilon || sum <= 0.)
		{
			std::cout << "Recovering composition from backup file: Fatal error!" << std::endl;
			std::cout << "Mesh point: " << j << " - Sum of mass fractions: " << sum << " - Error : " << sum - 1. << std::endl;
			OpenSMOKE::FatalErrorMessage("Wrong sum of mass fractions");
		}
		omega[iInert][j] += (1. - sum);
	}
}

void ReadFromBackupFile(const boost::filesystem::path path_file, std::vector<double>& x, std::vector<double>& T, std::vector<double>& P, std::vector<double>& m, std::vector< std::vector<double> >& omega, std::vector<std::string>& names_species)
{
	rapidxml::xml_document<> xml_main_input;
	std::vector<char> local_xml_input_string;
	OpenSMOKE::OpenInputFileXML(xml_main_input, local_xml_input_string, path_file);

	unsigned int index_T;
	unsigned int index_P;
	unsigned int index_MW;
	unsigned int number_of_massfractions_profiles;
	unsigned int number_of_additional_profiles;
	unsigned int index_mass_flow_rate = 7;				// [?]

	// Indices of T, P and MW
	{
		rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("t-p-mw");
		if (indices_node != 0)
		{
			std::stringstream values(indices_node->value());
			values >> index_T;
			values >> index_P;
			values >> index_MW;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the t-p-mw leaf");
	}

	// Additional
	{
		rapidxml::xml_node<>* additional_node = xml_main_input.first_node("opensmoke")->first_node("additional");
		if (additional_node != 0)
		{
			std::stringstream values(additional_node->value());
			values >> number_of_additional_profiles;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the additional leaf");
	}

	// Species (mass fractions)
	{
		rapidxml::xml_node<>* massfractions_node = xml_main_input.first_node("opensmoke")->first_node("mass-fractions");
		if (massfractions_node != 0)
		{
			std::stringstream values(massfractions_node->value());

			values >> number_of_massfractions_profiles;

			names_species.resize(number_of_massfractions_profiles);
			for (unsigned int j = 0; j<number_of_massfractions_profiles; j++)
			{
				std::string dummy;
				values >> dummy;
				names_species[j] = dummy;

				double dummy_double;
				double dummy_unsigned_int;
				values >> dummy_double;
				values >> dummy_unsigned_int;
			}
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the mass-fractions leaf");
	}

	// Read profiles
	std::vector< std::vector<double> > additional;
	{
		unsigned int number_of_abscissas;
		unsigned int number_of_ordinates;

		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("profiles-size");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				values >> number_of_abscissas;
				values >> number_of_ordinates;
			}
			else
				OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the profiles-size leaf");
		}

		omega.resize(number_of_massfractions_profiles);
		for (unsigned int j = 0; j<number_of_massfractions_profiles; j++)
			omega[j].resize(number_of_abscissas);

		additional.resize(number_of_additional_profiles);
		for (unsigned int j = 0; j<number_of_additional_profiles; j++)
			additional[j].resize(number_of_abscissas);

		rapidxml::xml_node<>* profiles_node = xml_main_input.first_node("opensmoke")->first_node("profiles");
		if (profiles_node != 0)
		{
			std::stringstream values(profiles_node->value());
			for (unsigned int i = 0; i<number_of_abscissas; i++)
			{
				for (unsigned int j = 0; j<number_of_additional_profiles; j++)
					values >> additional[j][i];
				for (unsigned int j = 0; j<number_of_massfractions_profiles; j++)
					values >> omega[j][i];
			}
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the profiles leaf");

		x = additional[0];
		std::transform(x.begin(), x.end(), x.begin(), std::bind1st(std::multiplies<double>(), 0.01));

		T = additional[index_T];
		P = additional[index_P];
		m = additional[index_mass_flow_rate];
	}
}
