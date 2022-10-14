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

namespace OpenSMOKE
{

	ExtendedPressureLogarithmicRateExpression::ExtendedPressureLogarithmicRateExpression()
	{
	}
	
	ExtendedPressureLogarithmicRateExpression::ExtendedPressureLogarithmicRateExpression(const ExtendedPressureLogarithmicRateExpression& orig)
	{
	}
	
	//ExtendedPressureLogarithmicRateExpression::~ExtendedPressureLogarithmicRateExpression()
	//{
	//}

	void ExtendedPressureLogarithmicRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedPressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Error:  " << message									<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
		exit(-1);
	}

	void ExtendedPressureLogarithmicRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedPressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message								<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
	}

	/* HERE HE IS ACTUALLY SETTING UP THE FORMALISM */
	void ExtendedPressureLogarithmicRateExpression::Setup(std::vector<std::string> coefficients)
	{
		/* It divides by 6 because each arrhenius expression is described by m*6 parameters (see notebook),
		where m is the number of arrhenius we use in total for the extended formalism of PLOG. 
		while -2 is due to 2 additional parameters, which are used to convert A and Ea that is whay they are not considered 
		
		So n is actually the number of Arrhenius expressions we use to describe the reaction overall*/
		const unsigned int n = (unsigned int)((coefficients.size() - 2) / 6);

					// A has the lenght of coefficients -2 
		             conversion_A_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 2]);
					// E has the lenght of coefficients -1 
		const double conversion_E  = boost::lexical_cast<double>(coefficients[coefficients.size() - 1]);

		// Resize species to a zero length vector or INITIALIZES?
		species_.resize(0);

		// Create list of species
		// the for spans the coefficients strings vector every 6 positions because it checks the only the species
		for (unsigned int i = 0; i < coefficients.size()-2; i += 6)
		{
			// The first species is MIX so he initializes found as false
			bool found = false;
			// loop over species_, if the species in coefficient is equal to a previously introduced species it breaks the iteration
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i] == species_[k])
				{
					found = true;
					break;
				}
			
			// this part add species in positions (i) and the index of species, which is always in position i+1
			if (found == false)
			{
				//add species string
				species_.push_back(coefficients[i]);
				//add species index, which is converted to int using boost
				species_indices_.push_back(boost::lexical_cast<int>(coefficients[i + 1]));
			}
		}

		// Populate local coefficients
		// HERE He local_coefficients, which is a vector of vector of doubles
		std::vector< std::vector<double> > local_coefficients(species_.size());

		// The first dimension is detacted from the number of species
		for (unsigned int i = 0; i < coefficients.size()-2; i += 6)
		{
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i] == species_[k])
				{
					// populates pressure
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 2]));
					// populates lnA
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 3]));
					// populates n
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 4]));
					// populates Beta
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 5]));
				}
		}

		// Check for data consistency //ciao
		{
			CheckForMixtureParameters();
		}

		// Memory allocation
		N_.resize(species_.size());
		lnA_.resize(species_.size());
		Beta_.resize(species_.size());
		E_over_R_.resize(species_.size());
		p_.resize(species_.size());
		lnp_.resize(species_.size());

		// Setup for single species
		for (unsigned int k = 0; k < species_.size(); k++)
			Setup(k, local_coefficients[k], conversion_E);
	}

	void ExtendedPressureLogarithmicRateExpression::Setup(unsigned int index, std::vector<double> coefficients, const double conversion_E)
	{
		// n is the number of arrhenius expressions that are used to describe one single species
		// in fact it is established by dividing coefficients by 4, which correspond to one vector in previous local_coefficients
		// COOL: THEY CAN ALSO HAVE DIFFERENT NUMBER OF PRESSURES
		N_[index] = (unsigned int)(coefficients.size()/4);
		
		// so we resize the object containing the parameters values into a vector of vector of doubles
		// first dimension being species
		// second dimension being N_[index]
		lnA_[index].resize(N_[index]);
		Beta_[index].resize(N_[index]);
		E_over_R_[index].resize(N_[index]);
		p_[index].resize(N_[index]);
		lnp_[index].resize(N_[index]);

		const double threshold = 1.e-32;
		
		// populating the objects by adding the values from coefficients to 
		unsigned count = 0;
		for(int i=0;i<N_[index];i++)
		{
			if (coefficients[count] < 0.){
				std::cout << "The parameter is " << coefficients[count] << std::endl;
				ErrorMessage("The PLOGMX/PLOGSP option can be used only with non-negative frequency factors!");
			}

			//pressure is converted from atm to Pascal, with a little error, because it should be done from bar?
			p_[index][i] = coefficients[count++]*101325.;
			lnp_[index][i] = std::log(p_[index][i]);
			lnA_[index][i] = std::log(std::max(coefficients[count++]*conversion_A_, threshold)) ;
			Beta_[index][i] = coefficients[count++];
			E_over_R_[index][i] = (coefficients[count++] * conversion_E) / PhysicalConstants::R_J_kmol;
		}

		// Checking if the values are provided in the correct order, as pressure has to be given in increasing order
		for(int i=1;i<N_[index];i++)
			if (p_[index][i] <= p_[index][i-1])
				ErrorMessage("The points on the pressure axis (PLOGMX/PLOGSP) are not provided in the correct order!");
	}

	double ExtendedPressureLogarithmicRateExpression::KineticConstant(const double T, const double P, const double cTot, const double* c)
	{
		double ctot_minus_cspecies = cTot;
		for (unsigned int k = 0; k < species_.size(); k++)
			if (species_indices_[k] != -1)	ctot_minus_cspecies -= c[species_indices_[k]];

		double kinetic_constant = 0.;
		for (unsigned int k = 0; k < species_.size(); k++)
		{
			if (species_indices_[k] != -1)
				kinetic_constant += c[species_indices_[k]] * KineticConstant(k, T, P);
			else
				kinetic_constant += ctot_minus_cspecies * KineticConstant(k, T, P);
		}
		kinetic_constant /= cTot;

		return kinetic_constant;
	}

	// these functions i don't think i need to know what they do for now
	// ciao
	double ExtendedPressureLogarithmicRateExpression::KineticConstant(unsigned int index, const double T, const double P)
	{


		if (P <= p_[index][0])
		{
			return std::exp(lnA_[index][0] + Beta_[index][0] * std::log(T) - E_over_R_[index][0] / T);
		}
		else if (P >= p_[index][N_[index] - 1])
		{
			return std::exp(lnA_[index][N_[index] - 1] + Beta_[index][N_[index] - 1] * std::log(T) - E_over_R_[index][N_[index] - 1] / T);
		}
		else
		{
			int i=0;
			for(i=0;i<N_[index]-1;i++)
				if (P<p_[index][i+1])
					break;
			
			double ln_kA = lnA_[index][i]+Beta_[index][i]*std::log(T)-E_over_R_[index][i]/T;
			double ln_kB = lnA_[index][i+1]+Beta_[index][i+1]*std::log(T)-E_over_R_[index][i+1]/T;
			return	std::exp( ln_kA+(ln_kB-ln_kA)*(std::log(P)-lnp_[index][i])/(lnp_[index][i+1]-lnp_[index][i]) );
		}
	}

	double ExtendedPressureLogarithmicRateExpression::KineticConstant_optismoke(std::string bath_gas, const double T, const double P)
	{
		
		int index = 0;		

		for (unsigned int k = 0; k < species_.size(); k++){
			
			if (species_[k] == bath_gas){
				
				index = k;
				break;
			}	
		}
			
		if (P <= p_[index][0])
		{
			return lnA_[index][0] + Beta_[index][0] * std::log(T) - E_over_R_[index][0] / T;
		}
		else if (P >= p_[index][N_[index] - 1])
		{
			return lnA_[index][N_[index] - 1] + Beta_[index][N_[index] - 1] * std::log(T) - E_over_R_[index][N_[index] - 1] / T;
		}
		else
		{
			int i=0;
			for(i=0;i<N_[index]-1;i++)
				if (P<p_[index][i+1])
					break;
			
			double ln_kA = lnA_[index][i]+Beta_[index][i]*std::log(T)-E_over_R_[index][i]/T;
			double ln_kB = lnA_[index][i+1]+Beta_[index][i+1]*std::log(T)-E_over_R_[index][i+1]/T;
			return	ln_kA+(ln_kB-ln_kA)*(std::log(P)-lnp_[index][i])/(lnp_[index][i+1]-lnp_[index][i]);
		}
	}


	void ExtendedPressureLogarithmicRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		// initializes this vector of strings called coefficients
		std::vector<std::string> coefficients;
		// initializes a double n
		double n;
		// reads fInput until the first white space and stores it in n as a double
		fInput >> n;
		// resizes coefficients vector to n elements
		coefficients.resize(int(n));
		// it separates all the string element wise 
		for (unsigned int i = 0; i<(unsigned int)(n); i++)
			fInput >> coefficients[i];

		Setup(coefficients);
	}

	void ExtendedPressureLogarithmicRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " ";
		fOutput << "Extended Pressure Logarithmic Interpolation" << std::endl;
		for (unsigned int k = 0; k < species_.size(); k++)
			WriteShortSummaryOnASCIIFile(k, fOutput);
	}

	void ExtendedPressureLogarithmicRateExpression::WriteShortSummaryOnASCIIFile(const unsigned int index, std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " "; 
		fOutput << "Bath species: " << species_[index] << std::endl;
		for (int j=0;j<N_[index];j++)
		{
			fOutput << std::setw(9)									<< " "; 
			fOutput << std::scientific << j+1						<< "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << p_[index][j]/101325.     << "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << std::exp(lnA_[index][j])/conversion_A_ << "\t";
			fOutput << std::setw(8) << std::setprecision(2) << std::fixed << std::right << Beta_[index][j];
			fOutput << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E_over_R_[index][j] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal << std::endl;
		}
	}

	void ExtendedPressureLogarithmicRateExpression::WriteCHEMKINOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput.unsetf(std::ios_base::floatfield);
		fOutput.precision(6);

		for (unsigned int k = 0; k < species_.size(); k++)
		{
			for (int j = 0; j < N_[k]; j++)
			{
				if (species_indices_[k] == -1)
				{
					fOutput << " PLOGMX / ";
				}
				else
				{
					fOutput << " PLOGSP / ";
					fOutput << species_[k] << "  ";
				}

				fOutput << std::showpoint << std::setw(12) << std::left << p_[k][j] / 101325.;
				fOutput << std::showpoint << std::setw(12) << std::left << std::exp(lnA_[k][j]) / conversion_A_;
				fOutput << std::showpoint << std::setw(12) << std::left << Beta_[k][j];
				fOutput << std::showpoint << std::setw(12) << std::left << E_over_R_[k][j] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal;;

				fOutput << "/" << std::endl;
			}
		}
	}

	// Check if the kinetic parameters for the mixture are provided
	void ExtendedPressureLogarithmicRateExpression::CheckForMixtureParameters()
	{
		// checks for all the species if at least one has the index -1
		// meaning if the kinetic parameters for MIX are provided, is yes it returns, if not "error message"
		for (unsigned int k = 0; k < species_.size(); k++)
			if (species_indices_[k] == -1)
				return;

		ErrorMessage("The kinetic parameters have to be provided for the mixture (through PLOGMX)");
	}

}
