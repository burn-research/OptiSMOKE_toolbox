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


#ifndef OpenSMOKE_ExtendedPressureLogarithmicRateExpression_H
#define	OpenSMOKE_ExtendedPressureLogarithmicRateExpression_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>

namespace OpenSMOKE
{
	//!  A class to manage a Pressure Dependence through logarithmic interpolation rate expression
	/*!
		 TODO
	*/

	class ExtendedPressureLogarithmicRateExpression
	{
		public:

			/**
			* Default constructor
			*/
			ExtendedPressureLogarithmicRateExpression();

			/**
			* Default copy constructor
			*/
			ExtendedPressureLogarithmicRateExpression(const ExtendedPressureLogarithmicRateExpression& orig);

			/**
			* Default destructor
			*/
			//virtual ~ExtendedPressureLogarithmicRateExpression();

			/**
			*@brief Prepares all the data to be used for evaluating the reaction rate
					THIS COEFFICIENT SHOULD BE THE VECTOR OF STRINGS CONTAINING THE SPECIES LIST OF THE MECHANISM */
			void Setup(std::vector<std::string> coefficients);

			/**
			*@brief Reads from a file the data about the reaction
			*/
			void ReadFromASCIIFile(std::istream& fInput);

			/**
			*@brief Writes a short summary on file
			*/
			void WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const;

			/**
			*@brief Evaluates the kinetic constant
			*/
			double KineticConstant(const double T, const double P, const double cTot, const double* c);

			/**
			*@brief Writes the reaction in CHEMKIN format
			*/
			void WriteCHEMKINOnASCIIFile(std::ostream& fOutput) const;


			//

			// RESTART FROM HERE
			/*double lnA_extendedPLOG(const unsigned int j, const unsigned int k); */


			//
			//std::vector<std::string> species_;				//!< list of bath species
			//std::vector<int> species_indices_;				//!< indices of bath species (0-based) (-1 means mixture)
			//std::vector<int> N_;							//!< number of points along the pressure axis					
			//std::vector< std::vector<double> > lnA_;		//!< logarithm (natural) of factor frequencies
			//std::vector< std::vector<double> > Beta_;		//!< temperature exponents
			//std::vector< std::vector<double> > E_over_R_;	//!< normalized activation energy [K]
			//std::vector< std::vector<double> > p_;			//!< list of pressure points [Pa]
			//std::vector< std::vector<double> > lnp_;		//!< logarithm of pressure points (for efficiency reasons only)

			std::vector<std::string> species() { return species_; }
			//int species_indices(const unsigned int j) { return species_indices_[j]; }
			std::vector<double> lnA(const unsigned int j) { return lnA_[j]; }
			std::vector<double> pressure(const unsigned int j) { return p_[j]; }
			std::vector<double> Beta(const unsigned int j) { return Beta_[j]; }
			std::vector<double> E_over_R(const unsigned int j) { return E_over_R_[j]; }


			void modify_A(int bath_gas_pos, double new_value, std::vector<double> nominal_values)
			{
				for(int i = 0; i < lnA_[bath_gas_pos].size(); i++ )
				{
					//std::cout<<"The lnA snominal value for the MX is 	"<<nominal_values[i]<<std::endl;
					lnA_[bath_gas_pos][i]= nominal_values[i]+new_value*std::log(10);
					//std::cout<<"The lnA new value for the MX is 			"<<lnA_[0][i]<<std::endl;
				}

			}

			void modify_Beta(int bath_gas_pos, double new_value, std::vector<double> nominal_values)
			{
				for(int i = 0; i < Beta_[bath_gas_pos].size(); i++ )
				{	
					//std::cout<<"The Beta nominal value for the MX is 	"<<nominal_values[i]<<std::endl;
					Beta_[bath_gas_pos][i]= nominal_values[i]+new_value;
					//std::cout<<"The Beta new value for the MX is 		"<<Beta_[0][i]<<std::endl;
				}

			}

			void modify_E_over_R(int bath_gas_pos, double new_value, std::vector<double> nominal_values)
			{
				for(int i = 0; i < E_over_R_[bath_gas_pos].size(); i++ )
				{	
					//std::cout<<"The E_over_R nominal value for the MX is 	"<<nominal_values[i]<<std::endl;
					E_over_R_[bath_gas_pos][i]= nominal_values[i]+new_value;
					//std::cout<<"The E_over_R new value for the MX is 		"<<E_over_R_[0][i]  <<std::endl;
				}

			}

			int get_bath_gas_position(std::string bath_gas){	

				int index = 0;

				for (unsigned int k = 0; k < species_.size(); k++){
					if (species_[k] == bath_gas){
						
						index = k;
						break;
					}
				}

				return index;
			}

                        /**
 *                         *@brief Evaluates the kinetic constant for a single bath gas!
 *                                                 */
                        double KineticConstant_optismoke(std::string bath_gas, const double T, const double P);

		private:

			/**
			*@brief Prepares all the data to be used for evaluating the reaction rate
			*/
			void Setup(unsigned int index, std::vector<double> coefficients, const double conversion_E);

			/**
			*@brief Evaluates the kinetic constant
			*/
			double KineticConstant(const unsigned int index, const double T, const double P);


			/**
			*@brief Writes a short summary on file
			*/
			void WriteShortSummaryOnASCIIFile(const unsigned int index, std::ostream& fOutput) const;

			/**
			*@brief Check if the kinetic parameters for the mixture are provided
			*/
			void CheckForMixtureParameters();

        private:
            
			/**
			*@brief Returns an error message
			*/
			void ErrorMessage(const std::string message);

			/**
			*@brief Returns a warning message
			*/
			void WarningMessage(const std::string message);

	private:

			double conversion_A_;							//!< conversion factor for frequency factor
			std::vector<std::string> species_;				//!< list of bath species
			std::vector<int> species_indices_;				//!< indices of bath species (0-based) (-1 means mixture)
			std::vector<int> N_;							//!< number of points along the pressure axis					
			std::vector< std::vector<double> > lnA_;		//!< logarithm (natural) of factor frequencies
			std::vector< std::vector<double> > Beta_;		//!< temperature exponents
			std::vector< std::vector<double> > E_over_R_;	//!< normalized activation energy [K]
			std::vector< std::vector<double> > p_;			//!< list of pressure points [Pa]
			std::vector< std::vector<double> > lnp_;		//!< logarithm of pressure points (for efficiency reasons only)
	};

}

#include "ExtendedPressureLogarithmicRateExpression.hpp"

#endif	/* OpenSMOKE_ExtendedPressureLogarithmicRateExpression_H */

