//ciao
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
|	License                                                               |
|                                                                         |
|   Copyright(C) 2018 Alberto Cuoci and Mattia Bissoli                    |
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
//ciao
#ifndef OpenSMOKE_IgnitionDelayTimes_Analyzer_H
#define	OpenSMOKE_IgnitionDelayTimes_Analyzer_H

#include "Grammar_IgnitionDelayTimes.h"

namespace OpenSMOKE
{

	//!  A class for estimating ignition delay times from temperature, pressure and composition temporal profiles
	/*!
	The purpose of this class is to estimate ignition delay times from temperature, pressure and composition
	temporal profiles in non-isothermal batch reactors (constant volume, constant pressure or user-defined volume)
	*/

	class IgnitionDelayTimes_Analyzer
	{
	public:

		/**
		*@brief Default constructor
		*/
		IgnitionDelayTimes_Analyzer();

		/**
		*@brief Setup the options from an external dictionary
		*@param dictionary external disctionary containing the options
		*@param thermodynamicsMap map containing the thermodynamic data
		*/
		template<typename Thermodynamics>
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMap);

		/**
		*@brief Returns true if the analyzer is currently active
		*/
		bool is_active() const { return is_active_; }

		/**
		*@brief Returns true if verbose option in turned on
		*/
		bool is_verbose() const { return is_verbose_; }

		/**
		*@brief Returns the idt (in s) based on the maximum temperature slope
		*/
		double temperature_slope_tau() const { return temperature_slope_max_tau_; }

		/**
		*@brief Returns the idt (in s) based on the temperature increase
		*/
		double temperature_increase_tau() const { return temperature_increase_tau_; }

		// !MF
		/**
		*@brief Returns the idt (in s) based on the maximum temperature
		*/
		double temperature_max_tau() const { return temperature_max_tau_; }

		/**
		*@brief Returns the idt (in s) based on the maximum pressure
		*/
		double pressure_max_tau() const { return pressure_max_tau_; }

		/**
		*@brief Returns the idt (in s) based on the maximum pressure slope
		*/
		double pressure_slope_tau() const { return pressure_slope_max_tau_; }

		/**
		*@brief Returns the idt (in s) based on the pressure increase
		*/
		double pressure_increase_tau() const { return pressure_increase_tau_; }

		/**
		*@brief Returns the idt (in s) based on the maximum species mole fraction
		*/
		std::vector<double> species_max_tau() const { return species_max_tau_; }

		/**
		*@brief Returns the idt (in s) based on the maximum species mole fraction slope
		*/
		std::vector<double> species_slope_tau() const { return species_slope_max_tau_; }

		/**
		*@brief Returns the species indices for the evaluation of idt
		*/
		std::vector<unsigned int> species_index() const { return species_index_; }

		/**
		*@brief Returns the idt (in s) based on the target of the absolute molefraction of species 
		*/
		std::vector<double> species_target_abs_tau() const { return species_target_abs_tau_; }
		
		/**
		*@brief Returns the idt (in s) based on the target of the relative molefraction of species 
		*/
		std::vector<double> species_target_rel_tau() const { return species_target_rel_tau_; }

		/**
		*@brief Returns the idt (in s) based on the target of the absolute concentration of species 
		*/
		std::vector<double> species_target_c_tau() const { return species_target_c_tau_; }
		
		/**
		*@brief Returns the idt (in s) based on the minimum incerception of species indicated by index
		*@param index of the species specified in the list in the OpenSMOKE++ input file
		*/

		double species_intercept_min_tau(int index) const { return species_intercept_min_tau_[index]-species_intercept_min_x_[index] / (-species_intercept_min_slope_[index]); }

		/**
		*@brief Returns the idt (in s) based on the maximum incerception of species indicated by index
		*@param index of the species specified in the list in the OpenSMOKE++ input file
		*/
		double species_intercept_max_tau(int index) const { return species_intercept_max_tau_[index]-species_intercept_max_x_[index]/species_intercept_max_slope_[index]; }
		
		/**
		*@brief Returns the list of indices of the species used in species_intercept_max_tau
		*/
		std::vector<unsigned int> species_intercept_max_index() const { return species_intercept_max_index_; }
		
		/**
		*@brief Returns the list of indices of the species used in species_intercept_min_tau
		*/
		std::vector<unsigned int> species_intercept_min_index() const { return species_intercept_min_index_; }
		// !MF

		/**
		*@brief Extracts the ignition delay time from a provided set of temperature, pressure and composition (mole fractions)
		*@param t the current time (in s)
		*@param T the current temperature (in K)
		*@param P the current pressure (in Pa)
		*@param x the current molar fractions
		*/
		void Analyze(const double t, const double T, const double P, const double* x);

		/**
		*@brief Prints the calculated ignition delay times on a file
		*@param file_name path to the output file
		*/
		void PrintOnFile(const boost::filesystem::path file_name);

		/**
		*@brief Prints the header line for parameteric analysis
		*@param counter index of column
		*@param fOutput the output file
		*/
		void PrintHeaderLine(unsigned int& counter, std::ostream& fOutput);

		/**
		*@brief Prints the ignition delay times on a file for parametric analysis
		*@param fOutputthe output file
		*/
		void Print(std::ostream& fOutput);

		/**
		*@brief Resets the ignition delay times
		*/
		void Reset();


	private:

		bool is_active_;								//!< true if the analyzer is active
		bool is_pressure_;								//!< true if the user want to estimate the ignition delay time based on the pressure slope
		bool is_temperature_;							//!< true if the user want to estimate the ignition delay time based on the temperature slope
		bool is_species_slope_;							//!< true if the user want to estimate the ignition delay time based on the species slopes

		double x_threshold_;							//!< minimum mole fraction for evaluation of i.d.t. on the basis of max mole fractions of species
		double time_minimum_;							//!< minimum time (in s) for starting the evalution of i.d.t.
		double time_minimum_interval_;					//!< minimum time interval (in s) for the evalution of i.d.t. based on slopes


		// Standard criteria: max and max-slope

		// List of species
		std::vector<unsigned int> species_index_;		//!< indices (0-index based) of selected species for evaluation of i.d.t.
		std::vector<std::string> species_name_;			//!< names of selected species for evaluation of i.d.t.

		// Maximum (on T, P, and species)
		double temperature_max_value_ = 0.;				//!< current max temperature [K]
		double temperature_max_tau_ = 0.;				//!< current values of i.d.t. based on max temperature (in s)
		double pressure_max_tau_ = 0.;					//!< current values of i.d.t. based on max pressure(in s)
		double pressure_max_value_ = 0.;				//!< current max pressure [Pa]
		std::vector<double> species_max_value_;			//!< current max values of mole fractions
		std::vector<double> species_max_tau_;			//!< current values of i.d.t. based on the species max mole fractions (in s)

		// Maximum slope (on T, P, and species)
		double temperature_slope_max_value_;			//!< current max value of temperature slope (in K/s)
		double temperature_slope_max_tau_;				//!< current value of i.d.t. based on the temperature slope (in s)
		double pressure_slope_max_value_;				//!< current max value of pressure slope (in Pa/s)
		double pressure_slope_max_tau_;					//!< current value of i.d.t. based on the pressure slope (in s)
		std::vector<double> species_slope_max_value_;	//!< current max values of mole fraction slopes (in 1/s)
		std::vector<double> species_slope_max_tau_;		//!< current values of i.d.t. based on the species slope (in s)


		// Additional criteria: targets (abs and rel), intercepts (min and max)

		// Species: target (absolute mole fraction)
		std::vector<double>					species_target_abs_value_;		//!< absolute mole fraction target values
		std::vector<double>					species_target_abs_tau_;		//!< ignition delay time (in s)
		std::vector<double>					species_target_abs_value_old_;	//!< previous mole fraction value (for interpolation purposes)
		std::vector<std::string>			species_target_abs_name_;		//!< species names
		std::vector<unsigned int>			species_target_abs_index_;		//!< indices of species (0-index based)

		// Species: target (relative mole fraction)
		std::vector<double>					species_target_rel_value_;			//!< relative (with respect to the maximum mole fraction) target value
		std::vector<double>					species_target_rel_tau_;			//!< ignition delay time (in s)
		std::vector<double>					species_target_rel_history_time_;	//!< list of times (in s)
		std::vector< std::vector<double> >	species_target_rel_history_;		//!< list of mole fraction profiles
		std::vector<std::string>			species_target_rel_name_;			//!< species names
		std::vector<unsigned int>			species_target_rel_index_;			//!< indices of species (0-index based)

		// Species: target (absolute concentration)
		std::vector<double>					species_target_c_value_;		//!< absolute mole fraction target values
		std::vector<double>					species_target_c_tau_;			//!< ignition delay time (in s)
		std::vector<double>					species_target_c_value_old_;	//!< previous mole fraction value (for interpolation purposes)
		std::vector<std::string>			species_target_c_name_;			//!< species names
		std::vector<unsigned int>			species_target_c_index_;		//!< indices of species (0-index based)

		// Species: min intercept
		std::vector<double>			species_intercept_min_tau_;		//!< ignition delay time (in s)
		std::vector<double>			species_intercept_min_x_;		//!< mole fraction where the slope is minimum
		std::vector<double>			species_intercept_min_slope_;	//!< minimum slope of mole fraction profile (in 1/s)
		std::vector<std::string>	species_intercept_min_name_;	//!< species names
		std::vector<unsigned int>	species_intercept_min_index_;	//!< indices of species (0-index based)
		std::vector<double>			species_intercept_min_x_old_;	//!< previous mole fraction value 

		// Species: max intercept
		std::vector<double>			species_intercept_max_tau_;		//!< ignition delay time (in s)
		std::vector<double>			species_intercept_max_x_;		//!< mole fraction where the slope is maximum
		std::vector<double>			species_intercept_max_slope_;	//!< maximum slope of mole fraction profile (in 1/s)
		std::vector<std::string>	species_intercept_max_name_;	//!< species names
		std::vector<unsigned int>	species_intercept_max_index_;	//!< indices of species (0-index based)
		std::vector<double>			species_intercept_max_x_old_;	//!< previous mole fraction value 

		// Temperature: relative increase (delta)
		double temperature_increase_;					//!< user-defined temperature increase for i.d.t.
		double temperature_increase_tau_;				//!< current value of i.d.t. based on the temperature increase
		double T0_;										//!< inlet/initial temperature (in K)

		// Pressure: relative increase (delta)
		double pressure_increase_;						//!< user-defined pressure increase for i.d.t.
		double pressure_increase_tau_;					//!< current value of i.d.t. based on the pressure increase
		double P0_;										//!< inlet/initial pressure (in Pa)

		// Old values (for derivatives)
		double tOld_;									//!< previous time (in s)
		double TOld_;									//!< previous temperature (in K)
		double POld_;									//!< previous pressure (in Pa)
		std::vector<double> xOld_;						//!< previous mole fractions (in s)


		// Rapid Compression Machine (RCM)
		bool is_rcm_;									//!< true if rapid compression machine analysis is required
		bool is_verbose_;								//!< if true, the regularized histories of temeprature are written on a file
		std::vector<double> history_time_;				//!< time history [s]
		std::vector<double> history_temperature_;		//!< temperature hystory [K]
		double filter_width_;							//!< filter width [s]
		double temperature_derivative_threshold_;		//!< threshold for temperature derivative [K/s]
		double regularization_dt_;						//!< regularization time interval [s]
		Eigen::MatrixXd reg_max_;						//!< ignition delay time detailes

	};
}

#include "IgnitionDelayTimes_Analyzer.hpp"


#endif	/* OpenSMOKE_IgnitionDelayTimes_Analyzer_H */
