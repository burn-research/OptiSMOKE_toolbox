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

bool compare_head(const Eigen::VectorXd& lhs, const Eigen::VectorXd& rhs)
{
	return lhs(0) > rhs(0);
}

Eigen::MatrixXd sorted_rows_by_head(Eigen::MatrixXd A)
{
	std::vector<Eigen::VectorXd> vec;
	for (int64_t i = 0; i < A.rows(); ++i)
		vec.push_back(A.row(i));

	std::sort(vec.begin(), vec.end(), &compare_head);

	for (int64_t i = 0; i < A.rows(); ++i)
		A.row(i) = vec[i];

	return A;
}
namespace OpenSMOKE
{
	IgnitionDelayTimes_Analyzer::IgnitionDelayTimes_Analyzer()
	{
		is_active_ = false;
		is_temperature_ = true;
		is_pressure_ = true;
		is_species_slope_ = true;

		x_threshold_ = 1.e-12;
		time_minimum_ = 1.e-12;
		time_minimum_interval_ = 1.e-12;

		temperature_increase_ = 0.;
		pressure_increase_ = 0.;

		is_rcm_ = false;
		filter_width_ = 1.e-4;						// filter width [s]
		temperature_derivative_threshold_ = 1.e3;	// threshold for temperature derivative [K/s]
		regularization_dt_ = 1.e-3;					// regularization time interval [s]
		is_verbose_ = false;						// if true, the regularized profiles are written on a file

		Reset();
	}

	void IgnitionDelayTimes_Analyzer::Reset()
	{
		tOld_ = 0.;
		TOld_ = 0.;
		POld_ = 0.;
		std::fill(xOld_.begin(), xOld_.end(), 0.);

		temperature_slope_max_value_ = 0.;
		pressure_slope_max_value_ = 0.;
		temperature_max_value_ = 0.;
		pressure_max_value_ = 0.;
		std::fill(species_max_value_.begin(), species_max_value_.end(), 0.);
		std::fill(species_slope_max_value_.begin(), species_slope_max_value_.end(), 0.);

		temperature_slope_max_tau_ = 0.;
		pressure_slope_max_tau_ = 0.;
		temperature_max_tau_ = 0.;
		pressure_max_tau_ = 0.;
		std::fill(species_max_tau_.begin(), species_max_tau_.end(), 0.);
		std::fill(species_slope_max_tau_.begin(), species_slope_max_tau_.end(), 0.);

		std::fill(species_target_abs_tau_.begin(), species_target_abs_tau_.end(), 0.);
		std::fill(species_target_rel_tau_.begin(), species_target_rel_tau_.end(), 0.);
		std::fill(species_target_abs_value_old_.begin(), species_target_abs_value_old_.end(), 0.);

		std::fill(species_target_c_tau_.begin(), species_target_c_tau_.end(), 0.);
		std::fill(species_target_c_value_old_.begin(), species_target_c_value_old_.end(), 0.);

		std::fill(species_intercept_max_tau_.begin(), species_intercept_max_tau_.end(), 0.);
		std::fill(species_intercept_max_x_.begin(), species_intercept_max_x_.end(), 0.);
		std::fill(species_intercept_max_slope_.begin(), species_intercept_max_slope_.end(), 0.);
		std::fill(species_intercept_max_x_old_.begin(), species_intercept_max_x_old_.end(), 0.);

		std::fill(species_intercept_min_tau_.begin(), species_intercept_min_tau_.end(), 0.);
		std::fill(species_intercept_min_x_.begin(), species_intercept_min_x_.end(), 0.);
		std::fill(species_intercept_min_slope_.begin(), species_intercept_min_slope_.end(), 0.);
		std::fill(species_intercept_min_x_old_.begin(), species_intercept_min_x_old_.end(), 0.);

		T0_ = 0.;
		temperature_increase_tau_ = 0.;
		P0_ = 0.;
		pressure_increase_tau_ = 0.;

		history_time_.resize(0);
		history_temperature_.resize(0);

		if (species_target_rel_index_.size() != 0)
		{
			species_target_rel_history_time_.resize(0);
			for (unsigned int i=0;i< species_target_rel_index_.size();i++)
				species_target_rel_history_[i].resize(0);
		}
	}

	template<typename Thermodynamics>
	void IgnitionDelayTimes_Analyzer::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML)
	{
		Grammar_IgnitionDelayTimes grammar_idts;
		dictionary.SetGrammar(grammar_idts);

		// Temperature
		if (dictionary.CheckOption("@Temperature") == true)
		{
			dictionary.ReadBool("@Temperature", is_temperature_);
			is_active_ = true;
		}

		// Pressure
		if (dictionary.CheckOption("@Pressure") == true)
		{
			dictionary.ReadBool("@Pressure", is_pressure_);
			is_active_ = true;
		}

		// Max species
		if (dictionary.CheckOption("@Species") == true)
		{
			dictionary.ReadOption("@Species", species_name_);

			for (unsigned int i = 0; i < species_name_.size(); i++)
				species_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_name_[i]) - 1);

			species_max_value_.resize(species_name_.size());
			species_max_tau_.resize(species_name_.size());
			species_slope_max_value_.resize(species_name_.size());
			species_slope_max_tau_.resize(species_name_.size());
			xOld_.resize(species_name_.size());

			is_active_ = true;
		}

		// Max species: intercept
		if (dictionary.CheckOption("@SpeciesMaxIntercept") == true)
		{
			dictionary.ReadOption("@SpeciesMaxIntercept", species_intercept_max_name_);

			for (unsigned int i = 0; i < species_intercept_max_name_.size(); i++)
				species_intercept_max_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_intercept_max_name_[i]) - 1);

			species_intercept_max_tau_.resize(species_intercept_max_name_.size());
			species_intercept_max_x_.resize(species_intercept_max_name_.size());
			species_intercept_max_slope_.resize(species_intercept_max_name_.size());
			species_intercept_max_x_old_.resize(species_intercept_max_name_.size());

			is_active_ = true;
		}

		// Min species: intercept
		if (dictionary.CheckOption("@SpeciesMinIntercept") == true)
		{
			dictionary.ReadOption("@SpeciesMinIntercept", species_intercept_min_name_);

			for (unsigned int i = 0; i < species_intercept_min_name_.size(); i++)
				species_intercept_min_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_intercept_min_name_[i]) - 1);

			species_intercept_min_tau_.resize(species_intercept_min_name_.size());
			species_intercept_min_x_.resize(species_intercept_min_name_.size());
			species_intercept_min_slope_.resize(species_intercept_min_name_.size());
			species_intercept_min_x_old_.resize(species_intercept_min_name_.size());
			
			is_active_ = true;
		}

		// Slope species
		if (dictionary.CheckOption("@SpeciesSlope") == true)
		{
			dictionary.ReadBool("@SpeciesSlope", is_species_slope_);
			is_active_ = true;
		}

		// Max species
		if (dictionary.CheckOption("@SpeciesThreshold") == true)
			dictionary.ReadDouble("@SpeciesThreshold", x_threshold_);

		// Minimum time for calculations
		if (dictionary.CheckOption("@MinimumTime") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@MinimumTime", value, units);
			if (units == "s")		  time_minimum_ = value;
			else if (units == "ms")   time_minimum_ = value / 1000.;
			else if (units == "min")  time_minimum_ = value * 60.;
			else if (units == "h")    time_minimum_ = value * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		// Minimum time interval for derivative calculations
		if (dictionary.CheckOption("@MinimumTimeInterval") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@MinimumTimeInterval", value, units);
			if (units == "s")		  time_minimum_interval_ = value;
			else if (units == "ms")   time_minimum_interval_ = value / 1000.;
			else if (units == "min")  time_minimum_interval_ = value * 60.;
			else if (units == "h")    time_minimum_interval_ = value * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		// Temperature
		if (dictionary.CheckOption("@TemperatureIncrease") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@TemperatureIncrease", value, units);
			if (units == "K")		 temperature_increase_ = value;
			else if (units == "C")   temperature_increase_ = value;
			else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
		}

		// Pressure
		if (dictionary.CheckOption("@PressureIncrease") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@PressureIncrease", value, units);
			if (units == "Pa")			pressure_increase_ = value;
			else if (units == "atm")	pressure_increase_ = value * 101325.;
			else if (units == "bar")	pressure_increase_ = value * 100000.;
			else OpenSMOKE::FatalErrorMessage("Unknown pressure units");
		}

		if (dictionary.CheckOption("@TargetMoleFractions") == true)
		{
			dictionary.ReadOption("@TargetMoleFractions", species_target_abs_name_, species_target_abs_value_);

			for (unsigned int i = 0; i < species_target_abs_name_.size(); i++)
				species_target_abs_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_target_abs_name_[i]) - 1);

			species_target_abs_tau_.resize(species_target_abs_name_.size());
			species_target_abs_value_old_.resize(species_target_abs_name_.size());

			is_active_ = true;
		}

		if (dictionary.CheckOption("@TargetRelativeMoleFractions") == true)
		{
			dictionary.ReadOption("@TargetRelativeMoleFractions", species_target_rel_name_, species_target_rel_value_);

			for (unsigned int i = 0; i < species_target_rel_name_.size(); i++)
				species_target_rel_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_target_rel_name_[i]) - 1);

			species_target_rel_tau_.resize(species_target_rel_name_.size());
			
			// History
			species_target_rel_history_.resize(species_target_rel_name_.size());
			species_target_rel_history_time_.reserve(3000);
			for (unsigned int i = 0; i < species_target_rel_name_.size(); i++)
				species_target_rel_history_[i].reserve(3000);

			is_active_ = true;
		}

		if (dictionary.CheckOption("@TargetConcentrations") == true)
		{
			std::vector<std::string> list_values;
			dictionary.ReadOption("@TargetConcentrations", list_values);

			if (list_values.size() % 2 == 0)
				OpenSMOKE::FatalErrorMessage("@TargetConcentrations option has a wrong number of elements.");

			const int n = (list_values.size() - 1) / 2;
			for (unsigned int i = 0; i < n; i++)
			{
				species_target_c_name_.push_back(list_values[i*2]);
				species_target_c_value_.push_back(boost::lexical_cast<double>(list_values[i*2+1]));
			}

			const std::string units = list_values[list_values.size() - 1];
			for (unsigned int i = 0; i < species_target_c_name_.size(); i++)
			{
				if (units == "kmol/m3")			species_target_c_value_[i] *= 1.e+0;
				else if (units == "mol/m3")		species_target_c_value_[i] *= 1.e-3;
				else if (units == "mol/dm3")	species_target_c_value_[i] *= 1.e+0;
				else if (units == "mol/cm3")	species_target_c_value_[i] *= 1.e+3;
				else OpenSMOKE::FatalErrorMessage("Unknown concentration units. Available units: kmol/m3 | mol/m3 | mol/cm3 | mol/dm3");
			}

			for (unsigned int i = 0; i < species_target_c_name_.size(); i++)
				species_target_c_index_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_target_c_name_[i]) - 1);

			species_target_c_tau_.resize(species_target_c_name_.size());
			species_target_c_value_old_.resize(species_target_c_name_.size());

			is_active_ = true;
		}

		// Rapid Compression Machine
		if (dictionary.CheckOption("@RapidCompressionMachine") == true)
		{
			dictionary.ReadBool("@RapidCompressionMachine", is_rcm_);

			history_time_.reserve(3000);
			history_temperature_.reserve(3000);

			if (dictionary.CheckOption("@FilterWidth") == true)
			{
				std::string units;
				double value;
				dictionary.ReadMeasure("@FilterWidth", value, units);
				if (units == "s")		  filter_width_ = value;
				else if (units == "ms")   filter_width_ = value / 1000.;
				else if (units == "min")  filter_width_ = value * 60.;
				else if (units == "h")    filter_width_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("@FilterWidth: unknown time units. Allowed units: s | ms | min | h");
			}

			if (dictionary.CheckOption("@RegularizationTimeInterval") == true)
			{
				std::string units;
				double value;
				dictionary.ReadMeasure("@RegularizationTimeInterval", value, units);
				if (units == "s")		  regularization_dt_ = value;
				else if (units == "ms")   regularization_dt_ = value / 1000.;
				else if (units == "min")  regularization_dt_ = value * 60.;
				else if (units == "h")    regularization_dt_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("@RegularizationTimeInterval: unknown time units. Allowed units: s | ms | min | h");
			}

			if (dictionary.CheckOption("@TemperatureDerivativeThreshold") == true)
			{
				std::string units;
				double value;
				dictionary.ReadMeasure("@TemperatureDerivativeThreshold", value, units);
				if (units == "K/s")		  temperature_derivative_threshold_ = value;
				else if (units == "C/s")  temperature_derivative_threshold_ = value;
				else if (units == "K/ms") temperature_derivative_threshold_ = value * 1000.;
				else if (units == "C/ms") temperature_derivative_threshold_ = value * 1000.;
				else OpenSMOKE::FatalErrorMessage("@TemperatureDerivativeThreshold: unknown time units. Allowed units: K/s | K/ms | C/s | C/ms");
			}

			if (dictionary.CheckOption("@Verbose") == true)
				dictionary.ReadBool("@Verbose", is_verbose_);
		}
	}

	void IgnitionDelayTimes_Analyzer::Analyze(const double t, const double T, const double P, const double* x)
	{
		// Initial temperature and pressure
		if (T0_ == 0.)	T0_ = T;
		if (P0_ == 0.)	P0_ = P;

		// Maximum criteria
		if (t>time_minimum_)
		{
			if (is_temperature_ == true)
			{
				if ( (T > temperature_max_value_) && (T - T0_) > 1.)
				{
					temperature_max_value_ = T;
					temperature_max_tau_ = t;
				}
			}

			if (is_pressure_ == true)
			{
				if ( (P > pressure_max_value_) && (P - P0_) > 10.)
				{
					pressure_max_value_ = P;
					pressure_max_tau_ = t;
				}
			}

			for (unsigned int i = 0; i < species_index_.size(); i++)
			{
				if ((x[species_index_[i]] > species_max_value_[i]) && (x[species_index_[i]]>x_threshold_))
				{
					species_max_value_[i] = x[species_index_[i]];
					species_max_tau_[i] = t;
				}
			}
		}

		// Criteria based on slopes and targets
		if ((t>time_minimum_) && (t - tOld_) > time_minimum_interval_)
		{
			if (is_temperature_ == true)
			{
				const double temperature_slope_current = (T - TOld_) / (t - tOld_);
				if (temperature_slope_current > temperature_slope_max_value_)
				{
					temperature_slope_max_value_ = temperature_slope_current;
					temperature_slope_max_tau_ = t;
				}
			}

			if (is_pressure_ == true)
			{
				const double pressure_slope_current = (P - POld_) / (t - tOld_);
				if (pressure_slope_current > pressure_slope_max_value_)
				{
					pressure_slope_max_value_ = pressure_slope_current;
					pressure_slope_max_tau_ = t;
				}
			}

			if (temperature_increase_ != 0. && temperature_increase_tau_ == 0.)
			{
				if (T > (T0_ + temperature_increase_))
					temperature_increase_tau_ = tOld_ + (T0_ + temperature_increase_ - TOld_)/ (T - TOld_)*(t - tOld_);
			}
			//ciao
			if (pressure_increase_ != 0. && pressure_increase_tau_ == 0.)
			{
				if (P > (P0_ + pressure_increase_))
					pressure_increase_tau_ = tOld_ + (P0_ + pressure_increase_ - POld_)/ (P - POld_)*(t - tOld_);
			}

			// Species
			if (is_species_slope_ == true)
			{
				for (unsigned int i = 0; i < species_index_.size(); i++)
				{
					const double species_slope_current = (x[species_index_[i]] - xOld_[i]) / (t - tOld_);

					if ((species_slope_current > species_slope_max_value_[i]) && (x[species_index_[i]] > x_threshold_))
					{
						species_slope_max_value_[i] = species_slope_current;
						species_slope_max_tau_[i] = t;
					}
				}
			}

			// Target species: mole fractions
			for (unsigned int i = 0; i < species_target_abs_index_.size(); i++)
			{
				if ((x[species_target_abs_index_[i]] > species_target_abs_value_[i]) && (species_target_abs_tau_[i] == 0.))
				{
					const double m = (x[species_target_abs_index_[i]] - species_target_abs_value_old_[i]) / (t - tOld_);
					const double q = species_target_abs_value_old_[i] - m * tOld_;
					species_target_abs_tau_[i] = (species_target_abs_value_[i] - q) / m;
				}
			}

			// Target species: mole fractions
			for (unsigned int i = 0; i < species_target_c_index_.size(); i++)
			{
				const double cTot = P / PhysicalConstants::R_J_kmol / T;	// total concentration (kmol/m3)

				if ((x[species_target_c_index_[i]]*cTot > species_target_c_value_[i]) && (species_target_c_tau_[i] == 0.))
				{
					const double m = (x[species_target_c_index_[i]]*cTot - species_target_c_value_old_[i]) / (t - tOld_);
					const double q = species_target_c_value_old_[i] - m * tOld_;
					species_target_c_tau_[i] = (species_target_c_value_[i] - q) / m;
				}
			}

			// Species intercept: max
			if (species_intercept_max_index_.size() != 0)
			{
				for (unsigned int i = 0; i < species_intercept_max_index_.size(); i++)
				{
					const double species_slope_current = (x[species_intercept_max_index_[i]] - species_intercept_max_x_old_[i]) / (t - tOld_);

					if ((species_slope_current > species_intercept_max_slope_[i]) && (x[species_intercept_max_index_[i]] > x_threshold_))
					{
						species_intercept_max_slope_[i] = species_slope_current;
						species_intercept_max_tau_[i] = t;
						species_intercept_max_x_[i] = x[species_intercept_max_index_[i]];
					}
				}
			}

			// Species intercept: min
			if (species_intercept_min_index_.size() != 0)
			{
				for (unsigned int i = 0; i < species_intercept_min_index_.size(); i++)
				{
					const double species_slope_current = -(x[species_intercept_min_index_[i]] - species_intercept_min_x_old_[i]) / (t - tOld_);

					if ((species_slope_current > species_intercept_min_slope_[i]) && (x[species_intercept_min_index_[i]] > x_threshold_))
					{
						species_intercept_min_slope_[i] = species_slope_current;
						species_intercept_min_tau_[i] = t;
						species_intercept_min_x_[i] = x[species_intercept_min_index_[i]];
					}
				}
			}
		}

		// Save data for rapid compression machine
		if (is_rcm_ == true)
		{
			const double min_t = 1.e-8;			// minimum time (below this there is no sense in looking at idts)
			const double min_deltat = 1.e-10;	// minimum time separation between two successive points
			const double min_T = 300.;			// minimum temperature (below this value we assume there is no ignition)

			if ( (t > min_t) && ((t - tOld_) > min_deltat) && (T>min_T) )
			{
				history_time_.push_back(t);
				history_temperature_.push_back(T);
			}
		}

		// Save data
		if (species_target_rel_index_.size() != 0)
		{
			const double min_t = 1.e-8;			// minimum time (below this there is no sense in looking at idts)
			const double min_deltat = 1.e-10;	// minimum time separation between two successive points
			const double min_T = 300.;			// minimum temperature (below this value we assume there is no ignition)

			if ((t > min_t) && ((t - tOld_) > min_deltat) && (T > min_T))
			{
				species_target_rel_history_time_.push_back(t);
				for (unsigned int i = 0; i < species_target_rel_index_.size(); i++)
					species_target_rel_history_[i].push_back(x[species_target_rel_index_[i]]);
			}
		}

		// Save previous values
		tOld_ = t;
		TOld_ = T;
		POld_ = P;
		for (unsigned int i = 0; i < species_index_.size(); i++)
			xOld_[i] = x[species_index_[i]];
		
		// Save previous values: target species
		for (unsigned int i = 0; i < species_target_abs_index_.size(); i++)
			species_target_abs_value_old_[i] = x[species_target_abs_index_[i]];

		// Save previous values: target species
		for (unsigned int i = 0; i < species_target_c_index_.size(); i++)
		{
			const double cTot = P / PhysicalConstants::R_J_kmol / T;
			species_target_c_value_old_[i] = x[species_target_c_index_[i]]*cTot;
		}

		// Save previous values: intercept
		for (unsigned int i = 0; i < species_intercept_max_index_.size(); i++)
			species_intercept_max_x_old_[i] = x[species_intercept_max_index_[i]];
		for (unsigned int i = 0; i < species_intercept_min_index_.size(); i++)
			species_intercept_min_x_old_[i] = x[species_intercept_min_index_[i]];
	}

	void IgnitionDelayTimes_Analyzer::PrintOnFile(const boost::filesystem::path file_name)
	{
		if (is_rcm_ == false)
		{
			std::ofstream fOutput(file_name.c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);

			fOutput << "------------------------------------------------------------" << std::endl;
			fOutput << "Criterion                     Tau[s]         Value          " << std::endl;
			fOutput << "------------------------------------------------------------" << std::endl;

			if (is_temperature_ == true)
			{
				fOutput << std::setw(30) << std::left << "T(max)";
				fOutput << std::setw(15) << std::left << temperature_max_tau_;
				fOutput << std::setw(15) << std::left << temperature_max_value_;
				fOutput << std::endl;

				fOutput << std::setw(30) << std::left << "T(slope)";
				fOutput << std::setw(15) << std::left << temperature_slope_max_tau_;
				fOutput << std::setw(15) << std::left << temperature_slope_max_value_;
				fOutput << std::endl;
			}

			if (is_pressure_ == true)
			{
				fOutput << std::setw(30) << std::left << "P(max)";
				fOutput << std::setw(15) << std::left << pressure_max_tau_;
				fOutput << std::setw(15) << std::left << pressure_max_value_;
				fOutput << std::endl;

				fOutput << std::setw(30) << std::left << "P(slope)";
				fOutput << std::setw(15) << std::left << pressure_slope_max_tau_;
				fOutput << std::setw(15) << std::left << pressure_slope_max_value_;
				fOutput << std::endl;
			}

			for (unsigned int i = 0; i < species_index_.size(); i++)
			{
				if (is_species_slope_ == true)
				{
					std::string label = species_name_[i] + "(slope)";
					fOutput << std::setw(30) << std::left << label;
					fOutput << std::setw(15) << std::left << species_slope_max_tau_[i];
					fOutput << std::setw(15) << std::left << species_slope_max_value_[i];
					fOutput << std::endl;
				}

				{
					std::string label = species_name_[i] + "(max)";
					fOutput << std::setw(30) << std::left << label;
					fOutput << std::setw(15) << std::left << species_max_tau_[i];
					fOutput << std::setw(15) << std::left << species_max_value_[i];
					fOutput << std::endl;
				}
			}

			for (unsigned int i = 0; i < species_target_abs_index_.size(); i++)
			{
				std::string label = species_target_abs_name_[i] + "(tgt-abs)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << species_target_abs_tau_[i];
				fOutput << std::setw(15) << std::left << species_target_abs_value_[i];
				fOutput << std::endl;
			}

			for (unsigned int i = 0; i < species_target_rel_index_.size(); i++)
			{
				// Find maximum value
				const double max_value = *std::max_element(species_target_rel_history_[i].begin(), species_target_rel_history_[i].end());
				const double target_value = species_target_rel_value_[i] * max_value;

				// Interpolation
				unsigned int index = 0;
				for (unsigned int j=1;j< species_target_rel_history_time_.size();j++)
					if (species_target_rel_history_[i][j]  >= target_value)
					{
						if (species_target_rel_history_time_[j] - species_target_rel_history_time_[j - 1] >= 1e-10)
						{
							const double m = (species_target_rel_history_[i][j] - species_target_rel_history_[i][j - 1]) /
								(species_target_rel_history_time_[j] - species_target_rel_history_time_[j - 1]);
							const double q = species_target_rel_history_[i][j - 1] - m * species_target_rel_history_time_[j - 1];
							species_target_rel_tau_[i] = (target_value - q) / m;
						}
						else
						{
							species_target_rel_tau_[i] = species_target_rel_history_time_[j];
						}

						break;
					}

				std::string label = species_target_rel_name_[i] + "(tgt-rel)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << species_target_rel_tau_[i];
				fOutput << std::setw(15) << std::left << species_target_rel_value_[i];
				fOutput << std::endl;
			}

			for (unsigned int i = 0; i < species_target_c_index_.size(); i++)
			{
				std::string label = species_target_c_name_[i] + "(tgt-c)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << species_target_c_tau_[i];
				fOutput << std::setw(15) << std::left << species_target_c_value_[i];
				fOutput << std::endl;
			}

			for (unsigned int i = 0; i < species_intercept_max_index_.size(); i++)
			{
				double tau = 0.;
				if (species_intercept_max_slope_[i] > 0.)
					tau = species_intercept_max_tau_[i] -
					species_intercept_max_x_[i] / species_intercept_max_slope_[i];

				std::string label = species_intercept_max_name_[i] + "(int-max)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << tau;
				fOutput << std::setw(15) << std::left << species_intercept_max_slope_[i];
				fOutput << std::endl;
			}

			for (unsigned int i = 0; i < species_intercept_min_index_.size(); i++)
			{
				double tau = 0.;
				if (species_intercept_min_slope_[i] > 0.)
					tau = species_intercept_min_tau_[i] -
					species_intercept_min_x_[i] / (-species_intercept_min_slope_[i]);

				std::string label = species_intercept_min_name_[i] + "(int-min)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << tau;
				fOutput << std::setw(15) << std::left << -species_intercept_min_slope_[i];
				fOutput << std::endl;
			}
			
			if (temperature_increase_ != 0.)
			{
				fOutput << std::setw(30) << std::left << "T(increase)";
				fOutput << std::setw(15) << std::left << temperature_increase_tau_;
				fOutput << std::setw(15) << std::left << T0_ + temperature_increase_;
				fOutput << std::endl;
			}

			if (pressure_increase_ != 0.)
			{
				fOutput << std::setw(30) << std::left << "P(increase)";
				fOutput << std::setw(15) << std::left << pressure_increase_tau_;
				fOutput << std::setw(15) << std::left << P0_ + pressure_increase_;
				fOutput << std::endl;
			}

			fOutput.close();
		}

		if (is_rcm_ == true)
		{
			// Preliminary filtering
			{
				std::cout << " * RCM: Preliminary filtering on history files: width=" << filter_width_ << " s" << std::endl;

				const unsigned int n = history_time_.size();
				std::vector<double> history_time_filtered_(n);
				std::vector<double> history_temperature_filtered_(n);

				OpenSMOKE::ApplyFilter(filter_width_, 5, n,
					history_time_filtered_.data(), history_temperature_filtered_.data(),
					history_time_.data(), history_temperature_.data());

				history_time_ = history_time_filtered_;
				history_temperature_ = history_temperature_filtered_;
			}
			
			// Regularization
			{
				std::cout << " * RCM: Regularization of temperature profile: dt=" << regularization_dt_ << " s" << std::endl;

				FixedProfile profile_temperature(history_time_.size(), history_time_.data(), history_temperature_.data());

				std::vector<double>reg_time;
				std::vector<double>reg_temperature;

				unsigned int count = 0;
				for (;;)
				{
					count++;
					const double t = regularization_dt_ * count;
					if (t >= history_time_[history_time_.size() - 1])
						break;

					reg_time.push_back(t);
					reg_temperature.push_back(profile_temperature.Interpolate(t));
				}

				std::vector<double>reg_dtemperature(reg_temperature.size());
				std::vector<double>reg_dlog10temperature(reg_temperature.size());
				std::fill(reg_dtemperature.begin(), reg_dtemperature.end(), 0.);
				std::fill(reg_dlog10temperature.begin(), reg_dlog10temperature.end(), 0.);
				for (unsigned int i = 1; i < reg_temperature.size() - 1; i++)
				{
					reg_dtemperature[i] = (reg_temperature[i + 1] - reg_temperature[i]) /
						(reg_time[i + 1] - reg_time[i]);

					if (reg_dtemperature[i] > 0.)
						reg_dlog10temperature[i] = std::log10(reg_dtemperature[i]);
					else
						reg_dlog10temperature[i] = -std::log10(-reg_dtemperature[i]);
				}

				// Search for maxima
				{
					std::vector<double> reg_max_time;
					std::vector<double> reg_max_temperature;
					std::vector<double> reg_max_dtemperature;
					std::vector<double> reg_max_dlog10temperature;

					const double max_dlog10temperature_threshold = std::log10(temperature_derivative_threshold_);

					const unsigned int n = reg_time.size();
					for (unsigned int i = 1; i < n - 1; i++)
					{
						if (reg_dlog10temperature[i] > reg_dlog10temperature[i - 1] && reg_dlog10temperature[i] > reg_dlog10temperature[i + 1])
							if (reg_dlog10temperature[i] > max_dlog10temperature_threshold)
							{
								reg_max_time.push_back(reg_time[i]);
								reg_max_temperature.push_back(reg_temperature[i]);
								reg_max_dtemperature.push_back(reg_dtemperature[i]);
								reg_max_dlog10temperature.push_back(reg_dlog10temperature[i]);
							}
					}

					Eigen::MatrixXd temp(reg_max_time.size(), 4);
					for (unsigned int i = 0; i < reg_max_time.size(); i++)
					{
						temp(i, 0) = reg_max_dlog10temperature[i];
						temp(i, 1) = reg_max_time[i];
						temp(i, 2) = reg_max_temperature[i];
						temp(i, 3) = reg_max_dtemperature[i];
					}

					reg_max_ = sorted_rows_by_head(temp);
				}

				if (is_verbose_ == true)
				{
					boost::filesystem::path file_name_rcm = file_name;
					file_name_rcm = file_name_rcm.replace_extension("hst");

					std::ofstream fOutput(file_name_rcm.c_str(), std::ios::out);
					fOutput.setf(std::ios::scientific);

					fOutput << "Time[s](1)      T[K](2)         dTdt[K/s](3)    dlog10Tdt[K/s](4)" << std::endl;
					for (unsigned int i = 0; i < reg_time.size(); i++)
					{
						fOutput << std::setw(16) << std::scientific << std::left << reg_time[i];
						fOutput << std::setw(16) << std::scientific << std::left << reg_temperature[i];
						fOutput << std::setw(16) << std::scientific << std::left << reg_dtemperature[i];
						fOutput << std::setw(16) << std::scientific << std::left << reg_dlog10temperature[i];
						fOutput << std::endl;
					}

					fOutput.close();
				}
			}

			// Write on files
			{
				std::ofstream fOutput(file_name.c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);

				fOutput << "-------------------------------------------------------------" << std::endl;
				fOutput << "Criterion      Tau[s]         T[K]           Value[K/s]      " << std::endl;
				fOutput << "-------------------------------------------------------------" << std::endl;

				for(unsigned int i=0;i<3;i++)
				{
					if (reg_max_.rows() >= i + 1)
					{
						std::stringstream index; index << i + 1;
						std::string label = "T" + index.str() + "(slope)";
						fOutput << std::setw(15) << std::left << label;
						fOutput << std::setw(15) << std::left << reg_max_(i, 1);
						fOutput << std::setw(15) << std::left << reg_max_(i, 2);
						fOutput << std::setw(15) << std::left << reg_max_(i, 3);
						fOutput << std::endl;
					}
				}

				fOutput.close();
			}
		}
	}

	void IgnitionDelayTimes_Analyzer::PrintHeaderLine(unsigned int& counter, std::ostream& fOutput)
	{
		if (is_rcm_ == false)
		{
			if (is_temperature_ == true)
			{
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_T(max)[s]", counter);
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_T(slope)[s]", counter);
			}

			if (is_pressure_ == true)
			{
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_P(max)[s]", counter);
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_P(slope)[s]", counter);
			}

			for (unsigned int i = 0; i < species_index_.size(); i++)
			{
				if (is_species_slope_ == true)
				{
					std::string label = "tau_" + species_name_[i] + "(slope)[s]";
					OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
				}

				{
					std::string label = "tau_" + species_name_[i] + "(max)[s]";
					OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
				}
			}

			for (unsigned int i = 0; i < species_target_abs_index_.size(); i++)
			{
				std::string label = "tau_" + species_target_abs_name_[i] + "(tgt-abs)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			for (unsigned int i = 0; i < species_target_rel_index_.size(); i++)
			{
				std::string label = "tau_" + species_target_rel_name_[i] + "(tgt-rel)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			for (unsigned int i = 0; i < species_target_c_index_.size(); i++)
			{
				std::string label = "tau_" + species_target_c_name_[i] + "(tgt-c)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			for (unsigned int i = 0; i < species_intercept_max_index_.size(); i++)
			{
				std::string label = "tau_" + species_intercept_max_name_[i] + "(int-max)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			for (unsigned int i = 0; i < species_intercept_min_index_.size(); i++)
			{
				std::string label = "tau_" + species_intercept_min_name_[i] + "(int-min)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			if (temperature_increase_ != 0.)
			{
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_T(delta)[s]", counter);
			}

			if (pressure_increase_ != 0.)
			{
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_P(delta)[s]", counter);
			}
		}
		else
		{
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau1_T[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "T1_T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "dTdt1_T[K/s]", counter);
			
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau2_T[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "T2_T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "dTdt2_T[K/s]", counter);
			
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau3_T[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "T3_T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "dTdt3_T[K/s]", counter);
		}
	}

	void IgnitionDelayTimes_Analyzer::Print(std::ostream& fOutput)
	{
		if (is_rcm_ == false)
		{
			if (is_temperature_ == true)
			{
				fOutput << std::setw(25) << std::left << temperature_max_tau_;
				fOutput << std::setw(25) << std::left << temperature_slope_max_tau_;
			}

			if (is_pressure_ == true)
			{
				fOutput << std::setw(25) << std::left << pressure_max_tau_;
				fOutput << std::setw(25) << std::left << pressure_slope_max_tau_;
			}

			for (unsigned int i = 0; i < species_index_.size(); i++)
			{
				if (is_species_slope_ == true)
					fOutput << std::setw(30) << std::left << species_slope_max_tau_[i];
				fOutput << std::setw(30) << std::left << species_max_tau_[i];
			}

			for (unsigned int i = 0; i < species_target_abs_index_.size(); i++)
			{
				fOutput << std::setw(30) << std::left << species_target_abs_tau_[i];
			}

			for (unsigned int i = 0; i < species_target_rel_index_.size(); i++)
			{
				fOutput << std::setw(30) << std::left << species_target_rel_tau_[i];
			}

			for (unsigned int i = 0; i < species_target_c_index_.size(); i++)
			{
				fOutput << std::setw(30) << std::left << species_target_c_tau_[i];
			}

			for (unsigned int i = 0; i < species_intercept_max_index_.size(); i++)
			{
				double tau = 0.;
				if (species_intercept_max_slope_[i] > 0.)
					tau = species_intercept_max_tau_[i] -
					species_intercept_max_x_[i] / species_intercept_max_slope_[i];
				
				fOutput << std::setw(30) << std::left << tau;
			}

			for (unsigned int i = 0; i < species_intercept_min_index_.size(); i++)
			{
				double tau = 0.;
				if (species_intercept_min_slope_[i] > 0.)
					tau = species_intercept_min_tau_[i] -
					species_intercept_min_x_[i] / (-species_intercept_min_slope_[i]);

				fOutput << std::setw(30) << std::left << tau;
			}

			if (temperature_increase_ != 0.)
				fOutput << std::setw(25) << std::left << temperature_increase_tau_;

			if (pressure_increase_ != 0.)
				fOutput << std::setw(25) << std::left << pressure_increase_tau_;
		}
		else
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if (reg_max_.rows() >= i + 1)
				{
					fOutput << std::setw(25) << std::left << reg_max_(i, 1);
					fOutput << std::setw(25) << std::left << reg_max_(i, 2);
					fOutput << std::setw(25) << std::left << reg_max_(i, 3);
				}
				else
				{
					fOutput << std::setw(25) << std::left << 0;
					fOutput << std::setw(25) << std::left << 0;
					fOutput << std::setw(25) << std::left << 0;
				}
			}
		}
	}
}
