/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	      |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/json.hpp>

namespace json = boost::json;

namespace OptiSMOKE{
    
    class DataManager
    {
    public:
        /// @brief Default constructor
        DataManager();

        /// @brief Default destructor
        ~DataManager();

        void ReadExperimentalData(std::vector<std::string> experimental_data_files);

        inline const std::vector<std::string>& dataset_names() const {return dataset_names_;};
        inline const std::vector<std::string>& solver_name() const {return solver_name_;};
        inline const std::vector<std::string>& QoI() const {return QoI_;};
        inline const std::vector<std::string>& QoI_target() const {return QoI_target_;};
        inline const std::vector<bool>& multiple_input() const {return multiple_input_;};
        inline const std::vector<std::vector<std::string>>& input_paths() const {return input_paths_;};
        inline const std::vector<std::vector<std::vector<std::string>>>& ordinates_label() const {return ordinates_label_;};
        inline const std::vector<std::vector<std::vector<std::string>>>& abscissae_label() const {return abscissae_label_;};
        inline const std::vector<std::vector<std::vector<std::string>>>& uncertainty_kind() const {return uncertainty_kind_;};
        inline const std::vector<std::vector<std::vector<double>>>& expdata_x() const {return expdata_x_;};
        inline const std::vector<std::vector<std::vector<double>>>& expdata_y() const {return expdata_y_;};
        inline const std::vector<std::vector<std::vector<double>>>& uncertainty() const {return uncertainty_;};
        inline const std::vector<bool>& save_simulations() const {return save_simulations_;};

    private:
        std::vector<std::string> dataset_names_;
        std::vector<std::string> solver_name_;
        std::vector<std::string> QoI_;
        std::vector<std::string> QoI_target_;
        std::vector<bool> multiple_input_;
        std::vector<bool> save_simulations_;

        std::vector<std::vector<std::string>> input_paths_;

        // This blocks here has to go three dimensions
        // since a datasets file can have multiple series
        // dimension one: number of files
        // dimension two: number of datasets whithin each file
        // dimension three: number point in each datasets
        std::vector<std::vector<std::vector<std::string>>> ordinates_label_;
        std::vector<std::vector<std::vector<std::string>>> abscissae_label_;
        std::vector<std::vector<std::vector<std::string>>> uncertainty_kind_;
        std::vector<std::vector<std::vector<double>>> expdata_x_;
        std::vector<std::vector<std::vector<double>>> expdata_y_;
        std::vector<std::vector<std::vector<double>>> uncertainty_;
        std::vector<std::vector<std::vector<double>>> standard_deviations_;

        // Default sigma for standard deviation if it is not present
        // inside the file
        double default_sigma = 2;

        void OrderData();

        void ComputeStandardDeviations();
    };
} // namespace OptiSMOKE

#include "DataManager.hpp"
#endif // DATAMANAGER_H