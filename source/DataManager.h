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

        void ReadExperimentalData(fs::path experimental_data_files);

    protected:
    private:
        std::string dataset_name;
        std::string solver_name;
        std::string QoI;
        std::string QoI_target;
        bool multiple_input;

        std::vector<double> expdata_x;
        std::vector<double> expdata_y;
        std::vector<double> uncertainty;
        std::string ordinates_label;
        std::string abscissae_label;
        std::string uncertainty_kind;
    };
} // namespace OptiSMOKE

#include "DataManager.hpp"
#endif // DATAMANAGER_H