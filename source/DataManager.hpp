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

namespace OptiSMOKE{

    DataManager::DataManager(){}

    DataManager::~DataManager(){}

    void DataManager::ReadExperimentalData(fs::path experimental_data){
        try
        {
            boost::property_tree::ptree ptree;
            boost::property_tree::read_json(experimental_data.c_str(), ptree);
        
            dataset_name = ptree.get<std::string>("name");
            solver_name = ptree.get<std::string>("solver");
            QoI = ptree.get<std::string>("QoI");
            QoI_target = ptree.get<std::string>("QoI_target");
            multiple_input = ptree.get<bool>("multiple_input");
            
            // ordinates_label = ptree.get<std::string>("data.ordinates_label");
            // abscissae_label = ptree.get<std::string>("data.abscissae_label");
            // uncertainty_kind = ptree.get<std::string>("data.uncertainty_kind");
            
            BOOST_FOREACH(boost::property_tree::ptree::value_type& node, 
                        ptree.get_child("data.abscissae"))
            {
                expdata_x.push_back(node.second.get_value<double>());
            }
            BOOST_FOREACH(boost::property_tree::ptree::value_type& node, 
                        ptree.get_child("data.ordinates"))
            {
                expdata_y.push_back(node.second.get_value<double>());
            }
            BOOST_FOREACH(boost::property_tree::ptree::value_type& node, 
                        ptree.get_child("data.uncertainty"))
            {
                uncertainty.push_back(node.second.get_value<double>());
            }
        }
        catch(const std::exception& e)
        {
            OptiSMOKE::FatalErrorMessage(e.what());
        }        
    }

} // namespace OptiSMOKE