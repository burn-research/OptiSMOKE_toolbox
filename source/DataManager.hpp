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

    void DataManager::ReadExperimentalData(std::vector<std::string> experimental_data_files){
        try
        {
            input_paths_.resize(experimental_data_files.size());
            ordinates_label_.resize(experimental_data_files.size());
            abscissae_label_.resize(experimental_data_files.size());
            uncertainty_kind_.resize(experimental_data_files.size());
        
            expdata_x_.resize(experimental_data_files.size());
            expdata_y_.resize(experimental_data_files.size());
            uncertainty_.resize(experimental_data_files.size());

            for(unsigned int i = 0; i < experimental_data_files.size(); i++){
                // Dont know if it is the best way to do it 
                // however still take this and will see
                boost::property_tree::ptree ptree;
                boost::property_tree::read_json(experimental_data_files[i].c_str(), ptree);
            
                dataset_names_.push_back(ptree.get<std::string>("name"));
                solver_name_.push_back(ptree.get<std::string>("solver"));
                QoI_.push_back(ptree.get<std::string>("QoI"));
                QoI_target_.push_back(ptree.get<std::string>("QoI_target"));
                multiple_input_.push_back(ptree.get<bool>("multiple_input"));

                BOOST_FOREACH(boost::property_tree::ptree::value_type &node, ptree.get_child("OS_Input_File")){
                    input_paths_[i].push_back(node.second.get_value<std::string>());            
                }

                BOOST_FOREACH(boost::property_tree::ptree::value_type &node, ptree.get_child("data")){
                    assert(node.first.empty());
                    abscissae_label_[i].push_back(node.second.get<std::string>("abscissae_label"));
                    ordinates_label_[i].push_back(node.second.get<std::string>("ordinates_label"));
                    uncertainty_kind_[i].push_back(node.second.get<std::string>("uncertainty_kind"));
                    
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("abscissae")){
                        assert(node2.first.empty());
                        expdata_x_[i].push_back(node2.second.get_value<double>());
                    }
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("ordinates")){
                        assert(node2.first.empty());
                        expdata_y_[i].push_back(node2.second.get_value<double>());
                    }
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("uncertainty")){
                        assert(node2.first.empty());
                        uncertainty_[i].push_back(node2.second.get_value<double>());
                    }
                }
            }
        }
        catch(const std::exception& e)
        {
            OptiSMOKE::FatalErrorMessage(e.what());
        }        
    }

} // namespace OptiSMOKE