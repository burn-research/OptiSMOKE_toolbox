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
                save_simulations_.push_back(ptree.get<bool>("save_simulations_data"));

                BOOST_FOREACH(boost::property_tree::ptree::value_type &node, ptree.get_child("OS_Input_File")){
                    assert(node.first.empty());
                    input_paths_[i].push_back(node.second.get_value<std::string>());
                }

                abscissae_label_[i].resize(ptree.get_child("data").size());
                ordinates_label_[i].resize(ptree.get_child("data").size());
                uncertainty_kind_[i].resize(ptree.get_child("data").size());
                expdata_x_[i].resize(ptree.get_child("data").size());
                expdata_y_[i].resize(ptree.get_child("data").size());
                uncertainty_[i].resize(ptree.get_child("data").size());

                unsigned int count = 0;
                BOOST_FOREACH(boost::property_tree::ptree::value_type &node, ptree.get_child("data")){
                    assert(node.first.empty());
                    abscissae_label_[i][count].push_back(node.second.get<std::string>("abscissae_label"));
                    ordinates_label_[i][count].push_back(node.second.get<std::string>("ordinates_label"));
                    
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("abscissae")){
                        assert(node2.first.empty());
                        expdata_x_[i][count].push_back(node2.second.get_value<double>());
                    }
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("ordinates")){
                        assert(node2.first.empty());
                        expdata_y_[i][count].push_back(node2.second.get_value<double>());
                    }
                    
                    boost::optional<std::string> uncertainty_node = node.second.get_optional<std::string>("uncertainty_kind");
                    if(uncertainty_node){
                        uncertainty_kind_[i][count].push_back(node.second.get<std::string>("uncertainty_kind"));
                        BOOST_FOREACH(boost::property_tree::ptree::value_type &node2, node.second.get_child("uncertainty")){
                            assert(node2.first.empty());
                            uncertainty_[i][count].push_back(node2.second.get_value<double>());
                        }
                    }
                    else {
                        // TODO: implement the standard deviation
                        // std::cout << "The uncertainty for the datasets ";
                        // std::cout << experimental_data_files[i].c_str();
                        // std::cout << " is not provided!" << std::endl;
                        uncertainty_kind_[i][count].push_back("relative");
                        for(unsigned int j=0; j < expdata_y_[i][count].size(); j++){
                            uncertainty_[i][count].push_back(0);
                        }

                    }
                    count += 1;
                }
            }
            
            OrderData();
        }
        catch(const std::exception& e)
        {
            OptiSMOKE::FatalErrorMessage(e.what());
        }        
    }

    void DataManager::OrderData(){
        // Number of experimental data files
        unsigned int num =  input_paths_.size(); 
        
        // tmp variables
        std::vector<std::string> dataset_names_tmp;
        std::vector<std::string> solver_name_tmp;
        std::vector<std::string> QoI_tmp;
        std::vector<std::string> QoI_target_tmp;
        std::vector<bool> multiple_input_tmp;
        std::vector<bool> save_simulations_tmp;
        std::vector<std::vector<std::string>> input_paths_tmp;

        std::vector<std::vector<std::vector<std::string>>> ordinates_label_tmp;
        std::vector<std::vector<std::vector<std::string>>> abscissae_label_tmp;
        std::vector<std::vector<std::vector<std::string>>> uncertainty_kind_tmp;
        std::vector<std::vector<std::vector<double>>> expdata_x_tmp;
        std::vector<std::vector<std::vector<double>>> expdata_y_tmp;
        std::vector<std::vector<std::vector<double>>> uncertainty_tmp;

        // Save position of the files
        std::vector<int> batch;
        std::vector<int> pfr;
        std::vector<int> psr;
        std::vector<int> premixed;
        std::vector<int> counterflow;

        for(unsigned int i = 0; i < num; i++){
            if(solver_name_[i] == "BatchReactor")
                batch.push_back(i);
            else if (solver_name_[i] == "PlugFlowReactor")
                pfr.push_back(i);
            else if (solver_name_[i] == "PerfectlyStirredReactor")
                psr.push_back(i);
            else if (solver_name_[i] == "PremixedLaminarFlame1D")
                premixed.push_back(i);
            else if (solver_name_[i] == "CounterFlowFlame1D")
                counterflow.push_back(i);
            else{
                std::string error_str;
                error_str.append("Invalid solver inside datasets ");
                error_str.append(dataset_names_[i]); 
                error_str.append(", available are:\n\t\t\t\t");
                error_str.append("Batchreactor | PlugFlowreactor | PerfectlyStirredReactor | PremixedLaminarFlame1D | CounterFlowFlame1D");
                OptiSMOKE::FatalErrorMessage(error_str);
            }
        }

        // This is code repetition 
        // TODO: fix this!!!
        if(batch.size() != 0){
            for(unsigned int i = 0; i<batch.size(); i++){
                unsigned int pos = batch[i];

                dataset_names_tmp.push_back(dataset_names_[pos]);
                solver_name_tmp.push_back(solver_name_[pos]);
                QoI_tmp.push_back(QoI_[pos]);
                QoI_target_tmp.push_back(QoI_target_[pos]);
                multiple_input_tmp.push_back(multiple_input_[pos]);
                input_paths_tmp.push_back(input_paths_[pos]);
                abscissae_label_tmp.push_back(abscissae_label_[pos]);
                ordinates_label_tmp.push_back(ordinates_label_[pos]);
                uncertainty_kind_tmp.push_back(uncertainty_kind_[pos]);
                expdata_x_tmp.push_back(expdata_x_[pos]);
                expdata_y_tmp.push_back(expdata_y_[pos]);
                uncertainty_tmp.push_back(uncertainty_[pos]);
                save_simulations_tmp.push_back(save_simulations_[pos]);
            }
        }
        
        if(pfr.size() != 0){
            for(unsigned int i = 0; i<pfr.size(); i++){
                unsigned int pos = pfr[i];
                
                dataset_names_tmp.push_back(dataset_names_[pos]);
                solver_name_tmp.push_back(solver_name_[pos]);
                QoI_tmp.push_back(QoI_[pos]);
                QoI_target_tmp.push_back(QoI_target_[pos]);
                multiple_input_tmp.push_back(multiple_input_[pos]);
                input_paths_tmp.push_back(input_paths_[pos]);
                abscissae_label_tmp.push_back(abscissae_label_[pos]);
                ordinates_label_tmp.push_back(ordinates_label_[pos]);
                uncertainty_kind_tmp.push_back(uncertainty_kind_[pos]);
                expdata_x_tmp.push_back(expdata_x_[pos]);
                expdata_y_tmp.push_back(expdata_y_[pos]);
                uncertainty_tmp.push_back(uncertainty_[pos]);
                save_simulations_tmp.push_back(save_simulations_[pos]);
            }
        }
        
        if(psr.size() != 0){
            for(unsigned int i = 0; i<psr.size(); i++){
                unsigned int pos = psr[i];
                
                dataset_names_tmp.push_back(dataset_names_[pos]);
                solver_name_tmp.push_back(solver_name_[pos]);
                QoI_tmp.push_back(QoI_[pos]);
                QoI_target_tmp.push_back(QoI_target_[pos]);
                multiple_input_tmp.push_back(multiple_input_[pos]);
                input_paths_tmp.push_back(input_paths_[pos]);
                abscissae_label_tmp.push_back(abscissae_label_[pos]);
                ordinates_label_tmp.push_back(ordinates_label_[pos]);
                uncertainty_kind_tmp.push_back(uncertainty_kind_[pos]);
                expdata_x_tmp.push_back(expdata_x_[pos]);
                expdata_y_tmp.push_back(expdata_y_[pos]);
                uncertainty_tmp.push_back(uncertainty_[pos]);
                save_simulations_tmp.push_back(save_simulations_[pos]);
            }
        }
        
        if(premixed.size() != 0){
            for(unsigned int i = 0; i<premixed.size(); i++){
                unsigned int pos = premixed[i];
                
                dataset_names_tmp.push_back(dataset_names_[pos]);
                solver_name_tmp.push_back(solver_name_[pos]);
                QoI_tmp.push_back(QoI_[pos]);
                QoI_target_tmp.push_back(QoI_target_[pos]);
                multiple_input_tmp.push_back(multiple_input_[pos]);
                input_paths_tmp.push_back(input_paths_[pos]);
                abscissae_label_tmp.push_back(abscissae_label_[pos]);
                ordinates_label_tmp.push_back(ordinates_label_[pos]);
                uncertainty_kind_tmp.push_back(uncertainty_kind_[pos]);
                expdata_x_tmp.push_back(expdata_x_[pos]);
                expdata_y_tmp.push_back(expdata_y_[pos]);
                uncertainty_tmp.push_back(uncertainty_[pos]);
                save_simulations_tmp.push_back(save_simulations_[pos]);
            }
        }

        if(counterflow.size() != 0){
            for(unsigned int i = 0; i<counterflow.size(); i++){
                unsigned int pos = counterflow[i];
                
                dataset_names_tmp.push_back(dataset_names_[pos]);
                solver_name_tmp.push_back(solver_name_[pos]);
                QoI_tmp.push_back(QoI_[pos]);
                QoI_target_tmp.push_back(QoI_target_[pos]);
                multiple_input_tmp.push_back(multiple_input_[pos]);
                input_paths_tmp.push_back(input_paths_[pos]);
                abscissae_label_tmp.push_back(abscissae_label_[pos]);
                ordinates_label_tmp.push_back(ordinates_label_[pos]);
                uncertainty_kind_tmp.push_back(uncertainty_kind_[pos]);
                expdata_x_tmp.push_back(expdata_x_[pos]);
                expdata_y_tmp.push_back(expdata_y_[pos]);
                uncertainty_tmp.push_back(uncertainty_[pos]);
                save_simulations_tmp.push_back(save_simulations_[pos]);
            }
        }
        
        dataset_names_ = dataset_names_tmp; 
        solver_name_ = solver_name_tmp;
        QoI_ = QoI_tmp;
        QoI_target_ = QoI_target_tmp;
        multiple_input_ = multiple_input_tmp;
        input_paths_ = input_paths_tmp;
        abscissae_label_ = abscissae_label_tmp;
        ordinates_label_ = ordinates_label_tmp;
        uncertainty_kind_ = uncertainty_kind_tmp;
        expdata_x_ = expdata_x_tmp;
        expdata_y_ = expdata_y_tmp;
        uncertainty_ = uncertainty_tmp;
        save_simulations_ = save_simulations_tmp;
    }

    void DataManager::ComputeStandardDeviations(){
        standard_deviations_.resize(expdata_y_.size());
        for(unsigned int i = 0; i < expdata_y_.size(); i++){
            standard_deviations_[i].resize(expdata_y_[i].size());
            for(unsigned int j = 0; j < expdata_y_[i].size(); j++){
                standard_deviations_[i][j].resize(expdata_y_[i][j].size());
                for(unsigned int k = 0; k < expdata_y_[i][j].size(); k++){
                    // standard_deviations_[i][j][k] = uncertainty_[i][j][k] * (expdata_y_[i][j][k]) / Sigma_vector[i];
                }
            }
        }
    }

} // namespace OptiSMOKE