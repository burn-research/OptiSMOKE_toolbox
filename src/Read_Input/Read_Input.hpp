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

namespace OpenSMOKE
{
	void Read_Input::ReadInfo(const char* plugin_input_file, bool print_out)
	{

		// Read data from main dictionary
		std::string input_file_name = plugin_input_file; 
		std::string main_dictionary_name = "OptimizationSetup";
		// Defines the grammar rules
		OpenSMOKE::Grammar_OptiSMOKEpp grammar_optismokepp;
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name).SetGrammar(grammar_optismokepp);
		
		// Read number of parallel threads for OpenMP
		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfThreads") == true)
		{
			dictionaries(main_dictionary_name).ReadInt("@NumberOfThreads", number_of_threads);
			if(print_out)
			{
				std::cout<<""<<std::endl;
				std::cout<<"Solving in parallel using "<<number_of_threads<<" threads."<<std::endl;
				std::cout<<""<<std::endl;
			}
			
		} else
		{
			number_of_threads=1;
		}
		// Read setting for penalty function (default: true)

		Debug_print_out = false;
		if (dictionaries(main_dictionary_name).CheckOption("@Debug") == true)
			dictionaries(main_dictionary_name).ReadBool("@Debug", Debug_print_out);

		Debug_Sim = false;
		if (dictionaries(main_dictionary_name).CheckOption("@DebugSim") == true)
			dictionaries(main_dictionary_name).ReadBool("@DebugSim", Debug_Sim);

		penalty_function = true;
		if (dictionaries(main_dictionary_name).CheckOption("@PenaltyFunction") == true)
			dictionaries(main_dictionary_name).ReadBool("@PenaltyFunction", penalty_function);

		UseBootStrap = false;
		
		if (dictionaries(main_dictionary_name).CheckOption("@UseBootStrap") == true)
                        dictionaries(main_dictionary_name).ReadBool("@UseBootStrap", UseBootStrap);
		CKI_File_Read = false;
		
		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfBootstrapVariations") == true)
                {
                        dictionaries(main_dictionary_name).ReadInt("@NumberOfBootstrapVariations", numberOfBootstrapVariations);
                }
		print_indexes = false;
		if (dictionaries(main_dictionary_name).CheckOption("@PrintIndexes") == true)
                        dictionaries(main_dictionary_name).ReadBool("@PrintIndexes", print_indexes);
		if (dictionaries(main_dictionary_name).CheckOption("@PrintSplines") == true)
                        dictionaries(main_dictionary_name).ReadBool("@PrintSplines", print_splines);
		if (dictionaries(main_dictionary_name).CheckOption("@PrintBootstrap") == true)
                        dictionaries(main_dictionary_name).ReadBool("@PrintBootstrap", print_bootstrap);
		
		// Read Nominal Kinetic scheme
		if (dictionaries(main_dictionary_name).CheckOption("@NominalKineticsFolder") == true)
		{
       			dictionaries(main_dictionary_name).ReadPath("@NominalKineticsFolder", path_nominal_kinetics_output);
       			OpenSMOKE::CheckKineticsFolder(path_nominal_kinetics_output);

       	} else
		{
			std::cout.setstate(std::ios_base::failbit); // Disable video output
			std::string name_of_rapid_kinetics_subdictionary;
			if (dictionaries(main_dictionary_name).CheckOption("@NominalKineticsPreProcessor") == true)
				dictionaries(main_dictionary_name).ReadDictionary("@NominalKineticsPreProcessor", name_of_rapid_kinetics_subdictionary);

			OpenSMOKE::Grammar_RapidKineticMechanism grammar_rapid_kinetics;
			dictionaries(name_of_rapid_kinetics_subdictionary).SetGrammar(grammar_rapid_kinetics);
			//boost::filesystem::path path_input_thermodynamics;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Thermodynamics") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Thermodynamics", path_input_thermodynamics);
		
			//boost::filesystem::path path_input_kinetics;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Kinetics") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Kinetics", path_input_kinetics);

			boost::filesystem::path path_input_transport;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Transport", path_input_transport);
		
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_nominal_kinetics_output);
			
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
			{
				OpenSMOKE::RapidKineticMechanismWithTransport(path_nominal_kinetics_output, path_input_transport.c_str(), path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
			} else
			{
				Read_transport = false;
				OpenSMOKE::RapidKineticMechanismWithoutTransport( path_nominal_kinetics_output, path_input_thermodynamics.c_str(), path_input_kinetics.c_str() );
			}
			std::cout.clear(); // Re-enable video output
			CKI_File_Read = true;
		}
		ReadNominalKinetics();

		// Read Kinetic scheme
		if (dictionaries(main_dictionary_name).CheckOption("@KineticsFolder") == true)
		{
     			dictionaries(main_dictionary_name).ReadPath("@KineticsFolder", path_kinetics_output);
       			OpenSMOKE::CheckKineticsFolder(path_kinetics_output);
       		} else
		{
			std::cout.setstate(std::ios_base::failbit); // Disable video output
			std::string name_of_rapid_kinetics_subdictionary;
			if (dictionaries(main_dictionary_name).CheckOption("@KineticsPreProcessor") == true)
				dictionaries(main_dictionary_name).ReadDictionary("@KineticsPreProcessor", name_of_rapid_kinetics_subdictionary);

			OpenSMOKE::Grammar_RapidKineticMechanism grammar_rapid_kinetics;
			dictionaries(name_of_rapid_kinetics_subdictionary).SetGrammar(grammar_rapid_kinetics);
		
			//boost::filesystem::path path_input_thermodynamics;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Thermodynamics") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Thermodynamics", path_input_thermodynamics);
		
			//boost::filesystem::path path_input_kinetics;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Kinetics") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Kinetics", path_input_kinetics);

			boost::filesystem::path path_input_transport;
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Transport", path_input_transport);
		
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
				dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);
			
			if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
			{
				OpenSMOKE::RapidKineticMechanismWithTransport(path_kinetics_output, path_input_transport.c_str(), path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
			} else
			{
				OpenSMOKE::RapidKineticMechanismWithoutTransport( path_kinetics_output, path_input_thermodynamics.c_str(), path_input_kinetics.c_str() );
			}
			std::cout.clear(); // Re-enable video output
			CKI_File_Read = true;
		}
		//std::cout<<"Before read kinetics"<<std::endl;
		//Fuction for reading kinetics
		ReadKinetics(print_out);

		// Read name of optimized mechanism 
		if (dictionaries(main_dictionary_name).CheckOption("@NameOfOptimizedKineticsFolder") == true)
		{
			dictionaries(main_dictionary_name).ReadString("@NameOfOptimizedKineticsFolder", path_folder);
		} else
		{
			path_folder = "Optimized_kinetics";
		}

		N_of_batch_datasets = 0;
		N_of_plugflow_datasets = 0;
		N_of_psr_datasets = 0;
		N_of_laminar_flame_datasets = 0;

		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfBatchDatasets") == true)
		{
			dictionaries(main_dictionary_name).ReadInt("@NumberOfBatchDatasets", N_of_batch_datasets);
		}
		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfPlugFlowDatasets") == true)
		{
			dictionaries(main_dictionary_name).ReadInt("@NumberOfPlugFlowDatasets", N_of_plugflow_datasets);
		}
		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfPerfectlyStirredDatasets") == true)
		{
			dictionaries(main_dictionary_name).ReadInt("@NumberOfPerfectlyStirredDatasets", N_of_psr_datasets);
		}
		if (dictionaries(main_dictionary_name).CheckOption("@NumberOfLaminarFlameDatasets") == true)
		{
			dictionaries(main_dictionary_name).ReadInt("@NumberOfLaminarFlameDatasets", N_of_laminar_flame_datasets);
		}

		// Read species of intereset
		if (dictionaries(main_dictionary_name).CheckOption("@SpeciesOfInterest") == true)
			dictionaries(main_dictionary_name).ReadOption("@SpeciesOfInterest", species_of_interest);
		
		if (dictionaries(main_dictionary_name).CheckOption("@ParametersBoundaries") == true)
        {
            dictionaries(main_dictionary_name).ReadString("@ParametersBoundaries", boundaries_method);
		}
		// Read Objective function type
		if (dictionaries(main_dictionary_name).CheckOption("@ObjectiveFunctionType") == true)
		{
			dictionaries(main_dictionary_name).ReadString("@ObjectiveFunctionType", Objective_Function);
			if(print_out)
			{
				std::cout<<""<<std::endl;
				std::cout<<"Using type "<<Objective_Function<<" for calculating the objective function."<<std::endl;
				std::cout<<""<<std::endl;
			}
		} else
		{
			if(print_out)
			{
				std::cout<<""<<std::endl;
				std::cout<<"No type specified for calculating the objective function. Using L2-norm."<<std::endl;	
				std::cout<<""<<std::endl;
			}
			Objective_Function = "L2-norm";
		}
		//std::cout<<"@ObjectiveFunctionType is "<<Objective_Function<<std::endl;

		if (dictionaries(main_dictionary_name).CheckOption("@SigmaExpDistribution") == true)
                {
                        dictionaries(main_dictionary_name).ReadInt("@SigmaExpDistribution", SigmaExpDistribution);
                        if(print_out)
                        {
                                std::cout<<"Using RE*y_exp =  "<<SigmaExpDistribution<<"sigma in normal distribution representing experiments."<<std::endl;
                        }
                } else
                {
			SigmaExpDistribution = 2;
                        if(print_out)
                        {
                                std::cout<<"Using RE*y_exp =  2sigma, as default. "<<std::endl;
                        }
                }
		
	
        if (dictionaries(main_dictionary_name).CheckOption("@AcceptedSigmaInKDistribution") == true)
        {
            dictionaries(main_dictionary_name).ReadInt("@AcceptedSigmaInKDistribution", AcceptedSigmaKDistribution);
            if(print_out)
        	{
                std::cout<<"Accept "<<AcceptedSigmaKDistribution<<"sigma in normal distribution representing the kinetic constant."<<std::endl;
            }
        } else {
			AcceptedSigmaKDistribution = 2;
		}
		// AB // Parameters distribution
		if (dictionaries(main_dictionary_name).CheckOption("@Parameters_Distribution") == true)
		{
			dictionaries(main_dictionary_name).ReadString("@Parameters_Distribution", distribution);
		}

		if (dictionaries(main_dictionary_name).CheckOption("@ListOfConstraints") == true)
		{
			udc_bool = true;
			std::string name_constraints_file;
			dictionaries(main_dictionary_name).ReadString("@ListOfConstraints", name_constraints_file);

			std::ifstream myfile(name_constraints_file);

			if (!myfile.is_open())
            		{
                		OpenSMOKE::FatalErrorMessage("Unable to open experimental file: " + name_constraints_file  + " correctly!");
            		}

			double temp_double;
			std::string temp_line;

            		while (std::getline(myfile,temp_line)){
				std::istringstream iss(temp_line);
				std::vector<double> tempv;
				while (iss >> temp_double) {
            				tempv.push_back(temp_double);
            		
        			}
				ud_constraints.push_back(tempv);
			}

		}

		// Read opensmoke input file names
		if (dictionaries(main_dictionary_name).CheckOption("@PathDatasetInputFiles") == true)
		{
			
			std::string path_dataset_input_files;
			// Read Path to files containing paths to input files x dataset
			dictionaries(main_dictionary_name).ReadString("@PathDatasetInputFiles", path_dataset_input_files);
			// vector of string containing all the paths
			std::vector<std::string> path_opensmoke_input_files;
			std::ifstream myfile(path_dataset_input_files);
                        std::string names;
                        while (myfile >> names)
			{
				path_opensmoke_input_files.push_back(names);
			}
			// list_of_opensmoke... is a vect of vect of str ...
			// in one dimensions it will have the size of the number of Datasets
			list_of_opensmoke_input_files.resize(path_opensmoke_input_files.size());
			for (int i=0; i < path_opensmoke_input_files.size(); i++)
			{
				// open the file corresponding to the directory of the input
				std::ifstream myfile_2(path_opensmoke_input_files[i]);
				std::string names_2;
				// populates the two by two vector of strings by pushing back all the names of the inputs corresponding to 
				// a certain dataset within its specific place in the first dimension of the vector.
				while (myfile_2 >> names_2) 
				{
						list_of_opensmoke_input_files[i].push_back(names_2);
				}
			}
		}
		else if (dictionaries(main_dictionary_name).CheckOption("@ListOfOpenSMOKEInputFiles") == true)
		{
			list_of_opensmoke_input_files.resize(1);
			dictionaries(main_dictionary_name).ReadOption("@ListOfOpenSMOKEInputFiles", list_of_opensmoke_input_files[0]);
		}
		else if (dictionaries(main_dictionary_name).CheckOption("@PathOpenSMOKEInputFiles") == true)
		{
		//	list_of_opensmoke_input_files.resize(1);
			std::string path_opensmoke_input_files;
			dictionaries(main_dictionary_name).ReadString("@PathOpenSMOKEInputFiles", path_opensmoke_input_files);

			std::ifstream myfile(path_opensmoke_input_files);
			std::string names;

			std::vector<std::string> a;
			while (myfile >> names) 
			{
				std::vector<std::string> a;
				a.push_back(names);
				//list_of_opensmoke_input_files[counter].push_back(names);
				list_of_opensmoke_input_files.push_back(a);
				
			}
		} else
		{
			OpenSMOKE::FatalErrorMessage("Either @PathDatasetInputFiles or @ListOfOpenSMOKEInputFiles or @PathOpenSMOKEInputFiles must be specified");		
		}	
		// Read experimental data file names
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfExperimentalDataFiles") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfExperimentalDataFiles", experimental_data_files);
		}
		else if (dictionaries(main_dictionary_name).CheckOption("@PathExperimentalDataFiles") == true)
		{
			std::string path_experimental_data_files;
			dictionaries(main_dictionary_name).ReadString("@PathExperimentalDataFiles", path_experimental_data_files);
			
			std::ifstream myfile(path_experimental_data_files);
			std::string names;
			//stores here the experimental_data_files names
			while (myfile >> names) 
			{
				experimental_data_files.push_back(names);
			}
		} else
		{
			OpenSMOKE::FatalErrorMessage("Either @ListOfExperimentalDataFiles or @PathExperimentalDataFiles must be specified");		
		}

		// HERE IT's very clear what it does, basically, it has 
		// Read Quantity of Interest, which can be 
		if (dictionaries(main_dictionary_name).CheckOption("@QuantityOfInterest") == true)
			dictionaries(main_dictionary_name).ReadOption("@QuantityOfInterest", QoI);
		
		// If only one Quantity of Interest has been specified, fill out the vector of QoI with that value.
		// if (QoI.size()==1)
		// {
		// 	QoI.resize(list_of_opensmoke_input_files.size());
		// 	for(int i=0; i< list_of_opensmoke_input_files.size(); i++)
		// 	{
		// 		QoI[i] = QoI[0];
		// 	}
		// }else


		

		// Check file extension of the experimental data files
		file_extension = boost::filesystem::extension(experimental_data_files[0]);
		if(print_out)
		{
			std::cout<<"File extention for exp file: "<<file_extension<<std::endl;
		}
		// HERE ADD -- CM
		// Use reading function according to file extension
		if (file_extension == ".xml")
		{
			ReadExpData_xml();
		} else if (file_extension ==".csv")
		{
			ReadExpData_csv();
		} else
		{		
			ReadExpData();
		}
		
		// Read tau calculation criteria
		if (dictionaries(main_dictionary_name).CheckOption("@CriteriaFortau") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@CriteriaFortau", tau_calc_type);
	
			if (tau_calc_type.size()==1)
			{
				if(print_out)
				{
					std::cout<<""<<std::endl;
					std::cout<<"Calculating ignition delay time based on definition \""<<tau_calc_type[0]<<"\" for all datasets "<<std::endl;
				}
			} else if (tau_calc_type.size()==Exp_data.size())
			{
				if(print_out)
				{
					std::cout<<""<<std::endl;
					for (int i=0; i<tau_calc_type.size(); i++)
					{
						std::cout<<"Calculating ignition delay time based on definition \""<<tau_calc_type[i]<<"\" for dataset "<< i+1 <<std::endl;
					}
				}
			} else
			{
				OpenSMOKE::FatalErrorMessage("Number of ignition delay time definitions [" + std::to_string(tau_calc_type.size()) + "] doesn't correspond to the number of experimental data sets [" + std::to_string(Exp_data.size()) + "]!");
			}
		

			// Create vector of tau calc strings (one for each simulation)
			Order_of_tau_calc.resize(QoI.size());
			if(tau_calc_type.size()>1)
			{
				//for (int i=0; i< list_of_opensmoke_input_files.size(); i++)
				for (int i=0; i< QoI.size(); i++)
				{
					if (QoI[i]=="tau")
					{
						Order_of_tau_calc[i].resize(list_of_opensmoke_input_files[i].size());
						for (int j=0; j < list_of_opensmoke_input_files[i].size(); j++)
						{
							Order_of_tau_calc[i][j] = tau_calc_type[i];
						}
					}
				}
			} else
			{
				for (int i=0; i< QoI.size(); i++)
                	        {
					if (QoI[i]=="tau")
					{
						Order_of_tau_calc[i].resize(list_of_opensmoke_input_files[i].size());
						for (int j=0; j < list_of_opensmoke_input_files[i].size(); j++)
						{
						Order_of_tau_calc[i][j] = tau_calc_type[0];
						}
					}
				}
			}
		}


		// EPLR - Which reactions for Optimization of A, n, Ea

		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_EPLR") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_EPLR", list_of_target_EPLR);

		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_BathGases_EPLR") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_BathGases_EPLR", list_of_bath_gases_EPLR);

		if (dictionaries(main_dictionary_name).CheckOption("@ListOfUncertaintyFactors_EPLR") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfUncertaintyFactors_EPLR", list_of_uncertainty_factors_EPLR);
		// EPLR - Which reactions for Optimization of A, n, Ea
		

		// Extended PLOG - Which reactions for Optimization of A, n, Ea
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions", list_of_target_extplog);

		if (dictionaries(main_dictionary_name).CheckOption("@ListOfUncertaintyFactors_ExtPLOG") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfUncertaintyFactors_ExtPLOG", list_of_uncertainty_factors_extplog);

		// Extended PLOG - Which reactions and which Third Bodies to Optimize
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions_TB") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions_TB", list_of_target_extended_plog_reactions);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Species") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Species", list_of_target_extended_plog_species);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMin_TBeff_ExtPLOG") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMin_TBeff_ExtPLOG", list_of_min_tb_extplog);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMax_TBeff_ExtPLOG") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMax_TBeff_ExtPLOG", list_of_max_tb_extplog);

		// Direct reactions - specific parameters to be optimized	
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_lnA") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_lnA", list_of_target_lnA);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_Beta") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_Beta", list_of_target_Beta);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_E_over_R") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_E_over_R", list_of_target_E_over_R);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_lnA_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_lnA_inf", list_of_target_lnA_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_Beta_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_Beta_inf", list_of_target_Beta_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_E_over_R_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_E_over_R_inf", list_of_target_E_over_R_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Reactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Reactions", list_of_target_thirdbody_reactions);
		
		// Classic PLOG - which reactions
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_classic_PLOG_Reactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_classic_PLOG_Reactions", list_of_target_classic_plog_reactions);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfUncertaintyFactors_classic_PLOG") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfUncertaintyFactors_classic_PLOG", list_of_uncertainty_factors_classic_plog);

		// List of third body species
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Species") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Species", list_of_target_thirdbody_species);
	
		// List of target reactions for uncertainty factors
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors", list_of_target_uncertainty_factors);
	
		// List of uncertainty factors
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfUncertaintyFactors") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfUncertaintyFactors", list_of_uncertainty_factors);
	
		// List of target inf reactions for uncertainty factors
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors_inf", list_of_target_uncertainty_factors_inf);
	
		// List of inf uncertainty factors
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfUncertaintyFactors_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfUncertaintyFactors_inf", list_of_uncertainty_factors_inf);

		// ICA for direct reactions
		if (dictionaries(main_dictionary_name).CheckOption("@ICA_DirectReactionsIndices") == true)
			dictionaries(main_dictionary_name).ReadOption("@ICA_DirectReactionsIndices", direct_reactions_indices_ica);
		// then based on how long they are @ListOfMinAbs_PCA and @ListOfMaxAbs_PCA I will encode a vector of zeros.
		if (dictionaries(main_dictionary_name).CheckOption("@ICA_AbsMin_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ICA_AbsMin_DirectReactions", abs_min_direct_reactions_ica);
		if (dictionaries(main_dictionary_name).CheckOption("@ICA_AbsMax_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ICA_AbsMax_DirectReactions", abs_max_direct_reactions_ica);
		if (dictionaries(main_dictionary_name).CheckOption("@ICA_AbsInit_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@ICA_AbsInit_DirectReactions", abs_init_direct_reactions_ica);
		
		for (int i =0; i < direct_reactions_indices_ica.size(); i++){
				std::cout<< "ICA reaction " << i+1 << " is " << direct_reactions_indices_ica[i]<<std::endl;
		}

		for (int i =0; i < abs_min_direct_reactions_ica.size(); i++){
				std::cout<< "ICA maximum value for variable " << i+1 << " is " << abs_max_direct_reactions_ica[i]   <<std::endl;
				std::cout<< "ICA initial value for variable " << i+1 << " is " << abs_init_direct_reactions_ica[i] <<std::endl;
				std::cout<< "ICA minimum value for variable " << i+1 << " is " << abs_min_direct_reactions_ica[i]  <<std::endl;
		}

		if (dictionaries(main_dictionary_name).CheckOption("@ICA_files_MixingMatrices") == true)
		{
			
			std::string path_ica_files;
			// Read Path to files containing paths to input files x dataset
			dictionaries(main_dictionary_name).ReadString("@ICA_files_MixingMatrices", path_ica_files);
			// Vector of string containing all the paths
			std::ifstream myfile(path_ica_files);

            std::string names;
            while (myfile >> names)
			{
				mixing_matrices_files_ica.push_back(names);
			}

			for (int i=0; i < mixing_matrices_files_ica.size(); i++)
			{
				std::cout<< "The ICA mixing matrix file indexed " << i+1 << " is " << mixing_matrices_files_ica[i] <<std::endl;
			}

			// Resizing the matrices object based on the number of files
			mixing_matrices_ica.resize(mixing_matrices_files_ica.size());

			for (int i=0; i < mixing_matrices_files_ica.size(); i++){
				
				// Open an ifstream for the i-th file and tries to read it, if not there, then a error is displayed
				std::ifstream myfile(mixing_matrices_files_ica[i]);
				if (!myfile.is_open())
				{
					OpenSMOKE::FatalErrorMessage("Unable to open the ica mixing matrix file: " + pca_files_direct_reactions[i] + " correctly! Verify the Path to file!");
				}

				// resize the lot for the i-th mixing matrix to  3 (one for each dimension).
				mixing_matrices_ica[i].resize(3);
				std::string line_2;

				int counter = 0;

				while (std::getline(myfile,line_2))
				{
					// get line by line
					std::stringstream ss_2(line_2);
					std::string token_2;
					// explore the line
					while (ss_2 >> token_2)
					{
						double quantity = boost::lexical_cast<double>(token_2);
						mixing_matrices_ica[i][counter].push_back(quantity);
					}

					counter++;
				}

			}

			// print the object for validation
			for (int i=0; i < mixing_matrices_ica.size(); i++){
				
				for (int j=0; j < mixing_matrices_ica[i][0].size(); j++){
					
					std::cout<< "The 0 ICA component for reaction "<< i << "is "<< mixing_matrices_ica[i][0][j]<<std::endl;
					std::cout<< "The 1 ICA component for reaction "<< i << "is "<< mixing_matrices_ica[i][1][j]<<std::endl;
					std::cout<< "The 2 ICA component for reaction "<< i << "is "<< mixing_matrices_ica[i][2][j]<<std::endl;
				}
			}

		}

		// AB // PCA for direct reactions
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_DirectReactions", pca_direct_reactions);
		// then based on how long they are @ListOfMinAbs_PCA and @ListOfMaxAbs_PCA I will encode a vector of zeros.
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MinAbs_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MinAbs_DirectReactions", pca_minabs_direct_reactions);
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MaxAbs_DirectReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MaxAbs_DirectReactions", pca_maxabs_direct_reactions);

		for (int i =0; i < pca_direct_reactions.size(); i++){
				std::cout<< "PCA reaction " << i+1 << " is " << pca_direct_reactions[i]<<std::endl;
		}

		for (int i =0; i < pca_maxabs_direct_reactions.size(); i++){
				std::cout<< "PCA maximum value for variable " << i+1 << " is " << pca_maxabs_direct_reactions[i]<<std::endl;
				std::cout<< "PCA minimum value for variable " << i+1 << " is " << pca_minabs_direct_reactions[i]<<std::endl;
		}

		if (dictionaries(main_dictionary_name).CheckOption("@PCA_files_DirectReactions") == true)
		{
			
			std::string path_pca_files;
			// Read Path to files containing paths to input files x dataset
			dictionaries(main_dictionary_name).ReadString("@PCA_files_DirectReactions", path_pca_files);
			// vector of string containing all the paths
			std::ifstream myfile(path_pca_files);
            std::string names;
            while (myfile >> names)
			{
				pca_files_direct_reactions.push_back(names);
			}

			for (int i=0; i < pca_files_direct_reactions.size(); i++)
			{
				std::cout<< "The PCA file indexed " << i+1 << " is " << pca_files_direct_reactions[i]<<std::endl;
			}

			pca_scaling_direct_reactions.resize(pca_files_direct_reactions.size());
			pca_eigenvectors_direct_reactions.resize(pca_files_direct_reactions.size());
			pca_centering_direct_reactions.resize(pca_files_direct_reactions.size());

			for (int i=0; i < pca_files_direct_reactions.size(); i++){

				std::ifstream myfile(pca_files_direct_reactions[i]);
				// alarm that rings when the experimental datafile is not there
				if (!myfile.is_open())
				{
					OpenSMOKE::FatalErrorMessage("Unable to open pca file: " + pca_files_direct_reactions[i] + " correctly!");
				}

				pca_eigenvectors_direct_reactions[i].resize(3);

				std::string line_2;
				int counter = 0;
				std::cout<<"It is going to start  reading the  body"<<std::endl;

				number_of_eigenvector_for_reaction.resize(pca_files_direct_reactions.size());
				while (std::getline(myfile,line_2))
				{
					// get line by line
					std::stringstream ss_2(line_2);
					std::string token_2;
					// explore the line
					while (ss_2 >> token_2)
					{
						double quantity = boost::lexical_cast<double>(token_2);

						if(counter == 0){
							pca_scaling_direct_reactions[i].push_back(quantity);
						} else if (counter == 1 || counter == 2 || counter == 3){
							pca_eigenvectors_direct_reactions[i][counter-1].push_back(quantity);

							if (counter == 3 && quantity ==0){
								
								number_of_eigenvector_for_reaction[i] = 2;
								std::cout<<"The number of retained eigenvectors for the reaction "<< i << " is "<<number_of_eigenvector_for_reaction[i]<<std::endl;
							} else {
								number_of_eigenvector_for_reaction[i] = 3;
								std::cout<<"The number of retained eigenvectors for the reaction "<< i << " is "<<number_of_eigenvector_for_reaction[i]<<std::endl;
							}
						} else if (counter == 4){
							pca_centering_direct_reactions[i].push_back(quantity);
						}
					}
					counter++;
				}

			}

			// for (int i=0; i < pca_scaling_direct_reactions.size(); i++){
				
			// 	for (int j=0; j < pca_scaling_direct_reactions[i].size(); j++){
					
			// 		std::cout<< "The PCA scaling factor for reaction " << i << " and variable "<< j << "is " << pca_scaling_direct_reactions[i][j]<<std::endl;
			// 	}
			// }

			// for (int i=0; i < pca_centering_direct_reactions.size(); i++){
				
			// 	for (int j=0; j < pca_centering_direct_reactions[i].size(); j++){
					
			// 		std::cout<< "The PCA centering factor for reaction " << i << " and variable "<< j << "is " << pca_centering_direct_reactions[i][j]<<std::endl;
			// 	}
			// }

			// for (int i=0; i < pca_eigenvectors_direct_reactions.size(); i++){
				
			// 	for (int j=0; j < pca_eigenvectors_direct_reactions[i][0].size(); j++){
					
			// 		std::cout<< "The 0 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_direct_reactions[i][0][j]<<std::endl;
			// 		std::cout<< "The 1 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_direct_reactions[i][1][j]<<std::endl;
			// 		std::cout<< "The 2 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_direct_reactions[i][2][j]<<std::endl;
			// 	}
			// }

		}

		// AB // PCA for inf pressure limit reactions
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_LimitPressureReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_LimitPressureReactions", pca_pinf_reactions);
		// then based on how long they are @ListOfMinAbs_PCA and @ListOfMaxAbs_PCA I will encode a vector of zeros.
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MinAbs_LimitPressureReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MinAbs_LimitPressureReactions", pca_minabs_pinf_reactions);
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MaxAbs_LimitPressureReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MaxAbs_LimitPressureReactions", pca_maxabs_pinf_reactions);

		for (int i =0; i < pca_pinf_reactions.size(); i++){
				std::cout<< "PCA reaction " << i+1 << " is " << pca_pinf_reactions[i]<<std::endl;
		}

		for (int i =0; i < pca_maxabs_pinf_reactions.size(); i++){
				std::cout<< "PCA maximum value for variable " << i+1 << " is " << pca_maxabs_pinf_reactions[i]<<std::endl;
				std::cout<< "PCA minimum value for variable " << i+1 << " is " << pca_minabs_pinf_reactions[i]<<std::endl;
		}

		if (dictionaries(main_dictionary_name).CheckOption("@PCA_files_LimitPressureReactions") == true)
		{
			
			std::string path_pca_files;
			// Read Path to files containing paths to input files x dataset
			dictionaries(main_dictionary_name).ReadString("@PCA_files_LimitPressureReactions", path_pca_files);
			// vector of string containing all the paths
			std::ifstream myfile(path_pca_files);
            std::string names;
            while (myfile >> names)
			{
				pca_files_pinf_reactions.push_back(names);
			}

			for (int i=0; i < pca_files_pinf_reactions.size(); i++)
			{
				std::cout<< "The PCA file indexed " << i+1 << " is " << pca_files_pinf_reactions[i]<<std::endl;
			}

			pca_scaling_pinf_reactions.resize(pca_files_pinf_reactions.size());
			pca_eigenvectors_pinf_reactions.resize(pca_files_pinf_reactions.size());
			pca_centering_pinf_reactions.resize(pca_files_pinf_reactions.size());

			for (int i=0; i < pca_files_pinf_reactions.size(); i++){

				std::ifstream myfile(pca_files_pinf_reactions[i]);
				// alarm that rings when the experimental datafile is not there
				if (!myfile.is_open())
				{
					OpenSMOKE::FatalErrorMessage("Unable to open pca file: " + pca_files_pinf_reactions[i] + " correctly!");
				}

				pca_eigenvectors_pinf_reactions[i].resize(3);

				std::string line_2;
				int counter = 0;
				std::cout<<"It is going to start  reading the  body"<<std::endl;
				while (std::getline(myfile,line_2))
				{
					// get line by line
					std::stringstream ss_2(line_2);
					std::string token_2;
					// explore the line
					while (ss_2 >> token_2)
					{
						double quantity = boost::lexical_cast<double>(token_2);

						if (counter == 0){
							pca_scaling_pinf_reactions[i].push_back(quantity);
						} else if (counter == 1 || counter == 2 || counter == 3){
							pca_eigenvectors_pinf_reactions[i][counter-1].push_back(quantity);	
						} else if (counter == 4){
							pca_centering_pinf_reactions[i].push_back(quantity);
						}
					}
					counter++;
				}

			}

			for (int i=0; i < pca_scaling_pinf_reactions.size(); i++){
				
				for (int j=0; j < pca_scaling_pinf_reactions[i].size(); j++){
					
					std::cout<< "The PCA scaling factor for reaction " << i << " and variable "<< j << "is " << pca_scaling_pinf_reactions[i][j]<<std::endl;
				}
			}

			for (int i=0; i < pca_centering_pinf_reactions.size(); i++){
				
				for (int j=0; j < pca_centering_pinf_reactions[i].size(); j++){
					
					std::cout<< "The PCA centering factor for reaction " << i << " and variable "<< j << "is " << pca_centering_pinf_reactions[i][j]<<std::endl;
				}
			}

			for (int i=0; i < pca_eigenvectors_pinf_reactions.size(); i++){
				
				for (int j=0; j < pca_eigenvectors_pinf_reactions[i][0].size(); j++){
					
					std::cout<< "The 0 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_pinf_reactions[i][0][j]<<std::endl;
					std::cout<< "The 1 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_pinf_reactions[i][1][j]<<std::endl;
					std::cout<< "The 2 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_pinf_reactions[i][2][j]<<std::endl;
				}
			}

		}


		// AB // PCA for PLOG reactions
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_PLOGReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_PLOGReactions", pca_plog_reactions);
		// then based on how long they are @ListOfMinAbs_PCA and @ListOfMaxAbs_PCA I will encode a vector of zeros.
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MinAbs_PLOGReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MinAbs_PLOGReactions", pca_minabs_plog_reactions);
		if (dictionaries(main_dictionary_name).CheckOption("@PCA_MaxAbs_PLOGReactions") == true)
			dictionaries(main_dictionary_name).ReadOption("@PCA_MaxAbs_PLOGReactions", pca_maxabs_plog_reactions);

		for (int i =0; i < pca_plog_reactions.size(); i++){
				std::cout<< "PCA reaction " << i+1 << " is " << pca_plog_reactions[i]<<std::endl;
		}

		for (int i =0; i < pca_maxabs_plog_reactions.size(); i++){
				std::cout<< "PCA maximum value for variable " << i+1 << " is " << pca_maxabs_plog_reactions[i]<<std::endl;
				std::cout<< "PCA minimum value for variable " << i+1 << " is " << pca_minabs_plog_reactions[i]<<std::endl;
		}

		if (dictionaries(main_dictionary_name).CheckOption("@PCA_files_PLOGReactions") == true)
		{
			
			std::string path_pca_files;
			// Read Path to files containing paths to input files x dataset
			dictionaries(main_dictionary_name).ReadString("@PCA_files_PLOGReactions", path_pca_files);
			// vector of string containing all the paths
			std::ifstream myfile(path_pca_files);
            std::string names;
            while (myfile >> names)
			{
				pca_files_plog_reactions.push_back(names);
			}

			for (int i=0; i < pca_files_plog_reactions.size(); i++)
			{
				std::cout<< "The PCA file indexed " << i+1 << " is " << pca_files_plog_reactions[i]<<std::endl;
			}

			pca_scaling_plog_reactions.resize(pca_files_plog_reactions.size());
			pca_eigenvectors_plog_reactions.resize(pca_files_plog_reactions.size());
			pca_centering_plog_reactions.resize(pca_files_plog_reactions.size());

			for (int i=0; i < pca_files_plog_reactions.size(); i++){

				std::ifstream myfile(pca_files_plog_reactions[i]);
				// alarm that rings when the experimental datafile is not there
				if (!myfile.is_open())
				{
					OpenSMOKE::FatalErrorMessage("Unable to open pca file: " + pca_files_plog_reactions[i] + " correctly!");
				}

				pca_eigenvectors_plog_reactions[i].resize(3);

				std::string line_2;
				int counter = 0;
				std::cout<<"It is going to start  reading the  body"<<std::endl;
				while (std::getline(myfile,line_2))
				{
					// get line by line
					std::stringstream ss_2(line_2);
					std::string token_2;
					// explore the line
					while (ss_2 >> token_2)
					{
						double quantity = boost::lexical_cast<double>(token_2);

						if(counter == 0){
							pca_scaling_plog_reactions[i].push_back(quantity);
						} else if (counter == 1 || counter == 2 || counter == 3){
							pca_eigenvectors_plog_reactions[i][counter-1].push_back(quantity);	
						} else if (counter == 4){
							pca_centering_plog_reactions[i].push_back(quantity);
						}
					}
					counter++;
				}

			}

			for (int i=0; i < pca_scaling_plog_reactions.size(); i++){
				
				for (int j=0; j < pca_scaling_plog_reactions[i].size(); j++){
					
					std::cout<< "The PCA scaling factor for reaction " << i << " and variable "<< j << "is " << pca_scaling_plog_reactions[i][j]<<std::endl;
				}
			}

			for (int i=0; i < pca_centering_plog_reactions.size(); i++){
				
				for (int j=0; j < pca_centering_plog_reactions[i].size(); j++){
					
					std::cout<< "The PCA centering factor for reaction " << i << " and variable "<< j << "is " << pca_centering_plog_reactions[i][j]<<std::endl;
				}
			}

			for (int i=0; i < pca_eigenvectors_plog_reactions.size(); i++){
				
				for (int j=0; j < pca_eigenvectors_plog_reactions[i][0].size(); j++){
					
					std::cout<< "The 0 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_plog_reactions[i][0][j]<<std::endl;
					std::cout<< "The 1 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_plog_reactions[i][1][j]<<std::endl;
					std::cout<< "The 2 PCA eigenvector component for reaction "<< i <<" and variable "<< j << "is "<< pca_eigenvectors_plog_reactions[i][2][j]<<std::endl;
				}
			}

		}


		// List of initial parameters (either specified or read from kinetics scheme)
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_lnA") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_lnA", list_of_initial_lnA);
		} else
		{
			for(int i=0; i < list_of_target_lnA.size(); i++)
			{
				list_of_initial_lnA.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML->A(list_of_target_lnA[i]-1))));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_Beta") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_Beta", list_of_initial_Beta);
		} else
		{
			for(int i=0; i < list_of_target_Beta.size(); i++)
			{
				list_of_initial_Beta.push_back(boost::lexical_cast<std::string>(kineticsMapXML->Beta(list_of_target_Beta[i]-1)));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_E_over_R") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_E_over_R", list_of_initial_E_over_R);
		} else
		{
			for(int i=0; i < list_of_target_E_over_R.size(); i++)
			{
				list_of_initial_E_over_R.push_back(boost::lexical_cast<std::string>(kineticsMapXML->E_over_R(list_of_target_E_over_R[i]-1)));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_lnA_inf") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_lnA_inf", list_of_initial_lnA_inf);
		} else
		{
			for(int i=0; i < list_of_target_lnA_inf.size(); i++)
			{
				int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_lnA_inf[i])-indices_of_falloff_reactions.begin();
				list_of_initial_lnA_inf.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction))));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_Beta_inf") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_Beta_inf", list_of_initial_Beta_inf);
		} else
		{
			for(int i=0; i < list_of_target_Beta_inf.size(); i++)
			{
				int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_Beta_inf[i])-indices_of_falloff_reactions.begin();
				list_of_initial_Beta_inf.push_back(boost::lexical_cast<std::string>(kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction)));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_E_over_R_inf") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_E_over_R_inf", list_of_initial_E_over_R_inf);
		} else
		{
			for(int i=0; i < list_of_target_E_over_R_inf.size(); i++)
			{
				int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_E_over_R_inf[i])-indices_of_falloff_reactions.begin();
				list_of_initial_E_over_R_inf.push_back(boost::lexical_cast<std::string>(kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction)));
			}
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfInitial_ThirdBody_Eff") == true)
		{
			dictionaries(main_dictionary_name).ReadOption("@ListOfInitial_ThirdBody_Eff", list_of_initial_thirdbody_eff);
		} else
		{
			for(int i=0; i< list_of_target_thirdbody_reactions.size(); i++)
			{
		        	int iSpecies = thermodynamicsMapXML->IndexOfSpecies(list_of_target_thirdbody_species[i]);
				list_of_initial_thirdbody_eff.push_back(boost::lexical_cast<std::string>(kineticsMapXML->ThirdBody(list_of_target_thirdbody_reactions[i]-1, iSpecies-1)));
			}
		}

		//std::cout << "Le classi iniziano qui!" << std::endl;
		if (dictionaries(main_dictionary_name).CheckOption("@ReactionsClasses") == true)
		{
			dictionaries(main_dictionary_name).ReadBool("@ReactionsClasses", Optimization4Classes);
		}
		if (dictionaries(main_dictionary_name).CheckOption("@ReactionsClassesDefinitions") == true)
		{
			dictionaries(main_dictionary_name).ReadPath("@ReactionsClassesDefinitions", ReactionClassesPath);
			ReadReactionClassesDefinition(ReactionClassesPath);
		}
		//std::cout << "Le classi finiscono qui!" << std::endl;
		// CREATE THE NOMINAL LIST OF EXTENDED PLOG HERE
		// Caclulate limits for the kinetic parameters
		ParameterLimits();

	
		// List of absolute maximum parameters
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_lnA") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_lnA", list_of_max_abs_lnA);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_Beta") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_Beta", list_of_max_abs_Beta);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_E_over_R") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_E_over_R", list_of_max_abs_E_over_R);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_lnA_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_lnA_inf", list_of_max_abs_lnA_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_Beta_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_Beta_inf", list_of_max_abs_Beta_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_E_over_R_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_E_over_R_inf", list_of_max_abs_E_over_R_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxAbs_ThirdBody_Eff") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxAbs_ThirdBody_Eff", list_of_max_abs_thirdbody_eff);
	
		// List of absolute minimum parameters
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_Beta") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_Beta", list_of_min_abs_Beta);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_E_over_R") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_E_over_R", list_of_min_abs_E_over_R);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_lnA_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_lnA_inf", list_of_min_abs_lnA_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_Beta_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_Beta_inf", list_of_min_abs_Beta_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_E_over_R_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_E_over_R_inf", list_of_min_abs_E_over_R_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinAbs_ThirdBody_Eff") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinAbs_ThirdBody_Eff", list_of_min_abs_thirdbody_eff);
	
		// List of relative maximum parameters
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_lnA") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_lnA", list_of_max_rel_lnA);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_Beta") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_Beta", list_of_max_rel_Beta);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_E_over_R") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_E_over_R", list_of_max_rel_E_over_R);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_lnA_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_lnA_inf", list_of_max_rel_lnA_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_Beta_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_Beta_inf", list_of_max_rel_Beta_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_E_over_R_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_E_over_R_inf", list_of_max_rel_E_over_R_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMaxRel_ThirdBody_Eff") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMaxRel_ThirdBody_Eff", list_of_max_rel_thirdbody_eff);
	
		// List of relative minimum parameters
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_lnA") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_lnA", list_of_min_rel_lnA);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_Beta") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_Beta", list_of_min_rel_Beta);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_E_over_R") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_E_over_R", list_of_min_rel_E_over_R);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_lnA_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_lnA_inf", list_of_min_rel_lnA_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_Beta_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_Beta_inf", list_of_min_rel_Beta_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_E_over_R_inf") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_E_over_R_inf", list_of_min_rel_E_over_R_inf);
		if (dictionaries(main_dictionary_name).CheckOption("@ListOfMinRel_ThirdBody_Eff") == true)
			dictionaries(main_dictionary_name).ReadOption("@ListOfMinRel_ThirdBody_Eff", list_of_min_rel_thirdbody_eff);
		// CM_dictionary //

		//  AB  // Initialization of default values for CM  //
		// Order of the basis functions, which refers to the order of the spline that will represent it 
		//which is degree + 1, or also the number of knots that will be used to describe the spline
		int m = 4;

		/* Degree of the basis functions */
		int g = 3;

		bool possibleNegativeOrdinates = true;
		/* Orders of magnitude of difference between the smallest and the largest
		possible value of the smoothing parameter lambda */
		double lambdaSearchInterval = 6;

		/* Number of steps in the for cycle for minimizing the smoothing parameter
		lambda */
		int numberOfStepsLambda = 13;

		/* Fraction of the range of a spline on the y-axis for determining which
		segments of the spline count as asymptotes. If the oscillations of the spline
		at one of its extremities are contained within a horizontal area with size
		determined by this value, the corresponding segment is identified as an
		asymptote */
		double fractionOfOrdinateRangeForAsymptoteIdentification = 0.005;

		/* Fraction of the range of a spline on the y-axis for determining which points
		count as well-defined maxima. In order to be considered a well-defined maximum,
		a point in a spline must not only have first derivative equal to 0 and negative
		second derivative, it must also be sufficiently distant from the two surrounding
		minima. The minimum admissible distance is determined using this variable */
		double fractionOfOrdinateRangeForMaximumIdentification = 0.025;

		/* Fraction of the range of the experimental data on the y-axis for determining
		whether a model can be considered a flat line with ordinate = 0 when compared to
		the experimental data. This is useful for identifying situations that are
		problematic for the calculation of the similarity indexes */
		double fractionOfExpHeightForZeroLineIdentification = 0.02;

		/* Minimum range on the x-axis of the models, before the addition of segments at
		the extremes, required for comparison with the experimental data, expressed as a
		fraction of the range of the experimental data on the x-axis */
		std::vector<double> fractionOfExpRangeForModelExtrapolation(1,0.5);

		/* Fraction of the range of the experimental data on the x-axis used to
		calculate the minimum shift possible */
		double fractionOfExpRangeForMinShift = 0.005;

		/* Fraction of the range of the experimental data on the x-axis used to shift
		the models around the position of perfect alignment between their well-defined
		maximum and the well-defined maximum of the experimental data */
		double fractionOfExpRangeForShiftAroundMaximum = 0.05;

		/* Specifies whether the program should attempt to line up the maxima of
		experimental data and models during the shift procedure */
		bool lineUpMaxima = true;

		/* Specifies whether the program should use the sum of the four dissimilarity
		indexes when calculating the shift, thus finding a single shift amount, or
		whether the program should consider each of the four indexes separately, thus
		obtaining four different values for the shift */
		bool useSumOfIndexesForAlignment = true;

		/* Number of trapezoids for the numerical calculation of the indexes */
		int numberOfTrapezoids = 99;

		/* Default value for the relative experimental error. Used in case relative
		errors are not provided along with the experimental data */
		double defaultRelativeError = 0.1;

		/* Specifies whether the program should use the similarity indexes to choose the
		best set of nodes for the experimental data spline */
		bool useIndexesToChooseExpNodes = false;

		/* Pascal's triangle */
		std::vector<std::vector<double>> pascalsTriangle;

		// TODO // Change HERE to have the CM dictionary //
		if (dictionaries(main_dictionary_name).CheckOption("@CurveMatchingOptions"))
		{
			std::string curve_matching_options;
			dictionaries(main_dictionary_name).ReadDictionary("@CurveMatchingOptions",curve_matching_options);

			// Defines the grammar rules
			OpenSMOKE::Grammar_CurveMatchingOptions grammar_curve_matching_options;
			dictionaries(curve_matching_options).SetGrammar(grammar_curve_matching_options);

			// Read options
			//if (dictionaries(curve_matching_options).CheckOption("@NumberOfBootstrapVariations") == true)
			//	dictionaries(curve_matching_options).ReadInt("@NumberOfBootstrapVariations", numberOfBootstrapVariations);

			if (dictionaries(curve_matching_options).CheckOption("@lineUpMaxima") == true)
				dictionaries(curve_matching_options).ReadBool("@lineUpMaxima", lineUpMaxima);

			if (dictionaries(curve_matching_options).CheckOption("@useSumOfIndexesForAlignment") == true)
				dictionaries(curve_matching_options).ReadBool("@useSumOfIndexesForAlignment", useSumOfIndexesForAlignment);

			if (dictionaries(curve_matching_options).CheckOption("@fractionOfExpRangeForModelExtrapolation") == true)
				dictionaries(curve_matching_options).ReadOption("@fractionOfExpRangeForModelExtrapolation", fractionOfExpRangeForModelExtrapolation);

		}
		// CM_dictionary //

		// Initialization of Dakota options and default values
		// Values for coliny_ea
		method = "coliny_ea";
		string_population_size = "50";
		fitness_type = "merit_function";
		mutation_type = "offset_normal";
		mutation_rate = "1.0";
		crossover_type = "two_point";
		crossover_rate = "0.0";
		replacement_type = "chc = 10";
		// Values for coliny_direct
		string_division = "major_dimension";
		max_boxsize_limit = "0.0";
		min_boxsize_limit = "1.0e-4";
		

		// Read dakota options from input file
		if (dictionaries(main_dictionary_name).CheckOption("@DakotaOptions"))
		{
			std::string dakota_options_dictionary;
			dictionaries(main_dictionary_name).ReadDictionary("@DakotaOptions",dakota_options_dictionary);

			// Defines the grammar rules
			OpenSMOKE::Grammar_DakotaOptions grammar_dakota_options;
			dictionaries(dakota_options_dictionary).SetGrammar(grammar_dakota_options);

			// Read options
			if (dictionaries(dakota_options_dictionary).CheckOption("@TabularDataFile") == true)
				dictionaries(dakota_options_dictionary).ReadString("@TabularDataFile", tabular_data_file);

			if (dictionaries(dakota_options_dictionary).CheckOption("@Method") == true)
				dictionaries(dakota_options_dictionary).ReadString("@Method", method);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MaxIterations") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MaxIterations", string_max_iterations);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MaxFunctionEvaluations") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MaxFunctionEvaluations", string_max_function_evaluations);

			if (dictionaries(dakota_options_dictionary).CheckOption("@ConvergenceTolerance") == true)
				dictionaries(dakota_options_dictionary).ReadString("@ConvergenceTolerance", string_convergence_tolerance);

			if (dictionaries(dakota_options_dictionary).CheckOption("@SolutionTarget") == true)
				dictionaries(dakota_options_dictionary).ReadString("@SolutionTarget", string_solution_target);

			if (dictionaries(dakota_options_dictionary).CheckOption("@Seed") == true)
				dictionaries(dakota_options_dictionary).ReadString("@Seed", string_seed);

			if (dictionaries(dakota_options_dictionary).CheckOption("@PopulationSize") == true)
				dictionaries(dakota_options_dictionary).ReadString("@PopulationSize", string_population_size);

			if (dictionaries(dakota_options_dictionary).CheckOption("@FitnessType") == true)
				dictionaries(dakota_options_dictionary).ReadString("@FitnessType", fitness_type);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MutationType") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MutationType", mutation_type);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MutationRate") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MutationRate", mutation_rate);

			if (dictionaries(dakota_options_dictionary).CheckOption("@CrossoverType") == true)
				dictionaries(dakota_options_dictionary).ReadString("@CrossoverType", crossover_type);

			if (dictionaries(dakota_options_dictionary).CheckOption("@CrossoverRate") == true)
				dictionaries(dakota_options_dictionary).ReadString("@CrossoverRate", crossover_rate);

			if (dictionaries(dakota_options_dictionary).CheckOption("@ReplacementType") == true)
				dictionaries(dakota_options_dictionary).ReadString("@ReplacementType", replacement_type);

			if (dictionaries(dakota_options_dictionary).CheckOption("@Division") == true)
				dictionaries(dakota_options_dictionary).ReadString("@Division", string_division);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MaxBoxsizeLimit") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MaxBoxsizeLimit", max_boxsize_limit);

			if (dictionaries(dakota_options_dictionary).CheckOption("@MinBoxsizeLimit") == true)
				dictionaries(dakota_options_dictionary).ReadString("@MinBoxsizeLimit", min_boxsize_limit);

			if (dictionaries(dakota_options_dictionary).CheckOption("@DiverseInput") == true)
				dictionaries(dakota_options_dictionary).ReadOption("@DiverseInput", diverse_dakota_input);
			
			gradient_option = false;	
 			if (dictionaries(dakota_options_dictionary).CheckOption("@Gradient") == true)
                                dictionaries(dakota_options_dictionary).ReadBool("@Gradient", gradient_option);
		}


		Initialize_the_pre_processor();
	}

	// ******************************************************** //
	//															//
	// 		Create Kinetics Pre-Processor for future use		//
	//															//
	// ******************************************************** //
	// Log file: print an empty log file


	void Read_Input::Initialize_the_pre_processor() {
		
		std::ofstream flog;
		{
			boost::filesystem::path file_name = "log";
			flog.open(file_name.c_str(), std::ios::out);
			flog.setf(std::ios::scientific);
		}

		std::cout.setstate(std::ios_base::failbit); // Disable video output

			//Initialize thermo reader
			typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;
			ThermoReader_CHEMKIN* thermoreader;
			
			//Read Thermo from file
			thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
			thermoreader->ReadFromFile(path_input_thermodynamics.string());

			// Read Kinetics File
			preprocessor_kinetics = new PreProcessorKinetics_CHEMKIN(flog);
			preprocessor_kinetics->ReadFromASCIIFile(path_input_kinetics.string());

			// Defining maps for species and their properties
			typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
			typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN> > PreProcessorSpecies_CHEMKIN_WithoutTransport;
			PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport;
			preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, *preprocessor_kinetics, flog);
			//check if there are errors in thermo/transport and kinetics
			CheckForFatalError(preprocessor_species_without_transport->Setup());
			CheckForFatalError(preprocessor_kinetics->ReadKineticsFromASCIIFile(preprocessor_species_without_transport->AtomicTable()));

		std::cout.clear(); // Re-enable video output

		delete thermoreader;
		delete preprocessor_species_without_transport;
	}
	
	// ******************************************************** //
	//															//
	// 		Create Kinetics Pre-Processor for future use		//
	//															//
	// ******************************************************** //

	// Function for reading Kinetics scheme
	void Read_Input::ReadKinetics(bool print_out)
	{
	
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenInputFileXML(doc, xml_string, path_kinetics_output / "kinetics.xml");
		double tStart_XML = OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new ThermodynamicsMap_CHEMKIN(doc,print_out);
		kineticsMapXML = new KineticsMap_CHEMKIN(*thermodynamicsMapXML,doc,print_out);
		if(!print_out)
		{
			std::cout.setstate(std::ios_base::failbit); // Disable video output
		}
		
		if (Read_transport == true)
			transportMapXML = new TransportPropertiesMap_CHEMKIN(doc);
		if(!print_out)
		{
			std::cout.clear(); // Re-enable video output
		}
		double tEnd_XML = OpenSMOKEGetCpuTime();
		if(print_out)
		{
			std::cout << "Time to read XML file: " << tEnd_XML - tStart_XML << std::endl;
		}

		// Extract third body fall off reaction information
		falloff_indices_of_thirdbody_species = kineticsMapXML->FallOffIndicesOfThirdbodySpecies();
		indices_of_falloff_reactions = kineticsMapXML->IndicesOfFalloffReactions();

		//read indices of extended plog reaction for optimization
		indices_of_extendedplogs_opt = kineticsMapXML->indices_of_extendedpressurelogopt_reactions();
		indices_of_extendedplogs     = kineticsMapXML->indices_of_extendedpressurelog_reactions();
		indices_of_classic_plogs     = kineticsMapXML->indices_of_pressurelog_reactions();

	}
	// Function for reading Nominal Kinetics scheme
	void Read_Input::ReadNominalKinetics()
	{
		std::cout.setstate(std::ios_base::failbit); // Disable video output
		rapidxml::xml_document<> doc_nominal;
		std::vector<char> xml_string_nominal;
		OpenInputFileXML(doc_nominal, xml_string_nominal, path_nominal_kinetics_output / "kinetics.xml");
		nominalthermodynamicsMapXML = new ThermodynamicsMap_CHEMKIN(doc_nominal,false);
		nominalkineticsMapXML = new KineticsMap_CHEMKIN(*nominalthermodynamicsMapXML,doc_nominal,false);

		if (Read_transport == true)
			nominaltransportMapXML = new TransportPropertiesMap_CHEMKIN(doc_nominal);

		nominal_indices_of_extendedplogs_opt = nominalkineticsMapXML->indices_of_extendedpressurelogopt_reactions();
		nominal_indices_of_extendedplogs 	 = nominalkineticsMapXML->indices_of_extendedpressurelog_reactions();
		nominal_indices_of_classic_plogs 	 = nominalkineticsMapXML->indices_of_pressurelog_reactions();
		
		std::cout.clear(); // Re-enable video output

	}
	// Function for reading experimental data from an xml file
	void Read_Input::ReadExpData_xml()
	{
		OpenSMOKE::FatalErrorMessage("Not able to read experimental data stored in xml format yet!");
	}
	// Function for reading experimental data from text files
    void Read_Input::ReadExpData()
    {
		// Resizing here the vectors to fill up while reading experimental data files.
		type_of_reactor.resize(experimental_data_files.size());
		QoI.resize(experimental_data_files.size());
		what_2_calc.resize(experimental_data_files.size());
		Sigma_vector.resize(experimental_data_files.size());
		type_KTP.resize(experimental_data_files.size());
		value_KTP.resize(experimental_data_files.size());

		// resize the first dimension of Expdata to be the same as the number of Datasets
		Exp_data.resize(experimental_data_files.size());
		standard_deviations.resize(experimental_data_files.size());

        // loop over 1D
        for (int i = 0; i < experimental_data_files.size(); i++) {
            // open the file in this buffer input stream
            std::ifstream myfile(experimental_data_files[i]);
            // alarm that rings when the experimental datafile is not there
            if (!myfile.is_open())
            {
                OpenSMOKE::FatalErrorMessage("Unable to open experimental file: " + experimental_data_files[i] + " correctly!");
            }
            
            //std::vector<std::string> Header_strings;
            std::string line;
            std::getline(myfile,line);
            std::stringstream ss(line);
            std::string token;

            int  n_of_target_species = 1;
            int counter = 0;
            while (ss >> token)
            {
				// std::cout<<"Token is " << token <<std::endl;
                if(counter == 0)
				{
                    type_of_reactor[i] = token;
                } 
				else if (counter == 1)
				{    
					QoI[i] = token;

					if (type_of_reactor[i] == "KTP")
					{
						++N_of_direct_measurement_dataset;

					} 
					else 
					{
						type_KTP[i]   = "NULL";
						value_KTP[i]  = 0;
					}
                } 
				else if (counter == 2)
				{
					// here the new keyword is introduced
					if (QoI[i] == "m_SP" || QoI[i] == "m_SP_time" || QoI[i] == "m_SP_out")
					{
						n_of_target_species = std::stoi(token);
						what_2_calc[i].resize(n_of_target_species);
						for (int j = 0; j < n_of_target_species; j++)
						{
							std::string temp_token;
							ss >> temp_token;
							what_2_calc[i][j] = temp_token;
							// std::cout<<"temp token is " << temp_token <<std::endl;
						}
					} 
					else
					{
						what_2_calc[i].resize(1);
						what_2_calc[i][0] = token;
					}
                } 
				else if (counter == 3)
				{
					if (type_of_reactor[i] == "KTP")
					{
						type_KTP[i]   = token;
					} 
					else
					{ 
						Sigma_vector[i] = token;
					}
                    
                } 
				else if (counter == 4)
				{
					if (type_of_reactor[i] == "KTP")
					{
						value_KTP[i]   = std::stod(token);
					}  
                }
                // Header_strings.push_back(token);
                counter = counter +1;
            }
            
			
			if (Sigma_vector[i].empty())
			{
				Sigma_vector[i] = std::to_string(default_sigma);
			}

            std::string line_2;
            std::string token_2;
            // here it will define the internal structure of experimental_data_files[i], which is a vector of vector
            
            std::vector<std::vector<std::string>> abscissa;
            std::vector<std::vector<std::string>> ordinate;
            std::vector<std::vector<std::string>> uncertainty;
            std::vector<std::vector<std::string>> t_eoc;
	    	// loop over until you can get line from myfile
            // int ext_counter = 0;

			abscissa.resize(n_of_target_species);
	    	ordinate.resize(n_of_target_species);
	    	uncertainty.resize(n_of_target_species);
	    	t_eoc.resize(n_of_target_species);
	    
            while (std::getline(myfile,line_2))
            {
				// get line by line
                std::stringstream ss_2(line_2);
                std::string token_2;
                
                int counter = 1;
	    		int counter_2 = 0;

				// explore the line
				while (ss_2 >> token_2)
                {
					if (QoI[i] == "m_SP" || QoI[i] == "m_SP_time" || QoI[i] == "m_SP_out")
					{
						// assign the current 1,2,3 
						if(counter == 1)
						{
							abscissa[counter_2].push_back(token_2);
						} 
						else if (counter == 2)
						{
							ordinate[counter_2].push_back(token_2);	
						} 
						else if (counter == 3)
						{
							uncertainty[counter_2].push_back(token_2);
						}
						// is counter a multiple of 3?
						if (counter % 3 == 0)
						{
							counter = 1;
							counter_2 = counter_2+1;
						} 
						else 
						{	
							//increase by one
							counter = counter+1;
						}
					} 
					else 
					{
						if(counter == 1)
						{
							abscissa[0].push_back(token_2);    
						} 
						else if (counter == 2)
						{
							ordinate[0].push_back(token_2);
						}
						else if (counter == 3){
							uncertainty[0].push_back(token_2); 
						} 
						else if (counter == 4)
						{
							t_eoc[0].push_back(token_2);
						}
						counter = counter +1;
					}
            	}
			}

			// ******************************************************** //
			//															//
			//				COMPUTE STANDARD DEVIATIONS					//
			//															//
			// ******************************************************** //

	    	standard_deviations[i].resize(n_of_target_species);
	    	for (int z = 0; z<n_of_target_species; z++)
	    	{
	    		for (int j = 0; j<abscissa[z].size(); j++)
	    		{	
					standard_deviations[i][z].push_back(std::stod(uncertainty[z][j])*std::stod(ordinate[z][j])/std::stod(Sigma_vector[i]));

					if (standard_deviations[i][z][j]<0)
					{
						OpenSMOKE::FatalErrorMessage("The standard deviation for dataset " + std::to_string(i) + " and point " + std::to_string(j) + "has a negative value!");
					}
	    		}
	    	}

			// ******************************************************** //
			//															//
			//					REPLACE 0s in Std Dev.					//
			//															//
			// ******************************************************** //

			for (int z = 0; z<abscissa.size(); z++) {
				
				std::vector<double> non_zero_standard_deviation;
				for (int j = 0; j<abscissa[z].size(); j++)
				{
					// First find all non-zero values
					if (standard_deviations[i][z][j]!=0)
					{
							non_zero_standard_deviation.push_back(standard_deviations[i][z][j]);
					}
				}

				// Then find minimum value
				double min_non_zero_standard_deveation = *std::min_element(non_zero_standard_deviation.begin(),non_zero_standard_deviation.end());
				
				for (int j = 0; j<abscissa[z].size(); j++) {
					// Then replace zero values with mininum non-zero value
					if (standard_deviations[i][z][j]==0) {
							standard_deviations[i][z][j]=min_non_zero_standard_deveation;
					}

					if(Debug_print_out) {
							std::cout<<"standard_deviations[i][z][j] = "<<standard_deviations[i][z][j]<<" i:"<<i<<" z:"<<z<<" j:"<<j<<std::endl;
					}
				}
			}

			// ******************************************************** //
			//															//
			//				CREATE FINAL DATA STRUCTURE					//
			//															//
			// ******************************************************** //
        	std::vector<std::vector<std::vector<std::string>>> Exp_data_string;
			Exp_data_string.resize(n_of_target_species);

			for (int z = 0; z<n_of_target_species; z++)
			{
				Exp_data_string[z].resize(4);
				for (int l = 0; l < abscissa[z].size(); l++)
				{
					Exp_data_string[z][0].push_back(abscissa[z][l]);
					Exp_data_string[z][1].push_back(ordinate[z][l]);
					Exp_data_string[z][2].push_back(uncertainty[z][l]);
					
					if (!t_eoc[z].empty()) {
						Exp_data_string[z][3].push_back(t_eoc[z][l]);
					} else {
						Exp_data_string[z][3].push_back("0.0");
					}
				}
			}
			// ******************************************************** //
			//				Translate string into double				//
			// ******************************************************** //
        	std::vector<std::vector<std::vector<double>>> Exp_data_temp;
        	Exp_data_temp.resize(Exp_data_string.size());

			for (int z=0; z<n_of_target_species; z++)
			{
					Exp_data_temp[z].resize(Exp_data_string[z].size());
						for (int k = 0; k < Exp_data_string[z].size(); k++)
						{
							for (int l = 0; l < Exp_data_string[z][k].size(); l++)
							{
									Exp_data_temp[z][k].push_back(std::stod(Exp_data_string[z][k][l]));
							}
						}
			}
			
			Exp_data[i].resize(Exp_data_temp.size());
			for (int z=0; z<n_of_target_species; z++)
			{
					Exp_data[i][z].resize(Exp_data_temp[z].size());
						for(int n=0; n < Exp_data_temp[z][0].size(); n++)
						{
							Exp_data[i][z][0].push_back(Exp_data_temp[z][0][n]);

							if (type_of_reactor[i] == "KTP"){
								
								Exp_data[i][z][1].push_back(std::log(Exp_data_temp[z][1][n]));
							} else {

								Exp_data[i][z][1].push_back(Exp_data_temp[z][1][n]);
							}

							
							Exp_data[i][z][2].push_back(Exp_data_temp[z][2][n]);
							Exp_data[i][z][3].push_back(Exp_data_temp[z][3][n]);

							if(Debug_print_out)
							{
								std::cout<<"Exp_data[i][z][0][j] = "<<Exp_data[i][z][0][n]<<" i:"<<i<<" z:"<<z<<" j:"<<n<<std::endl;
								std::cout<<"Exp_data[i][z][1][j] = "<<Exp_data[i][z][1][n]<<" i:"<<i<<" z:"<<z<<" j:"<<n<<std::endl;
								std::cout<<"Exp_data[i][z][2][j] = "<<Exp_data[i][z][2][n]<<" i:"<<i<<" z:"<<z<<" j:"<<n<<std::endl;
								std::cout<<"Exp_data[i][z][3][j] = "<<Exp_data[i][z][3][n]<<" i:"<<i<<" z:"<<z<<" j:"<<n<<std::endl;
							}
						}
			}
        	myfile.close();
    	}
	

    }


	// void Read_Input::BootStrapping_exp_data()
	// {
	// 	std::cout<<"now it's going to start bootstrap"<<std::endl;
	// 	bootstrapExp.resize(Exp_data.size());
	// 	std::cout<<"Number of bootstrap variations = "<<numberOfBootstrapVariations<<std::endl;
    // 		for (int a=0; a<Exp_data.size(); ++a)
	// 	{
	// 		indeed bootstrap variations are initially set to be equal to the experiments
	// 		bootstrapExp[a].resize(numberOfBootstrapVariations);

	// 		for (int b=0; b<numberOfBootstrapVariations; ++b)
	// 		{
    //     			bootstrapExp[a][b].resize(Exp_data[a][1].size());
	// 			for (int c=0; c <Exp_data[a][1].size(); ++c)
	// 			{
	// 				bootstrapExp[a][b][c]=Exp_data[a][1][c];
	// 			}
	// 		}		
	// 	}
    //             for (int k = 0; k < bootstrapExp.size(); k++)
    //             {
    //                    for (int l = 0; l < bootstrapExp[k].size(); l++)
    //                    {
    //                            for (int m = 0; m < bootstrapExp[k][l].size(); m++)
    //                            {
    //                                    std::cout<< "in dataset "<< k << "at row " << l <<"in column " << m <<" the value is "<< bootstrapExp[k][l][m] <<std::endl;
    //                           }
		
    //                    }
    //             }
	// 	std::cout<<"bootstrap object it's been initialized"<<std::endl;

	// 	for (int i=0; i<bootstrapExp.size(); ++i)
	// 	{
	// 	std::cout<< "size "<< bootstrapExp.size() <<std::endl;
    // 		for (int a=0; a<bootstrapExp[i][0].size(); ++a) 
	// 		{
	// 			std::cout<< "size 2 "<< bootstrapExp[i][0].size() <<std::endl;
    //     		Calculates the standard deviation for the single point
	// 			double stdDeviation = Exp_data[i][1][a]*Exp_data[i][2][a];
    //     		create an object of class std::default_random_engine, which generates pseudo-random numbers
    //     		std::default_random_engine generator;
    //     		initializes the seed of the random_engine_generator
    //     		generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    //     		initializes a normal distribution for this experimental point, with mean inputData[1][a] and standard deviation stdDeviation
    //     		std::normal_distribution<double> distribution(Exp_data[i][1][a],stdDeviation);
    //     		loop over the number of Bootstrap variations
    //     		for (int b=1; b<numberOfBootstrapVariations; ++b) 
	// 			{
    //         			generate a random number belonging to the distribution
    //         			double number = distribution(generator);
    //         			std::cout<< "random number "<< number <<std::endl;
	// 					if negative ordinates are not possible, then the minimum value is 0 for the random generated value.
    //        	 			this should only be a safety belt
    //         			if (number < 0 && possibleNegativeOrdinates == false)
	// 					{
	// 			 		std::cout<< "it enters the if" <<std::endl;
	// 			 			number = 0;
    //         			}
	// 			if logScale is false, returns the number
    //         			otherwise, return the log of the number to bootstrap.
	// 			if (logScale == false)
	// 			{
	// 				 std::cout<< "the value before "<< bootstrapExp[i][b][a] <<std::endl;
	// 				 bootstrapExp[i][b][a] = number;
	// 				 std::cout<< "the value after "<< bootstrapExp[i][b][a] <<std::endl;
	// 			}else
	// 			{
	// 				  bootstrapExp[i][b][a] = log(number);
	// 			}	
	// 		}
	// 	}
	// 	}
	// 	print out bootstrap
	// 	for (int k = 0; k < bootstrapExp.size(); k++)
	// 	{
	// 		for (int l = 0; l < bootstrapExp[k].size(); l++)
	// 		{
	// 			for (int m = 0; m < bootstrapExp[k][l].size(); m++)
	// 			{
	// 				std::cout<< "in dataset "<< k << "at row " << l <<"in column " << m <<" the value is "<< bootstrapExp[k][l][m] <<std::endl;
	// 			}
				
	// 		}
	// 	}
	// 	std::cout<<"it finishes bootstrap"<<std::endl;
	// 	splinesExp.resize(Exp_data.size());
	// 	for (int a=0; a<bootstrapExp.size(); ++a)
	// 	{
	// 		splinesExp[a].resize(numberOfBootstrapVariations);
	// 		for (int b=0; b<bootstrapExp[a].size(); ++b)
	// 		{
	// 			std::cout<<"it is going to spline"<<std::endl;
	// 			splinesExp[a][b].solve(Exp_data[a][0],bootstrapExp[a][b],0,0);
	// 			std::cout<<"it is going to remove asymptotes"<<std::endl;
	// 			splinesExp[a][b].removeAsymptotes();
	// 		}
	// 	}
	// 	std::cout<<"it finishes splines"<<std::endl;
	// }
	// Function for reading experimental data from csv files


	void Read_Input::ReadExpData_csv()
	{
		Exp_data.resize(experimental_data_files.size());
		for (int i = 0; i < experimental_data_files.size(); i++)
		{
			std::ifstream myfile(experimental_data_files[i]);
			if(!myfile.is_open())
			{
				OpenSMOKE::FatalErrorMessage("Unable to open experimental file: " + experimental_data_files[i] + " correctly!");	
			}
		        std::vector<std::vector<std::string> > Exp_data_string;
				std::vector<std::string> Header_strings;

        		std::string line;
        		std::getline(myfile,line);
		        std::stringstream ss(line);
		        std::string token;
		        while (std::getline(ss, token, ','))
		        {
                		Header_strings.push_back(token);
		        }

		        std::string line_2;
		        while (std::getline(myfile,line_2))
		        {
                		std::stringstream ss_2(line_2);
		                std::string token_2;
                		std::vector<std::string> Exp_data_row;
		                while (std::getline(ss_2, token_2, ',')) 
                		{
                                	Exp_data_row.push_back(token_2);
		                }
		               	Exp_data_string.push_back(Exp_data_row);
		        } 
			

			std::vector<std::vector<double>> Exp_data_temp;
			Exp_data_temp.resize(Exp_data_string.size());
			std::vector<int> column;
			double factor_unit;
			if (QoI[i]!="tau")
			{
		        	column.resize(species_of_interest.size());
		       		for (int j=0;j<species_of_interest.size();j++)
		      		{
		       	        	 for (int h=0;h < Header_strings.size();h++)
		       		         {
						if(!Header_strings[h].compare(0,species_of_interest[j].size()+2,species_of_interest[j] + "_x"))
			                        {
               		 		                column[j] = h;
               	                			break;
			                        }
               				 }
					if(!column[j])
					{
						OpenSMOKE::FatalErrorMessage("Unable to find species " + species_of_interest[j] + " in file " + experimental_data_files[i] + "!");	
					}
			        }
				for (int j =0; j < column.size(); j++)
				{
					for (int k = 0; k < Exp_data_string.size(); k++)
					{
						try
	                        		{
							Exp_data_temp[k].push_back(std::stod(Exp_data_string[k][column[j]]));
						} catch(...)
						{
							break;
						}	
					}
				}
				factor_unit = 1;
			} else
			{
		       	        for (int h=0;h < Header_strings.size();h++)
		       		{
					if(!Header_strings[h].compare(0,8,"IDT [us]"))
			                {
               		 			column.push_back(h);
						factor_unit = 0.000001;
               	                		break;
					} else if(!Header_strings[h].compare(0,8,"IDT [ms]"))
					{
						column.push_back(h);
						factor_unit = 0.001;
						break;
					} else if(!Header_strings[h].compare(0,7,"IDT [s]"))
					{
						column.push_back(h);
						factor_unit = 1;
						break;
					}
               			}
				for (int j =0; j < column.size(); j++)
				{
					for (int k = 0; k < Exp_data_string.size(); k++)
					{
						try
	                        		{
							Exp_data_temp[k].push_back(std::stod(Exp_data_string[k][column[j]]));
						} catch(...)
						{
							break;
						}
					} 
				}
			}
			int l=0;
			while(l<Exp_data_temp.size())
       			{
   				if(Exp_data_temp[l].empty())
		                {
        	                	Exp_data_temp.erase(Exp_data_temp.begin()+l);
        	       		} else
				{
					l++;
				}
	        	}
			for(int n=0; n < column.size(); n++)
			{
				for (int m=0; m < Exp_data_temp.size(); m++)
				{
					Exp_data[i].emplace_back();
					//Exp_data[i][m].push_back(factor_unit*Exp_data_temp[m][n]);
				}
			}
			myfile.close();
		}
	}

	// Function for calculating kinetic parameter limits
	void Read_Input::ParameterLimits()
	{

		double T_low = 300;
		double T_high = 2500;

		if(boundaries_method == "Furst")
		{
		// Initialize needed values
		std::vector<double> list_of_nominal_lnA_double;
		std::vector<double> list_of_nominal_Beta_double;
		std::vector<double> list_of_nominal_E_over_R_double;

		std::vector<double> list_of_min_abs_lnA_double;
		std::vector<double> list_of_max_abs_lnA_double;
		std::vector<double> list_of_min_abs_Beta_double;
		std::vector<double> list_of_max_abs_Beta_double;
		std::vector<double> list_of_min_abs_E_over_R_double;
		std::vector<double> list_of_max_abs_E_over_R_double;


		std::vector<double> kappa_lower_T_low;
		std::vector<double> kappa_upper_T_low;
		std::vector<double> kappa_lower_T_high;
		std::vector<double> kappa_upper_T_high;
		std::vector<double> Beta_1;
		std::vector<double> Beta_2;
		std::vector<double> lnA_1;
		std::vector<double> lnA_2;
		std::vector<double> E_over_R_1;
		std::vector<double> E_over_R_2;
		// Resize them
		list_of_nominal_lnA_double.resize(list_of_target_uncertainty_factors.size());
		list_of_nominal_Beta_double.resize(list_of_target_uncertainty_factors.size());
		list_of_nominal_E_over_R_double.resize(list_of_target_uncertainty_factors.size());

		list_of_min_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());
		list_of_max_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());
		list_of_min_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
		list_of_max_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
		list_of_min_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
		list_of_max_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());

		kappa_lower_T_low.resize(list_of_target_uncertainty_factors.size());
		kappa_upper_T_low.resize(list_of_target_uncertainty_factors.size());
		kappa_lower_T_high.resize(list_of_target_uncertainty_factors.size());
		kappa_upper_T_high.resize(list_of_target_uncertainty_factors.size());
		Beta_1.resize(list_of_target_uncertainty_factors.size());
		Beta_2.resize(list_of_target_uncertainty_factors.size());
		lnA_1.resize(list_of_target_uncertainty_factors.size());
		lnA_2.resize(list_of_target_uncertainty_factors.size());
		E_over_R_1.resize(list_of_target_uncertainty_factors.size());
		E_over_R_2.resize(list_of_target_uncertainty_factors.size());
		for (int i=0; i < list_of_target_uncertainty_factors.size(); i++)
		{
			// Nominal values of parameters
			list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML->A(list_of_target_uncertainty_factors[i]-1));
			list_of_nominal_Beta_double[i] = nominalkineticsMapXML->Beta(list_of_target_uncertainty_factors[i]-1);
			list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML->E_over_R(list_of_target_uncertainty_factors[i]-1);
			// Min and Max of lnA
			list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors[i]));
			list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i]));
			if (std::find(list_of_target_lnA.begin(),list_of_target_lnA.end(),list_of_target_uncertainty_factors[i]) != list_of_target_lnA.end())
			{
				list_of_min_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
				list_of_max_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
			}
			
			// Limiting values for the rate coefficient
			kappa_lower_T_low[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
			kappa_upper_T_low[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
			kappa_lower_T_high[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);
			kappa_upper_T_high[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);
			//if (!list_of_nominal_Beta_double[i]==0)
			//{
				// Calculating extreme values for Beta
				Beta_1[i] = (kappa_upper_T_low[i] - kappa_lower_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
				Beta_2[i] = (kappa_lower_T_low[i] - kappa_upper_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
				//Beta_1[i] = list_of_nominal_Beta_double[i]*(1+list_of_uncertainty_factors[i]*0.1);
				//Beta_2[i] = list_of_nominal_Beta_double[i]/(1+list_of_uncertainty_factors[i]*0.1);
				list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
				list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);
				if (std::find(list_of_target_Beta.begin(),list_of_target_Beta.end(),list_of_target_uncertainty_factors[i]) != list_of_target_Beta.end())
				{
					list_of_min_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
					list_of_max_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
				}
			//}
			//if (!list_of_nominal_E_over_R_double[i]==0)
			//{
				// Calculting extreame values of E_over_R
				lnA_1[i] = ( kappa_lower_T_high[i] - (T_low/T_high) * kappa_upper_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_1[i] = lnA_1[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_upper_T_low[i] * T_low;
				lnA_2[i] = ( kappa_upper_T_high[i] - (T_low/T_high) * kappa_lower_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_2[i] = lnA_2[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_lower_T_low[i] * T_low;
				list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
				list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);
				if (std::find(list_of_target_E_over_R.begin(),list_of_target_E_over_R.end(),list_of_target_uncertainty_factors[i]) != list_of_target_E_over_R.end())
				{
					list_of_min_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
					list_of_max_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
				}
			//}
		}

		// Initialize needed values
		std::vector<double> list_of_nominal_lnA_inf_double;
		std::vector<double> list_of_nominal_Beta_inf_double;
		std::vector<double> list_of_nominal_E_over_R_inf_double;

		std::vector<double> list_of_min_abs_lnA_inf_double;
		std::vector<double> list_of_max_abs_lnA_inf_double;
		std::vector<double> list_of_min_abs_Beta_inf_double;
		std::vector<double> list_of_max_abs_Beta_inf_double;
		std::vector<double> list_of_min_abs_E_over_R_inf_double;
		std::vector<double> list_of_max_abs_E_over_R_inf_double;

		std::vector<double> kappa_lower_T_low_inf;
		std::vector<double> kappa_upper_T_low_inf;
		std::vector<double> kappa_lower_T_high_inf;
		std::vector<double> kappa_upper_T_high_inf;
		std::vector<double> Beta_1_inf;
		std::vector<double> Beta_2_inf;
		std::vector<double> lnA_1_inf;
		std::vector<double> lnA_2_inf;
		std::vector<double> E_over_R_1_inf;
		std::vector<double> E_over_R_2_inf;

		// Resize them
		list_of_nominal_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_nominal_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_nominal_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

		list_of_min_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_min_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_min_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

		kappa_lower_T_low_inf.resize(list_of_target_uncertainty_factors_inf.size());
		kappa_upper_T_low_inf.resize(list_of_target_uncertainty_factors_inf.size());
		kappa_lower_T_high_inf.resize(list_of_target_uncertainty_factors_inf.size());
		kappa_upper_T_high_inf.resize(list_of_target_uncertainty_factors_inf.size());
		Beta_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		Beta_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
		lnA_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		lnA_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
		E_over_R_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		E_over_R_2_inf.resize(list_of_target_uncertainty_factors_inf.size());

		for (int i=0; i < list_of_target_uncertainty_factors_inf.size(); i++)
		{
			int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_uncertainty_factors_inf[i])-indices_of_falloff_reactions.begin();
			// Nominal values of inf parameters
			list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML->A_falloff_inf(pos_FallOff_Reaction));
			list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction);
			list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction);

			// Min and Max of lnA_inf
			list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors_inf[i]));
			list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,list_of_uncertainty_factors_inf[i]));
			if (std::find(list_of_target_lnA_inf.begin(),list_of_target_lnA_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_lnA_inf.end())
			{
				list_of_min_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
				list_of_max_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
			}
			
			// Limiting values for the rate coefficient
			kappa_lower_T_low_inf[i]  = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
			kappa_upper_T_low_inf[i]  = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
			kappa_lower_T_high_inf[i] = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);
			kappa_upper_T_high_inf[i] = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);

			//if (!list_of_nominal_Beta_inf_double[i]==0)
			//{
				// Calculating extreme values for Beta
				Beta_1_inf[i] = (kappa_upper_T_low_inf[i] - kappa_lower_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
				Beta_2_inf[i] = (kappa_lower_T_low_inf[i] - kappa_upper_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
				list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
				list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
				if (std::find(list_of_target_Beta_inf.begin(),list_of_target_Beta_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_Beta_inf.end())
				{
					list_of_min_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
					list_of_max_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
				}
			//}
			//if (!list_of_nominal_E_over_R_inf_double[i]==0)
			//{
				// Calculting extreame values of E_over_R
				lnA_1_inf[i] = ( kappa_lower_T_high_inf[i] - (T_low/T_high) * kappa_upper_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_1_inf[i] = lnA_1_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_upper_T_low_inf[i] * T_low;
				lnA_2_inf[i] = ( kappa_upper_T_high_inf[i] - (T_low/T_high) * kappa_lower_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_2_inf[i] = lnA_2_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_lower_T_low_inf[i] * T_low;
				list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
				list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
				if (std::find(list_of_target_E_over_R_inf.begin(),list_of_target_E_over_R_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_E_over_R_inf.end())
				{
					list_of_min_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
					list_of_max_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
				}
		
			//}
		
		}
		}
		if(boundaries_method == "Narrow")
		{
			std::vector<double> list_of_nominal_lnA_double;
			std::vector<double> list_of_nominal_Beta_double;
			std::vector<double> list_of_nominal_E_over_R_double;

			std::vector<double> list_of_min_abs_lnA_double;
			std::vector<double> list_of_max_abs_lnA_double;
			std::vector<double> list_of_min_abs_Beta_double;
			std::vector<double> list_of_max_abs_Beta_double;
			std::vector<double> list_of_min_abs_E_over_R_double;
			std::vector<double> list_of_max_abs_E_over_R_double;


			
			std::vector<double> Beta_1;
			std::vector<double> Beta_2;
			std::vector<double> lnA_1;
			std::vector<double> lnA_2;
			std::vector<double> E_over_R_1;
			std::vector<double> E_over_R_2;
			
			list_of_nominal_lnA_double.resize(list_of_target_uncertainty_factors.size());
			list_of_nominal_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_nominal_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
		
			list_of_min_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());
			list_of_min_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_min_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
		
			Beta_1.resize(list_of_target_uncertainty_factors.size());
			Beta_2.resize(list_of_target_uncertainty_factors.size());
			lnA_1.resize(list_of_target_uncertainty_factors.size());
			lnA_2.resize(list_of_target_uncertainty_factors.size());
			E_over_R_1.resize(list_of_target_uncertainty_factors.size());
			E_over_R_2.resize(list_of_target_uncertainty_factors.size());
			
			for (int i=0; i < list_of_target_uncertainty_factors.size(); i++)
			{
				list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML->A(list_of_target_uncertainty_factors[i]-1));
				list_of_nominal_Beta_double[i] = nominalkineticsMapXML->Beta(list_of_target_uncertainty_factors[i]-1);
				list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML->E_over_R(list_of_target_uncertainty_factors[i]-1);
				

				list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors[i]));
				list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i]));
				
				if (std::find(list_of_target_lnA.begin(),list_of_target_lnA.end(),list_of_target_uncertainty_factors[i]) != list_of_target_lnA.end())
				{	
					list_of_min_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
					list_of_max_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
				}

				Beta_1[i] = list_of_nominal_Beta_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				Beta_2[i] = list_of_nominal_Beta_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				
				list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
				list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);

				if (std::find(list_of_target_Beta.begin(),list_of_target_Beta.end(),list_of_target_uncertainty_factors[i]) != list_of_target_Beta.end())
				{
					list_of_min_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
					list_of_max_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
				}

				E_over_R_1[i] = list_of_nominal_E_over_R_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				E_over_R_2[i] = list_of_nominal_E_over_R_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				
				list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
				list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);
				
				if (std::find(list_of_target_E_over_R.begin(),list_of_target_E_over_R.end(),list_of_target_uncertainty_factors[i]) != list_of_target_E_over_R.end())
				{
					list_of_min_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
					list_of_max_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
				}
			}


		std::vector<double> list_of_nominal_lnA_inf_double;
		std::vector<double> list_of_nominal_Beta_inf_double;
		std::vector<double> list_of_nominal_E_over_R_inf_double;

		std::vector<double> list_of_min_abs_lnA_inf_double;
		std::vector<double> list_of_max_abs_lnA_inf_double;
		std::vector<double> list_of_min_abs_Beta_inf_double;
		std::vector<double> list_of_max_abs_Beta_inf_double;
		std::vector<double> list_of_min_abs_E_over_R_inf_double;
		std::vector<double> list_of_max_abs_E_over_R_inf_double;

		std::vector<double> Beta_1_inf;
		std::vector<double> Beta_2_inf;
		std::vector<double> lnA_1_inf;
		std::vector<double> lnA_2_inf;
		std::vector<double> E_over_R_1_inf;
		std::vector<double> E_over_R_2_inf;

		// Resize them
		list_of_nominal_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_nominal_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_nominal_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

		list_of_min_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_min_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_min_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
		list_of_max_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

		Beta_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		Beta_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
		lnA_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		lnA_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
		E_over_R_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
		E_over_R_2_inf.resize(list_of_target_uncertainty_factors_inf.size());

		for (int i=0; i < list_of_target_uncertainty_factors_inf.size(); i++)
		{
			int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_uncertainty_factors_inf[i])-indices_of_falloff_reactions.begin();
			// Nominal values of inf parameters
			list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML->A_falloff_inf(pos_FallOff_Reaction));
			list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction);
			list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction);

			// Min and Max of lnA_inf
			list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors_inf[i]));
			list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,list_of_uncertainty_factors_inf[i]));
			if (std::find(list_of_target_lnA_inf.begin(),list_of_target_lnA_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_lnA_inf.end())
			{
				list_of_min_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
				list_of_max_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
			}
			
			// Calculating extreme values for Beta
			Beta_1_inf[i] = list_of_nominal_Beta_inf_double[i]+std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) / std::log(T_high);
			Beta_2_inf[i] = list_of_nominal_Beta_inf_double[i]-std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) / std::log(T_high);;
			list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
			list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
			if (std::find(list_of_target_Beta_inf.begin(),list_of_target_Beta_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_Beta_inf.end())
			{
					list_of_min_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
					list_of_max_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
			}

			// Calculting extreame values of E_over_R
			E_over_R_1_inf[i] = list_of_nominal_E_over_R_inf_double[i]-std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) * T_low;
			E_over_R_2_inf[i] = list_of_nominal_E_over_R_inf_double[i]+std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) * T_low;
			list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
			list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
			if (std::find(list_of_target_E_over_R_inf.begin(),list_of_target_E_over_R_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_E_over_R_inf.end())
			{
				list_of_min_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
				list_of_max_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
			}
	
		}

		}

		if(boundaries_method == "Re-parametrization")
		{
			//here i added the A_scaled parameter
			std::cout<<"We are now withing Re-parametrization mode" << std::endl;
			std::vector<double> list_of_nominal_A_scaled_double;
			std::vector<double> list_of_min_abs_A_scaled_double;
			std::vector<double> list_of_max_abs_A_scaled_double;
			list_of_nominal_A_scaled_double.resize(list_of_target_uncertainty_factors.size());
			list_of_min_abs_A_scaled_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_A_scaled_double.resize(list_of_target_uncertainty_factors.size());
			std::cout<<"Everything has been created and re-sized correctly for A_scaled" << std::endl;

			std::vector<double> list_of_nominal_lnA_double;
			std::vector<double> list_of_nominal_Beta_double;
			std::vector<double> list_of_nominal_E_over_R_double;

			std::vector<double> list_of_min_abs_lnA_double;
			std::vector<double> list_of_max_abs_lnA_double;
			std::vector<double> list_of_min_abs_Beta_double;
			std::vector<double> list_of_max_abs_Beta_double;
			std::vector<double> list_of_min_abs_E_over_R_double;
			std::vector<double> list_of_max_abs_E_over_R_double;

			std::vector<double> Beta_1;
			std::vector<double> Beta_2;
			std::vector<double> lnA_1;
			std::vector<double> lnA_2;			
			std::vector<double> E_over_R_1;
			std::vector<double> E_over_R_2;
			
			list_of_nominal_lnA_double.resize(list_of_target_uncertainty_factors.size());
			list_of_nominal_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_nominal_E_over_R_double.resize(list_of_target_uncertainty_factors.size());

			list_of_min_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_lnA_double.resize(list_of_target_uncertainty_factors.size());

			list_of_min_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_Beta_double.resize(list_of_target_uncertainty_factors.size());
			list_of_min_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
			list_of_max_abs_E_over_R_double.resize(list_of_target_uncertainty_factors.size());
		
			Beta_1.resize(list_of_target_uncertainty_factors.size());
			Beta_2.resize(list_of_target_uncertainty_factors.size());
			lnA_1.resize(list_of_target_uncertainty_factors.size());
			lnA_2.resize(list_of_target_uncertainty_factors.size());
			E_over_R_1.resize(list_of_target_uncertainty_factors.size());
			E_over_R_2.resize(list_of_target_uncertainty_factors.size());
			
			for (int i=0; i < list_of_target_uncertainty_factors.size(); i++)
			{
				list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML->A(list_of_target_uncertainty_factors[i]-1));
				list_of_nominal_Beta_double[i] = nominalkineticsMapXML->Beta(list_of_target_uncertainty_factors[i]-1);
				list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML->E_over_R(list_of_target_uncertainty_factors[i]-1);

				//here i compute A_scaled_nom
				list_of_nominal_A_scaled_double[i] = std::exp(list_of_nominal_lnA_double[i]) * std::pow(1000,list_of_nominal_Beta_double[i]) * std::exp(-1*list_of_nominal_E_over_R_double[i]/1000);
				std::cout<<"The nominal value of A_scaled for reaction " << i <<"th is "<< list_of_nominal_A_scaled_double[i]<< std::endl;			

				list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors[i]));
				list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i]));
				
				if (std::find(list_of_target_lnA.begin(),list_of_target_lnA.end(),list_of_target_uncertainty_factors[i]) != list_of_target_lnA.end())
				{	
					list_of_min_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
					list_of_max_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
				}

				//Here I compute the Maximum/Minimum of A_scaled
				list_of_min_abs_A_scaled_double[i] = std::exp(list_of_min_abs_lnA_double[i]) * std::pow(1000,list_of_nominal_Beta_double[i]) * std::exp(-1*list_of_nominal_E_over_R_double[i]/1000);
				std::cout<<"The min value of A_scaled for reaction " << i <<"th is "<< list_of_min_abs_A_scaled_double[i]<< std::endl;	
				list_of_max_abs_A_scaled_double[i] = std::exp(list_of_max_abs_lnA_double[i]) * std::pow(1000,list_of_nominal_Beta_double[i]) * std::exp(-1*list_of_nominal_E_over_R_double[i]/1000);
				std::cout<<"The max value of A_scaled for reaction " << i <<"th is "<< list_of_max_abs_A_scaled_double[i]<< std::endl;	
				
				if (std::find(list_of_target_lnA.begin(),list_of_target_lnA.end(),list_of_target_uncertainty_factors[i]) != list_of_target_lnA.end())
				{	
					list_of_nom_abs_A_scaled.push_back(boost::lexical_cast<std::string>(list_of_nominal_A_scaled_double[i]));
					list_of_min_abs_A_scaled.push_back(boost::lexical_cast<std::string>(list_of_min_abs_A_scaled_double[i]));
					list_of_max_abs_A_scaled.push_back(boost::lexical_cast<std::string>(list_of_max_abs_A_scaled_double[i]));
				}

				Beta_1[i] = list_of_nominal_Beta_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				Beta_2[i] = list_of_nominal_Beta_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				
				list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
				list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);

				if (std::find(list_of_target_Beta.begin(),list_of_target_Beta.end(),list_of_target_uncertainty_factors[i]) != list_of_target_Beta.end())
				{
					list_of_min_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
					list_of_max_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
				}

				E_over_R_1[i] = list_of_nominal_E_over_R_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				E_over_R_2[i] = list_of_nominal_E_over_R_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				
				list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
				list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);
				
				if (std::find(list_of_target_E_over_R.begin(),list_of_target_E_over_R.end(),list_of_target_uncertainty_factors[i]) != list_of_target_E_over_R.end())
				{
					list_of_min_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
					list_of_max_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
				}
			}
		}


		// EPLR - Alpha, Beta, Eps 
		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{						
			list_of_nominal_lnA_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_EPLR.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_EPLR[i]));
			list_of_max_lnA_EPLR.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_EPLR[i]));
		}

		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{					
			list_of_nominal_ER_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_ER_EPLR.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i]))*T_low));
			list_of_max_ER_EPLR.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i]))*T_low));
		}

		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{
			list_of_nominal_Beta_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_EPLR.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i])) / std::log(T_high)));
			list_of_max_Beta_EPLR.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i])) / std::log(T_high)));
		}
		
		// EPLR - Alpha, Beta, Eps 

		// Extended PLOG - Alpha, Beta, Eps 
		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{						
			list_of_nominal_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_extplog[i]));
			list_of_max_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_extplog[i]));
		}

		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{					
			list_of_nominal_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_extplog[i]))*T_low));
			list_of_max_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_extplog[i]))*T_low));
		}

		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{
			list_of_nominal_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_extplog[i])) / std::log(T_high)));
			list_of_max_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_extplog[i])) / std::log(T_high)));
		}

		// Extended PLOGS - THIRD BODIES!
		int pos_extended_plog_reaction;
		int pos_extended_plog_species;
		for (int i=0; i < list_of_target_extended_plog_reactions.size(); i++ )
		{
			// find the position of the reaction having index list_of_target_extended_plog_reactions[i] within indices_of_extendedplogs
			pos_extended_plog_reaction = std::find(indices_of_extendedplogs_opt.begin(),indices_of_extendedplogs_opt.end(),list_of_target_extended_plog_reactions[i])-indices_of_extendedplogs_opt.begin();
			for (int k=0; k < nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species().size(); k++ )
			{		
				std::string temp_string = nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species()[k];
				if(!nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species()[k].compare(0,temp_string.length(),list_of_target_extended_plog_species[i]))
			    {
               		pos_extended_plog_species = k;
               	        break;
			        
               	}
			}
			// get the value of third body			
			list_of_nominal_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).ThirdBody(pos_extended_plog_species)[0]));
			list_of_min_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(list_of_min_tb_extplog[i]));
			list_of_max_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(list_of_max_tb_extplog[i]));
		}


 		// CLASSIC PLOG - Alpha, Beta, Eps
		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{						
			list_of_nominal_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_classic_plog[i]));
			list_of_max_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_classic_plog[i]));
		}

		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{			
			// the nominal random variable is 0, so that Eps_0 = Esp_0 + D is verified:		
			list_of_nominal_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			// The minimum and the maximum values of the random variable are then computed as follows:
			list_of_min_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i]))*T_low));
			list_of_max_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i]))*T_low));
		}

		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{
			list_of_nominal_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i])) / std::log(T_high)));
			list_of_max_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i])) / std::log(T_high)));
		}
	
		// std::cout << "Questa è la fine di parameters limit" << std::endl;
	}

	// Function for preparing Dakota input in string format
	void Read_Input::DakotaInputString()
	{	
		std::string param_name_string;
		std::string initial_values_string;
		std::string lower_bounds_string;
		std::string upper_bounds_string;
		std::string std_deviations_string;

		// if i live it like this direct reactions will always be needed to start, which is fine so far.
		if (pca_direct_reactions.size() > 0) {
								
								// add direct reactions
								int cnt = 0;
								int pc_idx = 1;

								for (int i=0; i<pca_maxabs_direct_reactions.size(); i++){
									
									std::string name_temp = "'PC" + std::to_string(pc_idx) + "_R" + std::to_string(pca_direct_reactions[cnt]) + "'";
									param_name_string+= name_temp + " ";

									pc_idx++;

									double rmd = (i+1)%number_of_eigenvector_for_reaction[0];

									if (rmd ==0){
										cnt++;
										pc_idx = 1;
									}

									initial_values_string += " 0 ";
									lower_bounds_string   += boost::lexical_cast<std::string>(pca_minabs_direct_reactions[i]) + " ";
									upper_bounds_string   += boost::lexical_cast<std::string>(pca_maxabs_direct_reactions[i]) + " ";
									std_deviations_string += " 0 ";
								}

								number_of_parameters = pca_maxabs_direct_reactions.size();

								// add P_inf reactions
								cnt = 0;
								pc_idx = 1;

								for (int i=0; i<pca_maxabs_pinf_reactions.size(); i++){
									
									std::string name_temp = "'PC" + std::to_string(pc_idx) + "_R" + std::to_string(pca_pinf_reactions[cnt]) + "_inf'";
									param_name_string+= name_temp + " ";

									pc_idx++;

									double rmd = (i+1)%3;

									if (rmd ==0){
										cnt++;
										pc_idx = 1;
									}

									initial_values_string += " 0 ";
									lower_bounds_string   += boost::lexical_cast<std::string>(pca_minabs_pinf_reactions[i]) + " ";
									upper_bounds_string   += boost::lexical_cast<std::string>(pca_maxabs_pinf_reactions[i]) + " ";
									std_deviations_string += " 0 ";
								}

								number_of_parameters = number_of_parameters + pca_maxabs_pinf_reactions.size();

								// add PLOG reactions
								cnt = 0;
								pc_idx = 1;

								for (int i=0; i<pca_maxabs_plog_reactions.size(); i++){
									
									std::string name_temp = "'PC" + std::to_string(pc_idx) + "_PLOG" + std::to_string(pca_plog_reactions[cnt]) + "'";
									param_name_string+= name_temp + " ";

									pc_idx++;

									double rmd = (i+1)%3;

									if (rmd ==0){
										cnt++;
										pc_idx = 1;
									}

									initial_values_string += " 0 ";
									lower_bounds_string   += boost::lexical_cast<std::string>(pca_minabs_plog_reactions[i]) + " ";
									upper_bounds_string   += boost::lexical_cast<std::string>(pca_maxabs_plog_reactions[i]) + " ";
									std_deviations_string += " 0 ";
								}

								// add extended PLOG parameters for Third Bodies
								std::vector<std::string> name_vec_ExtPlog;
								name_vec_ExtPlog.resize(list_of_target_extended_plog_reactions.size());
								// for NOW it does not make sense to me to give the possibility to set minimum and maximum values, as they are to many parameters.
								// But this is something i probably need to do
								for (int i=0; i< list_of_target_extended_plog_reactions.size(); i++)
								{
									name_vec_ExtPlog[i] = "'ExtPLOG_" + std::to_string(list_of_target_extended_plog_reactions[i]) + "_SP_" + list_of_target_extended_plog_species[i] + "'";
									param_name_string   += name_vec_ExtPlog[i] + " "; 

									//filling up the strings 
									initial_values_string += list_of_nominal_TB_ExtPLOG[i] + " ";
									lower_bounds_string   += list_of_min_TB_ExtPLOG[i]     + " ";
									upper_bounds_string   += list_of_max_TB_ExtPLOG[i]     + " ";
								}


								number_of_parameters = number_of_parameters + pca_maxabs_plog_reactions.size() + list_of_target_extended_plog_reactions.size();
		
		} else if (direct_reactions_indices_ica.size() > 0) {
			
			// add direct reactions
			int cnt = 0;
			int ic_idx = 1;

			for (int i=0; i<abs_max_direct_reactions_ica.size(); i++){
				
				std::string name_temp = "'IC" + std::to_string(ic_idx) + "_R" + std::to_string(direct_reactions_indices_ica[cnt]) + "'";
				param_name_string+= name_temp + " ";

				ic_idx++;

				double rmd = (i+1)%3;

				if (rmd ==0){
					cnt++;
					ic_idx = 1;
				}

				initial_values_string += boost::lexical_cast<std::string>(abs_init_direct_reactions_ica[i]) + " ";
				lower_bounds_string   += boost::lexical_cast<std::string>(abs_min_direct_reactions_ica[i])  + " ";
				upper_bounds_string   += boost::lexical_cast<std::string>(abs_max_direct_reactions_ica[i])  + " ";
				std_deviations_string += " 0 ";
			}

			number_of_parameters = abs_max_direct_reactions_ica.size();

		}else {

					if(boundaries_method == "Re-parametrization"){
						
								std::vector<std::string> name_vec_A_scaled;
								name_vec_A_scaled.resize(list_of_target_lnA.size());
								for (int i=0; i< list_of_target_lnA.size(); i++)
								{
									name_vec_A_scaled[i] = "'A_scaled_R" + std::to_string(list_of_target_lnA[i]) + "'";
									param_name_string+= name_vec_A_scaled[i] + " ";
									//list_of_nominal_A_scaled_double[i] 
									initial_values_string+= list_of_nom_abs_A_scaled[i] + " ";
									if (list_of_min_rel_lnA.size()>0)
									{
										//lower_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A(list_of_target_lnA[i]-1)))+std::log(list_of_min_rel_lnA[i])) + " ";
									} else
									{
										lower_bounds_string+= list_of_min_abs_A_scaled[i] + " ";
										// AB - Filling up the standard deviation string
										std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_nom_abs_A_scaled[i]) - std::stod(list_of_min_abs_A_scaled[i]))/3) + " ";
									}
									if (list_of_max_rel_lnA.size()>0)
									{
										upper_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A(list_of_target_lnA[i]-1)))+std::log(list_of_max_rel_lnA[i])) + " ";
									} else
									{
										upper_bounds_string+= list_of_max_abs_A_scaled[i] + " ";
									}
								}

					} else {
								std::vector<std::string> name_vec_lnA;
								name_vec_lnA.resize(list_of_target_lnA.size());
								for (int i=0; i< list_of_target_lnA.size(); i++)
								{
									name_vec_lnA[i] = "'lnA_R" + std::to_string(list_of_target_lnA[i]) + "'";
									param_name_string+= name_vec_lnA[i] + " ";
									initial_values_string+= list_of_initial_lnA[i] + " ";
									if (list_of_min_rel_lnA.size()>0)
									{
										std::cout << "Forse entro qui!" << std::endl;
										lower_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A(list_of_target_lnA[i]-1)))+std::log(list_of_min_rel_lnA[i])) + " ";
									} else
									{
										lower_bounds_string+= list_of_min_abs_lnA[i] + " ";
										// AB - Filling up the standard deviation string
										std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_lnA[i]) - std::stod(list_of_min_abs_lnA[i]))/3) + " ";
									}
									if (list_of_max_rel_lnA.size()>0)
									{
										upper_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A(list_of_target_lnA[i]-1)))+std::log(list_of_max_rel_lnA[i])) + " ";
									} else
									{
										upper_bounds_string+= list_of_max_abs_lnA[i] + " ";
									}
								}
								// std::cout << "Sto ciclo for non lo passo!" << std::endl;
					}

					std::vector<std::string> name_vec_lnA_inf;
					name_vec_lnA_inf.resize(list_of_target_lnA_inf.size());
					for (int i=0; i< list_of_target_lnA_inf.size(); i++)
					{
						name_vec_lnA_inf[i] = "'lnA_R" + std::to_string(list_of_target_lnA_inf[i]) + "_inf'";	
						param_name_string+=name_vec_lnA_inf[i] + " ";
						initial_values_string+= list_of_initial_lnA_inf[i] + " ";
						if (list_of_min_rel_lnA_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_lnA_inf[i])-indices_of_falloff_reactions.begin();
							lower_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction)))+std::log(list_of_min_rel_lnA_inf[i])) + " ";
						} else
						{
							lower_bounds_string+= list_of_min_abs_lnA_inf[i] + " ";
							// AB - Filling up the standard deviation string
							std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_lnA_inf[i]) - std::stod(list_of_min_abs_lnA_inf[i]))/3) + " ";
						}
						if (list_of_max_rel_lnA_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_lnA_inf[i])-indices_of_falloff_reactions.begin();
							upper_bounds_string+= boost::lexical_cast<std::string>((std::log(kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction)))+std::log(list_of_max_rel_lnA_inf[i])) + " ";
						} else
						{
							upper_bounds_string+= list_of_max_abs_lnA_inf[i] + " ";
						}
					}
				
					std::vector<std::string> name_vec_Beta;
					name_vec_Beta.resize(list_of_target_Beta.size());
					for (int i=0; i< list_of_target_Beta.size(); i++)
					{
						name_vec_Beta[i] = "'Beta_R" + std::to_string(list_of_target_Beta[i]) + "'";	
						param_name_string+=name_vec_Beta[i] + " ";
						initial_values_string+= list_of_initial_Beta[i] + " ";
						if (list_of_min_rel_Beta.size()>0)
						{
							lower_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->Beta(list_of_target_Beta[i]-1))*list_of_min_rel_Beta[i]) + " ";
						} else
						{
							lower_bounds_string+= list_of_min_abs_Beta[i] + " ";
							std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_Beta[i]) - std::stod(list_of_min_abs_Beta[i]))/3) + " ";
						}
						if (list_of_max_rel_Beta.size()>0)
						{
							upper_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->Beta(list_of_target_Beta[i]-1))*list_of_max_rel_Beta[i]) + " ";
						} else
						{
							upper_bounds_string+= list_of_max_abs_Beta[i] + " ";
						}
					}
					std::vector<std::string> name_vec_Beta_inf;
					name_vec_Beta_inf.resize(list_of_target_Beta_inf.size());
					for (int i=0; i<list_of_target_Beta_inf.size(); i++)
					{
						name_vec_Beta_inf[i] = "'Beta_R" + std::to_string(list_of_target_Beta_inf[i]) + "_inf'";	
						param_name_string+=name_vec_Beta_inf[i] + " ";
						initial_values_string+= list_of_initial_Beta_inf[i] + " ";
						if (list_of_min_rel_Beta_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_Beta_inf[i])-indices_of_falloff_reactions.begin();
							lower_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction))*list_of_min_rel_Beta_inf[i]) + " ";
						} else
						{
							lower_bounds_string+= list_of_min_abs_Beta_inf[i] + " ";
							std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_Beta_inf[i]) - std::stod(list_of_min_abs_Beta_inf[i]))/3) + " ";
						}
						if (list_of_max_rel_Beta_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_Beta_inf[i])-indices_of_falloff_reactions.begin();
							upper_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction))*list_of_max_rel_Beta_inf[i]) + " ";
						} else
						{
							upper_bounds_string+= list_of_max_abs_Beta_inf[i] + " ";
						}
					}
				
					std::vector<std::string> name_vec_E_over_R;
					name_vec_E_over_R.resize(list_of_target_E_over_R.size());
					for (int i=0; i< list_of_target_E_over_R.size(); i++)
					{
						name_vec_E_over_R[i] = "'E_over_R_R" + std::to_string(list_of_target_E_over_R[i]) + "'";
						param_name_string+=name_vec_E_over_R[i] + " ";
						initial_values_string+= list_of_initial_E_over_R[i] + " ";
						if (list_of_min_rel_E_over_R.size()>0)
						{
							lower_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->E_over_R(list_of_target_E_over_R[i]-1))*list_of_min_rel_E_over_R[i]) + " ";
						} else
						{
							lower_bounds_string+= list_of_min_abs_E_over_R[i] + " ";
							std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_E_over_R[i]) - std::stod(list_of_min_abs_E_over_R[i]))/3) + " ";
						}
						if (list_of_max_rel_E_over_R.size()>0)
						{
							upper_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->E_over_R(list_of_target_E_over_R[i]-1))*list_of_max_rel_E_over_R[i]) + " ";
						} else
						{
							upper_bounds_string+= list_of_max_abs_E_over_R[i] + " ";
						}
					}
				
					std::vector<std::string> name_vec_E_over_R_inf;
					name_vec_E_over_R_inf.resize(list_of_target_E_over_R_inf.size());
					for (int i=0; i< list_of_target_E_over_R_inf.size(); i++)
					{
						name_vec_E_over_R_inf[i] = "'E_over_R_R" + std::to_string(list_of_target_E_over_R_inf[i]) + "_inf'";	
						param_name_string+=name_vec_E_over_R_inf[i] + " ";
						initial_values_string+= list_of_initial_E_over_R_inf[i] + " ";
						if (list_of_min_rel_E_over_R_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_E_over_R_inf[i])-indices_of_falloff_reactions.begin();
							lower_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction))*list_of_min_rel_E_over_R_inf[i]) + " ";
						} else
						{
							lower_bounds_string+= list_of_min_abs_E_over_R_inf[i] + " ";
							std_deviations_string+= boost::lexical_cast<std::string>((std::stod(list_of_initial_E_over_R_inf[i]) - std::stod(list_of_min_abs_E_over_R_inf[i]))/3) + " ";
						}
						if (list_of_max_rel_E_over_R_inf.size()>0)
						{
							int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_E_over_R_inf[i])-indices_of_falloff_reactions.begin();
							upper_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction))*list_of_max_rel_E_over_R_inf[i]) + " ";
						} else
						{
							upper_bounds_string+= list_of_max_abs_E_over_R_inf[i] + " ";
						}
					}

					// third body efficiencies!
					std::vector<std::string> name_vec_thirdbody;
					name_vec_thirdbody.resize(list_of_target_thirdbody_reactions.size());
					for (int i=0; i< list_of_target_thirdbody_reactions.size(); i++)
					{
						name_vec_thirdbody[i] = "'M_R" + std::to_string(list_of_target_thirdbody_reactions[i]) + "_" + list_of_target_thirdbody_species[i] + "'";
						param_name_string+=name_vec_thirdbody[i] + " ";
						initial_values_string+= list_of_initial_thirdbody_eff[i] + " ";
						if (list_of_min_abs_thirdbody_eff.size()>0)
						{
							lower_bounds_string+= list_of_min_abs_thirdbody_eff[i] + " ";
						} else
						{
							int iSpecies = thermodynamicsMapXML->IndexOfSpecies(list_of_target_thirdbody_species[i]);
							lower_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->ThirdBody(list_of_target_thirdbody_reactions[i]-1, iSpecies-1))*list_of_min_rel_thirdbody_eff[i]) + " ";
							//std_deviations_string+= boost::lexical_cast<std::string>((boost::lexical_cast<std::double>(list_of_initial_E_over_R_inf[i]) - boost::lexical_cast<std::double>(list_of_min_abs_E_over_R_inf[i]))/3) + " ";
						}
						if (list_of_max_abs_thirdbody_eff.size()>0)
						{
							upper_bounds_string+= list_of_max_abs_thirdbody_eff[i] + " ";
						} else
						{
								int iSpecies = thermodynamicsMapXML->IndexOfSpecies(list_of_target_thirdbody_species[i]);
							upper_bounds_string+= boost::lexical_cast<std::string>((kineticsMapXML->ThirdBody(list_of_target_thirdbody_reactions[i]-1, iSpecies-1))*list_of_max_rel_thirdbody_eff[i]) + " ";
						}
					}

					// CLASSIC PLOG REACTIONS
					std::vector<std::string> name_vec_lnA_classic_plog;
					name_vec_lnA_classic_plog.resize(list_of_target_classic_plog_reactions.size());

					for (int i=0; i< list_of_target_classic_plog_reactions.size(); i++)
					{
						name_vec_lnA_classic_plog[i] = "'lnA_classic_PLOG_" + std::to_string(list_of_target_classic_plog_reactions[i]) + "'";
						param_name_string   += name_vec_lnA_classic_plog[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_lnA_classic_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_lnA_classic_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_lnA_classic_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>(list_of_uncertainty_factors_classic_plog[i]/3) + " ";

					}

					std::vector<std::string> name_vec_ER_classic_plog;
					name_vec_ER_classic_plog.resize(list_of_target_classic_plog_reactions.size());

					for (int i=0; i< list_of_target_classic_plog_reactions.size(); i++)
					{
						name_vec_ER_classic_plog[i] = "'E_over_R_classic_PLOG_" + std::to_string(list_of_target_classic_plog_reactions[i]) + "'";
						param_name_string   += name_vec_ER_classic_plog[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_ER_classic_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_ER_classic_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_ER_classic_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_ER_classic_plog_coefficients[i]) - std::stod(list_of_min_ER_classic_plog_coefficients[i]))/3) + " ";
						
					}
					
					std::vector<std::string> name_vec_Beta_classic_plog;
					name_vec_Beta_classic_plog.resize(list_of_target_classic_plog_reactions.size());

					for (int i=0; i< list_of_target_classic_plog_reactions.size(); i++)
					{
						name_vec_Beta_classic_plog[i] = "'Beta_classic_PLOG_" + std::to_string(list_of_target_classic_plog_reactions[i]) + "'";
						param_name_string   += name_vec_Beta_classic_plog[i] + " ";

						initial_values_string += list_of_nominal_Beta_classic_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_Beta_classic_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_Beta_classic_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_Beta_classic_plog_coefficients[i]) - std::stod(list_of_min_Beta_classic_plog_coefficients[i]))/3) + " ";
					}

					// EPLR
					std::vector<std::string> name_vec_lnA_EPLR;
					name_vec_lnA_EPLR.resize(list_of_target_EPLR.size());

					for (int i=0; i< list_of_target_EPLR.size(); i++)
					{
						name_vec_lnA_EPLR[i] = "'lnA_EPLR_R" + std::to_string(list_of_target_EPLR[i]) + "_bath_" + list_of_bath_gases_EPLR[i] + "'";
						param_name_string   += name_vec_lnA_EPLR[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_lnA_EPLR[i] + " ";
						lower_bounds_string   += list_of_min_lnA_EPLR[i]     + " ";
						upper_bounds_string   += list_of_max_lnA_EPLR[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>(list_of_uncertainty_factors_EPLR[i]/3) + " ";

					}

					std::vector<std::string> name_vec_ER_EPLR;
					name_vec_ER_EPLR.resize(list_of_target_EPLR.size());

					for (int i=0; i< list_of_target_EPLR.size(); i++)
					{
						name_vec_ER_EPLR[i] = "'E_over_R_EPLR_" + std::to_string(list_of_target_EPLR[i])  + "_bath_" + list_of_bath_gases_EPLR[i] + "'";
						param_name_string   += name_vec_ER_EPLR[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_ER_EPLR[i] + " ";
						lower_bounds_string   += list_of_min_ER_EPLR[i]     + " ";
						upper_bounds_string   += list_of_max_ER_EPLR[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_ER_EPLR[i]) - std::stod(list_of_min_ER_EPLR[i]))/3) + " ";
						
					}
					
					std::vector<std::string> name_vec_Beta_EPLR;
					name_vec_Beta_EPLR.resize(list_of_target_EPLR.size());

					for (int i=0; i< list_of_target_EPLR.size(); i++)
					{
						name_vec_Beta_EPLR[i] = "'Beta_EPLR_" + std::to_string(list_of_target_EPLR[i]) + "_bath_" + list_of_bath_gases_EPLR[i] + "'";
						param_name_string   += name_vec_Beta_EPLR[i] + " ";

						initial_values_string += list_of_nominal_Beta_EPLR[i] + " ";
						lower_bounds_string   += list_of_min_Beta_EPLR[i]     + " ";
						upper_bounds_string   += list_of_max_Beta_EPLR[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_Beta_EPLR[i]) - std::stod(list_of_min_Beta_EPLR[i]))/3) + " ";
					}

					// EXTENDED PLOG REACTIONS
					std::vector<std::string> name_vec_lnA_EXT_plog;
					name_vec_lnA_EXT_plog.resize(list_of_target_extplog.size());

					for (int i=0; i< list_of_target_extplog.size(); i++)
					{
						name_vec_lnA_EXT_plog[i] = "'lnA_ExtPLOG_" + std::to_string(list_of_target_extplog[i]) + "'";
						param_name_string   += name_vec_lnA_EXT_plog[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_lnA_ext_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_lnA_ext_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_lnA_ext_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>(list_of_uncertainty_factors_extplog[i]/3) + " ";

					}

					std::vector<std::string> name_vec_ER_EXT_plog;
					name_vec_ER_EXT_plog.resize(list_of_target_extplog.size());

					for (int i=0; i< list_of_target_extplog.size(); i++)
					{
						name_vec_ER_EXT_plog[i] = "'E_over_R_ExtPLOG_" + std::to_string(list_of_target_extplog[i]) + "'";
						param_name_string   += name_vec_ER_EXT_plog[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_ER_ext_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_ER_ext_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_ER_ext_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_ER_ext_plog_coefficients[i]) - std::stod(list_of_min_ER_ext_plog_coefficients[i]))/3) + " ";
						
					}
					
					std::vector<std::string> name_vec_Beta_EXT_plog;
					name_vec_Beta_EXT_plog.resize(list_of_target_extplog.size());

					for (int i=0; i< list_of_target_extplog.size(); i++)
					{
						name_vec_Beta_EXT_plog[i] = "'Beta_ExtPLOG_" + std::to_string(list_of_target_extplog[i]) + "'";
						param_name_string   += name_vec_Beta_EXT_plog[i] + " ";

						initial_values_string += list_of_nominal_Beta_ext_plog_coefficients[i] + " ";
						lower_bounds_string   += list_of_min_Beta_ext_plog_coefficients[i]     + " ";
						upper_bounds_string   += list_of_max_Beta_ext_plog_coefficients[i]     + " ";
						std_deviations_string += boost::lexical_cast<std::string>((std::stod(list_of_nominal_Beta_ext_plog_coefficients[i]) - std::stod(list_of_min_Beta_ext_plog_coefficients[i]))/3) + " ";
					}

					// add extended PLOG parameters for Third Bodies
					std::vector<std::string> name_vec_ExtPlog;
					name_vec_ExtPlog.resize(list_of_target_extended_plog_reactions.size());
					// for NOW it does not make sense to me to give the possibility to set minimum and maximum values, as they are to many parameters.
					// But this is something i probably need to do
					for (int i=0; i< list_of_target_extended_plog_reactions.size(); i++)
					{
						name_vec_ExtPlog[i] = "'ExtPLOG_" + std::to_string(list_of_target_extended_plog_reactions[i]) + "_SP_" + list_of_target_extended_plog_species[i] + "'";
						param_name_string   += name_vec_ExtPlog[i] + " "; 

						//filling up the strings 
						initial_values_string += list_of_nominal_TB_ExtPLOG[i] + " ";
						lower_bounds_string   += list_of_min_TB_ExtPLOG[i]     + " ";
						upper_bounds_string   += list_of_max_TB_ExtPLOG[i]     + " ";
					}

		// Calculate number of uncertain parameters.
		// AB // ADDED also lnA plog extended * 2.
		// Still missing place for Activation Energy.
		number_of_parameters = list_of_target_lnA.size() + list_of_target_Beta.size() + list_of_target_E_over_R.size() + list_of_target_lnA_inf.size() + list_of_target_Beta_inf.size() + list_of_target_E_over_R_inf.size() + list_of_target_thirdbody_reactions.size() + list_of_target_extended_plog_reactions.size() + list_of_target_extplog.size()*3 + list_of_target_classic_plog_reactions.size()*3 +list_of_target_EPLR.size()*3;

		}					
		//std::cout << "\nCiao tito questa è la roba che hai!" << std::endl;
		//std::cout << "Num param: " << number_of_parameters << std::endl;
		//std::cout << "Descriptors: " << param_name_string << std::endl;
		//std::cout << "initial_point: " << initial_values_string << std::endl;
		//std::cout << "LB: " << lower_bounds_string << std::endl;
		//std::cout << "UB: " << upper_bounds_string << std::endl;
		//
		  if(!tabular_data_file.empty())
		  {
			dakota_options_string = "     environment,"
		  	"\n      tabular_data";
		 	dakota_options_string.append("\n 		tabular_data_file '" + tabular_data_file + "'");
		  }
			dakota_options_string.append("\n	method,"); 
			dakota_options_string.append("\n 		" + method);
		  if(!string_max_iterations.empty())
		  {
		  	dakota_options_string.append("\n		  max_iterations = " + string_max_iterations);
		  }
		  if(!string_max_function_evaluations.empty())
		  {
		  	dakota_options_string.append("\n 		  max_function_evaluations = " + string_max_function_evaluations);
		  }
		  if(!string_convergence_tolerance.empty())
		  {
		  	dakota_options_string.append("\n		  convergence_tolerance = " + string_convergence_tolerance);
		  }
		  if(!string_solution_target.empty())
		  {
		  	dakota_options_string.append("\n 		  solution_target = " + string_solution_target);
		  }
		  if(!string_seed.empty())
		  {
		  	dakota_options_string.append("\n 		  seed = " + string_seed);
		  }
		  if(!diverse_dakota_input.empty())
		  {
			dakota_options_string.append("\n");
			for (int i = 0; i < diverse_dakota_input.size(); i++)
			{
				dakota_options_string.append( " " + diverse_dakota_input[i]);
			}
		  } else if(method == "coliny_ea")
	   	  {
		  	dakota_options_string.append( "\n 		  population_size = " + string_population_size);
			dakota_options_string.append( "\n		  fitness_type " + fitness_type);
			dakota_options_string.append( "\n		  mutation_type " + mutation_type);
			dakota_options_string.append( "\n		  mutation_rate " + mutation_rate);
			dakota_options_string.append( "\n		  crossover_type " + crossover_type);
			dakota_options_string.append( "\n		  crossover_rate " + crossover_rate);
			dakota_options_string.append( "\n		  replacement_type " + replacement_type);
		  } else if(method == "coliny_direct")
		  {
			dakota_options_string.append( "\n                 division " + string_division);
			dakota_options_string.append( "\n                 max_boxsize_limit " + max_boxsize_limit);
			dakota_options_string.append( "\n                 min_boxsize_limit " + min_boxsize_limit);
		  }
			dakota_options_string.append( "\n	variables,");
			if( distribution == "uniform"){
				dakota_options_string.append( "\n		continuous_design = " + std::to_string(number_of_parameters));
				dakota_options_string.append( "\n		  descriptors " + param_name_string);
				dakota_options_string.append( "\n 		  initial_point " + initial_values_string);
				dakota_options_string.append( "\n 	 	  lower_bounds " + lower_bounds_string);
				dakota_options_string.append( "\n 		  upper_bounds " + upper_bounds_string);
			} else if (distribution == "normal"){
				dakota_options_string.append( "\n               active uncertain " );
				dakota_options_string.append( "\n		normal_uncertain = " + std::to_string(number_of_parameters));
				dakota_options_string.append( "\n		  descriptors " + param_name_string);
				dakota_options_string.append( "\n 		  means " + initial_values_string);
				dakota_options_string.append( "\n 		  std_deviations " + std_deviations_string);
			}

			dakota_options_string.append( "\n	interface,");
			dakota_options_string.append( "\n		direct");
			dakota_options_string.append( "\n		  analysis_driver = 'plugin_opensmoke'");
			dakota_options_string.append( "\n	responses,");
			dakota_options_string.append( "\n		num_objective_functions = 1");
			
			// Options to use other optimization method (gradient-based)

			if (gradient_option == true){
				dakota_options_string.append( "\n               numerical_gradients");
				dakota_options_string.append( "\n               	method_source dakota");
				dakota_options_string.append( "\n               	interval_type forward");
				dakota_options_string.append( "\n               	fd_step_size = 1.e-5");
			}
			else{
				dakota_options_string.append( "\n		no_gradients");
			}

			dakota_options_string.append( "\n		no_hessians");
	}

	//void Read_Input::PrintFinalMechanism(std::vector<double> best_parameters)
	void Read_Input::PrintFinalMechanism() {
  		
		// if the "path to the folder" do not exist create it.
		if(!boost::filesystem::exists(path_folder))
		{
			OpenSMOKE::CreateDirectory(path_folder);	
		}

		//Set CHEMKIN file name and directory to it
		std::string name = "kinetics_OptiSMOKEpp.CKI";
        boost::filesystem::path filename = path_folder + "/" + name;
		//Open the fKinetics ofstream in base out mode
        std::ofstream fKinetics;
		fKinetics.open(filename.c_str(), std::ios::out);
		//
		fKinetics << "! Optimized mechanism created from OptiSMOKE++" << std::endl;
		fKinetics << "! Authors: BRITE (BURN) and CRECK groups" << std::endl;		

		// PRINT OUT: elements
        fKinetics << "ELEMENTS" << std::endl;
		for (unsigned int k = 0; k < thermodynamicsMapXML->elements().size(); k++)
		{
					fKinetics << thermodynamicsMapXML->elements()[k] << std::endl;
		}
		fKinetics << "END" << std::endl;
		fKinetics << std::endl;
		// PRINT OUT: species
		fKinetics << "SPECIES" << std::endl;
		unsigned int count = 0;	
		for (unsigned int k = 0; k < thermodynamicsMapXML->NumberOfSpecies(); k++)
		{
			count++;
					fKinetics << thermodynamicsMapXML->NamesOfSpecies()[k] << "  ";
					if (count % 6 == 0)
			{
				fKinetics << std::endl;
			}
		}
		fKinetics << std::endl;
		fKinetics << "END" << std::endl;
		fKinetics << std::endl;

		// PRINT OUT: reactions
	    fKinetics << "REACTIONS" << std::endl;

        for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
        {
			// initialize the stringstream for the reaction data, a string for the specific reaction
			std::stringstream reaction_data;
        	std::string reaction_string;
			// retrieve list of species from the thermodynamic  map 
			std::vector<std::string> list_species = thermodynamicsMapXML->NamesOfSpecies();
	        // Get the reaction string and erase white spaces in the string
        	preprocessor_kinetics->reactions()[k].GetReactionString(thermodynamicsMapXML->NamesOfSpecies(), reaction_string);
			boost::erase_all(reaction_string, " ");
			// Set precision to the stringstream
			reaction_data.precision(6);
			// if the reaction tag is not one of these special formats
			if(	preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_LINDEMANN_FALLOFF && 
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_LINDEMANN_CABR && 
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_TROE_FALLOFF && 
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_TROE_CABR && 
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_SRI_FALLOFF && 
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_SRI_CABR &&
				preprocessor_kinetics->reactions()[k].Tag() != PhysicalConstants::REACTION_EXTENDEDFALLOFF )
			{
				// Print normal reaction! or even dummy values for PLOG!
				reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML->A(k) / preprocessor_kinetics->reactions()[k].A_conversion();
				reaction_data.precision(6);
				reaction_data.width(12);
				reaction_data << std::fixed << std::right << kineticsMapXML->Beta(k);
				reaction_data.precision(6);
				reaction_data << std::setw(17) << std::right << kineticsMapXML->E_over_R(k) * PhysicalConstants::R_cal_mol << std::endl;
				// If the reaction is duplicate
				if(preprocessor_kinetics->reactions()[k].IsDuplicate() == true)
				{
							reaction_data << " DUPLICATE" << std::endl;
				}
				if(preprocessor_kinetics->reactions()[k].IsExplicitlyReversible() == true)
				{
					reaction_data << " REV /  ";
					reaction_data.precision(4);
					reaction_data << std::scientific << preprocessor_kinetics->reactions()[k].A_reversible() / preprocessor_kinetics->reactions()[k].Arev_conversion() << "  ";
					reaction_data.precision(3);
					reaction_data <<std::fixed <<  preprocessor_kinetics->reactions()[k].Beta_reversible() << "  ";
					reaction_data.precision(2);
					reaction_data << preprocessor_kinetics->reactions()[k].E_over_R_reversible() * PhysicalConstants::R_cal_mol;
					reaction_data << "  /" << std::endl;
				}
				// If it is PLOG
				if(preprocessor_kinetics->reactions()[k].IsPressureLog() == true)
				{
					int pos_classic_plog_reaction = std::find(indices_of_classic_plogs.begin(),indices_of_classic_plogs.end(),k+1)-indices_of_classic_plogs.begin();
					reaction_data.unsetf(std::ios_base::floatfield);
					reaction_data.precision(6);

					for(unsigned int l = 0; l < preprocessor_kinetics->reactions()[k].plog_coefficients().size() - 2; l++)
					{
						if(l % 4 == 0)
						{
							reaction_data << " PLOG /  ";
						}

						//reaction_data << std::showpoint << std::setw(16) << std::scientific << std::left << preprocessor_kinetics->reactions()[k].plog_coefficients()[l];
						reaction_data << std::showpoint << std::setw(16) << std::scientific << std::left << kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).PLOG_to_Print()[l];	
						if((l+1) % 4 == 0 || l == preprocessor_kinetics->reactions()[k].plog_coefficients().size() - 3)
						{
							reaction_data << "/" << std::endl;
						}
					}                    
				}
				// If it is ExtPLOGopt // Modified //
				if (preprocessor_kinetics->reactions()[k].IsExtendedPressureLogOpt() == true)
				{
					int pos_extendedplogs_opt = std::find(indices_of_extendedplogs_opt.begin(),indices_of_extendedplogs_opt.end(),k+1)-indices_of_extendedplogs_opt.begin();
					kineticsMapXML->extendedplogopt_reactions(pos_extendedplogs_opt).WriteCHEMKINOnASCIIFile(reaction_data);
				}


				// WORK HERE // 
				if (preprocessor_kinetics->reactions()[k].IsExtendedPressureLog() == true)
				{
					// Low-pressure parameters
					int pos_extendedplog = std::find(indices_of_extendedplogs.begin(),indices_of_extendedplogs.end(),k+1)-indices_of_extendedplogs.begin();
					kineticsMapXML->extendedpressurelog_reactions(pos_extendedplog).WriteCHEMKINOnASCIIFile(reaction_data);
				}
				// WORK HERE!!!

				if(preprocessor_kinetics->reactions()[k].IsJanevLanger() == true)
				{
					reaction_data.unsetf(std::ios_base::floatfield);
					reaction_data.precision(6);
					for(unsigned int l = 0; l < preprocessor_kinetics->reactions()[k].janev_langer_coefficients().size(); l++)
					{
						if(l % 5 == 0)
						{
							reaction_data << " JAN /  ";
						}

						reaction_data << std::showpoint << preprocessor_kinetics->reactions()[k].janev_langer_coefficients()[l]<< " ";
						if((l+1) % 5 == 0 || l == preprocessor_kinetics->reactions()[k].janev_langer_coefficients().size() - 1)
						{
							reaction_data << "/" << std::endl;
						}
					}
				}

				if(preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_THIRDBODY)
				{
					std::vector<double> third_body_efficiencies_ = preprocessor_kinetics->reactions()[k].third_body_efficiencies();
					for(unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
					{
						int third_body_index = preprocessor_kinetics->reactions()[k].third_body_indices()[j];
							reaction_data << list_species[third_body_index] << "/ ";
							reaction_data.precision(2);
							reaction_data << std::showpoint << std::fixed << std::left <<kineticsMapXML->ThirdBody(k, third_body_index) << "/ ";
					}    
					if(preprocessor_kinetics->reactions()[k].third_body_efficiencies().size() != 0)
					{
							reaction_data << std::endl;
					}	
				}
	
				if(preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_CHEBYSHEV)
				{
					reaction_data.unsetf(std::ios_base::floatfield);          
		
					if( preprocessor_kinetics->reactions()[k].chebyshev_temperature_limits()[0] != 300 || preprocessor_kinetics->reactions()[k].chebyshev_temperature_limits()[1] != 2500 )
					{
							reaction_data << " TCHEB/ ";
							reaction_data.precision(1);
							reaction_data << std::showpoint <<std::fixed << preprocessor_kinetics->reactions()[k].chebyshev_temperature_limits()[0] << " " << preprocessor_kinetics->reactions()[k].chebyshev_temperature_limits()[1];
							reaction_data << " /" << std::endl;
					}
	
					if(preprocessor_kinetics->reactions()[k].chebyshev_pressure_limits()[0] != 0.001 ||  preprocessor_kinetics->reactions()[k].chebyshev_pressure_limits()[1] != 100)
					{
							reaction_data << " PCHEB/ ";
							reaction_data.precision(4);
							reaction_data << std::showpoint <<std::fixed << std::left << preprocessor_kinetics->reactions()[k].chebyshev_pressure_limits()[0] << " " << preprocessor_kinetics->reactions()[k].chebyshev_pressure_limits()[1];
							reaction_data << " /" << std::endl;
					}
		
					unsigned int chebyshev_size =	boost::lexical_cast<unsigned int>(preprocessor_kinetics->reactions()[k].chebyshev_coefficients()[0]) * boost::lexical_cast<unsigned int>(preprocessor_kinetics->reactions()[k].chebyshev_coefficients()[1]);
			
					reaction_data.unsetf(std::ios_base::floatfield);
					reaction_data.precision(6);
					for(unsigned int l=0;l<chebyshev_size+2;l++)
					{                    
							if(l%6 == 0)
				{
								reaction_data << " CHEB/ ";
							}
							if(l < 2)
				{
								reaction_data << std::noshowpoint << preprocessor_kinetics->reactions()[k].chebyshev_coefficients()[l] << " ";
							}
							else
				{
								reaction_data << std::showpoint << preprocessor_kinetics->reactions()[k].chebyshev_coefficients()[l] << " ";
							}
							if((l+1)%6 == 0)
				{
								reaction_data << " /" << std::endl;
							}
							if(l == chebyshev_size + 1 && (l+1)%6 != 0 )
				{
								reaction_data << " /" << std::endl;
				}
					}
				}
			}
			else if (preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
			{
				// High-pressure kinetic parameters
				reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML->A_falloff_inf(k) / preprocessor_kinetics->reactions()[k].A_inf_conversion();
				reaction_data.precision(3);
				reaction_data.width(9);
				reaction_data << std::fixed << std::right << kineticsMapXML->Beta_falloff_inf(k);
				reaction_data.precision(2);
				reaction_data << std::setw(13) << std::right << kineticsMapXML->E_over_R_falloff_inf(k) * PhysicalConstants::R_cal_mol << std::endl;
				reaction_data.unsetf(std::ios_base::floatfield);
				// Low-pressure parameters
				OpenSMOKE::ExtendedFallOff extendedFallOff;
				extendedFallOff.Setup(preprocessor_kinetics->reactions()[k].extendedfalloff_coefficients());
				extendedFallOff.WriteCHEMKINOnASCIIFile(reaction_data);
				// Add third body efficiencies
				bool iThirdBody_ = false;
				std::vector<double> third_body_efficiencies_ = preprocessor_kinetics->reactions()[k].third_body_efficiencies();
				for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
				{
					int third_body_index = preprocessor_kinetics->reactions()[k].third_body_indices()[j];
					reaction_data << list_species[third_body_index] << "/ ";
					reaction_data.precision(2);
					reaction_data << std::fixed << std::showpoint << kineticsMapXML->ThirdBody(k, third_body_index) << "/  ";
					iThirdBody_ = true;
				}
				if (iThirdBody_ == true)
				{
					reaction_data << std::endl;
				}
			}
			else
    		{
				int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),k+1)-indices_of_falloff_reactions.begin();
            	reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction) / preprocessor_kinetics->reactions()[k].A_inf_conversion();;
            	reaction_data.precision(6);
            	reaction_data.width(12);
				reaction_data <<std::fixed << std::right << kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction);
        	    reaction_data.precision(4);
            	reaction_data << std::setw(15) << std::right << kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction) * PhysicalConstants::R_cal_mol << std::endl;

            	if(preprocessor_kinetics->reactions()[k].IsDuplicate() == true)
				{
                	reaction_data << " DUPLICATE" << std::endl;
	            }
        	    if(preprocessor_kinetics->reactions()[k].IsExplicitlyReversible() == true)
            	{
                	reaction_data << " REV /  ";
					reaction_data.precision(4);
					reaction_data << std::scientific << preprocessor_kinetics->reactions()[k].A_reversible() / preprocessor_kinetics->reactions()[k].Arev_conversion()<< "  ";
					reaction_data.precision(3);
					reaction_data <<std::fixed << preprocessor_kinetics->reactions()[k].Beta_reversible() << "  ";
					reaction_data.precision(2);
					reaction_data << preprocessor_kinetics->reactions()[k].E_over_R_reversible() * PhysicalConstants::R_cal_mol;
					reaction_data << "  /" << std::endl;
				}            
            	if(	preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_FALLOFF      ||
					preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ||
					preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_FALLOFF         )
				{
					reaction_data << " LOW/";
				}
				else if(preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_CABR ||
					preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_CABR)
				{
						reaction_data << " HIGH/";
				}

				// print out the parameter values with the desired formatting and precision // For LOW or HIGH
				reaction_data.width(15);
				reaction_data.precision(6);
				reaction_data << std::right << std::scientific << kineticsMapXML->A(k) / preprocessor_kinetics->reactions()[k].A_conversion();
				reaction_data.precision(6);
				reaction_data.width(14);
				reaction_data <<std::fixed << std::right << kineticsMapXML->Beta(k);
				reaction_data.precision(4);
				reaction_data << std::setw(16) << std::right << kineticsMapXML->E_over_R(k) * PhysicalConstants::R_cal_mol << " /" << std::endl;

				// prints out the TROE parameters for the blending function F
				if(		preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_FALLOFF || 
						preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_CABR      )
				{
					reaction_data << "TROE/";
					reaction_data.width(11);
					reaction_data.unsetf(std::ios_base::floatfield);
					for(unsigned int j = 0; j < preprocessor_kinetics->reactions()[k].troe().size(); j++)
					{
							reaction_data.precision(4);
							reaction_data << std::showpoint << "   " << preprocessor_kinetics->reactions()[k].troe()[j];
					}
					reaction_data << "/" << std::endl;                        
				}

				if(		preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ||
						preprocessor_kinetics->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_CABR  	 )
				{
					reaction_data << "SRI/ ";
					for(unsigned int j = 0; j < preprocessor_kinetics->reactions()[k].sri().size(); j++)
					{
							reaction_data.precision(4);
							reaction_data << std::showpoint << "  " << preprocessor_kinetics->reactions()[k].sri()[j];
					}
		
					reaction_data << "/" << std::endl;
				}
		
				bool iThirdBody_ = false;
				std::vector<double> third_body_efficiencies_ = preprocessor_kinetics->reactions()[k].third_body_efficiencies();
				for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
				{
						int third_body_index = preprocessor_kinetics->reactions()[k].third_body_indices()[j];
						reaction_data << list_species[third_body_index] << "/ ";
						reaction_data.precision(5);
						reaction_data << std::fixed << std::showpoint << kineticsMapXML->ThirdBody(k, third_body_index) << "/  ";
						iThirdBody_ = true;
				}
            	
            	if(iThirdBody_ == true)
				{
                				reaction_data << std::endl;  
				}
        	}
			// WEDONTUSETHISFORMAT
			if(preprocessor_kinetics->reactions()[k].IsFit1())
			{
					reaction_data.unsetf(std::ios_base::floatfield);
					reaction_data.precision(6);
					for(unsigned int l = 0; l < preprocessor_kinetics->reactions()[k].fit1_coefficients().size(); l++)
					{
						reaction_data << " FIT1 /  ";
						reaction_data << std::showpoint << preprocessor_kinetics->reactions()[k].fit1_coefficients()[l]<< " ";
					}
					reaction_data << "/" << std::endl;
			}
			// WEDONTUSETHISFORMAT
			if(preprocessor_kinetics->reactions()[k].IsLandauTeller())
			{
					reaction_data.unsetf(std::ios_base::floatfield);
					reaction_data.precision(3);
					for(unsigned int l = 0; l < preprocessor_kinetics->reactions()[k].landau_teller_coefficients().size(); l++)
					{
						reaction_data << " LT /  ";
						reaction_data << std::showpoint << preprocessor_kinetics->reactions()[k].landau_teller_coefficients()[l]<< " ";
					}
					reaction_data << "/" << std::endl;
			}
			// WEDONTUSETHISFORMAT
	        if(preprocessor_kinetics->reactions()[k].IsFORD())
			{
				reaction_data.unsetf(std::ios_base::floatfield);
				reaction_data.precision(4);
				std::vector<double> reactant_lambda_ = preprocessor_kinetics->reactions()[k].reactant_lambda();
				std::vector<unsigned int> reactant_lambda_indices_ = preprocessor_kinetics->reactions()[k].reactant_lambda_indices();
				
				for(unsigned int l = 0; l < reactant_lambda_.size(); l++)
				{
					reaction_data << " FORD /  ";
							int index = reactant_lambda_indices_[l];
	
					reaction_data << list_species[index] << "  ";
	
					reaction_data << std::showpoint <<std::fixed << reactant_lambda_[l];
					reaction_data << "/" << std::endl;
				}
			}
			// Finally print out the string stream to ofstream and the game is done for reaction k.    
			fKinetics << reaction_data.str();
			fKinetics << std::endl;
        }
		// close the reactions part
		fKinetics << "END" << std::endl;
		fKinetics << std::endl; 

		// close the ofstream
		fKinetics.close();

		std::cout<<""<<std::endl;
		std::cout<<"Wrote optimized kinetics to "<<path_folder<<"/"<<name<<std::endl;
		std::cout<<""<<std::endl;
		
	}
	
	void Read_Input::ReadReactionClassesDefinition(boost::filesystem::path ReactionClassFile)
	{
		/*
			File classes defnition description:
			- Riga 1: Nome della classe
			- Riga 2: Indici delle reazioni da ottimizzare
			- Riga 3: Uncertainty factors
			- Riga 4: Target names è ridondante in realtà
			- Riga 5: Scaling factor lnA
			- Riga 6: Active or not 0 or 1
			- Riga 7: Scaling factor Beta
			- Riga 8: Active or not 0 or 1
			- Riga 9: Scaling factor E_over_R
			- Riga 10: Active or not 0 or 1
		*/

		boost::filesystem::ifstream fileHandler(ReactionClassFile.c_str());
		std::string line;
		int numberOfLines = 0;
		std::vector<std::string> content;
		while (getline(fileHandler, line)) {
			numberOfLines++;
			content.push_back(line);
		}

		numberOfReactionClasses = numberOfLines / 10;

		std::vector<std::string> reactionClassName;
		std::vector<std::string> str_reaction_index;
		std::vector<std::string> str_unc;
		std::vector<std::string> str_target;
		std::vector<std::string> str_scaling_lnA;
		std::vector<std::string> str_which_lnA;
		std::vector<std::string> str_scaling_Beta;
		std::vector<std::string> str_which_Beta;
		std::vector<std::string> str_scaling_E_over_R;
		std::vector<std::string> str_which_E_over_R;

		for(int i = 0; i<numberOfReactionClasses; i++){
			reactionClassName.push_back(content[0 + i * 10]);
			str_reaction_index.push_back(content[1 + i * 10]);
			str_unc.push_back(content[2 + i * 10]);
			str_target.push_back(content[3 + i * 10]);
			str_scaling_lnA.push_back(content[4 + i * 10]);
			str_which_lnA.push_back(content[5 + i * 10]);
			str_scaling_Beta.push_back(content[6 + i * 10]);
			str_which_Beta.push_back(content[7 + i * 10]);
			str_scaling_E_over_R.push_back(content[8 + i * 10]);
			str_which_E_over_R.push_back(content[9 + i * 10]);
		}

		for(int i = 0; i < numberOfReactionClasses; i++){
			std::vector<std::string> tmp_idx; // reaction index
			std::vector<int> tmp_int_idx;
			std::vector<int> tmp_int_idx_lnA;
			std::vector<int> tmp_int_idx_Beta;
			std::vector<int> tmp_int_idx_E_over_R;

			std::vector<std::string> tmp_unc; // uncertainty factors
			std::vector<double> tmp_double_unc;

			std::vector<std::string> tmp_scaling_lnA; // scaling_lnA
			std::vector<double> tmp_double_scaling_lnA;

			std::vector<std::string> tmp_scaling_Beta; // scaling_Beta
			std::vector<double> tmp_double_scaling_Beta;

			std::vector<std::string> tmp_scaling_E_over_R; // scaling_E_over_R
			std::vector<double> tmp_double_scaling_E_over_R;

			std::vector<std::string> tmp_which_lnA; // which_lnA
			std::vector<int> tmp_int_which_lnA;

			std::vector<std::string> tmp_which_Beta; // which_Beta
			std::vector<int> tmp_int_which_Beta;

			std::vector<std::string> tmp_which_E_over_R; // which E_over_R
			std::vector<int> tmp_int_which_E_over_R;

			boost::split(tmp_idx, str_reaction_index[i], boost::is_any_of(" "));
			boost::split(tmp_unc, str_unc[i], boost::is_any_of(" "));
			boost::split(tmp_scaling_lnA, str_scaling_lnA[i], boost::is_any_of(" "));
			boost::split(tmp_scaling_Beta, str_scaling_Beta[i], boost::is_any_of(" "));
			boost::split(tmp_scaling_E_over_R, str_scaling_E_over_R[i], boost::is_any_of(" "));
			boost::split(tmp_which_lnA, str_which_lnA[i], boost::is_any_of(" "));
			boost::split(tmp_which_Beta, str_which_Beta[i], boost::is_any_of(" "));
			boost::split(tmp_which_E_over_R, str_which_E_over_R[i], boost::is_any_of(" "));

			// They have all the same size so just one loop
			for(int j = 0; j < tmp_idx.size(); j++){
				tmp_double_unc.push_back(std::stod(tmp_unc[j]));
				tmp_int_idx.push_back(std::stoi(tmp_idx[j]));
				if(std::stoi(tmp_which_lnA[j]) != 0){
					tmp_double_scaling_lnA.push_back(std::stod(tmp_scaling_lnA[j]));
					tmp_int_idx_lnA.push_back(std::stoi(tmp_idx[j]));
				}
				if(std::stoi(tmp_which_Beta[j]) != 0){
					tmp_double_scaling_Beta.push_back(std::stod(tmp_scaling_Beta[j]));
					tmp_int_idx_Beta.push_back(std::stoi(tmp_idx[j]));
				}
				if(std::stoi(tmp_which_E_over_R[j]) != 0){
					tmp_double_scaling_E_over_R.push_back(std::stod(tmp_scaling_E_over_R[j]));
					tmp_int_idx_E_over_R.push_back(std::stoi(tmp_idx[j]));
				}
			}

			matrixOfReactionIndex.push_back(tmp_int_idx);
			matrixOflnA.push_back(tmp_int_idx_lnA);
			matrixOfBeta.push_back(tmp_int_idx_Beta);
			matrixOfEoverR.push_back(tmp_int_idx_E_over_R);
			matrixOfUnceratintyFactors.push_back(tmp_double_unc);

			matrixOfscalinglnA.push_back(tmp_double_scaling_lnA);
			matrixOfscalingBeta.push_back(tmp_double_scaling_Beta);
			matrixOfscalingEoverR.push_back(tmp_double_scaling_E_over_R);

			tmp_idx.clear();
			tmp_int_idx_lnA.clear();
			tmp_int_idx_Beta.clear();
			tmp_int_idx_E_over_R.clear();
			tmp_unc.clear();
			tmp_double_unc.clear();
			tmp_scaling_lnA.clear();
			tmp_double_scaling_lnA.clear();
			tmp_scaling_Beta.clear();
			tmp_double_scaling_Beta.clear();
			tmp_scaling_E_over_R.clear();
			tmp_double_scaling_E_over_R.clear();
			tmp_which_lnA.clear();
			tmp_which_Beta.clear();
			tmp_which_E_over_R.clear();
			tmp_int_which_lnA.clear();
			tmp_int_which_Beta.clear();
			tmp_int_which_E_over_R.clear();
		}

		for(int i = 0; i<matrixOfReactionIndex.size(); i++){
			list_of_target_lnA.push_back(matrixOfReactionIndex[i][0]);
			list_of_target_Beta.push_back(matrixOfReactionIndex[i][0]);
			list_of_target_E_over_R.push_back(matrixOfReactionIndex[i][0]);
			list_of_target_uncertainty_factors.push_back(matrixOfReactionIndex[i][0]);

			list_of_initial_lnA.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML->A(matrixOfReactionIndex[i][0]-1))));
			list_of_initial_Beta.push_back(boost::lexical_cast<std::string>(kineticsMapXML->Beta(matrixOfReactionIndex[i][0]-1)));
			list_of_initial_E_over_R.push_back(boost::lexical_cast<std::string>(kineticsMapXML->E_over_R(matrixOfReactionIndex[i][0]-1)));
		}

		for(int i = 0; i<matrixOfUnceratintyFactors.size(); i++){
			list_of_uncertainty_factors.push_back(matrixOfUnceratintyFactors[i][0]);
		}
		std::cout << "Target lnA" << std::endl;
		for(int i = 0; i < list_of_target_lnA.size(); i++){
			std::cout << list_of_target_lnA[i] << std::endl;
		}

		std::cout << "Target Beta" << std::endl;
		for(int i = 0; i < list_of_target_Beta.size(); i++){
			std::cout << list_of_target_Beta[i] << std::endl;
		}

		std::cout << "Target E_over_R" << std::endl;
		for(int i = 0; i < list_of_target_E_over_R.size(); i++){
			std::cout << list_of_target_E_over_R[i] << std::endl;
		}

		std::cout << "Matrix of lnA" << std::endl;
		for(int i = 0; i < matrixOflnA[1].size(); i++){
			std::cout << matrixOflnA[1][i] << std::endl;
		}
		std::cout << "Matrix Beta" << std::endl;
		for(int i = 0; i < matrixOfBeta[1].size(); i++){
			std::cout << matrixOfBeta[1][i] << std::endl;
		}
		std::cout << "Matrix of E over R" << std::endl;
		for(int i = 0; i < matrixOfEoverR[1].size(); i++){
			std::cout << matrixOfEoverR[1][i] << std::endl;
		}
		std::cout << "Matrix of scaling lnA" << std::endl;
		for(int i = 0; i < matrixOfscalinglnA[1].size(); i++){
			std::cout << matrixOfscalinglnA[1][i] << std::endl;
		}
		std::cout << "Matrix of scaling Beta" << std::endl;
		for(int i = 0; i < matrixOfscalingBeta[1].size(); i++){
			std::cout << matrixOfscalingBeta[1][i] << std::endl;
		}
		std::cout << "Matrix of scaling E over R" << std::endl;
		for(int i = 0; i < matrixOfscalingEoverR[1].size(); i++){
			std::cout << matrixOfscalingEoverR[1][i] << std::endl;
		}
	}
}

