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

  InputManager::InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary) : dictionary_(dictionary){

    input_file_name_ = "input.dic";
    main_dictionary_ = "OptiSMOKEpp";
    output_folder_ = "Output";
    kinetics_folder_ = "kinetics";
    optimized_kinetics_folder_ = "Optimized_kinetics";

    iXml_ = false;
    iTransport_ = false;
  }

  InputManager::~InputManager(){}

  void InputManager::SetInputOptions(int argc, char* argv[]){

    // Up to isMaster_ = false; should go in the ctor
    // however let's stay simple for the moment
    rank_ = 0;
    nprocs_ = 1;

#ifdef OPTISMOKE_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
#endif

    if (rank_ == 0)
      isMaster_ = true;
    else
      isMaster_ = false;

    //Input Options
    {
      po::options_description desc("Allowed options");
      desc.add_options()
        ("help", "Help Message")
        ("input", po::value<std::string>(), "Input File Path (default: \"input.dic\")");

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      if (vm.count("help")){
        if (rank_ == 0)
          std::cout << desc << std::endl;

#ifdef OPTISMOKE_USE_MPI
        MPI_Finalize();
#endif

        exit(0);
      }

      if (vm.count("input")){
        input_file_name_ = vm["input"].as<std::string>();
      }
    }
  }

  void InputManager::ReadDictionary(){
    // It can be trivial however this is for future 
    // parallelization with MPI see:
    // https://github.com/astagni/DoctorSMOKEpp/blob/main/src/DataManager.hpp
    // Remember that this all goes under rank=0
    if(rank_ == 0)
      ReadMainDictionary();

#ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(rank_>0)
      ReadMainDictionary();

    if(rank_ == 0){
      if(!iXml_){
        if(!iTransport_){
          OpenSMOKE::RapidKineticMechanismWithoutTransport(
              output_folder_ / kinetics_data_.chemkin_output(),
              kinetics_data_.chemkin_thermodynamics(),
              kinetics_data_.chemkin_kinetics());
        }
        else{
          OpenSMOKE::RapidKineticMechanismWithTransport(
              output_folder_ / kinetics_data_.chemkin_output(),
              kinetics_data_.chemkin_transport(),
              kinetics_data_.chemkin_thermodynamics(),
              kinetics_data_.chemkin_kinetics());
        }
      }
    }

    CreateMaps();

#ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  void InputManager::ReadMainDictionary(){
    dictionary_.ReadDictionariesFromFile(input_file_name_);
    dictionary_(main_dictionary_).SetGrammar(main_grammar_);

    dictionary_(main_dictionary_).ReadPath("@OutputFolder", output_folder_);
    if(!fs::exists(output_folder_))
      fs::create_directories(output_folder_);

    if(dictionary_(main_dictionary_).CheckOption("@KineticsFolder")){
      iXml_ = true;
      dictionary_(main_dictionary_).ReadPath("@KineticsFolder", kinetics_folder_);
      if(!fs::exists(kinetics_folder_)){
        OptiSMOKE::FatalErrorMessage("The @KineticsFolder path does not exists!");
      }
      OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
    }
    else if(dictionary_(main_dictionary_).CheckOption("@KineticsPreProcessor")){
      dictionary_(main_dictionary_).ReadDictionary("@KineticsPreProcessor", preprocessor_dictionary_);
      kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary_, iTransport_);
    }
    else{
      OptiSMOKE::FatalErrorMessage("Please provide the kinetic mechanism through one of the following keywords: @KineticsFolder | @KineticsPreProcessor");
    }

    // path data set input files
    dictionary_(main_dictionary_).ReadOption("@ListOfExperimentalDataFiles", path_experimental_data_files_);

    // optimization libraries
    dictionary_(main_dictionary_).ReadString("@OptimizationLibrary", optimization_library_);

    // Dictionaries
    if(optimization_library_ == "dakota"){
      // Dakota options
      dictionary_(main_dictionary_).ReadDictionary("@DakotaOptions", dakota_dictionary_);
      dakota_options_.SetupFromDictionary(dictionary_, dakota_dictionary_);
    }
    else if(optimization_library_ == "nlopt"){
      // NLOPT options
      dictionary_(main_dictionary_).ReadDictionary("@NLOPTOptions", nlopt_dictionary_);
      nlopt_options_.SetupFromDictionary(dictionary_, nlopt_dictionary_);
    }
    else
      OptiSMOKE::FatalErrorMessage("Unknown optimization library. Available are: dakota | nlopt");

    if(dictionary_(main_dictionary_).CheckOption("@CurveMatchingOptions")){
      dictionary_(main_dictionary_).ReadDictionary("@CurveMatchingOptions", curvematching_dictionary_);
      curvematching_options_.SetupFromDictionary(dictionary_, curvematching_dictionary_);
    }

    // Optimization setup
    dictionary_(main_dictionary_).ReadDictionary("@OptimizationSetup", optimization_setup_dictionary_);
    optimization_setup_.SetupFromDictionary(dictionary_, optimization_setup_dictionary_);

    // Optimization target
    dictionary_(main_dictionary_).ReadDictionary("@OptimizationTarget", optimization_target_dictionary_);
    optimization_target_.SetupFromDictionary(dictionary_, optimization_target_dictionary_);

  }

  void InputManager::CreateMaps(){
    // This goes under kinetics map

    fs::path path_kinetics_output;
    if (!iXml_) // To be interpreted on-the-fly
      path_kinetics_output = output_folder_ / kinetics_data_.chemkin_output();
    else if (iXml_) // Already in XML format
      path_kinetics_output = kinetics_folder_;

    std::cout.setstate(std::ios_base::failbit); // Disable video output
    boost::property_tree::ptree ptree;
    boost::property_tree::read_xml( (path_kinetics_output / "kinetics.xml").string(), ptree );

    thermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
    kineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML_, ptree);
    if(iTransport_)
      transportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);

    boost::property_tree::ptree nominal_ptree;
    boost::property_tree::read_xml( (path_kinetics_output / "kinetics.xml").string(), nominal_ptree);

    nominalthermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(nominal_ptree);
    nominalkineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*nominalthermodynamicsMapXML_, nominal_ptree);
    if(iTransport_)
      nominaltransportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(nominal_ptree);
    std::cout.clear(); // Re-enable video output
  }

  void InputManager::SetUpNLOPT(){

    FromTargetToInitialParameter();

    ComputeBoundaries();

    TargetsPreliminaryOptions();

    parametric_file_name_ = output_folder_ / "optimization.out";

    std::vector<std::string> initial_values_str;
    std::vector<std::string> lb_str;
    std::vector<std::string> ub_str;

    boost::split(initial_values_str, initial_values_string_, boost::is_any_of(" "));
    boost::split(lb_str, lower_bounds_string_, boost::is_any_of(" "));
    boost::split(ub_str, upper_bounds_string_, boost::is_any_of(" "));
    boost::erase_all(param_name_string_, "'");
    boost::split(param_str_, param_name_string_, boost::is_any_of(" "));

    initial_values_str.pop_back();
    lb_str.pop_back();
    ub_str.pop_back();
    param_str_.pop_back();

    initial_values_.resize(initial_values_str.size());
    std::transform(initial_values_str.begin(), initial_values_str.end(), initial_values_.begin(), [](const std::string& str){return std::stod(str);});

    lb_.resize(lb_str.size());
    std::transform(lb_str.begin(), lb_str.end(), lb_.begin(), [](const std::string& str) {return std::stod(str);});

    ub_.resize(ub_str.size());
    std::transform(ub_str.begin(), ub_str.end(), ub_.begin(), [](const std::string& str) {return std::stod(str);});
  }

  void InputManager::DakotaInputString(){

    FromTargetToInitialParameter();

    ComputeBoundaries();

    TargetsPreliminaryOptions();

    dakota_input_string_ = " environment,"
      "\n  tabular_data";
    dakota_input_string_.append("\n   tabular_data_file '" + output_folder_.string() + "/" + dakota_options_.tabular_data_file() + "'");

    dakota_input_string_.append("\n method,"); 
    dakota_input_string_.append("\n  " + dakota_options_.method());
    dakota_input_string_.append("\n   max_iterations = " + dakota_options_.max_iterations());
    dakota_input_string_.append("\n   max_function_evaluations = " + dakota_options_.max_function_evaluations());
    dakota_input_string_.append("\n   convergence_tolerance = " + dakota_options_.convergence_tolerance());
    dakota_input_string_.append("\n   solution_target = " + dakota_options_.solution_target());
    dakota_input_string_.append("\n   seed = " + dakota_options_.seed());

    if(dakota_options_.diverse_input()){
      dakota_input_string_.append("\n");
      for (int i = 0; i < dakota_options_.diverse_dakota_input().size(); i++)
      {
        dakota_input_string_.append( " " + dakota_options_.diverse_dakota_input()[i]);
      }
    } 
    else if(dakota_options_.method() == "coliny_ea"){
      dakota_input_string_.append("\n   population_size = " + dakota_options_.population_size());
      dakota_input_string_.append("\n   fitness_type " + dakota_options_.fitness_type());
      dakota_input_string_.append("\n   mutation_type " + dakota_options_.mutation_type());
      dakota_input_string_.append("\n   mutation_rate " + dakota_options_.mutation_rate());
      dakota_input_string_.append("\n   crossover_type " + dakota_options_.crossover_type());
      dakota_input_string_.append("\n   crossover_rate " + dakota_options_.crossover_rate());
      dakota_input_string_.append("\n   replacement_type " + dakota_options_.replacement_type());
    } 
    else if(dakota_options_.method() == "coliny_direct"){
      dakota_input_string_.append("\n   division " + dakota_options_.division());
      dakota_input_string_.append("\n   max_boxsize_limit " + dakota_options_.max_boxsize_limit());
      dakota_input_string_.append("\n   min_boxsize_limit " + dakota_options_.min_boxsize_limit());
    }
    else{
      OptiSMOKE::FatalErrorMessage("Available methods currently implemented are coliny_ea | coliny_direct");
    }

    dakota_input_string_.append("\n variables,");

    if(optimization_setup_.parameter_distribution() == "uniform"){
      dakota_input_string_.append("\n  continuous_design = " + std::to_string(optimization_target_.number_of_parameters()));
      dakota_input_string_.append("\n   descriptors " + param_name_string_);
      dakota_input_string_.append("\n   initial_point " + initial_values_string_);
      dakota_input_string_.append("\n   lower_bounds " + lower_bounds_string_);
      dakota_input_string_.append("\n   upper_bounds " + upper_bounds_string_);
    }
    else if (optimization_setup_.parameter_distribution() == "normal"){
      dakota_input_string_.append("\n  active uncertain " );
      dakota_input_string_.append("\n  normal_uncertain = " + std::to_string(optimization_target_.number_of_parameters()));
      dakota_input_string_.append("\n   descriptors " + param_name_string_);
      dakota_input_string_.append("\n   means " + initial_values_string_);
      dakota_input_string_.append("\n   std_deviations " + std_deviations_string_);
    }// Da fare check su consistenza nell' input sul tipo di parameter boundary

    dakota_input_string_.append("\n interface,");
    dakota_input_string_.append("\n  direct");
    dakota_input_string_.append("\n  analysis_driver = 'opensmoke_plugin'");
    dakota_input_string_.append("\n responses,");
    dakota_input_string_.append("\n  num_objective_functions = 1");

    // Options to use other optimization method (e.g. gradient-based)
    // Qua forse va messa la possibilità di fare altri tipi di gradienti 
    // accordingly to dakota sicuro lo faccio ora non c'ho voglia
    if (dakota_options_.dakota_gradient() == true){
      dakota_input_string_.append("\n  numerical_gradients");
      dakota_input_string_.append("\n  method_source dakota");
      dakota_input_string_.append("\n  interval_type forward");
      dakota_input_string_.append("\n  fd_step_size = 1.e-5");
    }
    else{
      dakota_input_string_.append("\n  no_gradients");
    }

    dakota_input_string_.append("\n  no_hessians");
  }

  void InputManager::FromTargetToInitialParameter(){

    // lnA
    for(int i=0; i < optimization_target_.list_of_target_lnA().size(); i++)
      list_of_initial_lnA_.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i]-1))));

    // Beta
    for(int i=0; i < optimization_target_.list_of_target_Beta().size(); i++)
      list_of_initial_Beta_.push_back(boost::lexical_cast<std::string>(kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i]-1)));

    // E_over_R
    for(int i=0; i < optimization_target_.list_of_target_E_over_R().size(); i++)
      list_of_initial_E_over_R.push_back(boost::lexical_cast<std::string>(kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i]-1)));

    // lnA_inf
    std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
    for(int i=0; i < optimization_target_.list_of_target_lnA_inf().size(); i++){
      int pos_FallOff_Reaction = std::find(
          indices_of_falloff_reactions.begin(), 
          indices_of_falloff_reactions.end(), 
          optimization_target_.list_of_target_lnA_inf()[i]
          )-indices_of_falloff_reactions.begin();
      list_of_initial_lnA_inf_.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction))));
    }

    // Beta_inf
    for(int i=0; i < optimization_target_.list_of_target_Beta_inf().size(); i++){
      int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
          indices_of_falloff_reactions.end(),
          optimization_target_.list_of_target_Beta_inf()[i])-indices_of_falloff_reactions.begin();
      list_of_initial_Beta_inf_.push_back(boost::lexical_cast<std::string>(kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction)));
    }

    // E/R inf
    for(int i=0; i < optimization_target_.list_of_target_E_over_R_inf().size(); i++){
      int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
          indices_of_falloff_reactions.end(),
          optimization_target_.list_of_target_E_over_R_inf()[i])-indices_of_falloff_reactions.begin();
      list_of_initial_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction)));
    }

    for(int i=0; i< optimization_target_.list_of_target_thirdbody_reactions().size(); i++){
      int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
      list_of_initial_thirdbody_eff_.push_back(
          boost::lexical_cast<std::string>(
            kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i]-1, iSpecies-1)
            )
          );
    }
  }     

  void InputManager::ComputeBoundaries(){

    double T_low = 300;
    double T_high = 2500;

    // Initialize needed values at the specific size
    std::vector<double> list_of_nominal_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_nominal_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_nominal_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());

    std::vector<double> list_of_min_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_max_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_min_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_max_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_min_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> list_of_max_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());

    std::vector<double> kappa_lower_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> kappa_upper_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> kappa_lower_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> kappa_upper_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> Beta_1(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> Beta_2(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> lnA_1(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> lnA_2(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> E_over_R_1(optimization_target_.list_of_target_uncertainty_factors().size());
    std::vector<double> E_over_R_2(optimization_target_.list_of_target_uncertainty_factors().size());


    std::vector<double> list_of_nominal_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_nominal_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_nominal_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());

    std::vector<double> list_of_min_abs_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_max_abs_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_min_abs_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_max_abs_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_min_abs_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> list_of_max_abs_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());

    std::vector<double> kappa_lower_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> kappa_upper_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> kappa_lower_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> kappa_upper_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> Beta_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> Beta_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> lnA_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> lnA_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> E_over_R_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
    std::vector<double> E_over_R_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());

    if(optimization_setup_.parameter_boundaries() == "Furst"){
      for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++){
        // Nominal values of parameters
        list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i]-1));
        list_of_nominal_Beta_double[i] = nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
        list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i]-1);

        // Min and Max of lnA
        list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i] + std::log(std::pow(10, -optimization_target_.list_of_uncertainty_factors()[i]));
        list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i] + std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i]));
        if (std::find(optimization_target_.list_of_target_lnA().begin(),
              optimization_target_.list_of_target_lnA().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_lnA().end())
        {
          list_of_min_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
          list_of_max_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
        }

        // Limiting values for the rate coefficient
        kappa_lower_T_low[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
        kappa_upper_T_low[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
        kappa_lower_T_high[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);
        kappa_upper_T_high[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);

        // Calculating extreme values for Beta
        Beta_1[i] = (kappa_upper_T_low[i] - kappa_lower_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
        Beta_2[i] = (kappa_lower_T_low[i] - kappa_upper_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));

        list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
        list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);

        if (std::find(optimization_target_.list_of_target_Beta().begin(),
              optimization_target_.list_of_target_Beta().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_Beta().end())
        {
          list_of_min_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
          list_of_max_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
        }

        // Calculting extreame values of E_over_R
        lnA_1[i] = ( kappa_lower_T_high[i] - (T_low/T_high) * kappa_upper_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
        E_over_R_1[i] = lnA_1[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_upper_T_low[i] * T_low;
        lnA_2[i] = ( kappa_upper_T_high[i] - (T_low/T_high) * kappa_lower_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
        E_over_R_2[i] = lnA_2[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_lower_T_low[i] * T_low;    
        list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
        list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);

        if (std::find(optimization_target_.list_of_target_E_over_R().begin(),
              optimization_target_.list_of_target_E_over_R().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_E_over_R().end()){
          list_of_min_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
          list_of_max_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
        }
      }

      std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();

      for (unsigned int i=0; i < optimization_target_.list_of_target_uncertainty_factors_inf().size(); i++){

        int pos_FallOff_Reaction = std::find(
            indices_of_falloff_reactions.begin(), 
            indices_of_falloff_reactions.end(), 
            optimization_target_.list_of_target_uncertainty_factors_inf()[i]
            )-indices_of_falloff_reactions.begin();

        // Nominal values of inf parameters
        list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction));
        list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);
        list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction);

        // Min and Max of lnA_inf
        list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-optimization_target_.list_of_uncertainty_factors_inf()[i]));
        list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_inf()[i]));
        if (std::find(optimization_target_.list_of_target_lnA_inf().begin(), 
              optimization_target_.list_of_target_lnA_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_lnA_inf().end())
        {
          list_of_min_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
          list_of_max_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
        }

        // Limiting values for the rate coefficient
        kappa_lower_T_low_inf[i]  = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
        kappa_upper_T_low_inf[i]  = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
        kappa_lower_T_high_inf[i] = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);
        kappa_upper_T_high_inf[i] = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);

        Beta_1_inf[i] = (kappa_upper_T_low_inf[i] - kappa_lower_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
        Beta_2_inf[i] = (kappa_lower_T_low_inf[i] - kappa_upper_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
        list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
        list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
        if (std::find(optimization_target_.list_of_target_Beta_inf().begin(),
              optimization_target_.list_of_target_Beta_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_Beta_inf().end())
        {
          list_of_min_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
          list_of_max_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
        }

        // Calculting extreame values of E_over_R
        lnA_1_inf[i] = ( kappa_lower_T_high_inf[i] - (T_low/T_high) * kappa_upper_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
        E_over_R_1_inf[i] = lnA_1_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_upper_T_low_inf[i] * T_low;
        lnA_2_inf[i] = ( kappa_upper_T_high_inf[i] - (T_low/T_high) * kappa_lower_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
        E_over_R_2_inf[i] = lnA_2_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_lower_T_low_inf[i] * T_low;
        list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
        list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
        if (std::find(optimization_target_.list_of_target_E_over_R_inf().begin(),
              optimization_target_.list_of_target_E_over_R_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_E_over_R_inf().end())
        {
          list_of_min_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
          list_of_max_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
        }
      }
    }

    if(optimization_setup_.parameter_boundaries() == "Narrow"){

      for (unsigned int i=0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++){
        list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i]-1));
        list_of_nominal_Beta_double[i] = nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
        list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i]-1);


        list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10, -optimization_target_.list_of_uncertainty_factors()[i]));
        list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i]));

        if (std::find(optimization_target_.list_of_target_lnA().begin(),
              optimization_target_.list_of_target_lnA().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_lnA().end())
        {
          list_of_min_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
          list_of_max_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
        }

        Beta_1[i] = list_of_nominal_Beta_double[i]+std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors()[i])) / std::log(T_high);
        Beta_2[i] = list_of_nominal_Beta_double[i]-std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors()[i])) / std::log(T_high);

        list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
        list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);

        if (std::find(optimization_target_.list_of_target_Beta().begin(),
              optimization_target_.list_of_target_Beta().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_Beta().end())
        {
          list_of_min_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
          list_of_max_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
        }

        E_over_R_1[i] = list_of_nominal_E_over_R_double[i]-std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors()[i])) * T_low;
        E_over_R_2[i] = list_of_nominal_E_over_R_double[i]+std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors()[i])) * T_low;

        list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
        list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);

        if (std::find(optimization_target_.list_of_target_E_over_R().begin(),
              optimization_target_.list_of_target_E_over_R().end(),
              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_E_over_R().end())
        {
          list_of_min_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
          list_of_max_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
        }
      }

      std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
      for (unsigned int i=0; i < optimization_target_.list_of_target_uncertainty_factors_inf().size(); i++){
        int pos_FallOff_Reaction = std::find(
            indices_of_falloff_reactions.begin(),
            indices_of_falloff_reactions.end(),
            optimization_target_.list_of_target_uncertainty_factors_inf()[i]
            )-indices_of_falloff_reactions.begin();
        // Nominal values of inf parameters
        list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction));
        list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);
        list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction);

        // Min and Max of lnA_inf
        list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-optimization_target_.list_of_uncertainty_factors_inf()[i]));
        list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_inf()[i]));
        if (std::find(optimization_target_.list_of_target_lnA_inf().begin(),
              optimization_target_.list_of_target_lnA_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_lnA_inf().end())
        {
          list_of_min_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
          list_of_max_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
        }

        // Calculating extreme values for Beta
        Beta_1_inf[i] = list_of_nominal_Beta_inf_double[i]+std::log(std::pow(10,optimization_target_.list_of_target_uncertainty_factors_inf()[i])) / std::log(T_high);
        Beta_2_inf[i] = list_of_nominal_Beta_inf_double[i]-std::log(std::pow(10,optimization_target_.list_of_target_uncertainty_factors_inf()[i])) / std::log(T_high);;
        list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
        list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
        if (std::find(optimization_target_.list_of_target_Beta_inf().begin(),
              optimization_target_.list_of_target_Beta_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_Beta_inf().end())
        {
          list_of_min_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
          list_of_max_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
        }

        // Calculting extreame values of E_over_R
        E_over_R_1_inf[i] = list_of_nominal_E_over_R_inf_double[i]-std::log(std::pow(10,optimization_target_.list_of_target_uncertainty_factors_inf()[i])) * T_low;
        E_over_R_2_inf[i] = list_of_nominal_E_over_R_inf_double[i]+std::log(std::pow(10,optimization_target_.list_of_target_uncertainty_factors_inf()[i])) * T_low;
        list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
        list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
        if (std::find(optimization_target_.list_of_target_E_over_R_inf().begin(),
              optimization_target_.list_of_target_E_over_R_inf().end(),
              optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_E_over_R_inf().end()){
          list_of_min_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
          list_of_max_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
        }
      }
    }

    if(optimization_setup_.parameter_boundaries() == "Re-parametrization"){
      OptiSMOKE::ErrorMessage("Compute Boundaries", "Re-Implementation not refactored yet");
    }

    // CLASSIC PLOG - Alpha, Beta, Eps
    for (int i=0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++ ){					
      list_of_nominal_lnA_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
      list_of_min_lnA_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(-optimization_target_.list_of_uncertainty_factors_classic_plog()[i]));
      list_of_max_lnA_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>( optimization_target_.list_of_uncertainty_factors_classic_plog()[i]));
    }

    for (int i=0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++){			
      // the nominal random variable is 0, so that Eps_0 = Esp_0 + D is verified:		
      list_of_nominal_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
      // The minimum and the maximum values of the random variable are then computed as follows:
      list_of_min_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_classic_plog()[i]))*T_low));
      list_of_max_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_classic_plog()[i]))*T_low));
    }

    for (int i=0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++){
      list_of_nominal_Beta_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
      list_of_min_Beta_classic_plog_coefficients_.push_back(
          boost::lexical_cast<std::string>(
            -std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_classic_plog()[i])) / std::log(T_high)
            )
          );
      list_of_max_Beta_classic_plog_coefficients_.push_back(
          boost::lexical_cast<std::string>(
            +std::log(std::pow(10,optimization_target_.list_of_uncertainty_factors_classic_plog()[i])) / std::log(T_high)
            )
          );
    }
  }

  void InputManager::TargetsPreliminaryOptions(){

    name_vec_lnA.resize(optimization_target_.list_of_target_lnA().size());
    for (int i=0; i< optimization_target_.list_of_target_lnA().size(); i++){
      name_vec_lnA[i] = "'lnA_R" + std::to_string(optimization_target_.list_of_target_lnA()[i]) + "'";
      param_name_string_ += name_vec_lnA[i] + " ";
      initial_values_string_ += list_of_initial_lnA_[i] + " ";

      if (optimization_target_.list_of_min_rel_lnA().size()>0){
        lower_bounds_string_ += boost::lexical_cast<std::string>(
            (std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i]-1))) + std::log(optimization_target_.list_of_min_rel_lnA()[i])
            ) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_lnA_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>(
            (std::stod(list_of_initial_lnA_[i]) - std::stod(list_of_min_abs_lnA_[i]))/3
            ) + " ";
      }

      if (optimization_target_.list_of_max_rel_lnA().size()>0){
        upper_bounds_string_ += boost::lexical_cast<std::string>(
            (std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i]-1)))+std::log(optimization_target_.list_of_max_rel_lnA()[i])
            ) + " ";
      }	
      else{
        upper_bounds_string_ += list_of_max_abs_lnA_[i] + " ";
      }
    }

    name_vec_lnA_inf.resize(optimization_target_.list_of_target_lnA_inf().size());
    std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
    for (int i=0; i < optimization_target_.list_of_target_lnA_inf().size(); i++){
      name_vec_lnA_inf[i] = "'lnA_R" + std::to_string(optimization_target_.list_of_target_lnA_inf()[i]) + "_inf'";	
      param_name_string_ += name_vec_lnA_inf[i] + " ";
      initial_values_string_ += list_of_initial_lnA_inf_[i] + " ";
      if (optimization_target_.list_of_min_rel_lnA_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_lnA_inf()[i])-indices_of_falloff_reactions.begin();
        lower_bounds_string_ += boost::lexical_cast<std::string>((std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction)))+std::log(optimization_target_.list_of_min_rel_lnA_inf()[i])) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_lnA_inf_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_initial_lnA_inf_[i]) - std::stod(list_of_min_abs_lnA_inf_[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_rel_lnA_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_lnA_inf()[i])-indices_of_falloff_reactions.begin();
        upper_bounds_string_ += boost::lexical_cast<std::string>((std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction)))+std::log(optimization_target_.list_of_max_rel_lnA_inf()[i])) + " ";
      }
      else{
        upper_bounds_string_ += list_of_max_abs_lnA_inf_[i] + " ";
      }
    }

    name_vec_Beta.resize(optimization_target_.list_of_target_Beta().size());
    for (int i=0; i< optimization_target_.list_of_target_Beta().size(); i++){
      name_vec_Beta[i] = "'Beta_R" + std::to_string(optimization_target_.list_of_target_Beta()[i]) + "'";	
      param_name_string_ += name_vec_Beta[i] + " ";
      initial_values_string_ += list_of_initial_Beta_[i] + " ";

      if (optimization_target_.list_of_min_rel_Beta().size()>0){
        lower_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i]-1))*optimization_target_.list_of_min_rel_Beta()[i]) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_Beta_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_initial_Beta_[i]) - std::stod(list_of_min_abs_Beta_[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_rel_Beta().size()>0){
        upper_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i]-1))*optimization_target_.list_of_max_rel_Beta()[i]) + " ";
      }
      else{
        upper_bounds_string_ += list_of_max_abs_Beta_[i] + " ";
      }
    }

    name_vec_Beta_inf.resize(optimization_target_.list_of_target_Beta_inf().size());
    for (int i=0; i<optimization_target_.list_of_target_Beta_inf().size(); i++){
      name_vec_Beta_inf[i] = "'Beta_R" + std::to_string(optimization_target_.list_of_target_Beta_inf()[i]) + "_inf'";	
      param_name_string_ += name_vec_Beta_inf[i] + " ";
      initial_values_string_ += list_of_initial_Beta_inf_[i] + " ";
      if (optimization_target_.list_of_min_rel_Beta_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_Beta_inf()[i])-indices_of_falloff_reactions.begin();
        lower_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction))*optimization_target_.list_of_min_rel_Beta_inf()[i]) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_Beta_inf_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_initial_Beta_inf_[i]) - std::stod(list_of_min_abs_Beta_inf_[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_rel_Beta_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_Beta_inf()[i])-indices_of_falloff_reactions.begin();
        upper_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction))*optimization_target_.list_of_max_rel_Beta_inf()[i]) + " ";
      }
      else{
        upper_bounds_string_ += list_of_max_abs_Beta_inf_[i] + " ";
      }
    }

    name_vec_E_over_R.resize(optimization_target_.list_of_target_E_over_R().size());
    for (int i=0; i< optimization_target_.list_of_target_E_over_R().size(); i++){
      name_vec_E_over_R[i] = "'E_over_R_R" + std::to_string(optimization_target_.list_of_target_E_over_R()[i]) + "'";
      param_name_string_ += name_vec_E_over_R[i] + " ";
      initial_values_string_ += list_of_initial_E_over_R[i] + " ";
      if (optimization_target_.list_of_min_rel_E_over_R().size()>0){
        lower_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i]-1))*optimization_target_.list_of_min_rel_E_over_R()[i]) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_E_over_R_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_initial_E_over_R[i]) - std::stod(list_of_min_abs_E_over_R_[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_rel_E_over_R().size()>0){
        upper_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i]-1))*optimization_target_.list_of_max_rel_E_over_R()[i]) + " ";
      }
      else{
        upper_bounds_string_ += list_of_max_abs_E_over_R_[i] + " ";
      }
    }

    name_vec_E_over_R_inf.resize(optimization_target_.list_of_target_E_over_R_inf().size());
    for (int i=0; i < optimization_target_.list_of_target_E_over_R_inf().size(); i++){
      name_vec_E_over_R_inf[i] = "'E_over_R_R" + std::to_string(optimization_target_.list_of_target_E_over_R_inf()[i]) + "_inf'";	
      param_name_string_ += name_vec_E_over_R_inf[i] + " ";
      initial_values_string_ += list_of_initial_E_over_R_inf_ [i] + " ";
      if (optimization_target_.list_of_min_rel_E_over_R_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_E_over_R_inf()[i])-indices_of_falloff_reactions.begin();
        lower_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction))*optimization_target_.list_of_min_rel_E_over_R_inf()[i]) + " ";
      }
      else{
        lower_bounds_string_ += list_of_min_abs_E_over_R_inf_[i] + " ";
        std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_initial_E_over_R_inf_[i]) - std::stod(list_of_min_abs_E_over_R_inf_[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_rel_E_over_R_inf().size()>0){
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),optimization_target_.list_of_target_E_over_R_inf()[i])-indices_of_falloff_reactions.begin();
        upper_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction))*optimization_target_.list_of_max_rel_E_over_R_inf()[i]) + " ";
      }
      else{
        upper_bounds_string_ += list_of_max_abs_E_over_R_inf_[i] + " ";
      }
    }

    // third body efficiencies
    name_vec_thirdbody.resize(optimization_target_.list_of_target_thirdbody_reactions().size());
    for (int i=0; i< optimization_target_.list_of_target_thirdbody_reactions().size(); i++){
      name_vec_thirdbody[i] = "'M_R" + std::to_string(optimization_target_.list_of_target_thirdbody_reactions()[i]) + "_" + optimization_target_.list_of_target_thirdbody_species()[i] + "'";
      param_name_string_ += name_vec_thirdbody[i] + " ";
      initial_values_string_ += list_of_initial_thirdbody_eff_[i] + " ";
      if (optimization_target_.list_of_min_abs_thirdbody_eff().size()>0){
        lower_bounds_string_ += optimization_target_.list_of_min_abs_thirdbody_eff()[i] + " ";
      }
      else{
        int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
        lower_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i]-1, iSpecies-1))*optimization_target_.list_of_min_rel_thirdbody_eff()[i]) + " ";
        //std_deviations_string+= boost::lexical_cast<std::string>((boost::lexical_cast<std::double>(list_of_initial_E_over_R_inf[i]) - boost::lexical_cast<std::double>(list_of_min_abs_E_over_R_inf[i]))/3) + " ";
      }

      if (optimization_target_.list_of_max_abs_thirdbody_eff().size()>0){
        upper_bounds_string_ += optimization_target_.list_of_max_abs_thirdbody_eff()[i] + " ";
      }
      else{
        int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
        upper_bounds_string_ += boost::lexical_cast<std::string>((kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i]-1, iSpecies-1))*optimization_target_.list_of_max_rel_thirdbody_eff()[i]) + " ";
      }
    }

    // CLASSIC PLOG REACTIONS
    name_vec_lnA_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
    for (int i=0; i< optimization_target_.list_of_target_classic_plog_reactions().size(); i++){
      name_vec_lnA_classic_plog[i] = "'lnA_classic_PLOG_" + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i]) + "'";
      param_name_string_ += name_vec_lnA_classic_plog[i] + " "; 

      //filling up the strings 
      initial_values_string_ += list_of_nominal_lnA_classic_plog_coefficients_[i] + " ";
      lower_bounds_string_   += list_of_min_lnA_classic_plog_coefficients_[i] + " ";
      upper_bounds_string_   += list_of_max_lnA_classic_plog_coefficients_[i] + " ";
      std_deviations_string_ += boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_classic_plog()[i]/3) + " ";
    }

    name_vec_ER_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
    for (int i=0; i< optimization_target_.list_of_target_classic_plog_reactions().size(); i++){
      name_vec_ER_classic_plog[i] = "'E_over_R_classic_PLOG_" + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i]) + "'";
      param_name_string_ += name_vec_ER_classic_plog[i] + " "; 

      //filling up the strings 
      initial_values_string_ += list_of_nominal_ER_classic_plog_coefficients_[i] + " ";
      lower_bounds_string_   += list_of_min_ER_classic_plog_coefficients_[i]     + " ";
      upper_bounds_string_   += list_of_max_ER_classic_plog_coefficients_[i]     + " ";
      std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_nominal_ER_classic_plog_coefficients_[i]) - std::stod(list_of_min_ER_classic_plog_coefficients_[i]))/3) + " ";
    }

    name_vec_Beta_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
    for (int i=0; i< optimization_target_.list_of_target_classic_plog_reactions().size(); i++){
      name_vec_Beta_classic_plog[i] = "'Beta_classic_PLOG_" + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i]) + "'";
      param_name_string_   += name_vec_Beta_classic_plog[i] + " ";

      initial_values_string_ += list_of_nominal_Beta_classic_plog_coefficients_[i] + " ";
      lower_bounds_string_   += list_of_min_Beta_classic_plog_coefficients_[i]     + " ";
      upper_bounds_string_   += list_of_max_Beta_classic_plog_coefficients_[i]     + " ";
      std_deviations_string_ += boost::lexical_cast<std::string>((std::stod(list_of_nominal_Beta_classic_plog_coefficients_[i]) - std::stod(list_of_min_Beta_classic_plog_coefficients_[i]))/3) + " ";
    }
  }

  void InputManager::ReadExperimentalDataFiles()
  {
    // this entire function is useless remove it 
    // In principle this is useless but pay attention to MPI alloc
    data_manager_.ReadExperimentalData(path_experimental_data_files_);
    dataset_names_ = data_manager_.dataset_names();
    input_paths_ = data_manager_.input_paths();
    solver_name_ = data_manager_.solver_name();
    QoI_ = data_manager_.QoI();
    QoI_target_ = data_manager_.QoI_target();
    multiple_input_ = data_manager_.multiple_input();
    ordinates_label_ = data_manager_.ordinates_label();
    abscissae_label_ = data_manager_.abscissae_label();
    uncertainty_kind_ = data_manager_.uncertainty_kind();
    expdata_x_ = data_manager_.expdata_x();
    expdata_y_ = data_manager_.expdata_y();
    uncertainty_ = data_manager_.uncertainty();
    reactor_mode_ = data_manager_.reactor_mode();
  }

  void InputManager::loadDistributions()
  {
    int tot_simulations = 0;
    for(unsigned int i = 0; i < data_manager_.dataset_names().size(); i++){
      tot_simulations += data_manager_.input_paths()[i].size();
    }

    // Flattening vector of vectors of input files
    for(auto && v : data_manager_.input_paths()){
      flatten_inputs_.insert(flatten_inputs_.end(), v.begin(), v.end());
    }

    if (nprocs_ > tot_simulations)
    {
#ifdef OPTISMOKE_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI

      if (rank_ == 0)
      {
        std::cout << "Warning: " << nprocs_ << " processors were set for reduction, but ";
        std::cout << tot_simulations << " simulations are needed." << std::endl;
        std::cout << nprocs_ - tot_simulations << " processor(s) unused" << std::endl;
      }

#ifdef OPTISMOKE_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI
    }

#ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI

    int avereactors = tot_simulations / nprocs_;
    int extra = tot_simulations % nprocs_;

    // Step1 assign the number of simulations to each process
    n_local_sim_.resize(nprocs_);
    offset_.resize(nprocs_);

    offset_[0] = 0;
    for (unsigned int i = 0; i < nprocs_; i++)
    {
      n_local_sim_[i] = (i < extra) ? avereactors + 1 : avereactors;
      if (i != nprocs_ - 1)
        offset_[i + 1] = offset_[i] + n_local_sim_[i];
    }
    // std::cout << "Proc: " << rank_ << " number of local sim: " << n_local_sim_[rank_] << std::endl;

    local_index_.resize(n_local_sim_[rank_]);
    global_index_.resize(n_local_sim_[rank_]);
    file_index_.resize(n_local_sim_[rank_]);

#ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI

    // Global index è l'indice dell'input nel vettore di tutti gli input
    // devo creare un vettroie con anche l'indice del file degli esperimenti
    // File index mi da l'indice dell'expdata files
    for (unsigned int i = 0; i < n_local_sim_[rank_]; i++)
    {
      local_index_[i] = i;
      global_index_[i] = rank_ + i * nprocs_;
      file_index_[i] =  data_manager_.index_exp_data_file()[global_index_[i]];
    }

    for(unsigned int j=0; j < nprocs_; j++)
    {
#ifdef OPTISMOKE_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI

      if(j == rank_){
        std::cout << "rank: " << rank_ << std::endl;
        for(unsigned int k=0; k < local_index_.size(); k++)
        {
          std::cout << "local: " << local_index_[k] << " global: " << global_index_[k];
          std::cout << " file_id: " << file_index_[k] ;
          std::cout << " input file name: " << flatten_inputs_[global_index_[k]].c_str() << std::endl;
        }
      }

#ifdef OPTISMOKE_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI
    }

#ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif // OPTISMOKE_USE_MPI
       // Step2
       // Step3
       // Step4
       // Step5
       // Step6

       // if(isMaster_){
       // 	std::cout << "Total number of simulations: " << tot_simulations << std::endl;
       // 	std::cout << "Average number: " << avereactors << std::endl;
       // 	std::cout << "Extra: " << extra << std::endl;
       // }
  }

} // namespace OptiSMOKE
