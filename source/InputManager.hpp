/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|             ____        __  _ _____ __  _______  __ __ ______                     |
|            / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____            |
|           / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \           |
|          / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           |
|          \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/            |
|              /_/                                        /_/   /_/                 |
|                                                                                   |
| --------------------------------------------------------------------------------- |
|  Please refer to the copyright statement and license                              |
|  information at the end of this file.                                             |
| --------------------------------------------------------------------------------- |
|                                                                                   |
|           Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                   |
|                    Andrea Bertolino <andrea.bertolino@ulb.be>                     |
|                    Magnus Fürst     <magnus.furst@ulb.ac.be>                      |
|                                                                                   |
|             [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>          |
|                 Department of Chemistry, Materials and Chemical Engineering       |
|                 Politecnico di Milano                                             |
|                 P.zza Leonardo da Vinci 32, 20133 Milano                          |
|                                                                                   |
|             [2] BRITE Research Group <https://brite-research.be>                  |
|                 Brussels Institute for Thermal-fluid systems and clean Energy     |
|                 Avenue F.D. Rooseveltlaan 50                                      |
|                 Bruxelles 1050 Brussel                                            |
|                                                                                   |
\* ------------------------------------------------------------------------------- */

namespace OptiSMOKE {
InputManager::InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary) : dictionary_(dictionary) {
  input_file_name_ = "input.dic";
  main_dictionary_ = "OptiSMOKEpp";
  output_folder_ = "Output";
  kinetics_folder_ = "kinetics";
  optimized_kinetics_folder_ = "Optimized_kinetics";
}

InputManager::~InputManager() {}

void InputManager::SetInputOptions(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("help",
                     "Help Message")("input", po::value<std::string>(), "Input File Path (default: \"input.dic\")");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
  }

  if (vm.count("input")) {
    input_file_name_ = vm["input"].as<std::string>();
  }
}

void InputManager::ReadDictionary() {
  dictionary_.ReadDictionariesFromFile(input_file_name_);
  dictionary_(main_dictionary_).SetGrammar(main_grammar_);
  dictionary_(main_dictionary_).ReadPath("@OutputFolder", output_folder_);
  if (!fs::exists(output_folder_)) {
    fs::create_directories(output_folder_);
  }

  bool iXml = false;
  bool iTransport = false;
  OptiSMOKE::OptionsKinetics kinetics_data_;
  if (dictionary_(main_dictionary_).CheckOption("@KineticsFolder")) {
    iXml = true;
    dictionary_(main_dictionary_).ReadPath("@KineticsFolder", kinetics_folder_);
    if (!fs::exists(kinetics_folder_)) {
      OptiSMOKE::FatalErrorMessage("The @KineticsFolder path does not exists!");
    }
    OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
  } else if (dictionary_(main_dictionary_).CheckOption("@KineticsPreProcessor")) {
    std::string preprocessor_dictionary;
    dictionary_(main_dictionary_).ReadDictionary("@KineticsPreProcessor", preprocessor_dictionary);
    kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary);
    // TODO there is a bug I had not time to investigate further the following lines are a workaround
    if (kinetics_data_.iTransport() == true) {
      iTransport = true;
    }
  } else {
    OptiSMOKE::FatalErrorMessage(
        "Please provide the kinetic mechanism through one of the following keywords: @KineticsFolder | "
        "@KineticsPreProcessor");
  }

  dictionary_(main_dictionary_).ReadOption("@ListOfExperimentalDataFiles", path_experimental_data_files_);
  dictionary_(main_dictionary_).ReadString("@OptimizationLibrary", optimization_library_);

  if (optimization_library_ == "dakota") {
    std::string dakota_dictionary;
    dictionary_(main_dictionary_).ReadDictionary("@DakotaOptions", dakota_dictionary);
    dakota_options_.SetupFromDictionary(dictionary_, dakota_dictionary);
  } else if (optimization_library_ == "nlopt") {
    // dictionary_(main_dictionary_).ReadDictionary("@NLOPTOptions", nlopt_dictionary_);
    // nlopt_options_.SetupFromDictionary(dictionary_, nlopt_dictionary_);
  } else {
    OptiSMOKE::FatalErrorMessage("Unknown optimization library. Available are: dakota | nlopt");
  }

  if (dictionary_(main_dictionary_).CheckOption("@CurveMatchingOptions")) {
    std::string curvematching_dictionary;
    dictionary_(main_dictionary_).ReadDictionary("@CurveMatchingOptions", curvematching_dictionary);
    curvematching_options_.SetupFromDictionary(dictionary_, curvematching_dictionary);
  }

  {
    std::string optimization_setup_dictionary;
    dictionary_(main_dictionary_).ReadDictionary("@OptimizationSetup", optimization_setup_dictionary);
    optimization_setup_.SetupFromDictionary(dictionary_, optimization_setup_dictionary);

    std::string optimization_target_dictionary;
    dictionary_(main_dictionary_).ReadDictionary("@OptimizationTarget", optimization_target_dictionary);
    optimization_target_.SetupFromDictionary(dictionary_, optimization_target_dictionary);
  }


  if (!iXml) {
    if (!iTransport) {
      OpenSMOKE::RapidKineticMechanismWithoutTransport(output_folder_ / kinetics_data_.chemkin_output(),
                                                       kinetics_data_.chemkin_thermodynamics(),
                                                       kinetics_data_.chemkin_kinetics());
    } else {
      OpenSMOKE::RapidKineticMechanismWithTransport(output_folder_ / kinetics_data_.chemkin_output(),
                                                    kinetics_data_.chemkin_transport(),
                                                    kinetics_data_.chemkin_thermodynamics(),
                                                    kinetics_data_.chemkin_kinetics());
    }
  }

  fs::path path_kinetics_output;
  if (!iXml) {  // To be interpreted on-the-fly
    path_kinetics_output = output_folder_ / kinetics_data_.chemkin_output();
  } else {  // Mechanism already provided in XML format
    path_kinetics_output = kinetics_folder_;
  }

  std::cout.setstate(std::ios_base::failbit);  // Disable video output
  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml((path_kinetics_output / "kinetics.xml").string(), ptree);

  thermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
  kineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML_, ptree);
  if (iTransport) {
    transportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
  }

  boost::property_tree::ptree nominal_ptree;
  boost::property_tree::read_xml((path_kinetics_output / "kinetics.xml").string(), nominal_ptree);

  nominalthermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(nominal_ptree);
  nominalkineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*nominalthermodynamicsMapXML_, nominal_ptree);
  if (iTransport) {
    nominaltransportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(nominal_ptree);
  }
  std::cout.clear();  // Re-enable video output
}

// void InputManager::DakotaInputString() {
//   FromTargetToInitialParameter();
//
//   ComputeBoundaries();
//
//   TargetsPreliminaryOptions();
//
//   dakota_input_string_ =
//       " environment,"
//       "\n  tabular_data";
//   dakota_input_string_.append("\n   tabular_data_file '" + output_folder_.string() + "/"
//                               + dakota_options_.tabular_data_file() + "'");
//
//   dakota_input_string_.append("\n method,");
//   dakota_input_string_.append("\n  " + dakota_options_.method());
//   dakota_input_string_.append("\n   max_iterations = " + dakota_options_.max_iterations());
//   dakota_input_string_.append("\n   max_function_evaluations = " + dakota_options_.max_function_evaluations());
//   dakota_input_string_.append("\n   convergence_tolerance = " + dakota_options_.convergence_tolerance());
//   dakota_input_string_.append("\n   solution_target = " + dakota_options_.solution_target());
//   dakota_input_string_.append("\n   seed = " + dakota_options_.seed());
//
//   if (dakota_options_.diverse_input()) {
//     dakota_input_string_.append("\n");
//     for (int i = 0; i < dakota_options_.diverse_dakota_input().size(); i++) {
//       dakota_input_string_.append(" " + dakota_options_.diverse_dakota_input()[i]);
//     }
//   } else if (dakota_options_.method() == "coliny_ea") {
//     dakota_input_string_.append("\n   population_size = " + dakota_options_.population_size());
//     dakota_input_string_.append("\n   fitness_type " + dakota_options_.fitness_type());
//     dakota_input_string_.append("\n   mutation_type " + dakota_options_.mutation_type());
//     dakota_input_string_.append("\n   mutation_rate " + dakota_options_.mutation_rate());
//     dakota_input_string_.append("\n   crossover_type " + dakota_options_.crossover_type());
//     dakota_input_string_.append("\n   crossover_rate " + dakota_options_.crossover_rate());
//     dakota_input_string_.append("\n   replacement_type " + dakota_options_.replacement_type());
//   } else if (dakota_options_.method() == "coliny_direct") {
//     dakota_input_string_.append("\n   division " + dakota_options_.division());
//     dakota_input_string_.append("\n   max_boxsize_limit " + dakota_options_.max_boxsize_limit());
//     dakota_input_string_.append("\n   min_boxsize_limit " + dakota_options_.min_boxsize_limit());
//   } else {
//     OptiSMOKE::FatalErrorMessage("Available methods currently implemented are coliny_ea | coliny_direct");
//   }
//
//   dakota_input_string_.append("\n variables,");
//
//   if (optimization_setup_.parameter_distribution() == "uniform") {
//     dakota_input_string_.append("\n  continuous_design = "
//                                 + std::to_string(optimization_target_.number_of_parameters()));
//     dakota_input_string_.append("\n   descriptors " + param_name_string_);
//     dakota_input_string_.append("\n   initial_point " + initial_values_string_);
//     dakota_input_string_.append("\n   lower_bounds " + lower_bounds_string_);
//     dakota_input_string_.append("\n   upper_bounds " + upper_bounds_string_);
//   } else if (optimization_setup_.parameter_distribution() == "normal") {
//     dakota_input_string_.append("\n  active uncertain ");
//     dakota_input_string_.append("\n  normal_uncertain = "
//                                 + std::to_string(optimization_target_.number_of_parameters()));
//     dakota_input_string_.append("\n   descriptors " + param_name_string_);
//     dakota_input_string_.append("\n   means " + initial_values_string_);
//     dakota_input_string_.append("\n   std_deviations " + std_deviations_string_);
//   }  // Da fare check su consistenza nell' input sul tipo di parameter boundary
//
//   dakota_input_string_.append("\n interface,");
//   dakota_input_string_.append("\n  direct");
//   dakota_input_string_.append("\n  analysis_driver = 'opensmoke_plugin'");
//   dakota_input_string_.append("\n responses,");
//   dakota_input_string_.append("\n  num_objective_functions = 1");
//
//   // Options to use other optimization method (e.g. gradient-based)
//   // Qua forse va messa la possibilità di fare altri tipi di gradienti
//   // accordingly to dakota sicuro lo faccio ora non c'ho voglia
//   if (dakota_options_.dakota_gradient() == true) {
//     dakota_input_string_.append("\n  numerical_gradients");
//     dakota_input_string_.append("\n  method_source dakota");
//     dakota_input_string_.append("\n  interval_type forward");
//     dakota_input_string_.append("\n  fd_step_size = 1.e-5");
//   } else {
//     dakota_input_string_.append("\n  no_gradients");
//   }
//
//   dakota_input_string_.append("\n  no_hessians");
// }

// void InputManager::FromTargetToInitialParameter() {
//   // lnA
//   for (int i = 0; i < optimization_target_.list_of_target_lnA().size(); i++) {
//     list_of_initial_lnA_.push_back(boost::lexical_cast<std::string>(
//         std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i] - 1))));
//   }
//
//   // Beta
//   for (int i = 0; i < optimization_target_.list_of_target_Beta().size(); i++) {
//     list_of_initial_Beta_.push_back(
//         boost::lexical_cast<std::string>(kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i] - 1)));
//   }
//
//   // E_over_R
//   for (int i = 0; i < optimization_target_.list_of_target_E_over_R().size(); i++) {
//     list_of_initial_E_over_R.push_back(boost::lexical_cast<std::string>(
//         kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i] - 1)));
//   }
//
//   // lnA_inf
//   std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
//   for (int i = 0; i < optimization_target_.list_of_target_lnA_inf().size(); i++) {
//     int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                          indices_of_falloff_reactions.end(),
//                                          optimization_target_.list_of_target_lnA_inf()[i])
//                                - indices_of_falloff_reactions.begin();
//     list_of_initial_lnA_inf_.push_back(
//         boost::lexical_cast<std::string>(std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction))));
//   }
//
//   // Beta_inf
//   for (int i = 0; i < optimization_target_.list_of_target_Beta_inf().size(); i++) {
//     int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                          indices_of_falloff_reactions.end(),
//                                          optimization_target_.list_of_target_Beta_inf()[i])
//                                - indices_of_falloff_reactions.begin();
//     list_of_initial_Beta_inf_.push_back(
//         boost::lexical_cast<std::string>(kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction)));
//   }
//
//   // E/R inf
//   for (int i = 0; i < optimization_target_.list_of_target_E_over_R_inf().size(); i++) {
//     int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                          indices_of_falloff_reactions.end(),
//                                          optimization_target_.list_of_target_E_over_R_inf()[i])
//                                - indices_of_falloff_reactions.begin();
//     list_of_initial_E_over_R_inf_.push_back(
//         boost::lexical_cast<std::string>(kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction)));
//   }
//
//   // Third Body Efficiencies
//   for (int i = 0; i < optimization_target_.list_of_target_thirdbody_reactions().size(); i++) {
//     int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
//     list_of_initial_thirdbody_eff_.push_back(boost::lexical_cast<std::string>(
//         kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i] - 1, iSpecies - 1)));
//   }
//
//   // FORD
//   for (int i = 0; i < optimization_target_.list_of_ford().size(); i++) {
//     // These int are 0-based see the minus 1
//     int iReaction = optimization_target_.list_of_ford()[i] - 1;
//     int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_species_ford()[i]) - 1;
//     list_of_initial_ford_.push_back(boost::lexical_cast<std::string>(
//         kineticsMapXML_->stoichiometry().reactionorders_matrix_reactants().coeff(iReaction, iSpecies)));
//   }
//
//   // RORD
//   for (int i = 0; i < optimization_target_.list_of_rord().size(); i++) {
//     int iReaction = optimization_target_.list_of_rord()[i] - 1;
//     int iSpecies = thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_species_rord()[i]) - 1;
//     list_of_initial_rord_.push_back(boost::lexical_cast<std::string>(
//         kineticsMapXML_->stoichiometry().reactionorders_matrix_products().coeff(iReaction, iSpecies)));
//   }
// }

// void InputManager::ComputeBoundaries() {
//   double T_low = 300;
//   double T_high = 2500;
//
//   // Initialize needed values at the specific size
//   std::vector<double> list_of_nominal_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> list_of_nominal_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double>
//   list_of_nominal_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
//
//   std::vector<double> list_of_min_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> list_of_max_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> list_of_min_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> list_of_max_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double>
//   list_of_min_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double>
//   list_of_max_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
//
//   std::vector<double> kappa_lower_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> kappa_upper_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> kappa_lower_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> kappa_upper_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> Beta_1(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> Beta_2(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> lnA_1(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> lnA_2(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> E_over_R_1(optimization_target_.list_of_target_uncertainty_factors().size());
//   std::vector<double> E_over_R_2(optimization_target_.list_of_target_uncertainty_factors().size());
//
//   std::vector<double> list_of_nominal_lnA_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_nominal_Beta_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_nominal_E_over_R_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//
//   std::vector<double> list_of_min_abs_lnA_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_max_abs_lnA_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_min_abs_Beta_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_max_abs_Beta_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_min_abs_E_over_R_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> list_of_max_abs_E_over_R_inf_double(
//       optimization_target_.list_of_target_uncertainty_factors_inf().size());
//
//   std::vector<double> kappa_lower_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> kappa_upper_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> kappa_lower_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> kappa_upper_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> Beta_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> Beta_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> lnA_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> lnA_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> E_over_R_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//   std::vector<double> E_over_R_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
//
//   if (optimization_setup_.parameter_boundaries() == "Furst") {
//     for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++) {
//       // Nominal values of parameters
//       list_of_nominal_lnA_double[i] =
//           std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i] - 1));
//       list_of_nominal_Beta_double[i] =
//           nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i] - 1);
//       list_of_nominal_E_over_R_double[i] =
//           nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i] - 1);
//
//       // Min and Max of lnA
//       list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]
//                                       + std::log(std::pow(10,
//                                       -optimization_target_.list_of_uncertainty_factors()[i]));
//       list_of_max_abs_lnA_double[i] =
//           list_of_nominal_lnA_double[i] + std::log(std::pow(10,
//           optimization_target_.list_of_uncertainty_factors()[i]));
//       if (std::find(optimization_target_.list_of_target_lnA().begin(),
//                     optimization_target_.list_of_target_lnA().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_lnA().end()) {
//         list_of_min_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
//         list_of_max_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
//       }
//
//       // Limiting values for the rate coefficient
//       kappa_lower_T_low[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i] * std::log(T_low)
//                              - list_of_nominal_E_over_R_double[i] * std::pow(T_low, -1);
//       kappa_upper_T_low[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i] * std::log(T_low)
//                              - list_of_nominal_E_over_R_double[i] * std::pow(T_low, -1);
//       kappa_lower_T_high[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i] * std::log(T_high)
//                               - list_of_nominal_E_over_R_double[i] * std::pow(T_high, -1);
//       kappa_upper_T_high[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i] * std::log(T_high)
//                               - list_of_nominal_E_over_R_double[i] * std::pow(T_high, -1);
//
//       // Calculating extreme values for Beta
//       Beta_1[i] =
//           (kappa_upper_T_low[i] - kappa_lower_T_high[i] - list_of_nominal_E_over_R_double[i] * (1 / T_high - 1 /
//           T_low)) / (std::log(T_low) - std::log(T_high));
//       Beta_2[i] =
//           (kappa_lower_T_low[i] - kappa_upper_T_high[i] - list_of_nominal_E_over_R_double[i] * (1 / T_high - 1 /
//           T_low)) / (std::log(T_low) - std::log(T_high));
//
//       list_of_min_abs_Beta_double[i] = std::min(Beta_1[i], Beta_2[i]);
//       list_of_max_abs_Beta_double[i] = std::max(Beta_1[i], Beta_2[i]);
//
//       if (std::find(optimization_target_.list_of_target_Beta().begin(),
//                     optimization_target_.list_of_target_Beta().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_Beta().end()) {
//         list_of_min_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
//         list_of_max_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
//       }
//
//       // Calculting extreame values of E_over_R
//       lnA_1[i] = (kappa_lower_T_high[i] - (T_low / T_high) * kappa_upper_T_low[i]
//                   - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low / T_high) * std::log(T_low)))
//                  / (1 - (T_low / T_high));
//       E_over_R_1[i] =
//           lnA_1[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_upper_T_low[i] * T_low;
//       lnA_2[i] = (kappa_upper_T_high[i] - (T_low / T_high) * kappa_lower_T_low[i]
//                   - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low / T_high) * std::log(T_low)))
//                  / (1 - (T_low / T_high));
//       E_over_R_2[i] =
//           lnA_2[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_lower_T_low[i] * T_low;
//       list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i], E_over_R_2[i]);
//       list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i], E_over_R_2[i]);
//
//       if (std::find(optimization_target_.list_of_target_E_over_R().begin(),
//                     optimization_target_.list_of_target_E_over_R().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_E_over_R().end()) {
//         list_of_min_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
//         list_of_max_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
//       }
//     }
//
//     std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
//
//     for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors_inf().size(); i++) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//
//       // Nominal values of inf parameters
//       list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction));
//       list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);
//       list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction);
//
//       // Min and Max of lnA_inf
//       list_of_min_abs_lnA_inf_double[i] =
//           list_of_nominal_lnA_inf_double[i]
//           + std::log(std::pow(10, -optimization_target_.list_of_uncertainty_factors_inf()[i]));
//       list_of_max_abs_lnA_inf_double[i] =
//           list_of_nominal_lnA_inf_double[i]
//           + std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_inf()[i]));
//       if (std::find(optimization_target_.list_of_target_lnA_inf().begin(),
//                     optimization_target_.list_of_target_lnA_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_lnA_inf().end()) {
//         list_of_min_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
//         list_of_max_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
//       }
//
//       // Limiting values for the rate coefficient
//       kappa_lower_T_low_inf[i] = list_of_min_abs_lnA_inf_double[i]
//                                  + list_of_nominal_Beta_inf_double[i] * std::log(T_low)
//                                  - list_of_nominal_E_over_R_inf_double[i] * std::pow(T_low, -1);
//       kappa_upper_T_low_inf[i] = list_of_max_abs_lnA_inf_double[i]
//                                  + list_of_nominal_Beta_inf_double[i] * std::log(T_low)
//                                  - list_of_nominal_E_over_R_inf_double[i] * std::pow(T_low, -1);
//       kappa_lower_T_high_inf[i] = list_of_min_abs_lnA_inf_double[i]
//                                   + list_of_nominal_Beta_inf_double[i] * std::log(T_high)
//                                   - list_of_nominal_E_over_R_inf_double[i] * std::pow(T_high, -1);
//       kappa_upper_T_high_inf[i] = list_of_max_abs_lnA_inf_double[i]
//                                   + list_of_nominal_Beta_inf_double[i] * std::log(T_high)
//                                   - list_of_nominal_E_over_R_inf_double[i] * std::pow(T_high, -1);
//
//       Beta_1_inf[i] = (kappa_upper_T_low_inf[i] - kappa_lower_T_high_inf[i]
//                        - list_of_nominal_E_over_R_inf_double[i] * (1 / T_high - 1 / T_low))
//                       / (std::log(T_low) - std::log(T_high));
//       Beta_2_inf[i] = (kappa_lower_T_low_inf[i] - kappa_upper_T_high_inf[i]
//                        - list_of_nominal_E_over_R_inf_double[i] * (1 / T_high - 1 / T_low))
//                       / (std::log(T_low) - std::log(T_high));
//       list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i], Beta_2_inf[i]);
//       list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i], Beta_2_inf[i]);
//       if (std::find(optimization_target_.list_of_target_Beta_inf().begin(),
//                     optimization_target_.list_of_target_Beta_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_Beta_inf().end()) {
//         list_of_min_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
//         list_of_max_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
//       }
//
//       // Calculting extreame values of E_over_R
//       lnA_1_inf[i] = (kappa_lower_T_high_inf[i] - (T_low / T_high) * kappa_upper_T_low_inf[i]
//                       - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low / T_high) * std::log(T_low)))
//                      / (1 - (T_low / T_high));
//       E_over_R_1_inf[i] = lnA_1_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low)
//                           - kappa_upper_T_low_inf[i] * T_low;
//       lnA_2_inf[i] = (kappa_upper_T_high_inf[i] - (T_low / T_high) * kappa_lower_T_low_inf[i]
//                       - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low / T_high) * std::log(T_low)))
//                      / (1 - (T_low / T_high));
//       E_over_R_2_inf[i] = lnA_2_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low)
//                           - kappa_lower_T_low_inf[i] * T_low;
//       list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i], E_over_R_2_inf[i]);
//       list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i], E_over_R_2_inf[i]);
//       if (std::find(optimization_target_.list_of_target_E_over_R_inf().begin(),
//                     optimization_target_.list_of_target_E_over_R_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_E_over_R_inf().end()) {
//         list_of_min_abs_E_over_R_inf_.push_back(
//             boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
//         list_of_max_abs_E_over_R_inf_.push_back(
//             boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
//       }
//     }
//   }
//
//   if (optimization_setup_.parameter_boundaries() == "Narrow") {
//     for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++) {
//       list_of_nominal_lnA_double[i] =
//           std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i] - 1));
//       list_of_nominal_Beta_double[i] =
//           nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i] - 1);
//       list_of_nominal_E_over_R_double[i] =
//           nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i] - 1);
//
//       list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]
//                                       + std::log(std::pow(10,
//                                       -optimization_target_.list_of_uncertainty_factors()[i]));
//       list_of_max_abs_lnA_double[i] =
//           list_of_nominal_lnA_double[i] + std::log(std::pow(10,
//           optimization_target_.list_of_uncertainty_factors()[i]));
//
//       if (std::find(optimization_target_.list_of_target_lnA().begin(),
//                     optimization_target_.list_of_target_lnA().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_lnA().end()) {
//         list_of_min_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
//         list_of_max_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
//       }
//
//       Beta_1[i] = list_of_nominal_Beta_double[i]
//                   + std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i])) / std::log(T_high);
//       Beta_2[i] = list_of_nominal_Beta_double[i]
//                   - std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i])) / std::log(T_high);
//
//       list_of_min_abs_Beta_double[i] = std::min(Beta_1[i], Beta_2[i]);
//       list_of_max_abs_Beta_double[i] = std::max(Beta_1[i], Beta_2[i]);
//
//       if (std::find(optimization_target_.list_of_target_Beta().begin(),
//                     optimization_target_.list_of_target_Beta().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_Beta().end()) {
//         list_of_min_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
//         list_of_max_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
//       }
//
//       E_over_R_1[i] = list_of_nominal_E_over_R_double[i]
//                       - std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i])) * T_low;
//       E_over_R_2[i] = list_of_nominal_E_over_R_double[i]
//                       + std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors()[i])) * T_low;
//
//       list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i], E_over_R_2[i]);
//       list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i], E_over_R_2[i]);
//
//       if (std::find(optimization_target_.list_of_target_E_over_R().begin(),
//                     optimization_target_.list_of_target_E_over_R().end(),
//                     optimization_target_.list_of_target_uncertainty_factors()[i])
//           != optimization_target_.list_of_target_E_over_R().end()) {
//         list_of_min_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
//         list_of_max_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
//       }
//     }
//
//     std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
//     for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors_inf().size(); i++) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       // Nominal values of inf parameters
//       list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction));
//       list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);
//       list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction);
//
//       // Min and Max of lnA_inf
//       list_of_min_abs_lnA_inf_double[i] =
//           list_of_nominal_lnA_inf_double[i]
//           + std::log(std::pow(10, -optimization_target_.list_of_uncertainty_factors_inf()[i]));
//       list_of_max_abs_lnA_inf_double[i] =
//           list_of_nominal_lnA_inf_double[i]
//           + std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_inf()[i]));
//       if (std::find(optimization_target_.list_of_target_lnA_inf().begin(),
//                     optimization_target_.list_of_target_lnA_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_lnA_inf().end()) {
//         list_of_min_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
//         list_of_max_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
//       }
//
//       // Calculating extreme values for Beta
//       Beta_1_inf[i] =
//           list_of_nominal_Beta_inf_double[i]
//           + std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors_inf()[i])) /
//           std::log(T_high);
//       Beta_2_inf[i] =
//           list_of_nominal_Beta_inf_double[i]
//           - std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors_inf()[i])) /
//           std::log(T_high);
//       ;
//       list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i], Beta_2_inf[i]);
//       list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i], Beta_2_inf[i]);
//       if (std::find(optimization_target_.list_of_target_Beta_inf().begin(),
//                     optimization_target_.list_of_target_Beta_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_Beta_inf().end()) {
//         list_of_min_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
//         list_of_max_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
//       }
//
//       // Calculting extreame values of E_over_R
//       E_over_R_1_inf[i] =
//           list_of_nominal_E_over_R_inf_double[i]
//           - std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors_inf()[i])) * T_low;
//       E_over_R_2_inf[i] =
//           list_of_nominal_E_over_R_inf_double[i]
//           + std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors_inf()[i])) * T_low;
//
//       list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i], E_over_R_2_inf[i]);
//       list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i], E_over_R_2_inf[i]);
//       if (std::find(optimization_target_.list_of_target_E_over_R_inf().begin(),
//                     optimization_target_.list_of_target_E_over_R_inf().end(),
//                     optimization_target_.list_of_target_uncertainty_factors_inf()[i])
//           != optimization_target_.list_of_target_E_over_R_inf().end()) {
//         list_of_min_abs_E_over_R_inf_.push_back(
//             boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
//         list_of_max_abs_E_over_R_inf_.push_back(
//             boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
//       }
//     }
//   }
//
//   if (optimization_setup_.parameter_boundaries() == "Re-parametrization") {
//     OptiSMOKE::ErrorMessage("Compute Boundaries", "Re-Implementation not refactored yet");
//   }
//
//   // CLASSIC PLOG - Alpha, Beta, Eps
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     list_of_nominal_lnA_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
//     list_of_min_lnA_classic_plog_coefficients_.push_back(
//         boost::lexical_cast<std::string>(-optimization_target_.list_of_uncertainty_factors_classic_plog()[i]));
//     list_of_max_lnA_classic_plog_coefficients_.push_back(
//         boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_classic_plog()[i]));
//   }
//
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     // the nominal random variable is 0, so that Eps_0 = Esp_0 + D is verified:
//     list_of_nominal_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
//     // The minimum and the maximum values of the random variable are then computed as follows:
//     list_of_min_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(
//         -std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_classic_plog()[i])) * T_low));
//     list_of_max_ER_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(
//         +std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_classic_plog()[i])) * T_low));
//   }
//
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     list_of_nominal_Beta_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(0));
//     list_of_min_Beta_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(
//         -std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_classic_plog()[i]))
//         / std::log(T_high)));
//     list_of_max_Beta_classic_plog_coefficients_.push_back(boost::lexical_cast<std::string>(
//         +std::log(std::pow(10, optimization_target_.list_of_uncertainty_factors_classic_plog()[i]))
//         / std::log(T_high)));
//   }
// }

// void InputManager::TargetsPreliminaryOptions() {
//   name_vec_lnA.resize(optimization_target_.list_of_target_lnA().size());
//   for (int i = 0; i < optimization_target_.list_of_target_lnA().size(); i++) {
//     name_vec_lnA[i] = "'lnA_R" + std::to_string(optimization_target_.list_of_target_lnA()[i]) + "'";
//     param_name_string_ += name_vec_lnA[i] + " ";
//     initial_values_string_ += list_of_initial_lnA_[i] + " ";
//
//     if (optimization_target_.list_of_min_rel_lnA().size() > 0) {
//       lower_bounds_string_ += boost::lexical_cast<std::string>(
//                                   (std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i] - 1)))
//                                   + std::log(optimization_target_.list_of_min_rel_lnA()[i]))
//                               + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_lnA_[i] + " ";
//       std_deviations_string_ += boost::lexical_cast<std::string>(
//                                     (std::stod(list_of_initial_lnA_[i]) - std::stod(list_of_min_abs_lnA_[i])) / 3)
//                                 + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_lnA().size() > 0) {
//       upper_bounds_string_ += boost::lexical_cast<std::string>(
//                                   (std::log(kineticsMapXML_->A(optimization_target_.list_of_target_lnA()[i] - 1)))
//                                   + std::log(optimization_target_.list_of_max_rel_lnA()[i]))
//                               + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_lnA_[i] + " ";
//     }
//   }
//
//   name_vec_lnA_inf.resize(optimization_target_.list_of_target_lnA_inf().size());
//   std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();
//   for (int i = 0; i < optimization_target_.list_of_target_lnA_inf().size(); i++) {
//     name_vec_lnA_inf[i] = "'lnA_R" + std::to_string(optimization_target_.list_of_target_lnA_inf()[i]) + "_inf'";
//     param_name_string_ += name_vec_lnA_inf[i] + " ";
//     initial_values_string_ += list_of_initial_lnA_inf_[i] + " ";
//     if (optimization_target_.list_of_min_rel_lnA_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_lnA_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       lower_bounds_string_ +=
//           boost::lexical_cast<std::string>((std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction)))
//                                            + std::log(optimization_target_.list_of_min_rel_lnA_inf()[i]))
//           + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_lnA_inf_[i] + " ";
//       std_deviations_string_ +=
//           boost::lexical_cast<std::string>(
//               (std::stod(list_of_initial_lnA_inf_[i]) - std::stod(list_of_min_abs_lnA_inf_[i])) / 3)
//           + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_lnA_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_lnA_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       upper_bounds_string_ +=
//           boost::lexical_cast<std::string>((std::log(kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction)))
//                                            + std::log(optimization_target_.list_of_max_rel_lnA_inf()[i]))
//           + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_lnA_inf_[i] + " ";
//     }
//   }
//
//   name_vec_Beta.resize(optimization_target_.list_of_target_Beta().size());
//   for (int i = 0; i < optimization_target_.list_of_target_Beta().size(); i++) {
//     name_vec_Beta[i] = "'Beta_R" + std::to_string(optimization_target_.list_of_target_Beta()[i]) + "'";
//     param_name_string_ += name_vec_Beta[i] + " ";
//     initial_values_string_ += list_of_initial_Beta_[i] + " ";
//
//     if (optimization_target_.list_of_min_rel_Beta().size() > 0) {
//       lower_bounds_string_ +=
//           boost::lexical_cast<std::string>((kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i] - 1))
//                                            * optimization_target_.list_of_min_rel_Beta()[i])
//           + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_Beta_[i] + " ";
//       std_deviations_string_ += boost::lexical_cast<std::string>(
//                                     (std::stod(list_of_initial_Beta_[i]) - std::stod(list_of_min_abs_Beta_[i])) / 3)
//                                 + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_Beta().size() > 0) {
//       upper_bounds_string_ +=
//           boost::lexical_cast<std::string>((kineticsMapXML_->Beta(optimization_target_.list_of_target_Beta()[i] - 1))
//                                            * optimization_target_.list_of_max_rel_Beta()[i])
//           + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_Beta_[i] + " ";
//     }
//   }
//
//   name_vec_Beta_inf.resize(optimization_target_.list_of_target_Beta_inf().size());
//   for (int i = 0; i < optimization_target_.list_of_target_Beta_inf().size(); i++) {
//     name_vec_Beta_inf[i] = "'Beta_R" + std::to_string(optimization_target_.list_of_target_Beta_inf()[i]) + "_inf'";
//     param_name_string_ += name_vec_Beta_inf[i] + " ";
//     initial_values_string_ += list_of_initial_Beta_inf_[i] + " ";
//     if (optimization_target_.list_of_min_rel_Beta_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_Beta_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       lower_bounds_string_ +=
//       boost::lexical_cast<std::string>((kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction))
//                                                                * optimization_target_.list_of_min_rel_Beta_inf()[i])
//                               + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_Beta_inf_[i] + " ";
//       std_deviations_string_ +=
//           boost::lexical_cast<std::string>(
//               (std::stod(list_of_initial_Beta_inf_[i]) - std::stod(list_of_min_abs_Beta_inf_[i])) / 3)
//           + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_Beta_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_Beta_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       upper_bounds_string_ +=
//       boost::lexical_cast<std::string>((kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction))
//                                                                * optimization_target_.list_of_max_rel_Beta_inf()[i])
//                               + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_Beta_inf_[i] + " ";
//     }
//   }
//
//   name_vec_E_over_R.resize(optimization_target_.list_of_target_E_over_R().size());
//   for (int i = 0; i < optimization_target_.list_of_target_E_over_R().size(); i++) {
//     name_vec_E_over_R[i] = "'E_over_R_R" + std::to_string(optimization_target_.list_of_target_E_over_R()[i]) + "'";
//     param_name_string_ += name_vec_E_over_R[i] + " ";
//     initial_values_string_ += list_of_initial_E_over_R[i] + " ";
//     if (optimization_target_.list_of_min_rel_E_over_R().size() > 0) {
//       lower_bounds_string_ += boost::lexical_cast<std::string>(
//                                   (kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i] - 1))
//                                   * optimization_target_.list_of_min_rel_E_over_R()[i])
//                               + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_E_over_R_[i] + " ";
//       std_deviations_string_ +=
//           boost::lexical_cast<std::string>(
//               (std::stod(list_of_initial_E_over_R[i]) - std::stod(list_of_min_abs_E_over_R_[i])) / 3)
//           + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_E_over_R().size() > 0) {
//       upper_bounds_string_ += boost::lexical_cast<std::string>(
//                                   (kineticsMapXML_->E_over_R(optimization_target_.list_of_target_E_over_R()[i] - 1))
//                                   * optimization_target_.list_of_max_rel_E_over_R()[i])
//                               + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_E_over_R_[i] + " ";
//     }
//   }
//
//   name_vec_E_over_R_inf.resize(optimization_target_.list_of_target_E_over_R_inf().size());
//   for (int i = 0; i < optimization_target_.list_of_target_E_over_R_inf().size(); i++) {
//     name_vec_E_over_R_inf[i] =
//         "'E_over_R_R" + std::to_string(optimization_target_.list_of_target_E_over_R_inf()[i]) + "_inf'";
//     param_name_string_ += name_vec_E_over_R_inf[i] + " ";
//     initial_values_string_ += list_of_initial_E_over_R_inf_[i] + " ";
//     if (optimization_target_.list_of_min_rel_E_over_R_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_E_over_R_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       lower_bounds_string_ +=
//           boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction))
//                                            * optimization_target_.list_of_min_rel_E_over_R_inf()[i])
//           + " ";
//     } else {
//       lower_bounds_string_ += list_of_min_abs_E_over_R_inf_[i] + " ";
//       std_deviations_string_ +=
//           boost::lexical_cast<std::string>(
//               (std::stod(list_of_initial_E_over_R_inf_[i]) - std::stod(list_of_min_abs_E_over_R_inf_[i])) / 3)
//           + " ";
//     }
//
//     if (optimization_target_.list_of_max_rel_E_over_R_inf().size() > 0) {
//       int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),
//                                            indices_of_falloff_reactions.end(),
//                                            optimization_target_.list_of_target_E_over_R_inf()[i])
//                                  - indices_of_falloff_reactions.begin();
//       upper_bounds_string_ +=
//           boost::lexical_cast<std::string>((kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction))
//                                            * optimization_target_.list_of_max_rel_E_over_R_inf()[i])
//           + " ";
//     } else {
//       upper_bounds_string_ += list_of_max_abs_E_over_R_inf_[i] + " ";
//     }
//   }
//
//   // third body efficiencies
//   name_vec_thirdbody.resize(optimization_target_.list_of_target_thirdbody_reactions().size());
//   for (int i = 0; i < optimization_target_.list_of_target_thirdbody_reactions().size(); i++) {
//     name_vec_thirdbody[i] = "'M_R" + std::to_string(optimization_target_.list_of_target_thirdbody_reactions()[i]) +
//     "_"
//                             + optimization_target_.list_of_target_thirdbody_species()[i] + "'";
//     param_name_string_ += name_vec_thirdbody[i] + " ";
//     initial_values_string_ += list_of_initial_thirdbody_eff_[i] + " ";
//     if (optimization_target_.list_of_min_abs_thirdbody_eff().size() > 0) {
//       lower_bounds_string_ += optimization_target_.list_of_min_abs_thirdbody_eff()[i] + " ";
//     } else {
//       int iSpecies =
//       thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
//       lower_bounds_string_ +=
//           boost::lexical_cast<std::string>(
//               (kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i] - 1,
//                                           iSpecies - 1))
//               * optimization_target_.list_of_min_rel_thirdbody_eff()[i])
//           + " ";
//       // std_deviations_string+=
//       // boost::lexical_cast<std::string>((boost::lexical_cast<std::double>(list_of_initial_E_over_R_inf[i]) -
//       // boost::lexical_cast<std::double>(list_of_min_abs_E_over_R_inf[i]))/3) + " ";
//     }
//
//     if (optimization_target_.list_of_max_abs_thirdbody_eff().size() > 0) {
//       upper_bounds_string_ += optimization_target_.list_of_max_abs_thirdbody_eff()[i] + " ";
//     } else {
//       int iSpecies =
//       thermodynamicsMapXML_->IndexOfSpecies(optimization_target_.list_of_target_thirdbody_species()[i]);
//       upper_bounds_string_ +=
//           boost::lexical_cast<std::string>(
//               (kineticsMapXML_->ThirdBody(optimization_target_.list_of_target_thirdbody_reactions()[i] - 1,
//                                           iSpecies - 1))
//               * optimization_target_.list_of_max_rel_thirdbody_eff()[i])
//           + " ";
//     }
//   }
//
//   // CLASSIC PLOG REACTIONS
//   name_vec_lnA_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     name_vec_lnA_classic_plog[i] =
//         "'lnA_classic_PLOG_" + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i]) + "'";
//     param_name_string_ += name_vec_lnA_classic_plog[i] + " ";
//
//     // filling up the strings
//     initial_values_string_ += list_of_nominal_lnA_classic_plog_coefficients_[i] + " ";
//     lower_bounds_string_ += list_of_min_lnA_classic_plog_coefficients_[i] + " ";
//     upper_bounds_string_ += list_of_max_lnA_classic_plog_coefficients_[i] + " ";
//     std_deviations_string_ +=
//         boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_classic_plog()[i] / 3) + "
//         ";
//   }
//
//   name_vec_ER_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     name_vec_ER_classic_plog[i] = "'E_over_R_classic_PLOG_"
//                                   + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i])
//                                   + "'";
//     param_name_string_ += name_vec_ER_classic_plog[i] + " ";
//
//     // filling up the strings
//     initial_values_string_ += list_of_nominal_ER_classic_plog_coefficients_[i] + " ";
//     lower_bounds_string_ += list_of_min_ER_classic_plog_coefficients_[i] + " ";
//     upper_bounds_string_ += list_of_max_ER_classic_plog_coefficients_[i] + " ";
//     std_deviations_string_ +=
//         boost::lexical_cast<std::string>((std::stod(list_of_nominal_ER_classic_plog_coefficients_[i])
//                                           - std::stod(list_of_min_ER_classic_plog_coefficients_[i]))
//                                          / 3)
//         + " ";
//   }
//
//   name_vec_Beta_classic_plog.resize(optimization_target_.list_of_target_classic_plog_reactions().size());
//   for (int i = 0; i < optimization_target_.list_of_target_classic_plog_reactions().size(); i++) {
//     name_vec_Beta_classic_plog[i] =
//         "'Beta_classic_PLOG_" + std::to_string(optimization_target_.list_of_target_classic_plog_reactions()[i]) +
//         "'";
//     param_name_string_ += name_vec_Beta_classic_plog[i] + " ";
//
//     initial_values_string_ += list_of_nominal_Beta_classic_plog_coefficients_[i] + " ";
//     lower_bounds_string_ += list_of_min_Beta_classic_plog_coefficients_[i] + " ";
//     upper_bounds_string_ += list_of_max_Beta_classic_plog_coefficients_[i] + " ";
//     std_deviations_string_ +=
//         boost::lexical_cast<std::string>((std::stod(list_of_nominal_Beta_classic_plog_coefficients_[i])
//                                           - std::stod(list_of_min_Beta_classic_plog_coefficients_[i]))
//                                          / 3)
//         + " ";
//   }
//
//   // RPBMR REACTIONS
//   // name_vec_lnA_rpbmr.resize(optimization_target_.list_of_target_rpbmr_reactions().size());
//   // for (int i = 0; i < optimization_target_.list_of_target_rpbmr_reactions().size(); i++) {
//   //   name_vec_lnA_rpbmr[i] =
//   //       "'lnA_RPBMR_" + std::to_string(optimization_target_.list_of_target_rpbmr_reactions()[i]) + "'";
//   //   param_name_string_ += name_vec_lnA_rpbmr[i] + " ";
//   //
//   //   // filling up the strings
//   //   initial_values_string_ += list_of_nominal_lnA_rpbmr_coefficients_[i] + " ";
//   //   lower_bounds_string_ += list_of_min_lnA_rpbmr_coefficients_[i] + " ";
//   //   upper_bounds_string_ += list_of_max_lnA_rpbmr_coefficients_[i] + " ";
//   //   std_deviations_string_ +=
//   //       boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_rpbmr()[i] / 3) + " ";
//   // }
//   //
//   // name_vec_ER_rpbmr.resize(optimization_target_.list_of_target_rpbmr_reactions().size());
//   // for (int i = 0; i < optimization_target_.list_of_target_rpbmr_reactions().size(); i++) {
//   //   name_vec_ER_rpbmr[i] =
//   //       "'ER_RPBMR_" + std::to_string(optimization_target_.list_of_target_rpbmr_reactions()[i]) + "'";
//   //   param_name_string_ += name_vec_ER_rpbmr[i] + " ";
//   //
//   //   // filling up the strings
//   //   initial_values_string_ += list_of_nominal_E_over_R_rpbmr_coefficients_[i] + " ";
//   //   lower_bounds_string_ += list_of_min_E_over_R_rpbmr_coefficients_[i] + " ";
//   //   upper_bounds_string_ += list_of_max_E_over_R_rpbmr_coefficients_[i] + " ";
//   //   std_deviations_string_ +=
//   //       boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_rpbmr()[i] / 3) + " ";
//   // }
//   //
//   // name_vec_Beta_rpbmr.resize(optimization_target_.list_of_target_rpbmr_reactions().size());
//   // for (int i = 0; i < optimization_target_.list_of_target_rpbmr_reactions().size(); i++) {
//   //   name_vec_Beta_rpbmr[i] =
//   //       "'Beta_RPBMR_" + std::to_string(optimization_target_.list_of_target_rpbmr_reactions()[i]) + "'";
//   //   param_name_string_ += name_vec_Beta_rpbmr[i] + " ";
//   //
//   //   // filling up the strings
//   //   initial_values_string_ += list_of_nominal_Beta_rpbmr_coefficients_[i] + " ";
//   //   lower_bounds_string_ += list_of_min_Beta_rpbmr_coefficients_[i] + " ";
//   //   upper_bounds_string_ += list_of_max_Beta_rpbmr_coefficients_[i] + " ";
//   //   std_deviations_string_ +=
//   //       boost::lexical_cast<std::string>(optimization_target_.list_of_uncertainty_factors_rpbmr()[i] / 3) + " ";
//   // }
//
//   // FORD
//   name_vec_ford.resize(optimization_target_.list_of_ford().size());
//   for (int i = 0; i < optimization_target_.list_of_ford().size(); i++) {
//     name_vec_ford[i] = "'FORD_R" + std::to_string(optimization_target_.list_of_ford()[i]) + "_"
//                        + optimization_target_.list_of_species_ford()[i] + "'";
//     param_name_string_ += name_vec_ford[i] + " ";
//     initial_values_string_ += list_of_initial_ford_[i] + " ";
//     if (optimization_target_.list_of_min_abs_FORD().size() > 0) {
//       lower_bounds_string_ += boost::lexical_cast<std::string>(optimization_target_.list_of_min_abs_FORD()[i]) + " ";
//     } else {
//       // TODO IMPLEMENT RELATIVE CHANGES
//       OptiSMOKE::FatalErrorMessage("No relative changes implemented for FORD!");
//     }
//
//     if (optimization_target_.list_of_max_abs_FORD().size() > 0) {
//       upper_bounds_string_ += boost::lexical_cast<std::string>(optimization_target_.list_of_max_abs_FORD()[i]) + " ";
//     } else {
//       OptiSMOKE::FatalErrorMessage("No relative changes implementes for FORD!");
//     }
//   }
//
//   // RORD
//   name_vec_rord.resize(optimization_target_.list_of_rord().size());
//   for (int i = 0; i < optimization_target_.list_of_rord().size(); i++) {
//     name_vec_rord[i] = "'RORD_R" + std::to_string(optimization_target_.list_of_rord()[i]) + "_"
//                        + optimization_target_.list_of_species_rord()[i] + "'";
//     param_name_string_ += name_vec_rord[i] + " ";
//     initial_values_string_ += list_of_initial_rord_[i] + " ";
//     if (optimization_target_.list_of_min_abs_RORD().size() > 0) {
//       lower_bounds_string_ += boost::lexical_cast<std::string>(optimization_target_.list_of_min_abs_RORD()[i]) + " ";
//     } else {
//       // TODO IMPLEMENT RELATIVE CHANGES
//       OptiSMOKE::FatalErrorMessage("No relative changes implemented for RORD!");
//     }
//
//     if (optimization_target_.list_of_max_abs_RORD().size() > 0) {
//       upper_bounds_string_ += boost::lexical_cast<std::string>(optimization_target_.list_of_max_abs_RORD()[i]) + " ";
//     } else {
//       OptiSMOKE::FatalErrorMessage("No relative changes implementes for RORD!");
//     }
//   }
// }
}  // namespace OptiSMOKE


// void InputManager::SetUpNLOPT() {
//   FromTargetToInitialParameter();
//
//   ComputeBoundaries();
//
//   TargetsPreliminaryOptions();
//
//   parametric_file_name_ = output_folder_ / "optimization.out";
//
//   std::vector<std::string> initial_values_str;
//   std::vector<std::string> lb_str;
//   std::vector<std::string> ub_str;
//
//   boost::split(initial_values_str, initial_values_string_, boost::is_any_of(" "));
//   boost::split(lb_str, lower_bounds_string_, boost::is_any_of(" "));
//   boost::split(ub_str, upper_bounds_string_, boost::is_any_of(" "));
//   boost::erase_all(param_name_string_, "'");
//   boost::split(param_str_, param_name_string_, boost::is_any_of(" "));
//
//   initial_values_str.pop_back();
//   lb_str.pop_back();
//   ub_str.pop_back();
//   param_str_.pop_back();
//
//   initial_values_.resize(initial_values_str.size());
//   std::transform(initial_values_str.begin(),
//                  initial_values_str.end(),
//                  initial_values_.begin(),
//                  [](const std::string& str) { return std::stod(str); });
//
//   lb_.resize(lb_str.size());
//   std::transform(lb_str.begin(), lb_str.end(), lb_.begin(), [](const std::string& str) { return std::stod(str); });
//
//   ub_.resize(ub_str.size());
//   std::transform(ub_str.begin(), ub_str.end(), ub_.begin(), [](const std::string& str) { return std::stod(str); });
// }
/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli                                              |
|                                                                                   |
|   Permission is hereby granted, free of charge, to any person obtaining a copy    |
|   of this software and associated documentation files (the "Software"), to deal   |
|   in the Software without restriction, including without limitation the rights    |
|   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       |
|   copies of the Software, and to permit persons to whom the Software is           |
|   furnished to do so, subject to the following conditions:                        |
|                                                                                   |
|   The above copyright notice and this permission notice shall be included in all  |
|   copies or substantial portions of the Software.                                 |
|                                                                                   |
|   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      |
|   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        |
|   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     |
|   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          |
|   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   |
|   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   |
|   SOFTWARE.                                                                       |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
