/*-----------------------------------------------------------------------*\
|                                                                         |
|       ____            _  ______ __  __  ____  _  ________               |
|      / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|              |
|     | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _      |
|     | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_    |
|     | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|   |
|      \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|     |
|            | |                                                          |
|            |_|                                                          |
|                                                                         |
\*-----------------------------------------------------------------------*/

namespace OptiSMOKE {
DataManager::DataManager() {}

DataManager::~DataManager() {}

void DataManager::ReadExperimentalData(vector<string> &experimental_data_files) {
  input_paths_.resize(experimental_data_files.size());
  ordinates_label_.resize(experimental_data_files.size());
  abscissae_label_.resize(experimental_data_files.size());
  uncertainty_kind_.resize(experimental_data_files.size());

  expdata_x_.resize(experimental_data_files.size());
  expdata_y_.resize(experimental_data_files.size());
  uncertainty_.resize(experimental_data_files.size());

  try {
    for (unsigned int i = 0; i < input_paths_.size(); i++) {
      ptree ptree;
      read_json(experimental_data_files[i].c_str(), ptree);

      dataset_names_.push_back(ptree.get<string>("name"));
      solver_name_.push_back(ptree.get<string>("solver"));
      QoI_.push_back(ptree.get<string>("QoI"));
      QoI_target_.push_back(ptree.get<string>("QoI_target"));
      multiple_input_.push_back(ptree.get<bool>("multiple_input"));

      boost::optional<string> mode = ptree.get_optional<string>("reactor_mode");
      if (mode) {
        reactor_mode_.push_back(ptree.get<string>("reactor_mode"));
      } else {
        reactor_mode_.push_back("");
      }

      BOOST_FOREACH (ptree::value_type &node, ptree.get_child("OS_Input_File")) {
        assert(node.first.empty());
        input_paths_[i].push_back(node.second.get_value<string>());
      }

      abscissae_label_[i].resize(ptree.get_child("data").size());
      ordinates_label_[i].resize(ptree.get_child("data").size());
      uncertainty_kind_[i].resize(ptree.get_child("data").size());
      expdata_x_[i].resize(ptree.get_child("data").size());
      expdata_y_[i].resize(ptree.get_child("data").size());
      uncertainty_[i].resize(ptree.get_child("data").size());

      unsigned int count = 0;
      BOOST_FOREACH (ptree::value_type &node, ptree.get_child("data")) {
        assert(node.first.empty());
        abscissae_label_[i][count] = node.second.get<string>("abscissae_label");
        ordinates_label_[i][count] = node.second.get<string>("ordinates_label");

        BOOST_FOREACH (ptree::value_type &node2, node.second.get_child("abscissae")) {
          assert(node2.first.empty());
          expdata_x_[i][count].push_back(node2.second.get_value<double>());
        }

        BOOST_FOREACH (ptree::value_type &node2, node.second.get_child("ordinates")) {
          assert(node2.first.empty());
          expdata_y_[i][count].push_back(node2.second.get_value<double>());
        }

        boost::optional<string> uncertainty_node =
            node.second.get_optional<string>("uncertainty_kind");
        if (uncertainty_node) {
          uncertainty_kind_[i][count] = node.second.get<string>("uncertainty_kind");
          BOOST_FOREACH (ptree::value_type &node2,
                         node.second.get_child("uncertainty")) {
            assert(node2.first.empty());
            uncertainty_[i][count].push_back(node2.second.get_value<double>());
          }
        } else {
          // TODO: implement the standard deviation
          // std::cout << "The uncertainty for the datasets ";
          // std::cout << experimental_data_files[i].c_str();
          // std::cout << " is not provided!" << std::endl;
          uncertainty_kind_[i][count] = "relative";
          for (unsigned int j = 0; j < expdata_y_[i][count].size(); j++) {
            uncertainty_[i][count].push_back(0);
          }
        }
        count += 1;
      }
    }

    OrderData();
  } catch (const std::exception &e) { OptiSMOKE::FatalErrorMessage(e.what()); }
}

void DataManager::OrderData() {
  /// TODO: This is a very bad function and it needs a major improvement
  /// But since at the moment I don't have any time lessgo with it the problem
  /// arises from the entire class itself by the way
  // Number of experimental data files
  unsigned int num = input_paths_.size();

  // tmp variables
  vector<string> dataset_names_tmp;
  vector<string> solver_name_tmp;
  vector<string> reactor_mode_tmp;
  vector<string> QoI_tmp;
  vector<string> QoI_target_tmp;
  vector<bool> multiple_input_tmp;
  vector<vector<string>> input_paths_tmp;

  vector<vector<string>> ordinates_label_tmp;
  vector<vector<string>> abscissae_label_tmp;
  vector<vector<string>> uncertainty_kind_tmp;
  vector<vector<vector<double>>> expdata_x_tmp;
  vector<vector<vector<double>>> expdata_y_tmp;
  vector<vector<vector<double>>> uncertainty_tmp;

  // Save position of the files
  vector<int> batch;
  vector<int> pfr;
  vector<int> psr;
  vector<int> premixed;
  vector<int> counterflow;

  for (unsigned int i = 0; i < num; i++) {
    if (solver_name_[i] == "BatchReactor") {
      batch.push_back(i);
    } else if (solver_name_[i] == "PlugFlowReactor") {
      pfr.push_back(i);
    } else if (solver_name_[i] == "PerfectlyStirredReactor") {
      psr.push_back(i);
    } else if (solver_name_[i] == "PremixedLaminarFlame1D") {
      premixed.push_back(i);
    } else if (solver_name_[i] == "CounterFlowFlame1D") {
      counterflow.push_back(i);
    } else {
      string error_str;
      error_str.append("Invalid solver inside datasets ");
      error_str.append(dataset_names_[i]);
      error_str.append(", available are:\n\t\t\t\t");
      error_str.append(
          "Batchreactor | PlugFlowreactor | PerfectlyStirredReactor | "
          "PremixedLaminarFlame1D | "
          "CounterFlowFlame1D");
      OptiSMOKE::FatalErrorMessage(error_str);
    }
  }

  if (batch.size() != 0) {
    for (unsigned int i = 0; i < batch.size(); i++) {
      unsigned int pos = batch[i];
      dataset_names_tmp.push_back(dataset_names_[pos]);
      solver_name_tmp.push_back(solver_name_[pos]);
      reactor_mode_tmp.push_back(reactor_mode_[pos]);
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
    }
  }

  if (pfr.size() != 0) {
    for (unsigned int i = 0; i < pfr.size(); i++) {
      unsigned int pos = pfr[i];
      dataset_names_tmp.push_back(dataset_names_[pos]);
      solver_name_tmp.push_back(solver_name_[pos]);
      reactor_mode_tmp.push_back(reactor_mode_[pos]);
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
    }
  }

  if (psr.size() != 0) {
    for (unsigned int i = 0; i < psr.size(); i++) {
      unsigned int pos = psr[i];
      dataset_names_tmp.push_back(dataset_names_[pos]);
      solver_name_tmp.push_back(solver_name_[pos]);
      reactor_mode_tmp.push_back(reactor_mode_[pos]);
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
    }
  }

  if (premixed.size() != 0) {
    for (unsigned int i = 0; i < premixed.size(); i++) {
      unsigned int pos = premixed[i];
      dataset_names_tmp.push_back(dataset_names_[pos]);
      solver_name_tmp.push_back(solver_name_[pos]);
      reactor_mode_tmp.push_back(reactor_mode_[pos]);
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
    }
  }

  if (counterflow.size() != 0) {
    for (unsigned int i = 0; i < counterflow.size(); i++) {
      unsigned int pos = counterflow[i];
      dataset_names_tmp.push_back(dataset_names_[pos]);
      solver_name_tmp.push_back(solver_name_[pos]);
      reactor_mode_tmp.push_back(reactor_mode_[pos]);
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
    }
  }

  dataset_names_    = dataset_names_tmp;
  reactor_mode_     = reactor_mode_tmp;
  solver_name_      = solver_name_tmp;
  QoI_              = QoI_tmp;
  QoI_target_       = QoI_target_tmp;
  multiple_input_   = multiple_input_tmp;
  input_paths_      = input_paths_tmp;
  abscissae_label_  = abscissae_label_tmp;
  ordinates_label_  = ordinates_label_tmp;
  uncertainty_kind_ = uncertainty_kind_tmp;
  expdata_x_        = expdata_x_tmp;
  expdata_y_        = expdata_y_tmp;
  uncertainty_      = uncertainty_tmp;
}

void DataManager::ComputeStandardDeviations() {
  // Not yet implemented
  standard_deviations_.resize(expdata_y_.size());
  for (unsigned int i = 0; i < expdata_y_.size(); i++) {
    standard_deviations_[i].resize(expdata_y_[i].size());
    for (unsigned int j = 0; j < expdata_y_[i].size(); j++) {
      standard_deviations_[i][j].resize(expdata_y_[i][j].size());
      for (unsigned int k = 0; k < expdata_y_[i][j].size(); k++) {
        // standard_deviations_[i][j][k] = uncertainty_[i][j][k] * (expdata_y_[i][j][k])
        // / Sigma_vector[i];
      }
    }
  }
}
}  // namespace OptiSMOKE
