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

#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <boost/foreach.hpp>
#include <boost/json.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
using namespace boost::property_tree;

namespace OptiSMOKE {

class DataManager {
 public:
  DataManager();

  ~DataManager();

  void ReadExperimentalData(vector<string> &experimental_data_files);

  const vector<string> &dataset_names() const { return dataset_names_; };

  const vector<string> &solver_name() const { return solver_name_; };

  const vector<string> &QoI() const { return QoI_; };

  const vector<string> &QoI_target() const { return QoI_target_; };

  const vector<bool> &multiple_input() const { return multiple_input_; };

  const vector<vector<string>> &input_paths() const { return input_paths_; };

  const vector<vector<string>> &ordinates_label() const { return ordinates_label_; };

  const vector<vector<string>> &abscissae_label() const { return abscissae_label_; };

  const vector<vector<string>> &uncertainty_kind() const { return uncertainty_kind_; };

  const vector<vector<vector<double>>> &expdata_x() const { return expdata_x_; };

  const vector<vector<vector<double>>> &expdata_y() const { return expdata_y_; };

  const vector<vector<vector<double>>> &uncertainty() const { return uncertainty_; };

  const vector<string> &reactor_mode() const { return reactor_mode_; };

 private:
  vector<string> dataset_names_;
  vector<string> solver_name_;
  vector<string> QoI_;
  vector<string> QoI_target_;
  vector<bool> multiple_input_;
  vector<string> reactor_mode_;

  vector<vector<string>> input_paths_;

  // This blocks here has to go three dimensions
  // since a datasets file can have multiple series
  // dimension one: number of files
  // dimension two: number of datasets whithin each file
  // dimension three: number point in each datasets
  vector<vector<string>> ordinates_label_;
  vector<vector<string>> abscissae_label_;
  vector<vector<string>> uncertainty_kind_;
  vector<vector<vector<double>>> expdata_x_;
  vector<vector<vector<double>>> expdata_y_;
  vector<vector<vector<double>>> uncertainty_;
  vector<vector<vector<double>>> standard_deviations_;

  void OrderData();

  // Default sigma for standard deviation if it is not present
  // inside the file.
  const double default_sigma = 2;
  void ComputeStandardDeviations();
};
}  // namespace OptiSMOKE

#include "DataManager.hpp"
#endif  // DATAMANAGER_H
