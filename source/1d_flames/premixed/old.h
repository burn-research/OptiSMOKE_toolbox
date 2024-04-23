#ifndef OPTISMOKE_PREMIXEDLAMINARFLAME1D_H
#define OPTISMOKE_PREMIXEDLAMINARFLAME1D_H

// Utilities
#include <idealreactors/utilities/Utilities>
#include <maps/TransportPropertiesMap_CHEMKIN.h>
#include <utilities/ropa/OnTheFlyROPA.h>
#include <utilities/cema/OnTheFlyCEMA.h>
#include <utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h>
#include <utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h>
#include <utilities/soot/hmom/HMOM.h>
#include <idealreactors/utilities/Grammar_LewisNumbers.h>

// 1D grid
#include <premixedlaminarflame1d/utilities/Utilities.h>  // This is in premixed utilities.
#include <utilities/grids/adaptive/Grid1D.h>
#include <utilities/grids/adaptive/Grammar_Grid1D.h>
#include <utilities/grids/adaptive/Adapter_Grid1D.h>

// Laminar Flame 1D
#include <premixedlaminarflame1d/OpenSMOKE_PremixedLaminarFlame1D.h>

// Thermodynamic equilibrium class
#include <thermodynamicequilibrium/ThermodynamicEquilibrium.h>

// Multicomponent transport library
#include <maps/MulticomponentTransportLibrary.h>

#include "grammar/grammar_premixedlaminarflame.h"

namespace OptiSMOKE {
class PremixedLaminarFlame1D {
 public:
  void Setup(const std::string input_file_name_, OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
             OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML,
             OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);

  void Solve();

  const double& LFS() const { return LFS_; };

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;
  OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML_;
  OpenSMOKE::Grid1D* grid;
  OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D* flame_premixed;
  OpenSMOKE::FlammabilityLimits* flammability_limits;
  OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
  OpenSMOKE::MulticomponentTransportLibrary* mutlib;
  OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
  OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
  OpenSMOKE::HMOM* hmom;
  DaeSMOKE::DaeSolver_Parameters* dae_parameters;
  NlsSMOKE::NonLinearSolver_Parameters* nls_parameters;
  NlsSMOKE::FalseTransientSolver_Parameters* false_transient_parameters;

  std::vector<OpenSMOKE::OpenSMOKEVectorDouble> inlet_omega;
  std::vector<double> equivalence_ratios;
  std::vector<double> inlet_T;
  std::vector<double> P_Pa;
  double outlet_T;
  OpenSMOKE::OpenSMOKEVectorDouble outlet_omega;

  double LFS_;  // cm/s
};
}  // namespace OptiSMOKE

#include "PremixedLaminarFlame1D.hpp"
#endif  // OPTISMOKE_PREMIXEDLAMINARFLAME1D_H
