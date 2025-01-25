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
|                 Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano   |
|                                                                                   |
|             [2] BRITE Research Group <https://brite-research.be>                  |
|                 Brussels Institute for Thermal-fluid systems and clean Energy     |
|                 Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel              |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class GrammarPremixedLaminarFlame : public OpenSMOKE::OpenSMOKE_DictionaryGrammar {
 protected:
  virtual void DefineRules() {
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
                                                      OpenSMOKE::SINGLE_PATH,
                                                      "Name of the folder containing the kinetic scheme (XML Version)",
                                                      true,
                                                      "@KineticsPreProcessor",
                                                      "none",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Backup",
                                                      OpenSMOKE::SINGLE_PATH,
                                                      "Name of backup file (XML Version)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@DontUseBackupGrid",
        OpenSMOKE::SINGLE_BOOL,
        "If true, the user defined grid is used, instead of the grid corresponding to the backup file (default:false)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@KineticsPreProcessor",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Name of the dictionary containing the list of kinetic files to be interpreted",
        true,
        "@KineticsFolder",
        "none",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Type of simulation: BurnerStabilized | FlameSpeed",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Grid",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary defining the mesh",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@InletStream",
        OpenSMOKE::VECTOR_STRING,
        "Name of the dictionary/dictionaries defining the inlet gas composition, temperature and pressure",
        true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutletStream",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary defining the outlet gas composition, "
                                                      "temperature and pressure (this is only for initialization)",
                                                      false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletVelocity",
                                               OpenSMOKE::SINGLE_MEASURE,
                                               "Inlet velocity (first guess in case of flame speed calculations)",
                                               true,
                                               "@InletMassFlux",
                                               "none",
                                               "none"));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletMassFlux",
                                               OpenSMOKE::SINGLE_MEASURE,
                                               "Inlet mass flux (first guess in case of flame speed calculations)",
                                               true,
                                               "@InletVelocity",
                                               "none",
                                               "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@SensitivityAnalysis",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Dictionary containing additional options for solving the sensitivity analysis",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@DaeParameters",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Dictionary containing the numerical parameters for solving the stiff DAE system",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@NlsParameters",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Dictionary containing the numerical parameters for solving the NL system",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@FalseTransientParameters",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Dictionary containing the numerical parameters for solving the pseudo-transient phase",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output",
                                                      OpenSMOKE::SINGLE_PATH,
                                                      "Name of the folder where to write the output files",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@UseNlsSolver",
        OpenSMOKE::SINGLE_BOOL,
        "Use the NLS solver to solve the steady state problems, after the DAE solver (default: true)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@UseDaeSolver",
        OpenSMOKE::SINGLE_BOOL,
        "Use the Dae solver (instead of NLS) to solve the steady state problems (default: true)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@PolimiSoot",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Name of the dictionary defining the rules for analyzing soot calculated using the Polimi mechanism",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@HMOM",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Name of the dictionary defining the rules for applying the Hybrid Method of Moments (HMOM)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@OnTheFlyPostProcessing",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Dictionary specifying the details for carrying out the post-processing analyses (on the fly)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedTemperatureProfile",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of dictionary describing the temperature profile",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@FixedSpecificMassFlowRateProfile",
        OpenSMOKE::SINGLE_DICTIONARY,
        "Name of dictionary describing the specific (i.e. per unit area) mass flow rate profile",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Soret",
                                                      OpenSMOKE::SINGLE_BOOL,
                                                      "Add Soret effect (default: true)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@SoretKuoCorrelation",
        OpenSMOKE::SINGLE_BOOL,
        "Kuo's correlation used for calculating the thermal diffusion coefficients (default: false)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@FrozenMassDiffusivities",
        OpenSMOKE::SINGLE_BOOL,
        "Freeze mass diffuivity calculation during the Jacobian evaluation (default: false)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@SpeciesBundling",
        OpenSMOKE::SINGLE_DOUBLE,
        "Estimation of mass diffusion coefficients through the species bundling according the the specified maximum "
        "relative error (example: @SpeciesBundling 0.1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Radiation",
                                                      OpenSMOKE::SINGLE_BOOL,
                                                      "Radiative heat transfer (default: none)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature",
                                                      OpenSMOKE::SINGLE_MEASURE,
                                                      "Environment temperature (default: 298.15 K)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeGasMassFractions",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Derivative gas mass fractions (default: backward)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeGasTemperature",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Derivative gas temperature (default: backward)",
                                                      false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LewisNumbers",
                                               OpenSMOKE::SINGLE_DICTIONARY,
                                               "Name of dictionary containing the list of Lewis numbers of species",
                                               false,
                                               "none",
                                               "none",
                                               "@StefanMaxwell @MulticomponentTransport"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StefanMaxwell",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Stefan-Maxwell theory for diffusion fluxes",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "@LewisNumbers @MulticomponentTransport"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MulticomponentTransport",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Multicomponent transport based on the MuTLib",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "@StefanMaxwell @LewisNumbers"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@FixedOutletTemperature",
        OpenSMOKE::SINGLE_MEASURE,
        "Fixed outlet temperature (to be used only if @FixedSpecificMassFlowRateProfile is on)",
        false,
        "none",
        "@FixedSpecificMassFlowRateProfile",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallHeatExchangeCoefficient",
                                                      OpenSMOKE::SINGLE_MEASURE,
                                                      "Wall heat exchange coefficient (default: 0 W/m2/K)",
                                                      false,
                                                      "none",
                                                      "@InternalDiameter",
                                                      "@WallHeatNusseltNumber"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallHeatNusseltNumber",
                                                      OpenSMOKE::SINGLE_DOUBLE,
                                                      "Wall heat Nusselt number (default: 0)",
                                                      false,
                                                      "none",
                                                      "@InternalDiameter",
                                                      "@WallHeatExchangeCoefficient"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InternalDiameter",
                                                      OpenSMOKE::SINGLE_MEASURE,
                                                      "Internal diameter (default: 1 cm)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallTemperatureProfile",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of dictionary describing the wall temperature profile",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TaylorArisCorrection",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Taylor-Aris correction: none (default), mass-heat, mass, heat",
                                                      false,
                                                      "none",
                                                      "@InternalDiameter",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsModifier",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary defining the kinetics modifier",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WriteRHS",
                                                      OpenSMOKE::SINGLE_BOOL,
                                                      "Writes on a file the contributions to the RHS (default: false)",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FlammabilityLimits",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary defining the flammability limits options",
                                                      false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SimplifiedBoundaryFluxes",
                                               OpenSMOKE::SINGLE_BOOL,
                                               "The fluxes on the boundaries are calculated using an approximated "
                                               "formulation, more robust, but less accurate (default: false)",
                                               false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TestingOption",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "To easily make code experiments without recompiling",
                                                      false));
  }
};
}  // namespace OptiSMOKE

/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli, Andrea Bertolino, Magnus Fürst              |
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
