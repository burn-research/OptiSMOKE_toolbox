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
|        Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                      |
|                 Andrea Bertolino <andrea.bertolino@ulb.be>                        |
|                 Magnus Fürst     <magnus.furst@ulb.ac.be>                         |
|                                                                                   |
|          [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>             |
|              Department of Chemistry, Materials and Chemical Engineering          |
|              Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano      |
|                                                                                   |
|          [2] BRITE Research Group <https://brite-research.be>                     |
|              Brussels Institute for Thermal-fluid systems and clean Energy        |
|              Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                 |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class GrammarOptismoke : public OpenSMOKE::OpenSMOKE_DictionaryGrammar {
 protected:
  virtual void DefineRules() {
    // clang-format off
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DakotaOptions",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary with options for Dakota",
                                                      false,
                                                      "none",
                                                      "@OptimizationLibrary",
                                                      "@NLOPTOptions"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NLOPTOptions",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary with options for NLOPT",
                                                      false,
                                                      "none",
                                                      "@OptimizationLibrary",
                                                      "@DakotaOptions"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CurveMatchingOptions",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary with options for setting Curve Matching as the objective function",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationSetup",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary to set the options for the optimization",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationTarget",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary to set the targets for the optimization",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
                                                      OpenSMOKE::SINGLE_DICTIONARY,
                                                      "Name of the dictionary containing the list of kinetic files to be interpreted",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputFolder",
                                                      OpenSMOKE::SINGLE_PATH,
                                                      "Name of the folder where to write all the output file of the program",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfDataFiles",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "List of Experimental data files",
                                                      true));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationLibrary",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Name of the library containing the optimization routines",
                                                      true));
    // clang-format on
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
