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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it> |
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

#ifndef OPTIONS_KINETICS_H
#define OPTIONS_KINETICS_H

#include "idealreactors/utilities/Grammar_RapidKineticMechanism.h"

namespace OptiSMOKE {
class options_kinetics {
 public:
  options_kinetics();
  ~options_kinetics();

  void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                           std::string dictionary_name);

  // Access function

  const fs::path& chemkin_kinetics() const { return chemkin_kinetics_; };
  const fs::path& chemkin_thermodynamics() const { return chemkin_thermodynamics_; };
  const fs::path& chemkin_transport() const { return chemkin_transport_; };
  const fs::path& chemkin_output() const { return chemkin_output_; };
  const bool& iTransport() const {return iTransport_;};

 private:
  OpenSMOKE::Grammar_RapidKineticMechanism kinetics_grammar_;

  fs::path chemkin_kinetics_;
  fs::path chemkin_thermodynamics_;
  fs::path chemkin_transport_;
  fs::path chemkin_output_;

  bool iTransport_;
};
}  // namespace OptiSMOKE

#include "options_kinetics.hpp"

#endif  // OPTIONS_KINETICS_H
