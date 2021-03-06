/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_SurfaceOnTheFlyROPA_H
#define	OpenSMOKE_SurfaceOnTheFlyROPA_H

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"

namespace OpenSMOKE
{
	class SurfaceOnTheFlyROPA
	{
	public:

		SurfaceOnTheFlyROPA(	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsMap);

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, boost::filesystem::path path_kinetics_output);

		bool is_active() { return is_active_; }

		void Analyze(	std::ofstream& fOut, const int current_step, const double t, const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEVectorDouble& omega,
						const OpenSMOKE::OpenSMOKEVectorDouble& Z, const OpenSMOKE::OpenSMOKEVectorDouble& a, const OpenSMOKE::OpenSMOKEVectorDouble& gamma);

		void WriteHead(std::ofstream& fOut, const std::string reactor_type);

		void SetThreshold(const double threshold) { threshold_contributions_ = threshold; }
		void SetCompactMode(const bool compact_mode) { compact_mode_ = compact_mode; }
		void SetMergeForwardAndBackwardReactions(const bool flag) { merge_reaction_rates_ = flag; }
		void SetSpecies(const std::vector<std::string> list_of_species) { list_species_ = list_of_species; }
		void SetActive(const bool flag) { is_active_ = flag; }

		void Setup(const boost::filesystem::path& path_kinetics_output);

	private:

	//	enum writing_policy_type { FIXED_TIME_STEPS_, LIST_CONVERSIONS, LIST_TIMES };
	//	bool CheckForROPA(const int current_step, const double time, OpenSMOKE::OpenSMOKEVectorDouble& conversions);

	private:

		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
		OpenSMOKE::KineticsMap_Surface_CHEMKIN&			kineticsMap_;			//!< kinetic map

		std::vector<std::string> reaction_names_;
		std::vector<std::string> list_species_;

		double threshold_contributions_;
		double threshold_formation_rate_;

		bool merge_reaction_rates_;
		bool compact_mode_;

		bool is_active_;
		OpenSMOKE::OpenSMOKEVectorDouble R_;
		OpenSMOKE::OpenSMOKEVectorDouble P_;
		OpenSMOKE::OpenSMOKEVectorDouble D_;
		OpenSMOKE::OpenSMOKEVectorDouble rf_;
		OpenSMOKE::OpenSMOKEVectorDouble rb_;
	};
}

#include "SurfaceOnTheFlyROPA.hpp"

#endif	/* OpenSMOKE_SurfaceOnTheFlyROPA_H */

