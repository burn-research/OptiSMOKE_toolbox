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

#ifndef OpenSMOKE_BatchReactor_VolumeProfile_H
#define	OpenSMOKE_BatchReactor_VolumeProfile_H

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	class BatchReactor_VolumeProfile
	{
		public:
	
			enum type_of_law { BatchReactor_VolumeProfile_PressureCoefficient, BatchReactor_VolumeProfile_VolumeHistory, BatchReactor_VolumeProfile_FromPressureProfile };	
			
			BatchReactor_VolumeProfile();
			
			void SetProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& y);
			void SetPressureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& y);
			void SetPressureCoefficient(const double pressureCoefficient);
			
			double pressureCoefficient() const { return pressureCoefficient_; }
			double operator () (const double x) const;

			double type() const { return type_; }
		
		private:
		
			unsigned int n_;

			type_of_law type_;
			
			OpenSMOKE::OpenSMOKEVectorDouble x_;
			OpenSMOKE::OpenSMOKEVectorDouble y_;
			OpenSMOKE::OpenSMOKEVectorDouble m_;
			
			double x0_;
			double xf_;
			double y0_;
			double yf_;	

			double pressureCoefficient_;
	};
}

#include "BatchReactor_VolumeProfile.hpp"

#endif	/* OpenSMOKE_BatchReactor_VolumeProfile_H */

