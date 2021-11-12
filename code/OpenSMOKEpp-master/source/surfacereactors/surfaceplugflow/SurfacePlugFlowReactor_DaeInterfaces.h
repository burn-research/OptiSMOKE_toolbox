/*-----------------------------------------------------------------------*\
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
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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

#ifndef OpenSMOKE_SurfacePlugFlowReactor_DaeInterfaces_H
#define	OpenSMOKE_SurfacePlugFlowReactor_DaeInterfaces_H

#include "math/OpenSMOKEVector.h"
#include "math/external-dae-solvers/OpenSMOKE_DaeSystemObject.h"

#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#endif

#if OPENSMOKE_USE_DASPK == 1
#include "math/external-dae-solvers/OpenSMOKE_DASPK_DAE_Interface.h"
#include "math/external-dae-solvers/OpenSMOKE_DASPK_DAE.h"
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
#include "math/external-dae-solvers/OpenSMOKE_IDA_Sundials_Interface.h"
#include "math/external-dae-solvers/OpenSMOKE_IDA_Sundials.h"
#endif

#include "math/native-dae-solvers/MultiValueSolver"
#include "SurfacePlugFlowReactor.h"

namespace OpenSMOKE
{
	class DAESystem_OpenSMOKE_SurfacePlugFlowReactor
	{
	public:

		DAESystem_OpenSMOKE_SurfacePlugFlowReactor() {};

		void SetReactor(SurfacePlugFlowReactor* reactor)
		{
			reactor_ = reactor;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
			OpenSMOKE::ChangeDimensions(ne_, &y_, true);
			OpenSMOKE::ChangeDimensions(ne_, &dy_, false);
		}

		virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
		{
			y_.CopyFrom(Y.data());
			reactor_->DaeEquations(t, y_, dy_);
			dy_.CopyTo(DY.data());
		}

		virtual void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) {};

		void Print(const double t, const Eigen::VectorXd &Y)
		{
			y_.CopyFrom(Y.data());
			reactor_->DaePrint(t, y_);
		}

	private:

		SurfacePlugFlowReactor* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};

	#if OPENSMOKE_USE_BZZMATH == 1
	class DAESystem_BzzDae_SurfacePlugFlowReactor : public BzzDaeSystemObject 
	{
	public:

		DAESystem_BzzDae_SurfacePlugFlowReactor(	SurfacePlugFlowReactor& plugflow) :
		plugflow_(plugflow)
		{
			ChangeDimensions(plugflow_.NumberOfDaeEquations(), &y_, true);
			ChangeDimensions(plugflow_.NumberOfDaeEquations(), &dy_, false);
		}
	
		virtual void GetSystemFunctions(BzzVector &Y, double t, BzzVector &DY)
		{
			y_.CopyFrom(Y.GetHandle());
			plugflow_.DaeEquations(t, y_, dy_);
			dy_.CopyTo(DY.GetHandle());
		}

		void MyPrint(BzzVector &Y, double t)
		{
			y_.CopyFrom(Y.GetHandle());
			plugflow_.DaePrint(t, y_);
		}

		virtual void ObjectBzzPrint(void) {};

	private:

		SurfacePlugFlowReactor& plugflow_;
		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};
	#endif
	
	#if OPENSMOKE_USE_DASPK == 1
	class DAESystem_DASPK_SurfacePlugFlowReactor : public OpenSMOKE::OpenSMOKE_DaeSystemObject
	{
		SurfacePlugFlowReactor* plugflow_;

		DEFINE_DAESOLVERINTERFACE_DASPK_DAE(DAESystem_DASPK_SurfacePlugFlowReactor)

		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetSurfacePlugFlowReactor(SurfacePlugFlowReactor* plugflow)
		{
			plugflow_ = plugflow;
			ChangeDimensions(plugflow_->NumberOfDaeEquations(), &y_, true);
			ChangeDimensions(plugflow_->NumberOfDaeEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = plugflow_->DaeEquations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = plugflow_->DaePrint(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_DAESOLVERINTERFACE_DASPK_DAE(DAESystem_DASPK_SurfacePlugFlowReactor)
	#endif
	
	#if OPENSMOKE_USE_SUNDIALS == 1
	class DAESystem_IDA_SurfacePlugFlowReactor : public OpenSMOKE::OpenSMOKE_DaeSystemObject
	{
		DEFINE_DAESOLVERINTERFACE_IDA_Sundials(DAESystem_IDA_SurfacePlugFlowReactor)

		SurfacePlugFlowReactor* plugflow_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetSurfacePlugFlowReactor(SurfacePlugFlowReactor* plugflow)
		{
			plugflow_ = plugflow;
			ChangeDimensions(plugflow_->NumberOfDaeEquations(), &y_,  true);
			ChangeDimensions(plugflow_->NumberOfDaeEquations(), &dy_, true);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = plugflow_->DaeEquations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = plugflow_->DaePrint(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_DAESOLVERINTERFACE_IDA_Sundials(DAESystem_IDA_SurfacePlugFlowReactor)
	#endif
}

#endif	// OpenSMOKE_SurfacePlugFlowReactor_DaeInterfaces_H
