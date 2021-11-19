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

#ifndef OpenSMOKE_BatchReactor_OdeInterfaces_H
#define	OpenSMOKE_BatchReactor_OdeInterfaces_H

#include "math/OpenSMOKEVector.h"
#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"

#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#endif

#if OPENSMOKE_USE_DVODE == 1
#include "math/external-ode-solvers/OpenSMOKE_DVODE_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DVODE.h"
#endif

#if OPENSMOKE_USE_ODEPACK == 1
#include "math/external-ode-solvers/OpenSMOKE_DLSODE_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODE.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODA_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DLSODA.h"
#endif

#if OPENSMOKE_USE_DASPK == 1
#include "math/external-ode-solvers/OpenSMOKE_DASPK_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_DASPK.h"
#endif

#if OPENSMOKE_USE_RADAU == 1
#include "math/external-ode-solvers/OpenSMOKE_RADAU_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_RADAU.h"
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
#include "math/external-ode-solvers/OpenSMOKE_CVODE_Sundials_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_CVODE_Sundials.h"
#endif

#if OPENSMOKE_USE_MEBDF == 1
#include "math/external-ode-solvers/OpenSMOKE_MEBDF_Interface.h"
#include "math/external-ode-solvers/OpenSMOKE_MEBDF.h"
#endif

#include "BatchReactor.h"

namespace OpenSMOKE
{
	class ODESystem_OpenSMOKE_BatchReactor
	{
	public:

		ODESystem_OpenSMOKE_BatchReactor() {};

		void SetReactor(BatchReactor* reactor)
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
			reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(DY.data());
		}

		void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J)
		{
			y_.CopyFrom(Y.data());
			reactor_->DenseAnalyticalJacobian(t, y_, J);
		};

		void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::SparseMatrix<double> &J)
		{
			y_.CopyFrom(Y.data());
			reactor_->SparseAnalyticalJacobian(t, y_, J);
		};

		void Print(const double t, const Eigen::VectorXd &Y)
		{
			y_.CopyFrom(Y.data());
			reactor_->Print(t, y_);
		}

	private:

		BatchReactor* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};

	#if OPENSMOKE_USE_BZZMATH == 1
	class ODESystem_BzzOde_BatchReactor : public BzzOdeSystemObject 
	{
	public:

		ODESystem_BzzOde_BatchReactor(	BatchReactor& batch) :
		batch_(batch)
		{
			ChangeDimensions(batch_.NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_.NumberOfEquations(), &dy_, false);
		}
	
		virtual void GetSystemFunctions(BzzVector &Y, double t, BzzVector &DY)
		{
			y_.CopyFrom(Y.GetHandle());
			batch_.Equations(t, y_, dy_);
			dy_.CopyTo(DY.GetHandle());
		}

		void MyPrint(BzzVector &Y, double t)
		{
			y_.CopyFrom(Y.GetHandle());
			batch_.Print(t, y_);
		}

		virtual void ObjectBzzPrint(void) {};

	private:

		BatchReactor& batch_;
		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};
	#endif

	
	#if OPENSMOKE_USE_DVODE == 1
	class ODESystem_DVODE_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DVODE(ODESystem_DVODE_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DVODE(ODESystem_DVODE_BatchReactor)
	#endif

	#if OPENSMOKE_USE_ODEPACK == 1
	class ODESystem_DLSODE_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODE(ODESystem_DLSODE_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODE(ODESystem_DLSODE_BatchReactor)

	class ODESystem_DLSODA_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODA(ODESystem_DLSODA_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODA(ODESystem_DLSODA_BatchReactor)
	#endif

	#if OPENSMOKE_USE_DASPK == 1
	class ODESystem_DASPK_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		BatchReactor* batch_;

		DEFINE_ODESOLVERINTERFACE_DASPK(ODESystem_DASPK_BatchReactor)

		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DASPK(ODESystem_DASPK_BatchReactor)
	#endif

	#if OPENSMOKE_USE_RADAU == 1
	class ODESystem_RADAU5_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_RADAU(ODESystem_RADAU5_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_RADAU(ODESystem_RADAU5_BatchReactor)
	#endif

	#if OPENSMOKE_USE_SUNDIALS == 1
	class ODESystem_CVODE_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_CVODE_Sundials(ODESystem_CVODE_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_CVODE_Sundials(ODESystem_CVODE_BatchReactor)
	#endif

	#if OPENSMOKE_USE_MEBDF == 1
	class ODESystem_MEBDF_BatchReactor : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_MEBDF(ODESystem_MEBDF_BatchReactor)

		BatchReactor* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactor(BatchReactor* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
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
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_MEBDF(ODESystem_MEBDF_BatchReactor)
	#endif
}

#endif	// OpenSMOKE_BatchReactor_OdeInterfaces_H
