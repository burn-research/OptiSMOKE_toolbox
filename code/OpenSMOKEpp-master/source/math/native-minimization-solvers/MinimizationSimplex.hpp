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
|   Copyright(C) 2018  Alberto Cuoci                                      |
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

namespace OpenSMOKE
{

	const double MinimizationSimplex::ETA_REL_SIMPLEX = 1.e-1;
	const double MinimizationSimplex::ETA_ABS_SIMPLEX = .01;
	const double MinimizationSimplex::BETA = 0.5;
	const double MinimizationSimplex::GAMMA = 2.;
	const int    MinimizationSimplex::MAX_DEFAULT_ITERATIONS = 100000000;
	const double MinimizationSimplex::X_TOL_REL = 1.e-13;
	const double MinimizationSimplex::X_TOL_ABS = 1.e-25;
	const double MinimizationSimplex::F_TOL_REL = 1.e-13;
	const double MinimizationSimplex::F_TOL_ABS = 1.e-25;
	const double MinimizationSimplex::BIG_NUMBER = 1.e+30;

	bool MinimizationSimplex::unfeasible = false;

	double MachEpsDouble()
	{
		double macheps = 1.;
		double eps = 2.;
		while (eps != 1.)
		{
			macheps /= 2.;
			eps = 1. + macheps; // in double rounded form
		}
		return macheps * 2.;
	}

	void SortAndTrackIndicesIncreasing(std::vector<unsigned int>& indices, std::vector<double>& values)
	{
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort(begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] < values[b]; });
		std::sort(values.begin(), values.end());
		/*
		BzzVectorInt s(indices.size());
		BzzVector aux(indices.size());
		for (unsigned int i = 1; i <= indices.size(); i++)
		aux[i] = values[i - 1];
		Sort(&aux, &s);
		for (unsigned int i = 1; i <= indices.size(); i++)
		{
		values[i - 1] = aux[i];
		indices[i - 1] = s[i] - 1;
		}
		*/
	}

	void MinimizationSimplex::InitialSetup(const Eigen::VectorXd &x00)
	{
		yTol_ = -BIG_NUMBER;
		stop_ = 0;
		xMin_ = xi_ = x0_ = x00;
		iter_ = 0;
		numVariables_ = x0_.size();
		numVertices_ = numVariables_ + 1;
		v_ = new Eigen::VectorXd[numVertices_];
		f_.resize(numVertices_);
		aux_.resize(numVertices_);
		sorted_.resize(numVertices_);
		vB_.resize(numVariables_);

		v_[numVertices_ - 1] = xi_;
		f_(numVertices_ - 1) = fi_;
		h_.resize(numVariables_);
		xTolRel_ = X_TOL_REL;
		xTolAbs_ = X_TOL_ABS;
		fTolRel_ = F_TOL_REL;
		fTolAbs_ = F_TOL_ABS;
		control_ = -1;
		stop_ = 0;
		methodStatus_ = START;
	}

	void MinimizationSimplex::SetTolF(const double yt)
	{
		control_ = -1;
		yTol_ = yt;
	}

	void MinimizationSimplex::SetTolAbsF(const double tolAbs)
	{
		control_ = -1;
		fTolAbs_ = tolAbs;
	}

	void MinimizationSimplex::SetTolRelF(const double tolRelF)
	{
		control_ = -1;
		if (tolRelF > 10.*MachEpsDouble())
			fTolRel_ = tolRelF;
		else
			fTolRel_ = 10.*MachEpsDouble();
	}

	void MinimizationSimplex::SetTolAbsX(const double tolAbsX)
	{
		control_ = -1;
		xTolAbs_ = tolAbsX;
	}

	void MinimizationSimplex::SetTolRelX(const double tolRelX)
	{
		control_ = -1;
		if (tolRelX > 10. * MachEpsDouble())
			xTolRel_ = tolRelX;
		else
			xTolRel_ = 10. * MachEpsDouble();
	}

	char MinimizationSimplex::operator()()
	{
		return (*this)(0);
	}


	double MinimizationSimplex::MinFunction(const Eigen::VectorXd &x)
	{
		double F;
		int constraint = 0;

		switch (functionSimplexType_)
		{

		case TWO_IV1_IV2:

			zi_ = z0_;
			zi_(iVar1_ - 1) = x(0);
			zi_(iVar2_ - 1) = x(1);
			for (unsigned int i = 0; i < numVariables_; i++)
			{
				if (x(i) < xL_(i))
				{
					F = f0_ + 10. * (xL_(i) - x(i));
					constraint = 1;
				}
				if (x(i) > xU_(i))
				{
					F = f0_ + 10. * (x(i) - xU_(i));
					constraint = 1;
				}
			}

			if (constraint == 0)
				F = ptrFun(zi_);

			return F;

		case UNCONSTRAINED_MULTI:

			for (unsigned int i = 0; i < numVariables_; i++)
			{
				if (x(i) < xL_(i))
				{
					F = f0_ + 10. * (xL_(i) - x(i));
					constraint = 1;
				}
				if (x(i) > xU_(i))
				{
					F = f0_ + 10. * (x(i) - xU_(i));
					constraint = 1;
				}
			}

			if (constraint == 0)
				F = ptrFun(x);

			return F;
		}

		return 0.;
	}

	void MinimizationSimplex::Start(void)
	{
		double fh;

		for (unsigned int i = 0; i < numVariables_; i++)
			h_(i) = ETA_REL_SIMPLEX * std::fabs(xi_(i)) + ETA_ABS_SIMPLEX;

		for (unsigned int i = 0; i < numVariables_; i++)
		{
			x_ = xi_;
			x_(i) += h_(i);
			if (x_(i) > xU_(i))
			{
				x_(i) = xi_(i) - h_(i);
				if (x_(i) < xL_(i))
					x_(i) = xi_(i) + 0.01 * h_(i);
			}

			unfeasible = false;
			iterTotal_++;
			fh = MinFunction(x_);
			if (unfeasible == true)
			{
				x_(i) = xi_(i);
				x_(i) -= h_(i);
				if (x_(i) < xL_(i))
				{
					x_(i) = xi_(i) + h_(i);
					if (x_(i) > xU_(i))
						x_(i) = xi_(i) - .01 * h_(i);
				}
				unfeasible = false;
				iterTotal_++;
				fh = MinFunction(x_);

				if (unfeasible == true)
				{
					fh = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}
			}
			else
			{
				iterTotalFeasible_++;
				iter_++;
			}
			v_[i] = x_;
			f_(i) = fh;
		}

		for (unsigned int i = 0; i < aux_.size(); i++)
			aux_[i] = f_[i];
		SortAndTrackIndicesIncreasing(sorted_, aux_);

		fMin_ = f_(sorted_[0]);
		xMin_ = v_[sorted_[0]];

		if (fMin_ < yTol_)
		{
			control_ = 3;
			methodStatus_ = EXIT;
			return;
		}
		methodStatus_ = SIMPLEX_ACTION;
	}

	void MinimizationSimplex::Baricenter(void)
	{
		vB_.setZero();
		for (unsigned int i = 0; i < numVariables_; i++)
			vB_ += v_[sorted_[i]];

		vB_ /= double(numVariables_);
	}

	void MinimizationSimplex::SimplexAction(void)
	{
		double fh;

		if (fMin_ > 0. && f_.maxCoeff() < 1.e-20 && iterTotal_ > 2 * numVariables_ + 100)
		{
			std::cout << "WARNING Message: Robust " << iterSimplexCollapsed_ << " " << iterTotal_ << " " << fMin_ << std::endl;
			control_ = 2;
			methodStatus_ = EXIT;
			return;
		}
		if (fMin_ < yTol_)
		{
			control_ = 3;
			methodStatus_ = EXIT;
			return;
		}

		vP_ = v_[sorted_[numVertices_ - 1]];
		fP_ = f_(sorted_[numVertices_ - 1]);
		vA_ = v_[sorted_[numVariables_ - 1]];
		fA_ = f_(sorted_[numVariables_ - 1]);

		Baricenter();
		dvPvB_ = vB_ - vP_;
		normdbx_ = dvPvB_.lpNorm<2>() / double(numVariables_);
		normB_ = vB_.lpNorm<2>();

		double tol = xTolAbs_ + xTolRel_ * normB_;
		if (normdbx_ <= tol)
		{
			controlCollapsed_++;
			dfMinMax_ = std::fabs(fMin_ - f_(sorted_[numVariables_ - 1]));
			if (iterSimplexCollapsed_ == 0)
			{
				xCollapsed_ = xMin_;
				fCollapsed_ = fMin_;

			primoCollasso:

				dxa_ = v_[sorted_[numVariables_ - 1]] - v_[sorted_[0]];
				dxa_ *= 0.5;
				v_[sorted_[numVariables_ - 1]] = v_[sorted_[0]] + dxa_;

				unfeasible = false;
				iterTotal_++;
				f_(sorted_[numVariables_ - 1]) = MinFunction(v_[sorted_[numVariables_ - 1]]);

				if (unfeasible == true)
				{
					f_(sorted_[numVariables_ - 1]) = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}

				vC_ = v_[sorted_[0]] - dxa_;

				unfeasible = false;
				iterTotal_++;
				fC_ = MinFunction(vC_);
				if (unfeasible == true)
				{
					fC_ = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}

				if (fC_ < f_(sorted_[numVariables_ - 1]))
				{
					v_[sorted_[numVariables_ - 1]] = vC_;
					f_(sorted_[numVariables_ - 1]) = fC_;
				}

				iterSimplexCollapsed_ = 1;

				return;
			}
			else if (iterSimplexCollapsed_ == 1)
			{

				dxa_ = xMin_ - xCollapsed_;
				normdbx_ = dxa_.lpNorm<2>() / double(numVariables_);

				if (normdbx_ > 10. * tol && fMin_ < fCollapsed_ - 1.0001 * std::fabs(fMin_))
				{
					fCollapsed_ = fMin_;

					goto primoCollasso;
				}

				for (unsigned int i = 0; i < numVariables_; i++)
				{
					if ((i + 1) % 2 == 0)
					{
						v_[sorted_[numVariables_ - 1]](i) = v_[sorted_[0]](i) *
							(1. + .01 * double(controlCollapsed_));
						v_[sorted_[numVertices_ - 1]](i) = v_[sorted_[0]](i);
					}
					else
					{
						v_[sorted_[numVariables_ - 1]](i) = v_[sorted_[0]](i);
						v_[sorted_[numVertices_ - 1]](i) = v_[sorted_[0]](i) *
							(1. + .01 * double(controlCollapsed_));
					}
				}

				unfeasible = false;
				iterTotal_++;

				f_(sorted_[numVariables_ - 1]) = MinFunction(v_[sorted_[numVariables_ - 1]]);

				if (unfeasible == true)
				{
					f_(sorted_[numVariables_ - 1]) = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}

				unfeasible = false;
				iterTotal_++;

				f_(sorted_[numVertices_ - 1]) = MinFunction(v_[sorted_[numVertices_ - 1]]);

				if (unfeasible == true)
				{
					f_(sorted_[numVertices_ - 1]) = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}

				for (unsigned int i = 0; i < aux_.size(); i++)
					aux_[i] = f_[i];
				SortAndTrackIndicesIncreasing(sorted_, aux_);

				if (fMin_ < yTol_)
				{
					control_ = 3;
					methodStatus_ = EXIT;
					return;
				}

				iterSimplexCollapsed_ = 2;
				return;
			}
			else if (iterSimplexCollapsed_ >= 2)
			{
				dxa_ = xMin_ - xCollapsed_;
				normdbx_ = dxa_.lpNorm<2>() / double(numVariables_);

				if (normdbx_ > 100. * tol && fMin_ < fCollapsed_ - 1.0001 * std::fabs(fMin_))
				{
					fCollapsed_ = fMin_;
					goto primoCollasso;
				}

				if (dfMinMax_ < std::fabs(fMin_) * fTolRel_ + fTolAbs_ || iterSimplexCollapsed_ > 2)
				{
					control_ = 4;
					methodStatus_ = EXIT;
					return;
				}

				iterSimplexCollapsed_++;
				double pp = .01 / double(controlCollapsed_);
				v_[numVertices_ - 1] = xMin_;
				f_(numVertices_ - 1) = fMin_;

				for (unsigned int i = 0; i < numVariables_; i++)
				{
					h_(i) = pp * (1. + std::fabs(xMin_(i)));
					x_ = xMin_;
					x_(i) += h_(i);

					unfeasible = false;
					iterTotal_++;

					fh = MinFunction(x_);

					if (unfeasible == true)
					{
						x_(i) = xMin_(i) - h_(i);
						unfeasible = false;
						iterTotal_++;
						fh = MinFunction(x_);

						if (unfeasible == true)
						{
							fh = BIG_NUMBER;
						}
						else
						{
							iterTotalFeasible_++;
							iter_++;
						}
					}
					else
					{
						iterTotalFeasible_++;
						iter_++;
					}

					v_[i] = x_;
					f_(i) = fh;
				}

				for (unsigned int i = 0; i < aux_.size(); i++)
					aux_[i] = f_[i];
				SortAndTrackIndicesIncreasing(sorted_, aux_);

				if (fMin_ < yTol_)
				{
					control_ = 3;
					methodStatus_ = EXIT;
					return;
				}

				return;
			}
		}

		vR_ = vB_ + dvPvB_;

		unfeasible = false;
		iterTotal_++;
		fR_ = MinFunction(vR_);
		if (unfeasible == true)
		{
			fR_ = BIG_NUMBER;
		}
		else
		{
			iterTotalFeasible_++;
			iter_++;
		}

		if (fR_ < f_(sorted_[0]))
		{
			dxa_ = dvPvB_ * GAMMA;
			vE_ = vB_ + dxa_;

			unfeasible = false;
			iterTotal_++;
			fE_ = MinFunction(vE_);

			if (unfeasible == true)
			{
				fE_ = BIG_NUMBER;
			}
			else
			{
				iterTotalFeasible_++;
				iter_++;
			}

			if (fE_ <= fR_)
			{
				v_[sorted_[numVertices_ - 1]] = vE_;
				f_(sorted_[numVertices_ - 1]) = fE_;
			}
			else
			{
				v_[sorted_[numVertices_ - 1]] = vR_;
				f_(sorted_[numVertices_ - 1]) = fR_;
			}
		}
		else if (fR_ < fA_)
		{
			v_[sorted_[numVertices_ - 1]] = vR_;
			f_(sorted_[numVertices_ - 1]) = fR_;
		}
		else
		{
			dxa_ = dvPvB_ * BETA;

			if (fR_ < fP_)
			{
				v_[sorted_[numVertices_ - 1]] = vB_ + dxa_;

				unfeasible = false;
				iterTotal_++;

				f_(sorted_[numVertices_ - 1]) = MinFunction(v_[sorted_[numVertices_ - 1]]);

				if (unfeasible == true)
				{
					f_(sorted_[numVertices_ - 1]) = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}
			}
			else
			{
				v_[sorted_[numVertices_ - 1]] = vB_ - dxa_;

				unfeasible = false;
				iterTotal_++;

				f_(sorted_[numVertices_ - 1]) = MinFunction(v_[sorted_[numVertices_ - 1]]);

				if (unfeasible == true)
				{
					f_(sorted_[numVertices_ - 1]) = BIG_NUMBER;
				}
				else
				{
					iterTotalFeasible_++;
					iter_++;
				}
			}
		}

		for (unsigned int i = 0; i < aux_.size(); i++)
			aux_[i] = f_[i];
		SortAndTrackIndicesIncreasing(sorted_, aux_);

		fMin_ = f_(sorted_[0]);
		xMin_ = v_[sorted_[0]];

		if (fMin_ < yTol_)
		{
			control_ = 3;
			methodStatus_ = EXIT;
			return;
		}
	}

	MinimizationSimplex::MinimizationSimplex(void)
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		numVertices_ = numVariables_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		fMin_ = BIG_NUMBER;
		ptrFun = 0;
	}

	MinimizationSimplex::MinimizationSimplex
	(const MinimizationSimplex &rval)
	{
		ErrorMessage("Copy initializer is not available");
	}

	MinimizationSimplex::MinimizationSimplex
	(const Eigen::VectorXd &x00, double(*ptr)(const Eigen::VectorXd &x))
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;

		unfeasible = false;
		xL_ = x00;
		xU_ = x00;
		xL_.setConstant(-BIG_NUMBER);
		xU_.setConstant(BIG_NUMBER);
		f0_ = 0.;
		numVariables_ = x00.size();
		fi_ = MinFunction(x00); iterTotal_++;

		if (unfeasible == true)
			ErrorMessage("The starting point must be feasible");

		iterTotalFeasible_++;

		fMin_ = f0_ = fi_;
		InitialSetup(x00);
		xL_ = x00;
		xU_ = x00;
		xL_.setConstant(-BIG_NUMBER);
		xU_.setConstant(BIG_NUMBER);
	}

	void MinimizationSimplex::operator()
		(const Eigen::VectorXd &x00, double(*ptr)(const Eigen::VectorXd &x))
	{
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		if (ptrFun != 0)
			if (numVertices_ != 0)
				delete[] v_;

		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		unfeasible = false;
		fi_ = MinFunction(x00); iterTotal_++;

		if (unfeasible == true)
			ErrorMessage("The starting point must be feasible");

		iterTotalFeasible_++;

		fMin_ = f0_ = fi_;
		InitialSetup(x00);
		xL_ = x00;
		xU_ = x00;
		xL_.setConstant(-BIG_NUMBER);
		xU_.setConstant(BIG_NUMBER);
	}

	MinimizationSimplex::MinimizationSimplex
	(const Eigen::VectorXd &x00, const double ff, double(*ptr)(const Eigen::VectorXd &x))
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;

		unfeasible = false;
		fMin_ = f0_ = fi_ = ff;
		InitialSetup(x00);
		xL_ = x00;
		xU_ = x00;
		xL_.setConstant(-BIG_NUMBER);
		xU_.setConstant(BIG_NUMBER);
	}

	void MinimizationSimplex::operator()
		(const Eigen::VectorXd &x00, const double ff, double(*ptr)(const Eigen::VectorXd &x))
	{
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		if (ptrFun != 0)
			if (numVertices_ != 0)
				delete[] v_;

		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		unfeasible = true;
		fMin_ = f0_ = fi_ = ff;
		InitialSetup(x00);
		xL_ = x00;
		xU_ = x00;
		xL_.setConstant(-BIG_NUMBER);
		xU_.setConstant(BIG_NUMBER);
	}

	MinimizationSimplex::MinimizationSimplex
	(const Eigen::VectorXd &x00, double(*ptr)(const Eigen::VectorXd &x),
		const Eigen::VectorXd &xxL, const Eigen::VectorXd &xxU)
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		int constraint;
		constraint = 0;
		for (unsigned int i = 0; i < x00.size(); i++)
		{
			if (x00(i) < xxL(i))
			{
				constraint = 1;
				break;
			}
			if (x00(i) > xxU(i))
			{
				constraint = 1;
				break;
			}
		}
		InitialSetup(x00);
		xL_ = xxL;
		xU_ = xxU;

		if (constraint != 0)
			ErrorMessage("The starting point must be feasible");

		unfeasible = false;
		fi_ = MinFunction(x00);
		iterTotal_++;

		if (unfeasible == true)
			ErrorMessage("The starting point must be feasible");

		iterTotalFeasible_++;
		fMin_ = f0_ = fi_;
	}

	void MinimizationSimplex::operator()
		(const Eigen::VectorXd &x00, double(*ptr)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xxL, const Eigen::VectorXd &xxU)
	{
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		if (ptrFun != 0)
			if (numVertices_ != 0)
				delete[] v_;

		ptrFun = ptr;
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		int constraint;
		constraint = 0;
		for (unsigned i = 0; i < x00.size(); i++)
		{
			if (x00(i) < xxL(i))
			{
				constraint = 1;
				break;
			}
			if (x00(i) > xxU(i))
			{
				constraint = 1;
				break;
			}
		}

		if (constraint != 0)
			ErrorMessage("The starting point must be feasible");

		unfeasible = false;

		fi_ = MinFunction(x00);
		iterTotal_++;

		if (unfeasible == true)
			ErrorMessage("The starting point must be feasible");

		iterTotalFeasible_++;

		fMin_ = f0_ = fi_;
		InitialSetup(x00);
		xL_ = xxL;
		xU_ = xxU;
	}

	MinimizationSimplex::MinimizationSimplex
	(const Eigen::VectorXd &x00, const double ff, double(*ptr)(const Eigen::VectorXd& x),
		const Eigen::VectorXd &xxL, const Eigen::VectorXd &xxU)
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		fMin_ = f0_ = fi_ = ff;
		ptrFun = ptr;
		InitialSetup(x00);
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		xL_ = xxL;
		xU_ = xxU;
	}

	void MinimizationSimplex::operator()
		(const Eigen::VectorXd &x00, const double ff, double(*ptr)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xxL, const Eigen::VectorXd &xxU)
	{
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		if (ptrFun != 0)
			if (numVertices_ != 0)
				delete[] v_;

		fMin_ = f0_ = fi_ = ff;
		ptrFun = ptr;
		InitialSetup(x00);
		functionSimplexType_ = UNCONSTRAINED_MULTI;
		xL_ = xxL;
		xU_ = xxU;
	}

	MinimizationSimplex::MinimizationSimplex
	(const Eigen::VectorXd &xxS, const double ffS, double(*ptr)(const Eigen::VectorXd &x),
		const int iv1, const double ttL1, const double ttU1,
		const int iv2, const double ttL2, const double ttU2)
	{
		yTol_ = -BIG_NUMBER;
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		fMin_ = f0_ = fi_ = ffS;
		ptrFun = ptr;
		z0_ = xxS;
		iVar1_ = iv1;
		iVar2_ = iv2;
		Eigen::VectorXd x00(2);
		x00(0) = xxS(iv1 - 1);
		x00(1) = xxS(iv2 - 1);
		InitialSetup(x00);
		functionSimplexType_ = TWO_IV1_IV2;
		xL_.resize(2);
		xU_.resize(2);
		xL_(0) = ttL1;
		xL_(1) = ttL2;
		xU_(0) = ttU1;
		xU_(1) = ttU2;
	}

	void MinimizationSimplex::operator ()
		(const Eigen::VectorXd &xxS, const double ffS, double(*ptr)(const Eigen::VectorXd &x),
			const int iv1, const double ttL1, const double ttU1,
			const int iv2, const double ttL2, const double ttU2)
	{
		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		iterTotal_ = 0;
		iterTotalFeasible_ = 0;
		if (ptrFun != 0)
			if (numVertices_ != 0)
				delete[] v_;

		fMin_ = f0_ = fi_ = ffS;
		ptrFun = ptr;
		z0_ = xxS;
		iVar1_ = iv1;
		iVar2_ = iv2;
		Eigen::VectorXd x00(2);
		x00(0) = xxS(iv1 - 1);
		x00(1) = xxS(iv2 - 1);
		InitialSetup(x00);
		functionSimplexType_ = TWO_IV1_IV2;
		xL_.resize(2);
		xU_.resize(2);
		xL_(0) = ttL1;
		xL_(1) = ttL2;
		xU_(0) = ttU1;
		xU_(1) = ttU2;
	}

	void MinimizationSimplex::operator()
		(const Eigen::VectorXd &x00)
	{
		if (ptrFun == 0)
			ErrorMessage("Error in Initialization: objective function not appropriate");

		if (x00.size() != numVariables_)
			ErrorMessage("The new starting point dimensions do not match previous ones");

		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		xMin_ = xi_ = x0_ = x00;
		iter_ = 0;
		unfeasible = false;
		iterTotal_++;
		fi_ = MinFunction(x0_);

		if (unfeasible == true)
			ErrorMessage("The starting point must be feasible");

		iterTotalFeasible_++;

		fMin_ = f0_ = fi_;
		v_[numVertices_ - 1] = xi_;
		f_(numVertices_ - 1) = fi_;
		control_ = -1;
		stop_ = 0;
		methodStatus_ = START;
	}

	void MinimizationSimplex::operator() (const Eigen::VectorXd &x00, const double ff)
	{
		if (ptrFun == 0)
			ErrorMessage("Error in initialization: objective function not appropriate");

		if (x00.size() != numVariables_)
			ErrorMessage("The new starting point dimensions do not match previous ones");

		controlCollapsed_ = 0;
		iterSimplexCollapsed_ = 0;
		xMin_ = xi_ = x0_ = x00;
		iter_ = 0;
		fMin_ = f0_ = fi_ = ff;
		v_[numVertices_ - 1] = xi_;
		f_(numVertices_ - 1) = fi_;
		control_ = -1;
		stop_ = 0;
		methodStatus_ = START;
	}

	void MinimizationSimplex::Restart(const Eigen::VectorXd &x00)
	{
		iterSimplexCollapsed_ = 0;
		controlCollapsed_ = 0;
		(*this)(x00);
	}

	MinimizationSimplex::~MinimizationSimplex(void)
	{
		if (numVertices_ != 0)
			delete[] v_;
	}

	void MinimizationSimplex::FinalSummary(std::ostream& out)
	{
		out << "Starting point" << std::endl;
		for (unsigned int i = 0; i < x0_.size(); i++)
			out << x0_(i) << std::endl;
		out << std::endl;

		out << "Objective Function value in starting point" << std::endl;
		out << f0_ << std::endl;
		out << std::endl;

		out << "Number of iterations in the last call" << std::endl;
		out << iter_ << std::endl;
		out << std::endl;

		out << "Total number of iterations (feasible and unfeasible)" << std::endl;
		out << iterTotal_ << std::endl;
		out << std::endl;

		out << "Total number of feasible iterations" << std::endl;
		out << iterTotalFeasible_ << std::endl;
		out << std::endl;

		out << "Objective Function value in final point" << std::endl;
		out << fMin_ << std::endl;
		out << std::endl;

		out << "Optimal point" << std::endl;
		for (unsigned int i = 0; i < xMin_.size(); i++)
			out << xMin_(i) << std::endl;
		out << std::endl;

		out << "Function values in all final vertices" << std::endl;
		for (unsigned int i = 0; i < f_.size(); i++)
			out << f_(i) << std::endl;
		out << std::endl;

		out << "Exit status" << std::endl;
		if (control_ < 0)
		{
			out << " * Search failed before reaching the solution" << std::endl;
			if (control_ == -3)
				out << "   Unfeasible points" << std::endl;
		}
		else if (control_ == 0)
		{
			out << " * Found problems" << std::endl;
		}
		else
		{
			out << " * Search successfully reached the solution" << std::endl;
			if (control_ == 1)
				out << "   Maximum number of iterations" << std::endl;
			if (control_ == 2)
				out << "   Satisfied test on dx" << std::endl;
			if (control_ == 3)
				out << "   Satisfied fSolution <= fTol" << std::endl;
			if (control_ == 4)
				out << "   Satisfied test on the function values for simplex" << std::endl;
		}
		out << std::endl;
	}

	char MinimizationSimplex::operator ()(unsigned int ni)
	{
		if (ni <= 0)
			ni = MAX_DEFAULT_ITERATIONS;

		numIterations_ = ni;
		iter_ = 0;

		char ok;
		if (stop_ == 1)
		{
			control_ = -2;
			std::cout << "WARNING message: MinimizationSimplex used with stop_ == 1" << std::endl;
			return control_;
		}

		try
		{
			ok = MinimumSolve();
		}
		catch (...)
		{
			std::cout << "WARNING message: MinimizationSimplex exception handling" << std::endl;
			stop_ = 1;
			control_ = -2;
			throw;
		}

		return ok;
	}

	char MinimizationSimplex::MinimumSolve()
	{
		while (1)
		{
			if (iter_ >= numIterations_)
			{
				control_ = 1;
				break;
			}
			if (methodStatus_ == EXIT)
				break;
			switch (methodStatus_)
			{

			case EXIT:
				ErrorMessage("EXIT status detected");
				break;

			case START:
				Start();
				break;

			case SIMPLEX_ACTION:
				SimplexAction();
				break;
			}
		}
		return control_;
	}

	double MinimizationSimplex::GetSolution(Eigen::VectorXd& xx)
	{
		xx = xMin_;
		return fMin_;
	}

	void MinimizationSimplex::GetXSolution(Eigen::VectorXd& xS)
	{
		xS = xMin_;
	}

	void MinimizationSimplex::ErrorMessage(std::string const message)
	{
		std::cout << "MinimizationSimplex fatal error message" << std::endl;
		std::cout << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(-1);
	}
}
