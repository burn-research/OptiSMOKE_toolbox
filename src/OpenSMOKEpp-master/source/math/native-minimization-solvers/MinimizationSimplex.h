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

#ifndef OpenSMOKEpp_MinimizationSimplex_H
#define OpenSMOKEpp_MinimizationSimplex_H

#include <Eigen/Dense>
#include <vector>

//!  A class for multidimensional minimization based on the Simplex method
/*!
The purpose of this class is to carry out multidimensional minimization
based on the Simplex method
This class is derived from the BzzMinimizationRobust class available in the
BzzMath library (G. Buzzi-Ferraris, Metodi numerici e software in C++, Pearson Education Italia, 1998)
*/

namespace OpenSMOKE
{
	class MinimizationSimplex
	{

	public:

		/**
		*@brief Default constructor
		*/
		MinimizationSimplex();

		/**
		*@brief Copy constructor
		*@param rval object to be copied
		*/
		MinimizationSimplex(const MinimizationSimplex &rval);

		/**
		*@brief Default constructor
		*@param x0 first-guess solution (starting point)
		*@param f  pointer to function to be minimized
		*/
		MinimizationSimplex(const Eigen::VectorXd &x0, double(*f)(const Eigen::VectorXd &x));

		/**
		*@brief Default constructor
		*@param x0 first-guess solution (starting point)
		*@param f0 function value at first-guess solution
		*@param f  pointer to function to be minimized
		*/
		MinimizationSimplex(const Eigen::VectorXd &x0, const double f0, double(*f)(const Eigen::VectorXd &x));

		/**
		*@brief Default constructor
		*@param x0 first-guess solution (starting point)
		*@param f  pointer to function to be minimized
		*@param xL minimum constraints
		*@param xU minimum constraints
		*/
		MinimizationSimplex(const Eigen::VectorXd &x0, double(*f)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xL, const Eigen::VectorXd &xU);

		/**
		*@brief Default constructor
		*@param x0 first-guess solution (starting point)
		*@param f  pointer to function to be minimized
		*@param f0 function value at first-guess solution
		*@param xL minimum constraints
		*@param xU minimum constraints
		*/
		MinimizationSimplex(const Eigen::VectorXd &x0, const double f0, double(*f)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xL, const Eigen::VectorXd &xU);

		/**
		*@brief Default constructor (special case, 2 parameters)
		*@param xS first-guess solution (starting point)
		*@param fS function value at first-guess solution
		*@param f  pointer to function to be minimized
		*@param iv1 index of first variable (starting from 1)
		*@param tL1 minimum constraint on first variable
		*@param tU1 maximum constraint on first variable
		*@param iv2 index of second variable (starting from 1)
		*@param tL2 minimum constraint on first variable
		*@param tU2 maximum constraint on second variable
		*/
		MinimizationSimplex(const Eigen::VectorXd &xS, const double fS, double(*f)(const Eigen::VectorXd &x),
			const int iv1, const double tL1, const double tU1,
			const int iv2, const double tL2, const double tU2);

		/**
		*@brief Default destructor
		*/
		~MinimizationSimplex();

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*@param f  pointer to function to be minimized
		*/
		void operator () (const Eigen::VectorXd& x0, double(*f)(const Eigen::VectorXd &x));

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*@param f0 function value at first-guess solution
		*@param f  pointer to function to be minimized
		*/
		void operator () (const Eigen::VectorXd& x0, const double f0, double(*ptr)(const Eigen::VectorXd &x));

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*@param f  pointer to function to be minimized
		*@param xL minimum constraints
		*@param xU minimum constraints
		*/
		void operator () (const Eigen::VectorXd &x0, double(*f)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xL, const Eigen::VectorXd&xU);

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*@param f0 function value at first-guess solution
		*@param f  pointer to function to be minimized
		*@param xL minimum constraints
		*@param xU minimum constraints
		*/
		void operator () (const Eigen::VectorXd &x0, const double f0, double(*ptr)(const Eigen::VectorXd &x),
			const Eigen::VectorXd &xL, const Eigen::VectorXd &xU);

		/**
		*@brief Default initializer
		*@param xS first-guess solution (starting point)
		*@param fS function value at first-guess solution
		*@param f  pointer to function to be minimized
		*@param iv1 index of first variable (starting from 1)
		*@param tL1 minimum constraint on first variable
		*@param tU1 maximum constraint on first variable
		*@param iv2 index of second variable (starting from 1)
		*@param tL2 minimum constraint on first variable
		*@param tU2 maximum constraint on second variable
		*/
		void operator () (const Eigen::VectorXd &xS, const double fS, double(*f)(const Eigen::VectorXd &x),
			const int iv1, const double tL1, const double tU1,
			const int iv2, const double tL2, const double tU2);

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*@param f0 unction value at first-guess solution
		*/
		void operator () (const Eigen::VectorXd &x0, const double f0);

		/**
		*@brief Default initializer
		*@param x0 first-guess solution (starting point)
		*/
		void operator ()(const Eigen::VectorXd &x0);

		/**
		*@brief Restart the solution process
		*@param x0 first-guess solution (starting point)
		*/
		void Restart(const Eigen::VectorXd &x0);

		/**
		*@brief Recover the solution and returns the minimum value of objective function
		*@param x minimum solution
		*@return the minimum objective function
		*/
		double GetSolution(Eigen::VectorXd& xMin);

		/**
		*@brief Recover the solution
		*@param x minimum solution
		*/
		void GetXSolution(Eigen::VectorXd& xMin);

		/**
		*@brief Returns the minimum value of objective function
		*@return the minimum objective function
		*/
		inline double GetBzzMinimumF() const { return fMin_; }

		/**
		*@brief Returns the total number of iterations
		*@return the total number of iterations
		*/
		inline int TotalIterations() const { return iterTotal_; }

		/**
		*@brief Returns the total number of feasible iterations
		*@return the total number of feasible iterations
		*/
		inline int TotalFeasibleIterations() const { return iterTotalFeasible_; }

		/**
		*@brief Sets the tolerance
		*@param yt the tolerance
		*/
		void SetTolF(const double yt = -BIG_NUMBER);

		/**
		*@brief Sets the absolute tolerance on objective function (default 1e-25)
		*@param tolAbs the absolute tolerance on objective function
		*/
		void SetTolAbsF(const double tolAbs);

		/**
		*@brief Sets the relative tolerance on objective function (default 1e-13)
		*@param tolRel the relative tolerance on objective function
		*/
		void SetTolRelF(const double tolRelF);

		/**
		*@brief Sets the absolute tolerance on optimization unknowns (default 1e-25)
		*@param tolAbs the absolute tolerance on optimization unknowns
		*/
		void SetTolAbsX(const double tolAbsX);

		/**
		*@brief Sets the relative tolerance on optimization unknowns (default 1e-13)
		*@param tolAbs the relative tolerance on optimization unknowns
		*/
		void SetTolRelX(const double tolRelX);

		/**
		*@brief Calculates the solution
		*@return a flag indicating if the solution was successfully reached
		*/
		char operator ()();

		/**
		*@brief Calculates the solution
		*@param ni the maximum number of iterations
		*@return a flag indicating if the solution was successfully reached
		*/
		char operator ()(const unsigned int ni);

		/**
		*@brief Prints a final summary
		*@param out the stream were to write the final summary
		*/
		void FinalSummary(std::ostream& out);


	public:

		static bool unfeasible;		//!< boolean variable to specify unfeasible set of unknowns


	private:

		/**
		*@brief Objective function (external function)
		*@param x current values of unknown parameters
		*/
		double(*ptrFun)(const Eigen::VectorXd& x);

		/**
		*@brief Setup operations
		*@param x0 starting value
		*/
		void InitialSetup(const Eigen::VectorXd& x0);

		/**
		*@brief Start the minimization algorithm
		*/
		void Start();

		/**
		*@brief Internal calculation of baricenter
		*/
		void Baricenter();

		/**
		*@brief Step in minimization algorithm
		*/
		double MinFunction(const Eigen::VectorXd& x);

		/**
		*@brief Simplex algorithm
		*/
		void SimplexAction();

		/**
		*@brief Minimization algorithm
		*/
		char MinimumSolve();

		/**
		*@brief Fatal error message
		*@param message to be printed on the screen
		*/
		void ErrorMessage(const std::string message);


	private:

		enum MethodStatusSimplexPlus { EXIT, START, SIMPLEX_ACTION } methodStatus_;
		enum FunctionTypeSimplex { TWO_IV1_IV2, UNCONSTRAINED_MULTI, CONSTRAINED_MULTI } functionSimplexType_;

		unsigned int	numVariables_;
		unsigned int	numVertices_;
		unsigned int	niter_;
		unsigned int	iter_;
		unsigned int	iterTotal_;
		unsigned int	iterTotalFeasible_;
		unsigned int	numIterations_;
		unsigned int	iterSimplexCollapsed_;

		char control_;
		char controlCollapsed_;
		char stop_;

		Eigen::VectorXd	x0_;	//!< starting unknowns
		Eigen::VectorXd	xi_;	//!< current unknowns
		Eigen::VectorXd	xMin_;	//!< minimum unknowns
		Eigen::VectorXd	x_;		//!< generic unknowns

		Eigen::VectorXd	xCollapsed_;
		Eigen::VectorXd	h_;
		Eigen::VectorXd	*v_;		//!< vector of vertices coordinates
		Eigen::VectorXd	f_;			//!< function value in vertices
		Eigen::VectorXd	d_;			//!< control callapse
		Eigen::VectorXd	vB_;		//!< base
		Eigen::VectorXd	vR_;		//!< reflection
		Eigen::VectorXd	vE_;		//!< expansion
		Eigen::VectorXd	vC_;		//!< contraction
		Eigen::VectorXd	vA_;		//!<
		Eigen::VectorXd	vP_;		//!< projected vertex 
		Eigen::VectorXd	dvPvB_;		//!< 
		Eigen::VectorXd	dxa_;		//!< 
		Eigen::VectorXd	xL_;		//!< minimum constraints
		Eigen::VectorXd	xU_;		//!< maximum constraints

		Eigen::VectorXd	z0_;		//!< 
		Eigen::VectorXd zi_;		//!< 

		int iVar1_;					//!< 
		int iVar2_;					//!< 

		double	fi_;				//!< objective funtion current value
		double	f0_;				//!< objective funtion starting value
		double	fMin_;				//!< objective funtion minimum value

		double	fCollapsed_;		//!< 
		double	dfMinMax_;			//!< 
		double	normdbx_;			//!< 
		double	normB_;				//!< 

		double	xTolRel_;			//!< relative tolerance on unknowns
		double	xTolAbs_;			//!< absolute tolerance on unknowns
		double	fTolRel_;			//!< relative tolerance on unknowns
		double	fTolAbs_;			//!< absolute tolerance on unknowns

		double	fB_;				//!< base
		double	fR_;				//!< reflection
		double	fE_;				//!< expansion
		double	fA_;				//!< 
		double	fP_;				//!< 
		double	fC_;				//!<  contraction
		double  yTol_;				//!< 

		std::vector<unsigned int>	sorted_;	//!< 
		std::vector<double>			aux_;		//!< 

	private:

		const static double ETA_REL_SIMPLEX;
		const static double ETA_ABS_SIMPLEX;
		const static double BETA;
		const static double GAMMA;
		const static int    MAX_DEFAULT_ITERATIONS;
		const static double X_TOL_REL;
		const static double X_TOL_ABS;
		const static double F_TOL_REL;
		const static double F_TOL_ABS;
		const static double BIG_NUMBER;
	};
}

#include "MinimizationSimplex.hpp"

#endif // OpenSMOKEpp_MinimizationSimplex_H
