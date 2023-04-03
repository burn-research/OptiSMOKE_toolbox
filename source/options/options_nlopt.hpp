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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	      |
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

namespace OptiSMOKE
{
    options_nlopt::options_nlopt(){
		global_ = true;
		local_ = false;
		max_iterations_ = 1000;
		max_function_evaluations_ = 1000;
		convergence_tolerance_ = 1e-6;		
		solution_target_ = 1e-4;

		variant_ = "NONE";
    }

    void options_nlopt::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, 
											std::string dictionary_name)
	{

        dictionary_manager(dictionary_name).SetGrammar(nlopt_options_grammar_);


		dictionary_manager(dictionary_name).ReadString("@Algorithm", algorithm_);

		dictionary_manager(dictionary_name).ReadString("@Variant", variant_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@Global"))
			dictionary_manager(dictionary_name).ReadBool("@Global", global_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@Local"))
			dictionary_manager(dictionary_name).ReadBool("@Local", local_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@MaxIterations"))
			dictionary_manager(dictionary_name).ReadInt("@MaxIterations", max_iterations_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@MaxFunctionEvaluations"))
			dictionary_manager(dictionary_name).ReadInt("@MaxFunctionEvaluations", max_function_evaluations_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@ConvergenceTolerance"))
			dictionary_manager(dictionary_name).ReadDouble("@ConvergenceTolerance", convergence_tolerance_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@SolutionTarget"))
			dictionary_manager(dictionary_name).ReadDouble("@SolutionTarget", solution_target_);
	}

	void options_nlopt::SetupAlgorithm(){

		if (algorithm_ == "DIRECT"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::GN_DIRECT;
			else if (variant_ == "L")
				algo_int_ = nlopt::GN_DIRECT_L;
			else if (variant_ == "L_RAND")
				algo_int_ = nlopt::GN_DIRECT_L_RAND;
			else if (variant_ == "NOSCAL")
				algo_int_ = nlopt::GN_DIRECT_NOSCAL;
			else if (variant_ == "L_NOSCAL")
				algo_int_ = nlopt::GN_DIRECT_L_NOSCAL;
			else if (variant_ == "L_RAND_NOSCAL")
				algo_int_ = nlopt::GN_DIRECT_L_RAND_NOSCAL;
			else if (variant_ == "ORIG")
				algo_int_ = nlopt::GN_ORIG_DIRECT;
			else if (variant_ == "L_ORIG")
				algo_int_ = nlopt::GN_ORIG_DIRECT_L;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for DIRECT algorithms: \
						NONE | L | L_RAND | NOSCAL | L_NOSCAL | L_RAND_NOSCAL | ORIG | L_ORIG");
		}
		else if (algorithm_ == "CRS"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::GN_CRS2_LM;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for CRS algorithms: NONE ");
		}
		else if (algorithm_ == "MLSL"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::G_MLSL;
			else if (variant_ == "LDS")
				algo_int_ = nlopt::G_MLSL_LDS;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for MLSL algorithms: NONE | LDS");
		}
		else if (algorithm_ == "STOGO"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::GD_STOGO;
			else if (variant_ == "RAND")
				algo_int_ = nlopt::GD_STOGO_RAND;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for STOGO algorithms: NONE | RAND");
		}
		else if (algorithm_ == "ISRES"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::GN_ISRES;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for ISRES algorithms: NONE");
		}
		else if (algorithm_ == "ESCH"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::GN_ESCH;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for ESCH algorithms: NONE");
		}
		else if (algorithm_ == "COBYLA"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_COBYLA;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for COBYLA algorithms: NONE");
		}
		else if (algorithm_ == "BOBYQA"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_BOBYQA;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for BOBYQA algorithms: NONE");
		}
		else if (algorithm_ == "NEWUOA"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_NEWUOA_BOUND;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for NEWUOA algorithms: NONE");
		}
		else if (algorithm_ == "PRAXIS"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_PRAXIS;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for PRAXIS algorithms: NONE");
		}
		else if (algorithm_ == "NELDERMEAD"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_NELDERMEAD;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for NELDERMEAD algorithms: NONE");
		}
		else if (algorithm_ == "SBPLX"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LN_SBPLX;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for SBPLX algorithms: NONE");
		}
		else if (algorithm_ == "SLSQP"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LD_SLSQP;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for SLSQP algorithms: NONE");
		}
		else if (algorithm_ == "LBFGS"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LD_LBFGS;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for LBFGS algorithms: NONE");
		}
		else if (algorithm_ == "TNEWTON_PRECOND"){
			if (variant_ == "NONE" || variant_ == "")
				algo_int_ = nlopt::LD_TNEWTON_PRECOND_RESTART;
			else if (variant_ == "NO_RESTART")
				algo_int_ = nlopt::LD_TNEWTON_PRECOND;
			else if (variant_ == "NO_PRECOND")
				algo_int_ = nlopt::LD_TNEWTON_RESTART;
			else if (variant_ == "NO_RESTART_NO_PRECOND")
				algo_int_ = nlopt::LD_TNEWTON;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: NONE | NO_RESTART | NO_PRECOND | NO_RESTART_NO_PRECOND");
		}
		else if (algorithm_ == "SLM_VAR"){
			if (variant_ == "VAR2")
				algo_int_ = nlopt::LD_VAR2;
			else if (variant_ == "VAR1")
				algo_int_ = nlopt::LD_VAR1;
			else
				OptiSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: VAR2 | VAR1");
		}
		else
			OptiSMOKE::FatalErrorMessage("Unknown algorithm. Avaliable are: DIRECT | ");
	}	
}