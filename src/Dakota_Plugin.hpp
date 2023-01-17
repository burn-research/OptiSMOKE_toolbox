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
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst                                    |
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

namespace SIM {

void OpenSMOKEDirectApplicInterface::Setup_Plugin()
{
	// Debug parameter: If true will print kinetic parameters before and after changing them / 
	// Print out the simulation values (only first value if more than one species considered)
	Debug_ChangeParam = false;

	// Setup evaluation number
	eval_nr=0;

	// Extract data from input file
	// It Calls the ReadInfo function into the ObjectInput2, so now all the variables stored there i have within this object.
	ObjectInput2.ReadInfo(plugin_input_file_,true);
	// Adjust number of threads if parametric analysis is required
	#if defined(_OPENMP)
	{
		// Indicates that the number of threads available in subsequent parallel region can be adjusted by the run time.
		// A value that indicates if the number of threads available in subsequent parallel region can be adjusted by 
		// the runtime. If nonzero, the runtime can adjust the number of threads, if zero, the runtime will not 
		// dynamically adjust the number of threads.
		omp_set_dynamic(0);
		// Set the number of threads
		omp_set_num_threads(ObjectInput2.number_of_threads);
	}
	#endif

	// Initialize opensmoke reactor objects
	// Number of opensmoke input files
	nInputFiles = ObjectInput2.list_of_opensmoke_input_files.size();
	// type of reactor is a variable that stores the type of reactor only if i use one type of reactor! 
	// So even if i do want to change, here it's still fine, because it just resizes the objects of 
	// Initialize opensmoke reactor objects if different reactors are used
	

	std::cout<<"N of Batches " << ObjectInput2.N_of_batch_datasets    <<std::endl;
	std::cout<<"N of PFR     " << ObjectInput2.N_of_plugflow_datasets <<std::endl;
	std::cout<<"N of PSR     " << ObjectInput2.N_of_psr_datasets      <<std::endl;
	std::cout<<"N of Flames  " << ObjectInput2.N_of_laminar_flame_datasets      <<std::endl;
	std::cout<<"N of KTP     " << ObjectInput2.N_of_direct_measurement_dataset      <<std::endl;
	std::cout<<"N of inputs  " << ObjectInput2.list_of_opensmoke_input_files.size()      <<std::endl;
	

	int offset=0;
	if (ObjectInput2.N_of_batch_datasets > 0)
	{	
		batch_reactors.resize(ObjectInput2.N_of_batch_datasets);

		for (int i = 0; i < ObjectInput2.N_of_batch_datasets; i++)
		{
			batch_reactors[i] = new OpenSMOKE::BatchReactor_Plugin[ObjectInput2.list_of_opensmoke_input_files[i+offset].size()];			
		}
	} 
	offset += ObjectInput2.N_of_batch_datasets;
	
	if (ObjectInput2.N_of_plugflow_datasets > 0)
	{
		plugflow_reactors.resize(ObjectInput2.N_of_plugflow_datasets);
		//std::cout<<"Size of PFR dataset     " << plugflow_reactors.size()      <<std::endl;
		for (int i = 0; i < ObjectInput2.N_of_plugflow_datasets; i++)
		{
			plugflow_reactors[i] = new OpenSMOKE::PlugFlowReactor_Plugin[ObjectInput2.list_of_opensmoke_input_files[i+offset].size()];
			//std::cout<<"Size of set PFR in dataset" << i <<" is  " << ObjectInput2.list_of_opensmoke_input_files[i+offset].size()      <<std::endl;
		}
	} 
	offset += ObjectInput2.N_of_plugflow_datasets;
	//std::cout<<"offset is " << offset      <<std::endl;
	if (ObjectInput2.N_of_psr_datasets > 0)
	{
		psr_reactors.resize(ObjectInput2.N_of_psr_datasets);
		//std::cout<<"Size of PFR dataset     " << psr_reactors.size()      <<std::endl;
		for (int i = 0; i < ObjectInput2.N_of_psr_datasets; i++)
		{
			psr_reactors[i] = new OpenSMOKE::PerfectlyStirredReactor_Plugin[ObjectInput2.list_of_opensmoke_input_files[i+offset].size()];
		        //std::cout<<"Number of set PSR in dataset" << i <<" is     " << ObjectInput2.list_of_opensmoke_input_files[i+offset].size()      <<std::endl;
		}
	}
	offset += ObjectInput2.N_of_psr_datasets;

	//TODO// remove this, is useless.
	if (ObjectInput2.N_of_laminar_flame_datasets > 0)
	{
		//laminar_flames.resize(ObjectInput2.N_of_laminar_flame_datasets);
		//std::cout<<"Size of Laminar Flame dataset     " << laminar_flames.size()      <<std::endl;
		//for (int i = 0; i < ObjectInput2.N_of_laminar_flame_datasets; i++)
		//{
		//	laminar_flames[i].resize(ObjectInput2.list_of_opensmoke_input_files[i+offset].size());
		        //std::cout<<"Number of set Laminar Flames in dataset" << i <<" is     " << ObjectInput2.list_of_opensmoke_input_files[i+offset].size()      <<std::endl;
		//}
	}


	// ******************************************************** //
	//															//
	//				CREATE CONSTRAINTS FOR REACTIONs			//
	//															//
	// ******************************************************** //

	// Initializes vector of uncertainty factors for Direct reactions
	std::vector<double> f_factors;
	f_factors.resize(ObjectInput2.nominalkineticsMapXML->NumberOfReactions());

	for (int i=0; i< ObjectInput2.list_of_uncertainty_factors.size(); i++)
	{
		f_factors[ObjectInput2.list_of_target_uncertainty_factors[i]-1] = ObjectInput2.list_of_uncertainty_factors[i];
	}
	// Get nominal parameter values for Direct reactions 
	std::vector<double> lnA_0;
	std::vector<double> Beta_0;
	std::vector<double> E_over_R_0;	
	lnA_0.resize(ObjectInput2.nominalkineticsMapXML->NumberOfReactions());
	Beta_0.resize(ObjectInput2.nominalkineticsMapXML->NumberOfReactions());
	E_over_R_0.resize(ObjectInput2.nominalkineticsMapXML->NumberOfReactions());

	for (int i = 0; i < ObjectInput2.nominalkineticsMapXML->NumberOfReactions(); i++)
	{
		lnA_0[i] = std::log(ObjectInput2.nominalkineticsMapXML->A(i));
		Beta_0[i] = ObjectInput2.nominalkineticsMapXML->Beta(i);
		E_over_R_0[i] = ObjectInput2.nominalkineticsMapXML->E_over_R(i);
	}
	// Initializes vector of uncertainty factors for P_inf limit in TROE reactions
	std::vector<double> f_factors_inf;
	f_factors_inf.resize(ObjectInput2.nominalkineticsMapXML->NumberOfFallOffReactions());
	for (int i=0; i< ObjectInput2.list_of_uncertainty_factors_inf.size(); i++)
	{
		int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),ObjectInput2.list_of_target_uncertainty_factors_inf[i])-ObjectInput2.indices_of_falloff_reactions.begin();
		f_factors_inf[pos_FallOff_Reaction] = ObjectInput2.list_of_uncertainty_factors_inf[i];
	}
	std::vector<double> lnA_0_inf;
	std::vector<double> Beta_0_inf;
	std::vector<double> E_over_R_0_inf;
	lnA_0_inf.resize(ObjectInput2.nominalkineticsMapXML->NumberOfFallOffReactions());
	Beta_0_inf.resize(ObjectInput2.nominalkineticsMapXML->NumberOfFallOffReactions());
	E_over_R_0_inf.resize(ObjectInput2.nominalkineticsMapXML->NumberOfFallOffReactions());
	// Get nominal parameter values for P_inf limit in TROE reactions
	for (int i = 0; i < ObjectInput2.nominalkineticsMapXML->NumberOfFallOffReactions(); i++)
	{
		lnA_0_inf[i] = std::log(ObjectInput2.nominalkineticsMapXML->A_falloff_inf(i));
		Beta_0_inf[i] = ObjectInput2.nominalkineticsMapXML->Beta_falloff_inf(i);
		E_over_R_0_inf[i] = ObjectInput2.nominalkineticsMapXML->E_over_R_falloff_inf(i);
	}
	// COMPUTE NOMINAL, MAX, and MIN reaction rates for all target reactions	
	std::vector<std::vector<double>> k_0;
	std::vector<double> T_span = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
	k_0.resize(ObjectInput2.list_of_target_uncertainty_factors.size());
	k_upper.resize(ObjectInput2.list_of_target_uncertainty_factors.size());
	k_lower.resize(ObjectInput2.list_of_target_uncertainty_factors.size());

	// COMPUTE MAX/MIN reaction rates according to temperature dependent or constant uncertainty factors for Direct reactions//
	if (ObjectInput2.udc_bool == true)
	{
		//std::cout<< "It enters the part for user-defined constaints" <<std::endl;
		for (int j=0; j < ObjectInput2.list_of_target_uncertainty_factors.size(); j++)
		{
			k_0[j].resize(T_span.size());
			k_upper[j].resize(T_span.size());
			k_lower[j].resize(T_span.size());
			for (int i=0; i < T_span.size(); i++)
			{
				k_0[j][i] = std::exp(lnA_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) * std::pow(T_span[i],Beta_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) * std::exp((-1*E_over_R_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1])/T_span[i]);
				k_upper[j][i] = k_0[j][i] * std::pow(10,ObjectInput2.ud_constraints[j][i]) * ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
				k_lower[j][i] = k_0[j][i] * std::pow(10,-1*ObjectInput2.ud_constraints[j][i]) / ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
				std::cout<< "The value of f for the constraint j:" << j << " i:" << i << " is	" << ObjectInput2.ud_constraints[j][i] <<std::endl;
			}
		}
	}
	else
	{
		for (int j=0; j < ObjectInput2.list_of_target_uncertainty_factors.size(); j++)
		{
			k_0[j].resize(T_span.size());
			k_upper[j].resize(T_span.size());
			k_lower[j].resize(T_span.size());
			for (int i=0; i < T_span.size(); i++)
			{
				k_0[j][i] = std::exp(lnA_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) * std::pow(T_span[i],Beta_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) * std::exp((-1*E_over_R_0[ObjectInput2.list_of_target_uncertainty_factors[j]-1])/T_span[i]);
				k_upper[j][i] = k_0[j][i] * std::pow(10,f_factors[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) * ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
				k_lower[j][i] = k_0[j][i] * std::pow(10,-1*f_factors[ObjectInput2.list_of_target_uncertainty_factors[j]-1]) / ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
			}

			for (int i=0; i < T_span.size(); i++)
			{
				std::cout<< "The value of f for the lower constraint j:" << j << " i:" << i << " is     " << k_lower[j][i] <<std::endl;
			}

            for (int i=0; i < T_span.size(); i++)
            {
                std::cout<< "The value of f for the upper constraint j:" << j << " i:" << i << " is     " << k_upper[j][i] <<std::endl;
            }
		}
	}

	// COMPUTE MAX/MIN reaction rates according to constant uncertainty factors for P_inf limit in TROE reactions //
	std::vector<std::vector<double>> k_0_inf;
	k_0_inf.resize(ObjectInput2.list_of_target_uncertainty_factors.size());
	k_upper_inf.resize(ObjectInput2.list_of_target_uncertainty_factors.size());
	k_lower_inf.resize(ObjectInput2.list_of_target_uncertainty_factors.size());
	for (int j=0; j <ObjectInput2.list_of_target_uncertainty_factors_inf.size(); j++)
	{
		k_0_inf[j].resize(T_span.size());
		k_upper_inf[j].resize(T_span.size());
		k_lower_inf[j].resize(T_span.size());
		int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),ObjectInput2.list_of_target_uncertainty_factors_inf[j])-ObjectInput2.indices_of_falloff_reactions.begin();
		for (int i=0; i < T_span.size(); i++)
		{
			k_0_inf[j][i] = std::exp(lnA_0_inf[pos_FallOff_Reaction]) * std::pow(T_span[i],Beta_0_inf[pos_FallOff_Reaction]) * std::exp((-1*E_over_R_0_inf[pos_FallOff_Reaction])/T_span[i]);
			k_upper_inf[j][i] = k_0_inf[j][i] * std::pow(10,f_factors_inf[pos_FallOff_Reaction])* (double)ObjectInput2.AcceptedSigmaKDistribution/2;
			k_lower_inf[j][i] = k_0_inf[j][i] * std::pow(10,-1*f_factors_inf[pos_FallOff_Reaction])/ ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
		}
	}

	// COMPUTE MAX/MIN reaction rates according to constant uncertainty factors for PLOG reactions //
	std::vector<std::vector<std::vector<double>>> k_classic_plog;
	//resize first dimensions according to the number of classic plog we are interested in for optimisation
    k_classic_plog.resize(ObjectInput2.list_of_uncertainty_factors_classic_plog.size());
	k_UB_classic_plog.resize(ObjectInput2.list_of_uncertainty_factors_classic_plog.size());
	k_LB_classic_plog.resize(ObjectInput2.list_of_uncertainty_factors_classic_plog.size());
	
	for (int j=0; j <ObjectInput2.list_of_uncertainty_factors_classic_plog.size(); j++)
	{
		//finding the position of the target reaction in correspondent OS++ object
		int pos_classic_plog_reaction = std::find(ObjectInput2.indices_of_classic_plogs.begin(),ObjectInput2.indices_of_classic_plogs.end(),ObjectInput2.list_of_target_classic_plog_reactions[j])-ObjectInput2.indices_of_classic_plogs.begin();

		k_classic_plog[j].resize(ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
		k_UB_classic_plog[j].resize(ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
		k_LB_classic_plog[j].resize(ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());

		for (int k=0; k < ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++)
		{
			double A_CP = std::exp(ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k]);
			double n_CP = ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).Beta()[k];
			double E_over_R_CP = ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[k];

			k_classic_plog[j][k].resize(T_span.size());
			k_UB_classic_plog[j][k].resize(T_span.size());
			k_LB_classic_plog[j][k].resize(T_span.size());

			for (int i=0; i < T_span.size(); i++)
			{
				k_classic_plog[j][k][i]    = A_CP * std::pow(T_span[i],n_CP) * std::exp((-1*E_over_R_CP)/T_span[i]);
				k_UB_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10,ObjectInput2.list_of_uncertainty_factors_classic_plog[j])* (double)ObjectInput2.AcceptedSigmaKDistribution/2;
				k_LB_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10,-1*ObjectInput2.list_of_uncertainty_factors_classic_plog[j])/ ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
			}
		}		
	}

	// COMPUTE MAX/MIN reaction rates according to constant uncertainty factors for PLOG reactions //
	std::vector<std::vector<std::vector<double>>> k_ext_plog;
	//resize first dimensions according to the number of classic plog we are interested in for optimisation
    k_ext_plog.resize(ObjectInput2.list_of_target_extplog.size());
	k_UB_ext_plog.resize(ObjectInput2.list_of_target_extplog.size());
	k_LB_ext_plog.resize(ObjectInput2.list_of_target_extplog.size());


	for (int j=0; j <ObjectInput2.list_of_target_extplog.size(); j++)
	{
		int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),ObjectInput2.list_of_target_extplog[j])-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();


		k_ext_plog[j].resize(ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size());
		k_UB_ext_plog[j].resize(ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size());
		k_LB_ext_plog[j].resize(ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size());

		for (int k=0; k < ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size(); k++)
		{
			double A_EP = std::exp(ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0)[k]);
			double n_EP = ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0)[k];
			double E_over_R_EP = ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0)[k];

			k_ext_plog[j][k].resize(T_span.size());
			k_UB_ext_plog[j][k].resize(T_span.size());
			k_LB_ext_plog[j][k].resize(T_span.size());

			for (int i=0; i < T_span.size(); i++)
			{
				k_ext_plog[j][k][i]    = A_EP * std::pow(T_span[i],n_EP) * std::exp((-1*E_over_R_EP)/T_span[i]);
				k_UB_ext_plog[j][k][i] = k_ext_plog[j][k][i] * std::pow(10,ObjectInput2.list_of_uncertainty_factors_extplog[j])* (double)ObjectInput2.AcceptedSigmaKDistribution/2;
				k_LB_ext_plog[j][k][i] = k_ext_plog[j][k][i] * std::pow(10,-1*ObjectInput2.list_of_uncertainty_factors_extplog[j])/ ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
			}
		}		
	}

	// COMPUTE MAX/MIN for EPLR  reactions //

	//std::cout<< "OptiSMOKE arrives at the constraints" <<std::endl;
	std::vector<std::vector<std::vector<double>>> k_EPLR;
    k_EPLR.resize(ObjectInput2.list_of_target_EPLR.size());
	k_UB_EPLR.resize(ObjectInput2.list_of_target_EPLR.size());
	k_LB_EPLR.resize(ObjectInput2.list_of_target_EPLR.size());

	//std::cout<< "OptiSMOKE:1" <<std::endl;

	for (int j=0; j <ObjectInput2.list_of_target_EPLR.size(); j++)
	{
		// Find position of the reaction indexed with "list_of_target_EPLR[j]"
		int pos_EPLR = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),ObjectInput2.list_of_target_EPLR[j])-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();
		// Find position of the specific bath gas within the reaction indexed with "list_of_target_EPLR[j]"
		int pos_BathGas_EPLR = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).get_bath_gas_position(ObjectInput2.list_of_bath_gases_EPLR[j]);


		//std::cout<< "pos_EPLR "			<< pos_EPLR <<std::endl;
		//std::cout<< "pos_BathGas_EPLR "	<< pos_BathGas_EPLR <<std::endl;


		//std::cout<< "OptiSMOKE:2" <<std::endl;

		k_EPLR[j].resize(ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size());
		k_UB_EPLR[j].resize(ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size());
		k_LB_EPLR[j].resize(ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size());

		//std::cout<< "OptiSMOKE:3" <<std::endl;
	
		for (int k=0; k < ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size(); k++)
		{

			//std::cout<< "OptiSMOKE:4" <<std::endl;	
			double A_EP = std::exp(ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR)[k]);
			double n_EP = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).Beta(pos_BathGas_EPLR)[k];
			double E_over_R_EP = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).E_over_R(pos_BathGas_EPLR)[k];

			k_EPLR[j][k].resize(T_span.size());
			k_UB_EPLR[j][k].resize(T_span.size());
			k_LB_EPLR[j][k].resize(T_span.size());

			//std::cout<< "OptiSMOKE:5" <<std::endl;
			for (int i=0; i < T_span.size(); i++)
			{
				k_EPLR[j][k][i]    = A_EP * std::pow(T_span[i],n_EP) * std::exp((-1*E_over_R_EP)/T_span[i]);
				k_UB_EPLR[j][k][i] = k_EPLR[j][k][i] * std::pow(10,ObjectInput2.list_of_uncertainty_factors_EPLR[j]) * (double)ObjectInput2.AcceptedSigmaKDistribution/2;
				k_LB_EPLR[j][k][i] = k_EPLR[j][k][i] * std::pow(10,-1*ObjectInput2.list_of_uncertainty_factors_EPLR[j]) / ((double)ObjectInput2.AcceptedSigmaKDistribution/2);
			}
			//std::cout<< "OptiSMOKE:6" <<std::endl;
		}

		//std::cout<< "OptiSMOKE:" <<std::endl;		
	}

	// ******************************************************** //
	//															//
	//				CREATE SPLINES FOR EXPERIMENTS				//
	//															//
	// ******************************************************** //
	if (ObjectInput2.Objective_Function == "CurveMatching")
	{

		
			// resize with the  number of datasets
			splinesExp.resize(ObjectInput2.Exp_data.size());

			if(ObjectInput2.UseBootStrap == false)
			{

				for (int a=0; a<ObjectInput2.Exp_data.size(); ++a)
				{	
					
					// resize with the number of target species
					splinesExp[a].resize(ObjectInput2.Exp_data[a].size());

					//std::cout<<"The dataset where i have the block is "<< a+1 <<std::endl;

					
					for (int b = 0; b<ObjectInput2.Exp_data[a].size(); ++b) {

						

						splinesExp[a][b].resize(1);

						if(ObjectInput2.QoI[a] == "IDT")
						{
							// Initialize 1 vector to get the ordinates
							std::vector<double> temporary_vector;			
							temporary_vector.resize(ObjectInput2.Exp_data[a][b][1].size());
							// 
							for (int z=0; z<ObjectInput2.Exp_data[a][b][1].size(); ++z)
							{
									temporary_vector[z] = std::log(ObjectInput2.Exp_data[a][b][1][z]);
							}
							// 						
							splinesExp[a][b][0].solve(ObjectInput2.Exp_data[a][b][0],temporary_vector,0,0,ObjectInput2.print_splines);

						} else {


							//
							//std::cout<<"Inside the dataset "<< a << " the problematic species is "<< b+1<<std::endl;

							splinesExp[a][b][0].solve(ObjectInput2.Exp_data[a][b][0],ObjectInput2.Exp_data[a][b][1],0,0,ObjectInput2.print_splines);
							
						}
						//std::cout<<"it is going to remove asymptotes"<<std::endl;
						splinesExp[a][b][0].removeAsymptotes();
						
					}

				}

			} else {

				
				BootStrapping_exp_data(ObjectInput2.Exp_data);
                		
				for (int a=0; a<ObjectInput2.Exp_data.size(); ++a) {
					

					splinesExp[a].resize(ObjectInput2.Exp_data[a].size());

					for (int c=0; c<ObjectInput2.Exp_data[a].size(); ++c) {

						splinesExp[a][c].resize(ObjectInput2.numberOfBootstrapVariations);

						for (int b=0; b< ObjectInput2.numberOfBootstrapVariations; ++b) {
							
							

							if(ObjectInput2.QoI[a] == "IDT") {

								std::vector<double> temporary_vector;
								temporary_vector.resize(bootstrapExp[a][c][b].size());

								for (int z=0; z<temporary_vector.size(); ++z)
								{
										temporary_vector[z] = std::log(bootstrapExp[a][c][b][z]);
								}
																
								splinesExp[a][c][b].solve(ObjectInput2.Exp_data[a][c][0],temporary_vector,0,0,ObjectInput2.print_splines);

							} else {
								
								splinesExp[a][c][b].solve(ObjectInput2.Exp_data[a][c][0],bootstrapExp[a][c][b],0,0,ObjectInput2.print_splines);
							}

							//std::cout<<"it is going to remove asymptotes"<<std::endl;
							splinesExp[a][c][b].removeAsymptotes();
						}
					}
                }
			}	
	}
	//std::cout<< "It has initialized the empty splines for data correctly"<<std::endl;
}



int OpenSMOKEDirectApplicInterface::derived_map_ac(const Dakota::String& ac_name)
{
#ifdef MPI_DEBUG
  Cout << "analysis server " << analysisServerId << " invoking " << ac_name
       << " within SIM::OpenSMOKEDirectApplicInterface." << std::endl;
#endif // MPI_DEBUG

  if (multiProcAnalysisFlag) 
  {
    Cerr << "Error: plugin serial direct fn does not support multiprocessor "
	 << "analyses." << std::endl;
    Dakota::abort_handler(-1);
  }
  int fail_code = 0;
  if (ac_name == "plugin_opensmoke") 
  {
      fail_code = opensmoke_interface(xC, directFnASV[0], fnVals[0]);//, cnVals);
  }
  else 
  {
    Cerr << ac_name << " is not available as an analysis within "
         << "SIM::OpenSMOKEDirectApplicInterface." << std::endl;
    Dakota::abort_handler(Dakota::INTERFACE_ERROR);
  }

  // Failure capturing
  if (fail_code) 
  {
    std::string err_msg("Error evaluating plugin analysis_driver ");
    err_msg += ac_name;
    throw Dakota::FunctionEvalFailure(err_msg);
  }

  return 0;
}

void OpenSMOKEDirectApplicInterface::
wait_local_evaluations(Dakota::PRPQueue& prp_queue)
{
  if (multiProcAnalysisFlag) 
  {
    Cerr << "Error: plugin serial direct fn does not support multiprocessor "
	 << "analyses." << std::endl;
    Dakota::abort_handler(-1);
  }
  for (Dakota::PRPQueueIter prp_iter = prp_queue.begin();
       prp_iter != prp_queue.end(); prp_iter++) 
  {
    // For each job in the processing queue, evaluate the response
    int fn_eval_id = prp_iter->eval_id();
    const Dakota::Variables& vars = prp_iter->variables();
    const Dakota::ActiveSet& set  = prp_iter->active_set();
    Dakota::Response         resp = prp_iter->response(); // shared rep
    if (outputLevel > Dakota::SILENT_OUTPUT)
      Cout << "OpenSMOKEDirectApplicInterface:: evaluating function evaluation "
	   << fn_eval_id << " in batch mode." << std::endl;
      short asv = set.request_vector()[0];
      Dakota::Real& fn_val = resp.function_value_view(0);
      opensmoke_interface(vars.continuous_variables(), asv, fn_val);

    // indicate completion of job to ApplicationInterface schedulers
    completionSet.insert(fn_eval_id);
  }
}

// Convert the c_vars (vector of suggested values) from DAKOTA to actual kinetic parameters value and update the kinetic map
int OpenSMOKEDirectApplicInterface::opensmoke_interface(const Dakota::RealVector& c_vars, short asv, Dakota::Real& fn_val)
{
	// Evaluation counter	
	eval_nr++;
	// Offset based on c_vars vector
	unsigned int count = 0;
	// Change kinetic parameters

	if (ObjectInput2.pca_direct_reactions.size() > 0)
	{
				std::cout<<"The PCA reconstruction will be performed on sampled variables"<<std::endl;

				std::vector<std::vector<double>> direct_matrix_samples;
				direct_matrix_samples.resize(ObjectInput2.pca_direct_reactions.size()); 

				int var_cnt = 0;
				int cnt = 0;

				for (int i=0; i<ObjectInput2.pca_minabs_direct_reactions.size(); i++){

					std::cout<<"Counter is "<< cnt << std::endl;
					direct_matrix_samples[cnt].push_back(c_vars[var_cnt]);

					double rmd = (i+1)%ObjectInput2.number_of_eigenvector_for_reaction[0];

					std::cout<<"Reminder is "<< rmd << std::endl;

					if (rmd ==0){
						cnt++;
					}

					var_cnt++;	
				}

				std::vector<std::vector<double>> pinf_matrix_samples;
				pinf_matrix_samples.resize(ObjectInput2.pca_pinf_reactions.size());

				cnt = 0;
				if (ObjectInput2.pca_pinf_reactions.size()!=0){

					for (int i=0; i<ObjectInput2.pca_minabs_pinf_reactions.size(); i++){

						std::cout<<"Counter is "<< cnt << std::endl;
						pinf_matrix_samples[cnt].push_back(c_vars[var_cnt]);

						double rmd = (i+1)%3;

						std::cout<<"Reminder is "<< rmd << std::endl;

						if (rmd ==0){
							cnt++;
						}	
					}

					var_cnt++;
				}

				std::vector<std::vector<double>> plog_matrix_samples;
				plog_matrix_samples.resize(ObjectInput2.pca_plog_reactions.size()); 

				cnt = 0;
				if (ObjectInput2.pca_plog_reactions.size()!=0){

					for (int i=0; i<ObjectInput2.pca_minabs_plog_reactions.size(); i++){

						std::cout<<"Counter is "<< cnt << std::endl;
						plog_matrix_samples[cnt].push_back(c_vars[var_cnt]);

						double rmd = (i+1)%3;

						std::cout<<"Reminder is "<< rmd << std::endl;

						if (rmd ==0){
							cnt++;
						}	
					}

					var_cnt++;
				}

				// for (int i=0; i<z_matrix_samples.size(); i++){
				// 		for (int j=0; j<z_matrix_samples[i].size(); j++){
				// 			std::cout<<"The value of sampled Z for reaction "<< i << " and variable "<< j << " is "<< z_matrix_samples[i][j]<<std::endl;
				// 		}
				// }

				//REPLACEMENT OF DIRECT PARAMETERS
				std::vector<double> temp_vector;
				std::vector<std::vector<double>> reconstructed_direct_parameters;
				reconstructed_direct_parameters.resize(ObjectInput2.pca_direct_reactions.size());
				std::cout<<"Here "<<std::endl;
				for (int i=0; i<ObjectInput2.pca_direct_reactions.size(); i++){
					
					temp_vector = reconstuction_pca(i,direct_matrix_samples[i], ObjectInput2.pca_eigenvectors_direct_reactions, ObjectInput2.pca_scaling_direct_reactions, ObjectInput2.pca_centering_direct_reactions);
					for (int j=0; j<temp_vector.size(); j++){
					
						reconstructed_direct_parameters[i].push_back(temp_vector[j]);
						std::cout<<"The value of reconstructed parameter "<< j << " in reaction "<< i << " is "<< temp_vector[j]<<std::endl;
					}					
				}
				
				for (int i=0; i<ObjectInput2.pca_direct_reactions.size(); i++){
					
					ChangelnA(ObjectInput2.pca_direct_reactions[i], reconstructed_direct_parameters[i][0]);
					ChangeBeta(ObjectInput2.pca_direct_reactions[i], reconstructed_direct_parameters[i][1]);
					ChangeE_over_R(ObjectInput2.pca_direct_reactions[i], reconstructed_direct_parameters[i][2]);
				}

				//REPLACEMENT OF P_INF PARAMETERS
				std::vector<std::vector<double>> reconstructed_pinf_parameters;
				reconstructed_pinf_parameters.resize(ObjectInput2.pca_pinf_reactions.size());
				std::cout<<"Here "<<std::endl;

				for (int i=0; i<ObjectInput2.pca_pinf_reactions.size(); i++){
					
					temp_vector = reconstuction_pca(i,pinf_matrix_samples[i], ObjectInput2.pca_eigenvectors_pinf_reactions, ObjectInput2.pca_scaling_pinf_reactions, ObjectInput2.pca_centering_pinf_reactions);

					for (int j=0; j<temp_vector.size(); j++){
					
						reconstructed_pinf_parameters[i].push_back(temp_vector[j]);
						std::cout<<"The value of reconstructed parameter "<< j << " in reaction "<< i << " is "<< temp_vector[j]<<std::endl;
					}					
				}
				
				for (int i=0; i<ObjectInput2.pca_pinf_reactions.size(); i++){
					
					ChangelnA_inf(ObjectInput2.pca_pinf_reactions[i], reconstructed_pinf_parameters[i][0]);
					ChangeBeta_inf(ObjectInput2.pca_pinf_reactions[i], reconstructed_pinf_parameters[i][1]);
					ChangeE_over_R_inf(ObjectInput2.pca_pinf_reactions[i], reconstructed_pinf_parameters[i][2]);
				}


				//REPLACEMENT OF PLOG PARAMETERS
				std::vector<std::vector<double>> reconstructed_plog_parameters;
				reconstructed_plog_parameters.resize(ObjectInput2.pca_plog_reactions.size());
				std::cout<<"Here "<<std::endl;
				for (int i=0; i<ObjectInput2.pca_plog_reactions.size(); i++){
					
					temp_vector = reconstuction_pca(i,plog_matrix_samples[i], ObjectInput2.pca_eigenvectors_plog_reactions, ObjectInput2.pca_scaling_plog_reactions, ObjectInput2.pca_centering_plog_reactions);

					for (int j=0; j<temp_vector.size(); j++){
					
						reconstructed_plog_parameters[i].push_back(temp_vector[j]);
						std::cout<<"The value of reconstructed parameter "<< j << " in reaction "<< i << " is "<< temp_vector[j]<<std::endl;
					}					
				}
				
				for (int i=0; i<ObjectInput2.pca_plog_reactions.size(); i++){
					
					ChangelnA_classic_PLOG(ObjectInput2.pca_plog_reactions[i], std::exp(reconstructed_plog_parameters[i][0]));
					Change_Beta_classic_PLOG(ObjectInput2.pca_plog_reactions[i], reconstructed_plog_parameters[i][1]);
					Change_ER_classic_PLOG(ObjectInput2.pca_plog_reactions[i], reconstructed_plog_parameters[i][2]);
				}
			
	} 
	else if (ObjectInput2.direct_reactions_indices_ica.size() > 0){

				std::cout<<"The ICA mixing process is going to be performed"<<std::endl;

				std::vector<std::vector<double>> direct_matrix_samples;
				direct_matrix_samples.resize(ObjectInput2.direct_reactions_indices_ica.size()); 

				int var_cnt = 0;
				int cnt = 0;
				// organize the c_vars vector from Dakota
				for (int i=0; i<ObjectInput2.abs_min_direct_reactions_ica.size(); i++){
					
					direct_matrix_samples[cnt].push_back(c_vars[var_cnt]);

					double rmd = (i+1)%3;

					if (rmd ==0){
						cnt++;
					}

					var_cnt++;	
				}

				//REPLACEMENT OF DIRECT PARAMETERS
				std::vector<double> temp_vector;
				std::vector<std::vector<double>> mixed_direct_parameters;
				mixed_direct_parameters.resize(ObjectInput2.direct_reactions_indices_ica.size());

				for (int i=0; i<ObjectInput2.direct_reactions_indices_ica.size(); i++){
					
					temp_vector = ica_mixing(i,direct_matrix_samples[i], ObjectInput2.mixing_matrices_ica[i]);

					for (int j=0; j<temp_vector.size(); j++){
					
						mixed_direct_parameters[i].push_back(temp_vector[j]);
						std::cout<<"The value of reconstructed parameter "<< j << " in reaction "<< i << " is "<< temp_vector[j]<<std::endl;
					}					
				}
				
				for (int i=0; i<ObjectInput2.direct_reactions_indices_ica.size(); i++){
					
					ChangelnA(ObjectInput2.direct_reactions_indices_ica[i], mixed_direct_parameters[i][0]);
					ChangeBeta(ObjectInput2.direct_reactions_indices_ica[i], mixed_direct_parameters[i][1]);
					ChangeE_over_R(ObjectInput2.direct_reactions_indices_ica[i], mixed_direct_parameters[i][2]);
				}
			
	} 
	else if (ObjectInput2.boundaries_method == "Re-parametrization") {
				std::cout<<"YOU ARE WITHING RE-PARAMETRIZATION"<<std::endl;
				// STORE THE variables from c_vars first
				std::vector<double> store_A_scaled(ObjectInput2.list_of_target_lnA.size());
				std::vector<double> store_Beta(ObjectInput2.list_of_target_Beta.size());
				std::vector<double> store_E_over_R(ObjectInput2.list_of_target_E_over_R.size());
				// lnA
				if (ObjectInput2.list_of_target_lnA.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_lnA.size(); k++)
					{
						store_A_scaled[k] = c_vars[count];
						count=count+1;
					}
				}
				// Beta
				if (ObjectInput2.list_of_target_Beta.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_Beta.size(); k++)
					{
						store_Beta[k] = c_vars[count];
						count=count+1;
					}
				}
				// E_over_R
				if (ObjectInput2.list_of_target_E_over_R.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_E_over_R.size(); k++)
					{
						store_E_over_R[k] = c_vars[count];
						count=count+1;
						
					}
				}
				// TRANSLATE A_scaled into A and change all the variables into the mechanism
				if (ObjectInput2.list_of_target_lnA.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_lnA.size(); k++)
					{
						double A_from_rescale = std::log(store_A_scaled[k]/(std::pow(1000,store_Beta[k])*std::exp(-1*store_E_over_R[k]/1000)));
						std::cout<< "The value of scaled A is 					"<<store_A_scaled[k]<<std::endl;
						std::cout<< "The value of Beta is     					"<<store_Beta[k]<<std::endl;
						std::cout<< "The value of E_over_R is 					"<<store_E_over_R[k]<<std::endl;
						std::cout<< "The value of A after being re-scaled is 	"<<A_from_rescale<<std::endl;
						std::cout<< "	"<<std::endl;

						ChangelnA(ObjectInput2.list_of_target_lnA[k], A_from_rescale);
					}
				}
				// Beta
				if (ObjectInput2.list_of_target_Beta.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_Beta.size(); k++)
					{
						ChangeBeta(ObjectInput2.list_of_target_Beta[k], store_Beta[k]);
					}
				}
				// E_over_R
				if (ObjectInput2.list_of_target_E_over_R.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_E_over_R.size(); k++)
					{
						ChangeE_over_R(ObjectInput2.list_of_target_E_over_R[k], store_E_over_R[k]);	
					}
				}
	} 
	else {
		// lnA
		if (ObjectInput2.list_of_target_lnA.size()!=0)
		{
			if(ObjectInput2.Optimization4Classes == true){
				for(unsigned int k = 0; k < ObjectInput2.numberOfReactionClasses; k++){
					for(unsigned int z = 0; z < ObjectInput2.matrixOflnA[z].size(); z++){
						int target_reaction = ObjectInput2.matrixOflnA[k][z];
						if(eval_nr == 1){
							ChangelnA(target_reaction, std::log(ObjectInput2.kineticsMapXML -> A(target_reaction - 1)));
						}
						else{
							if(z == 0){
								ChangelnA(target_reaction, c_vars[count]);
							}
							else{
								double scaling_lnA = ObjectInput2.matrixOfscalinglnA[k][z];
								ChangelnA(target_reaction, std::log(scaling_lnA * ObjectInput2.nominalkineticsMapXML->A(ObjectInput2.matrixOfReactionIndex[k][z]-1))); // check
							}
						}
						count = count + 1;
					}
				}
			}
			else{
				for (unsigned int k = 0; k < ObjectInput2.list_of_target_lnA.size(); k++){
					ChangelnA(ObjectInput2.list_of_target_lnA[k], c_vars[count]);
					count=count+1;
				}
			}
		}
		// lnA_inf
		if (ObjectInput2.list_of_target_lnA_inf.size()!=0){
			for (unsigned int k = 0; k < ObjectInput2.list_of_target_lnA_inf.size(); k++){
				ChangelnA_inf(ObjectInput2.list_of_target_lnA_inf[k], c_vars[count]);
				count=count+1;
			}
		}
		// Beta
		if (ObjectInput2.list_of_target_Beta.size()!=0){
			if(ObjectInput2.Optimization4Classes == true){
				for(unsigned int k = 0; k < ObjectInput2.numberOfReactionClasses; k++){
					for(unsigned int z = 0; z < ObjectInput2.matrixOfBeta[z].size(); z++){
						int target_reaction = ObjectInput2.matrixOfBeta[k][z];
						if(eval_nr == 1){
							ChangeBeta(target_reaction, ObjectInput2.kineticsMapXML -> Beta(target_reaction - 1));
						}
						else{
							if(z == 0){
								ChangeBeta(target_reaction, c_vars[count]);
							}
							else{
								double scaling_Beta = ObjectInput2.matrixOfscalingBeta[k][z];
								ChangeBeta(target_reaction, scaling_Beta + ObjectInput2.nominalkineticsMapXML-> Beta(ObjectInput2.matrixOfReactionIndex[k][z]-1)); // check
							}
						}
						count = count + 1;
					}
				}
			}
			else{
				for (unsigned int k = 0; k < ObjectInput2.list_of_target_Beta.size(); k++){
					ChangeBeta(ObjectInput2.list_of_target_Beta[k], c_vars[count]);
					count=count+1;
				}
			}
		}
		// Beta_inf
		if (ObjectInput2.list_of_target_Beta_inf.size()!=0){
			for (unsigned int k = 0; k < ObjectInput2.list_of_target_Beta_inf.size(); k++){
				ChangeBeta_inf(ObjectInput2.list_of_target_Beta_inf[k], c_vars[count]);
				count=count+1;
			}
		}
		// E_over_R
		if (ObjectInput2.list_of_target_E_over_R.size()!=0){
			if(ObjectInput2.Optimization4Classes == true){
				for(unsigned int k = 0; k < ObjectInput2.numberOfReactionClasses; k++){
					for(unsigned int z = 0; z < ObjectInput2.matrixOfEoverR[z].size(); z++){
						int target_reaction = ObjectInput2.matrixOfEoverR[k][z];
						if(eval_nr == 1){
							ChangeE_over_R(target_reaction, ObjectInput2.kineticsMapXML -> E_over_R(target_reaction - 1));
						}
						else{
							if(z == 0){
								ChangeE_over_R(target_reaction, c_vars[count]);
							}
							else{
								double scaling_E_over_R = ObjectInput2.matrixOfscalingEoverR[k][z];
								ChangeE_over_R(target_reaction, scaling_E_over_R + ObjectInput2.nominalkineticsMapXML-> E_over_R(ObjectInput2.matrixOfReactionIndex[k][z]-1)); // check
							}
						}
						count = count + 1;
					}
				}
			}
			else{
				for (unsigned int k = 0; k < ObjectInput2.list_of_target_E_over_R.size(); k++){
					ChangeE_over_R(ObjectInput2.list_of_target_E_over_R[k], c_vars[count]);
					count=count+1;
				}
			}
		}
				// E_over_R_inf
				if (ObjectInput2.list_of_target_E_over_R_inf.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_E_over_R_inf.size(); k++)
					{
						ChangeE_over_R_inf(ObjectInput2.list_of_target_E_over_R_inf[k], c_vars[count]);
						count=count+1;
					}
				}
				// Third body efficiencies
				if (ObjectInput2.list_of_target_thirdbody_reactions.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_thirdbody_reactions.size(); k++)
					{
						ChangeThirdBody_Eff(ObjectInput2.list_of_target_thirdbody_reactions[k], ObjectInput2.list_of_target_thirdbody_species[k], c_vars[count]);
						count=count+1;
					}
				}

				// CLASSIC PLOG
				if (ObjectInput2.list_of_target_classic_plog_reactions.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_classic_plog_reactions.size(); k++)
					{	
							//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
							//std::cout<<"The new value is " << std::to_string(c_vars[count]) <<std::endl;
							ChangelnA_classic_PLOG(ObjectInput2.list_of_target_classic_plog_reactions[k],c_vars[count]);
							count=count+1;
							//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}

				if (ObjectInput2.list_of_target_classic_plog_reactions.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_classic_plog_reactions.size(); k++)
					{	
							//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
							//std::cout<<"The new value is " << std::to_string(c_vars[count]) <<std::endl;
							Change_ER_classic_PLOG(ObjectInput2.list_of_target_classic_plog_reactions[k],c_vars[count]);
							count=count+1;
							//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}
				
				if (ObjectInput2.list_of_target_classic_plog_reactions.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_classic_plog_reactions.size(); k++)
					{
						Change_Beta_classic_PLOG(ObjectInput2.list_of_target_classic_plog_reactions[k],c_vars[count]);
						count=count+1;
					}
				}

				// EPLR
				if (ObjectInput2.list_of_target_EPLR.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_EPLR.size(); k++)
					{	
						
						ChangelnA_EPLR(ObjectInput2.list_of_target_EPLR[k], c_vars[count], ObjectInput2.list_of_bath_gases_EPLR[k]);
						count=count+1;
					}
				}

				if (ObjectInput2.list_of_target_EPLR.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_EPLR.size(); k++)
					{	
						Change_ER_EPLR(ObjectInput2.list_of_target_EPLR[k], c_vars[count], ObjectInput2.list_of_bath_gases_EPLR[k]);
						count=count+1;
					}
				}
				
				if (ObjectInput2.list_of_target_EPLR.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_EPLR.size(); k++)
					{
						//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
						Change_Beta_EPLR(ObjectInput2.list_of_target_EPLR[k], c_vars[count], ObjectInput2.list_of_bath_gases_EPLR[k]);
						count=count+1;
						//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}

				// EXTENDED PLOG
				if (ObjectInput2.list_of_target_extplog.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_extplog.size(); k++)
					{	
							//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
							//std::cout<<"The new value is " << std::to_string(c_vars[count]) <<std::endl;
							ChangelnA_ExtPLOG(ObjectInput2.list_of_target_extplog[k], c_vars[count], ObjectInput2.list_of_target_extended_plog_reactions, ObjectInput2.list_of_target_extended_plog_species);
							count=count+1;
							//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}

				if (ObjectInput2.list_of_target_extplog.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_extplog.size(); k++)
					{	
							//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
							//std::cout<<"The new value is " << std::to_string(c_vars[count]) <<std::endl;
							Change_ER_ExtPLOG(ObjectInput2.list_of_target_extplog[k],c_vars[count], ObjectInput2.list_of_target_extended_plog_reactions, ObjectInput2.list_of_target_extended_plog_species);
							count=count+1;
							//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}
				
				if (ObjectInput2.list_of_target_extplog.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_extplog.size(); k++)
					{
						Change_Beta_ExtPLOG(ObjectInput2.list_of_target_extplog[k],c_vars[count], ObjectInput2.list_of_target_extended_plog_reactions, ObjectInput2.list_of_target_extended_plog_species);
						count=count+1;
					}
				}			

				// COME ON MAN // THIS IS IT
				if (ObjectInput2.list_of_target_extended_plog_reactions.size()!=0)
				{
					for (unsigned int k = 0; k < ObjectInput2.list_of_target_extended_plog_reactions.size(); k++)
					{	
							//std::cout<<"The Value of count before is " << std::to_string(count) <<std::endl;
							//std::cout<<"The new value is " << std::to_string(c_vars[count]) <<std::endl;
							Change_ExtPLOG_TB(ObjectInput2.list_of_target_extended_plog_reactions[k], ObjectInput2.list_of_target_extended_plog_species[k], c_vars[count]);
							count=count+1;
							//std::cout<<"The Value of count after is " << std::to_string(count) <<std::endl;
					}
				}
	}
	//std::cout<<"All the parameters have been replaced" <<std::endl;

	// Print out if penalty function is turned off
	if(!ObjectInput2.penalty_function)
	{
		std::cout<<"Penalty function is turned off!"<<std::endl;
	}

	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// 						EXECUTE SIMULATIONS 							//
	// 		Here, all simulations are executed if and only if no penalty	//
	// 		are detected. Then, based on different criteria, simulations	//
	// 		are performed on flames, batch, plug flow and perfectly			//
	// 		stirred reactors												//
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	violated_uncertainty = false;
	if (ObjectInput2.penalty_function && ObjectInput2.list_of_target_uncertainty_factors.size() !=0)
	{
		//fn_val == obj_function
		violated_uncertainty = Check_k(fn_val);
	} 

	if (violated_uncertainty == false)
	{
		// initialise a 2D vector for Sim_values
		std::vector<std::vector<std::vector<double>>> Sim_values;
		Sim_values.resize(ObjectInput2.list_of_opensmoke_input_files.size());

		//std::cout<<"The list of OS input file is "<< ObjectInput2.list_of_opensmoke_input_files.size() <<std::endl;

		// loop over the number of datasets in the considered database
		for( int i = 0; i < ObjectInput2.list_of_opensmoke_input_files.size(); i++)
		{
			// SOLVING BatchReactors
			Sim_values[i].resize(ObjectInput2.what_2_calc[i].size());

			if (ObjectInput2.type_of_reactor[i] == "Batch" || ObjectInput2.type_of_reactor[i] == "RCM") {
								

				if (ObjectInput2.QoI[i] == "IDT") {	

					Sim_values[i][0].resize(ObjectInput2.Exp_data[i][0][1].size());
					
					std::cout<<"Solving Batch with "<< ObjectInput2.what_2_calc[i][0] <<" IDT definition for dataset "<<i+1<<std::endl;

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						batch_reactors[i][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						batch_reactors[i][m].Update_and_Solve_Batch(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values[i][0][m] = batch_reactors[i][m].Solve_tau(ObjectInput2.what_2_calc[i][0])*std::pow(10,6);
						
						if (ObjectInput2.type_of_reactor[i] == "RCM"){
							Sim_values[i][0][m] = Sim_values[i][0][m] - ObjectInput2.Exp_data[i][0][3][m];
						} 
						
						if(ObjectInput2.Debug_Sim) {	
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;		
						}
					}							
				}

				if (ObjectInput2.QoI[i] == "Species"){

					Sim_values[i][0].resize(ObjectInput2.Exp_data[i][0][1].size());

					std::cout<<"Solving Batch with "<< ObjectInput2.what_2_calc[i][0] <<" as target species for dataset "<<i+1<<std::endl;

					batch_reactors[i][0].Setup(ObjectInput2.list_of_opensmoke_input_files[i][0], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
					batch_reactors[i][0].Update_and_Solve_Batch(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
					Sim_values[i][0] = batch_reactors[i][0].Solve_Species(ObjectInput2.Exp_data[i][0][0], ObjectInput2.what_2_calc[i][0]);


					if(ObjectInput2.Debug_Sim)
					{
						for (int m = 0; m < Sim_values[i][0].size(); m++)
						{
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;	
						}
					}			
				}

				if (ObjectInput2.QoI[i] == "m_SP") {

					std::cout<<"Solving Batch with "<< ObjectInput2.what_2_calc[i].size() <<" species for dataset "<<i+1<<std::endl;
					batch_reactors[i][0].Setup(ObjectInput2.list_of_opensmoke_input_files[i][0], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
					batch_reactors[i][0].Update_and_Solve_Batch(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);


					//////////////////
					Sim_values[i] = batch_reactors[i][0].Solve_Multiple_Species(ObjectInput2.Exp_data[i],ObjectInput2.what_2_calc[i]);
					//////////////////

					if(ObjectInput2.Debug_Sim){

						for(int w=0;w<Sim_values[i].size();w++){

							for(int v=0;v<Sim_values[i][w].size();v++){

								std::cout<<"x["<< i << "][" << w <<  "][" << v << "] = " << Sim_values[i][w][v]  << std::endl;		
							}
						}
					}	
				}
			}
			
			// SOLVING PlugFlowReactors
			if (ObjectInput2.type_of_reactor[i] == "PFR") {
				
				Sim_values[i].resize(ObjectInput2.what_2_calc[i].size());
				//Sim_values[i].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());

				if (ObjectInput2.QoI[i] == "IDT")
				{

					Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());
					std::cout<<"Solving PFR with "<< ObjectInput2.what_2_calc[i][0] <<" IDT definition for dataset "<<i+1<<std::endl;

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						plugflow_reactors[i-batch_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						plugflow_reactors[i-batch_reactors.size()][m].Update_and_Solve_PFR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values[i][0][m] = plugflow_reactors[i-batch_reactors.size()][m].Solve_tau(ObjectInput2.what_2_calc[i][0])*std::pow(10,6);
						
						if(ObjectInput2.Debug_Sim)
						{	
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;	
						}
					}
				}

				if (ObjectInput2.QoI[i] == "Species")
				{
					Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());
					std::cout<<"Solving PFR with "<< ObjectInput2.what_2_calc[i][0] <<" as target species for dataset "<<i+1<<std::endl;

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						plugflow_reactors[i-batch_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						plugflow_reactors[i-batch_reactors.size()][m].Update_and_Solve_PFR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values[i][0][m] = plugflow_reactors[i-batch_reactors.size()][m].Solve_Species(ObjectInput2.what_2_calc[i][0]);

						if(ObjectInput2.Debug_Sim)
						{
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;	
						}
					}		
				} 

				if (ObjectInput2.QoI[i] == "m_SP_time")
				{
					std::cout<<"Solving PFR with "<< ObjectInput2.what_2_calc[i].size()<<" species for dataset "<<i+1<<std::endl;
					//for (int m = 0; m < ObjectInput_.ListOfOpensmokeInputFiles()[i].size(); m++)
					//{
						plugflow_reactors[i-batch_reactors.size()][0].Setup(ObjectInput2.list_of_opensmoke_input_files[i][0], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						plugflow_reactors[i-batch_reactors.size()][0].Update_and_Solve_PFR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values[i] = plugflow_reactors[i-batch_reactors.size()][0].Solve_Multipl_Species_time_profile(ObjectInput2.Exp_data[i],ObjectInput2.what_2_calc[i]);

						if(ObjectInput2.Debug_Sim)
						{
							for(int w=0;w<Sim_values[i].size();w++)
							{
								for(int v=0;v<Sim_values[i][w].size();v++)
								{
									std::cout<<"x["<< i << "][" << w <<  "][" << v << "] = " << Sim_values[i][w][v]  << std::endl;		
								}
							}
						}
					//}		

				}
				
				if (ObjectInput2.QoI[i] == "m_SP_out")
				{
					std::cout<<"Solving PFR with "<< ObjectInput2.what_2_calc[i].size()<<" species for dataset "<<i+1<<std::endl;

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						// perform simulation for the single PFR!!!
						plugflow_reactors[i-batch_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						plugflow_reactors[i-batch_reactors.size()][m].Update_and_Solve_PFR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						// extract vector of compositions for species of interest!!!
						std::vector<double> temp_res = plugflow_reactors[i-batch_reactors.size()][m].Solve_Multipl_Species_outlet(ObjectInput2.what_2_calc[i]);

						// arrange the results in the variables of the simulations results!!!
						for (int t = 0; t < ObjectInput2.what_2_calc[i].size(); t++){
								Sim_values[i][t].push_back(temp_res[t]);
						}
						// print debugging if requested
						if(ObjectInput2.Debug_Sim)
						{
							for(int w=0;w<Sim_values[i].size();w++)
							{
								for(int v=0;v<Sim_values[i][w].size();v++)
								{
									std::cout<<"x["<< i << "][" << w <<  "][" << v << "] = " << Sim_values[i][w][v]  << std::endl;		
								}
							}
						}
					}		

				}

				if (ObjectInput2.QoI[i] == "X_in_T")
				{
					Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());
					std::cout << "Solving PFR with conversion of: " << 
						ObjectInput2.what_2_calc[i][0] << " for dataset " << i+1 << std::endl;
					for(int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						plugflow_reactors[i-batch_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], 
																			ObjectInput2.thermodynamicsMapXML, 
																			ObjectInput2.kineticsMapXML);

						plugflow_reactors[i-batch_reactors.size()][m].Update_and_Solve_PFR(ObjectInput2.thermodynamicsMapXML, 
																						ObjectInput2.kineticsMapXML);

						Sim_values[i][0][m] = plugflow_reactors[i-batch_reactors.size()][m].Solve_Outlet_Conversion(ObjectInput2.what_2_calc[i][0]);
						
						if(ObjectInput2.Debug_Sim)
						{
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;		
						}
					}
				}
			}
				// SOLVING PerfectlyStirredReactors
			if (ObjectInput2.type_of_reactor[i] == "PSR") {

				Sim_values[i].resize(ObjectInput2.what_2_calc[i].size());	
				std::cout<<"Solving PSR with "<< ObjectInput2.what_2_calc[i][0] <<" as target species for dataset "<<i+1<<std::endl;

				if (ObjectInput2.QoI[i] == "Species") {
					
					Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++) {
						psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Update_and_Solve_PSR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values[i][0][m]=psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Solve_Species(ObjectInput2.what_2_calc[i][0]);

						if(ObjectInput2.Debug_Sim) {
								
							std::cout<<"x["<< i << "][0][" << m << "] = " << Sim_values[i][0][m]  << std::endl;
						}
					}	
				} 

				if (ObjectInput2.QoI[i] == "m_SP")
				{
					std::cout<<"Solving PSR with "<< ObjectInput2.what_2_calc[i].size()<<" species for dataset "<<i+1<<std::endl;
					Sim_values[i].resize(ObjectInput2.what_2_calc[i].size());

					for (int z = 0; z < ObjectInput2.what_2_calc[i].size(); z++)
					{
						Sim_values[i][z].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());
					}

					std::vector<double> Sim_values_temp;
					Sim_values_temp.resize(ObjectInput2.what_2_calc[i].size());

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Setup(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Update_and_Solve_PSR(ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML);
						Sim_values_temp = psr_reactors[i-batch_reactors.size()-plugflow_reactors.size()][m].Solve_Multipl_Species(ObjectInput2.what_2_calc[i]);

						for (int z = 0; z < ObjectInput2.what_2_calc[i].size(); z++)
						{
							if (ObjectInput2.what_2_calc[i][z]=="DeltaTemp"){
								Sim_values[i][z][m] = Sim_values_temp[z]-ObjectInput2.Exp_data[i][z][0][m];
							} else {
								Sim_values[i][z][m] = Sim_values_temp[z];
							}
						}
					}
					if (ObjectInput2.Debug_Sim)
												{		
						for (int z = 0; z < ObjectInput2.what_2_calc[i].size(); z++)
						{
							for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
														{
								std::cout<<"x["<< i << "][" << z <<  "][" << m << "] = " << Sim_values[i][z][m]  << std::endl;	
							}
						}
					}
				}
			}

			if (ObjectInput2.type_of_reactor[i] == "LaminarFlame")
			{	
				Sim_values[i].resize(ObjectInput2.what_2_calc[i].size());
				std::cout<<"Solving Laminar Flame with "<< ObjectInput2.what_2_calc[i][0] <<" for dataset "<<i+1<<std::endl;

				if (ObjectInput2.QoI[i] == "LFS")
				{
					Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());

					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						//std::shared_ptr<OpenSMOKE::Flame1D_Plugin> ciao_dario = std::make_shared<OpenSMOKE::Flame1D_Plugin>();
						OpenSMOKE::Flame1D_Plugin* laminar_flames = new OpenSMOKE::Flame1D_Plugin;
						std::vector<double> empty_vector;
						
						Sim_values[i][0][m]=laminar_flames->Setup_and_Solve(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML, ObjectInput2.transportMapXML, "LFS", empty_vector)[0];
						
						if(ObjectInput2.Debug_Sim)
						{
							std::cout<<"LFS[i][m] = " << Sim_values[i][0][m] << " i = "<< i+1 <<" m = "<<m<<std::endl;	
						}
						//delete laminar_flames;
					}	
				}

				if (ObjectInput2.QoI[i] == "BSF")
				{
					//Sim_values[i][0].resize(ObjectInput2.list_of_opensmoke_input_files[i].size());
					for (int m = 0; m < ObjectInput2.list_of_opensmoke_input_files[i].size(); m++)
					{
						//std::shared_ptr<OpenSMOKE::Flame1D_Plugin> ciao_dario = std::make_shared<OpenSMOKE::Flame1D_Plugin>();
						OpenSMOKE::Flame1D_Plugin* laminar_flames = new OpenSMOKE::Flame1D_Plugin;

						Sim_values[i][0]=laminar_flames->Setup_and_Solve(ObjectInput2.list_of_opensmoke_input_files[i][m], ObjectInput2.thermodynamicsMapXML, ObjectInput2.kineticsMapXML, ObjectInput2.transportMapXML, ObjectInput2.what_2_calc[i][0], ObjectInput2.Exp_data[i][0][0]);
						
						if(ObjectInput2.Debug_Sim)
						{
							for (int k = 0; k < Sim_values[i][0].size(); k++){
								std::cout<<"Molar fraction of "<<ObjectInput2.what_2_calc[i][0]<<"["<< i <<"]["<< m <<"]["<< k <<"] = "<< Sim_values[i][m][k] <<std::endl;	
							}
						}
						//delete laminar_flames;
					}	
				}
			}

			if (ObjectInput2.type_of_reactor[i] == "KTP") {	


			Sim_values[i][0].resize(ObjectInput2.Exp_data[i][0][0].size());

				if (ObjectInput2.preprocessor_kinetics->reactions()[std::stoi(ObjectInput2.QoI[i])-1].IsExtendedPressureLog() == true)
				{

					int pos_EPLR = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),std::stoi(ObjectInput2.QoI[i]))-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();

					

					if (ObjectInput2.type_KTP[i] == "ConstP"){

						for (int m = 0; m<ObjectInput2.Exp_data[i][0][0].size(); m++){

							Sim_values[i][0][m] = ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).KineticConstant_optismoke(ObjectInput2.what_2_calc[i][0], ObjectInput2.Exp_data[i][0][0][m], ObjectInput2.value_KTP[i]*101325);

							if(ObjectInput2.Debug_Sim)
							{
								std::cout<<"KTP[i][m] = " << Sim_values[i][0][m] << " i = "<< i <<" m = "<<m<<std::endl;	
							}
						
						}
						
					} else if (ObjectInput2.type_KTP[i] == "ConstT"){

						//std::cout<< "The size of Exp_data[i][0][0] is "<< ObjectInput2.Exp_data[i][0][0].size() <<std::endl;

						for (int m = 0; m<ObjectInput2.Exp_data[i][0][0].size(); m++){
							
							Sim_values[i][0][m] = ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).KineticConstant_optismoke(ObjectInput2.what_2_calc[i][0], ObjectInput2.value_KTP[i], ObjectInput2.Exp_data[i][0][0][m]*101325);

							if(ObjectInput2.Debug_Sim)
							{
								std::cout<<"KTP[i][m] = " << Sim_values[i][0][m] << " i = "<< i <<" m = "<<m<<std::endl;	
							}
						}
					}

				} else {
					
					for (int m = 0; m<Sim_values[i][0].size(); m++){

						//compute the log of k from 
						Sim_values[i][0][m] = std::log(ObjectInput2.kineticsMapXML->A(std::stoi(ObjectInput2.QoI[i])-1) / ObjectInput2.preprocessor_kinetics->reactions()[std::stoi(ObjectInput2.QoI[i])-1].A_conversion()) +  ObjectInput2.kineticsMapXML->Beta(std::stoi(ObjectInput2.QoI[i])-1) * std::log(ObjectInput2.Exp_data[i][0][0][m]) - ObjectInput2.kineticsMapXML->E_over_R(std::stoi(ObjectInput2.QoI[i])-1) / ObjectInput2.Exp_data[i][0][0][m];
					
					if(ObjectInput2.Debug_Sim)
					{
						std::cout<<"KTP[i][m] = " << Sim_values[i][0][m] << " i = "<< i <<" m = "<<m<<std::endl;	
					}
					
					}

				}
				
			}		
		}

	  	fn_val = Calculate_ObjFunc(Sim_values);
		Sim_values.clear();
	}

	if (eval_nr==1)
	{
		prev_fn_val=fn_val;
	} else if (prev_fn_val>fn_val)
	{
		prev_fn_val=fn_val;
		// Overwrite best mechanism CKI file
  		if (ObjectInput2.CKI_File_Read)
  		{
  			/*std::vector<double> best_param_vec;
		  	for(int i = 0; i < c_vars.length(); i++)
  			{
				best_param_vec.push_back(c_vars[i]);
  			}*/
			//double tStart_CKI = OpenSMOKE::OpenSMOKEGetCpuTime();
  			ObjectInput2.PrintFinalMechanism();
			//double tEnd_CKI = OpenSMOKE::OpenSMOKEGetCpuTime();
			//std::cout << "Time to write CKI file: " << tEnd_CKI - tStart_CKI << std::endl;
		}
	}
  	return 0;
}
// Function for controlling that k and k_inf values for new set of parameters are within the uncertainty limits. 
// If not, the penalty function is activated and the objective fucntion value is set to 1e7.
bool OpenSMOKEDirectApplicInterface::Check_k(Dakota::Real& fn_val)
{
		std::vector<double> T_span = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
		for (int j=0; j < ObjectInput2.list_of_target_uncertainty_factors.size(); j++)
		{
			//std::cout<<"Checking the constrained direct reaction number " << ObjectInput2.list_of_target_uncertainty_factors[j]<<std::endl;

			std::vector<double> k_check;
			k_check.resize(T_span.size());
			for (int i=0; i < T_span.size(); i++)
			{
				// Calculating the k value for the new set of kinetic parameters for a temperature span of 300-3000 K
				k_check[i] = ObjectInput2.kineticsMapXML->A(ObjectInput2.list_of_target_uncertainty_factors[j]-1) * std::pow(T_span[i],ObjectInput2.kineticsMapXML->Beta(ObjectInput2.list_of_target_uncertainty_factors[j]-1)) * std::exp((-1*ObjectInput2.kineticsMapXML->E_over_R(ObjectInput2.list_of_target_uncertainty_factors[j]-1))/T_span[i]);
				// If at one temperature the k value is either below the lower bound or above the upper bound, 
				// forcefully set the objective function value to 1e7 and print out for which reaction the violation occured 
				if ((k_check[i]<=k_lower[j][i]) || (k_check[i]>=k_upper[j][i]))
				{
					if (ObjectInput2.Objective_Function == "CurveMatching")
					{
						fn_val = 1;
					}else {
						fn_val = 10000000;
					}
					std::cout<<"Violation for reaction "<<ObjectInput2.list_of_target_uncertainty_factors[j]<<std::endl;
					return true;
				}
			}

			//std::cout<<"Checked the constrained direct reaction number " << ObjectInput2.list_of_target_uncertainty_factors[j]<<std::endl;
		}

		for (int j=0; j <ObjectInput2.list_of_target_uncertainty_factors_inf.size(); j++)
		{
			//std::cout<<"Checking the constrained P inf reaction number " << ObjectInput2.list_of_target_uncertainty_factors_inf[j]<<std::endl;
			std::vector<double> k_check_inf;
			k_check_inf.resize(T_span.size());
			// Finding position of the fall off reaction
			int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),ObjectInput2.list_of_target_uncertainty_factors_inf[j])-ObjectInput2.indices_of_falloff_reactions.begin();
			for (int i=0; i < T_span.size(); i++)
			{
				// Calculating the k_inf value for the new set of kinetic parameters for a temperature span of 300-3000 K
				k_check_inf[i] = ObjectInput2.kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction) * std::pow(T_span[i],ObjectInput2.kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction)) * std::exp((-1*ObjectInput2.kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction))/T_span[i]);
				// If at one temperature the k_inf value is either below the lower bound or above the upper bound, 
				// forcefully set the objective function value to 1e7 and print out for which reaction the violation occured 
				if ((k_check_inf[i]<=k_lower_inf[j][i]) || (k_check_inf[i]>=k_upper_inf[j][i]))
				{
					if (ObjectInput2.Objective_Function == "CurveMatching")
					{
							fn_val = 1;
					}else {
							fn_val = 10000000;
					}
					std::cout<<"Violation for reaction "<<ObjectInput2.list_of_target_uncertainty_factors_inf[j]<<" (inf) " <<std::endl;
					return true;
				}
			}

			//std::cout<<"Checked the constrained P inf reaction number " << ObjectInput2.list_of_target_uncertainty_factors_inf[j]<<std::endl;
		}

		for (int j=0; j <ObjectInput2.list_of_uncertainty_factors_classic_plog.size(); j++)
		{
			//std::cout<<"Checking the constrained PLOG number " << ObjectInput2.list_of_uncertainty_factors_classic_plog[j]<<std::endl;
			int pos_classic_plog_reaction = std::find(ObjectInput2.indices_of_classic_plogs.begin(),ObjectInput2.indices_of_classic_plogs.end(),ObjectInput2.list_of_target_classic_plog_reactions[j])-ObjectInput2.indices_of_classic_plogs.begin();

			for (int k=0; k < ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++)
			{
				double A_CP_trial = std::exp(ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k]);
				double n_CP_trial = ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).Beta()[k];
				double E_over_R_CP_trial = ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[k];

				std::vector<double> k_check_CP;
				k_check_CP.resize(T_span.size());

				for (int i=0; i < T_span.size(); i++)
				{
					k_check_CP[i]   = A_CP_trial * std::pow(T_span[i],n_CP_trial) * std::exp((-1*E_over_R_CP_trial)/T_span[i]);

					if ((k_check_CP[i]<=k_LB_classic_plog[j][k][i]) || (k_check_CP[i]>=k_UB_classic_plog[j][k][i]))
					{
						if (ObjectInput2.Objective_Function == "CurveMatching")
						{
								fn_val = 1;
						}else {
								fn_val = 10000000;
						}
						std::cout<<"Violation for PLOG reaction "<<ObjectInput2.list_of_target_classic_plog_reactions[j] <<std::endl;
						return true;
					}
				}
			}

			//std::cout<<"Checked the constrained PLOG number " << ObjectInput2.list_of_target_classic_plog_reactions[j]<<std::endl;
		}

		for (int j=0; j <ObjectInput2.list_of_uncertainty_factors_extplog.size(); j++)
		{
			//std::cout<<"Checking the constrained ExtPLOG number " << ObjectInput2.list_of_uncertainty_factors_extplog[j]<<std::endl;
			int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),ObjectInput2.list_of_target_extplog[j])-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();

			for (int k=0; k < ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size(); k++)
			{
				double A_EP_trial = std::exp(ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0)[k]);
				double n_EP_trial = ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0)[k];
				double E_over_R_EP_trial = ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0)[k];

				std::vector<double> k_check_CP;
				k_check_CP.resize(T_span.size());

				for (int i=0; i < T_span.size(); i++)
				{
					k_check_CP[i]   = A_EP_trial * std::pow(T_span[i],n_EP_trial) * std::exp((-1*E_over_R_EP_trial)/T_span[i]);

					if ((k_check_CP[i]<=k_LB_ext_plog[j][k][i]) || (k_check_CP[i]>=k_UB_ext_plog[j][k][i]))
					{
						if (ObjectInput2.Objective_Function == "CurveMatching")
						{
								fn_val = 1;
						}else {
								fn_val = 10000000;
						}
						std::cout<<"Violation for ExtPLOG reaction "<<ObjectInput2.list_of_target_extplog[j] <<std::endl;
						return true;
					}
				}
			}

			//std::cout<<"Checked the constrained ExtPLOG number " << ObjectInput2.list_of_target_extplog[j]<<std::endl;
		}

		for (int j=0; j <ObjectInput2.list_of_target_EPLR.size(); j++)
		{
			//std::cout<<"Checking the constrained EPLR number " << ObjectInput2.list_of_target_EPLR[j] <<std::endl;

			int pos_EPLR          = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),ObjectInput2.list_of_target_EPLR[j])-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();
			int pos_BathGas_EPLR  = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).get_bath_gas_position(ObjectInput2.list_of_bath_gases_EPLR[j]);
			
			for (int k=0; k < ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size(); k++)
			{
				double A_EP_trial = std::exp(ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR)[k]);
				double n_EP_trial = ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).Beta(pos_BathGas_EPLR)[k];
				double E_over_R_EP_trial = ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).E_over_R(pos_BathGas_EPLR)[k];

				std::vector<double> k_check_CP;
				k_check_CP.resize(T_span.size());

				for (int i=0; i < T_span.size(); i++)
				{
					k_check_CP[i]   = A_EP_trial * std::pow(T_span[i],n_EP_trial) * std::exp((-1*E_over_R_EP_trial)/T_span[i]);

					if ((k_check_CP[i]<=k_LB_EPLR[j][k][i]) || (k_check_CP[i]>=k_UB_EPLR[j][k][i]))
					{
						if (ObjectInput2.Objective_Function == "CurveMatching")
						{
								fn_val = 1;
						}else {
								fn_val = 10000000;
						}
						std::cout<<"Violation for ExtPLOG reaction "<<ObjectInput2.list_of_target_EPLR[j] <<std::endl;
						return true;
					}
				}
			}

			//std::cout<<"Checked the constrained ExtPLOG number " << ObjectInput2.list_of_target_EPLR[j]<<std::endl;
		}

	return false;
}

std::vector <double> OpenSMOKEDirectApplicInterface::reconstuction_pca(int i, std::vector<double> &z_parameters, std::vector<std::vector<std::vector<double>>> &eigenvectors, std::vector<std::vector<double>> &scaling, std::vector<std::vector<double>> &centering)
{
		// ObjectInput2.pca_eigenvectors_direct_reactions

		// projection
		std::vector<double> sum_variable(3,0);

		for (int t=0; t<eigenvectors[i][0].size(); t++){

			for (int j=0; j<z_parameters.size(); j++){

					sum_variable[t] = sum_variable[t] + z_parameters[j]* eigenvectors[i][j][t];
			}
		}

		for (int z=0 ; z<sum_variable.size(); z++){

			std::cout<<"For reaction"<< i <<", the value of projected sample in position "<< z << " is "<< sum_variable[z]<<std::endl;
		}

		// unscaling
		std::vector<double> unscaled_parameters(3,0);
		for (int t=0; t<sum_variable.size(); t++){

			unscaled_parameters[t] = sum_variable[t]*scaling[i][t];
		}

		for (int z=0 ; z<unscaled_parameters.size(); z++){

			std::cout<<"For reaction"<< i <<", the value of unscaled sample in position "<< z << " is "<< unscaled_parameters[z] <<std::endl;
		}

		// decentering
		std::vector<double> decentered_parameters(3,0);
		for (int t=0; t<unscaled_parameters.size(); t++){

			decentered_parameters[t] = unscaled_parameters[t]+centering[i][t];
		}

		for (int z=0 ; z<decentered_parameters.size(); z++){

			std::cout<<"For reaction"<< i <<", the value of projected sample in position "<< z << " is "<< decentered_parameters[z] <<std::endl;
		}

		return decentered_parameters;
}

std::vector <double> OpenSMOKEDirectApplicInterface::ica_mixing(int i, std::vector<double> &unmixed_parameters, std::vector<std::vector<double>> &mixing_matrix)
{
		// mixing
		std::vector<double> mixed_parameters(3,0);

		for (int t=0; t<mixing_matrix.size(); t++){

			for (int j=0; j<mixing_matrix[t].size(); j++){

					mixed_parameters[t] = mixed_parameters[t] + unmixed_parameters[j]* mixing_matrix[t][j];
			}
		}

		for (int z=0 ; z<mixed_parameters.size(); z++){

			std::cout<<"For reaction"<< i <<", the value of projected sample in position "<< z << " is "<< mixed_parameters[z]<<std::endl;
		}

		return mixed_parameters;
}

// Function for changing lnA in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangelnA(int i, double lnA_new)
{
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change lnA[i-1] = " << std::log(ObjectInput2.kineticsMapXML->A(i-1))<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_A(i-1, std::exp(lnA_new));
			std::cout<<"After change lnA[i-1]  = " << std::log(ObjectInput2.kineticsMapXML->A(i-1))<< " i = "<< i <<std::endl;	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_A(i-1, std::exp(lnA_new));
		}
}
    
// Function for changing Beta in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangeBeta(int i, double Beta_new)
{
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change Beta[i-1] = " << ObjectInput2.kineticsMapXML->Beta(i-1)<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_Beta(i-1, Beta_new);
			std::cout<<"After change Beta[i-1]  = " << ObjectInput2.kineticsMapXML->Beta(i-1)<< " i = "<< i <<std::endl;	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_Beta(i-1, Beta_new);
		}
}
	
// Function for changing E_over_R in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangeE_over_R(int i, double E_over_R_new)
{
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change E_over_R[i-1] = " << ObjectInput2.kineticsMapXML->E_over_R(i-1)<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_E_over_R(i-1, E_over_R_new);
			std::cout<<"After change E_over_R[i-1]  = " << ObjectInput2.kineticsMapXML->E_over_R(i-1)<< " i = "<< i <<std::endl;	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_E_over_R(i-1, E_over_R_new);
		}
}

// Function for changing lnA_inf in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangelnA_inf(int i, double lnA_inf_new)
{
		// Finding position of the fall off reaction
		int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),i)-ObjectInput2.indices_of_falloff_reactions.begin();
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change lnA_inf__[pos_FallOff_Reaction] = " << std::log(ObjectInput2.kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction))<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_A_falloff_inf(pos_FallOff_Reaction, std::exp(lnA_inf_new));
			std::cout<<"After change lnA_inf__[pos_FallOff_Reaction] = " << std::log(ObjectInput2.kineticsMapXML->A_falloff_inf(pos_FallOff_Reaction))<< " i = "<< i <<std::endl;	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_A_falloff_inf(pos_FallOff_Reaction, std::exp(lnA_inf_new));
		}
}
    
// Function for changing Beta_inf in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangeBeta_inf(int i, double Beta_inf_new)
{
		// Finding position of the fall off reaction
		int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),i)-ObjectInput2.indices_of_falloff_reactions.begin();
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change Beta_inf__[pos_FallOff_Reaction] = " << ObjectInput2.kineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction)<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_Beta_falloff_inf(pos_FallOff_Reaction, Beta_inf_new);
			std::cout<<"After change Beta_inf__[pos_FallOff_Reaction]  = " << ObjectInput2.kineticsMapXML->Beta(pos_FallOff_Reaction)<< " i = "<< i <<std::endl;	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_Beta_falloff_inf(pos_FallOff_Reaction, Beta_inf_new);
		}
}
	
// Function for changing E_over_R_inf in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangeE_over_R_inf(int i, double E_over_R_inf_new)
{
		// Finding position of the fall off reaction
		int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),i)-ObjectInput2.indices_of_falloff_reactions.begin();
		// Chaning the parameter value and print out parameter value before and after the change if Debug_ChangeParam is true
		if (Debug_ChangeParam == true)
		{
			std::cout<<"Before change E_over_R_inf__[pos_FallOff_Reaction] = " << ObjectInput2.kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction)<< " i = "<< i <<std::endl;	
			ObjectInput2.kineticsMapXML->Set_E_over_R_falloff_inf(pos_FallOff_Reaction, E_over_R_inf_new);
			std::cout<<"After change E_over_R_inf__[pos_FallOff_Reaction]  = " << ObjectInput2.kineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction)<< " i = "<< i <<std::endl;
	
		} else	
		{
			ObjectInput2.kineticsMapXML->Set_E_over_R_falloff_inf(pos_FallOff_Reaction, E_over_R_inf_new);
		}
}

// Function for changing third body efficiencies in the kinetics.
void OpenSMOKEDirectApplicInterface::ChangeThirdBody_Eff(int i, std::string Species, double ThirdBody_Eff_new)
{
		// Finding position of the fall off reaction
		//int pos_FallOff_Reaction = std::find(ObjectInput2.indices_of_falloff_reactions.begin(),ObjectInput2.indices_of_falloff_reactions.end(),i)-ObjectInput2.indices_of_falloff_reactions.begin();
		// Finding the index of the third body species
        	int iSpecies = ObjectInput2.thermodynamicsMapXML->IndexOfSpecies(Species);
		// Finding position of the third body species to be changed
		//int pos_ThirdBody_Species = std::find(ObjectInput2.falloff_indices_of_thirdbody_species[pos_FallOff_Reaction].begin(),ObjectInput2.falloff_indices_of_thirdbody_species[pos_FallOff_Reaction].end(),iSpecies)-ObjectInput2.falloff_indices_of_thirdbody_species[pos_FallOff_Reaction].begin();
		// Changing the value of the thirdbody species
		if (Debug_ChangeParam ==true)
		{
			std::cout<<"Before change M[i][Species] = "<<ObjectInput2.kineticsMapXML->ThirdBody(i-1, iSpecies-1)<<" i = "<< i <<" "<<Species<<std::endl;
			ObjectInput2.kineticsMapXML->Set_ThirdBody(i-1, iSpecies-1, ThirdBody_Eff_new);
			std::cout<<"After change M[i][Species]  = "<<ObjectInput2.kineticsMapXML->ThirdBody(i-1, iSpecies-1)<<" i = "<< i <<" "<<Species<<std::endl;
		} else
		{
			ObjectInput2.kineticsMapXML->Set_ThirdBody(i-1, iSpecies-1, ThirdBody_Eff_new);
		}
}

void OpenSMOKEDirectApplicInterface::Change_ExtPLOG_TB(int r_index, std::string Species, double new_TB)
{
		int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();


		int pos_extended_plog_species = Extract_SpeciesPos_ExtPLOG(pos_extended_plog_reaction, Species);

		if (Debug_ChangeParam ==true)
		{
			ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_TB(pos_extended_plog_species, new_TB);
			for (int i=0; i < ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).ThirdBody(pos_extended_plog_species).size(); i++)
			{
				std::cout<<"After change ThirdBody["<< r_index <<"]["<<Species<<"]["<< i <<"] = "<<ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).ThirdBody(pos_extended_plog_species)[i] <<" "<<std::endl;
			}
		} else
		{
			ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_TB(pos_extended_plog_species, new_TB);
		}
}  

// CLASSIC PLOGs
void OpenSMOKEDirectApplicInterface::ChangelnA_classic_PLOG(int r_index, double lnA_plog_NEW)
{
		int pos_classic_plog_reaction = std::find(ObjectInput2.indices_of_classic_plogs.begin(),ObjectInput2.indices_of_classic_plogs.end(),r_index)-ObjectInput2.indices_of_classic_plogs.begin();
		
		if (Debug_ChangeParam ==true)
		{
			ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_lnA(lnA_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA()) ;
			
			for (int i=0; i < ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); i++)
			{
				std::cout<<"After change lnA_classic_PLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA()[i] << " " <<std::endl;
			}
		} else
		{
			ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_lnA(lnA_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA()) ;
		}
		std::cout<<"All the lnA_classic_PLOG have been replaced" <<std::endl;
}

void OpenSMOKEDirectApplicInterface::Change_ER_classic_PLOG(int r_index, double E_over_R_plog_NEW)
{
		// Finding position of the string Species within the reaction indexed as r_index
		//std::cout<<"The reaction is " << std::to_string(r_index) <<std::endl;
		//std::cout<<"The Species is " << Species <<std::endl;
		//std::cout<<"The pressure index is " << std::to_string(p_index) <<std::endl;

		//ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions()
		int pos_classic_plog_reaction = std::find(ObjectInput2.indices_of_classic_plogs.begin(),ObjectInput2.indices_of_classic_plogs.end(),r_index)-ObjectInput2.indices_of_classic_plogs.begin();
		//std::cout<<"The position of the reaction within the structure of extended plog is " << std::to_string(pos_extended_plog_reaction) <<std::endl;

		if (Debug_ChangeParam ==true)
		{
			ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_E_over_R(E_over_R_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()) ;
			
			for (int i=0; i < ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); i++)
			{
				std::cout<<"After change ER_classic_PLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[i] << " " <<std::endl;
			}
		} else
		{
			ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_E_over_R(E_over_R_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()) ;
		}
		//std::cout<<"All the lnA_classic_PLOG have been replaced" <<std::endl;
} 

void OpenSMOKEDirectApplicInterface::Change_Beta_classic_PLOG(int r_index, double Beta_plog_NEW)
{
	int pos_classic_plog_reaction = std::find(ObjectInput2.indices_of_classic_plogs.begin(),ObjectInput2.indices_of_classic_plogs.end(),r_index)-ObjectInput2.indices_of_classic_plogs.begin();

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_Beta(Beta_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).Beta()) ;
		
		for (int i=0; i < ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); i++)
		{
			 std::cout<<"After change Beta_classic_PLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).Beta()[i] << " " <<std::endl;
		}
	} else{
		ObjectInput2.kineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).modify_Beta(Beta_plog_NEW, ObjectInput2.nominalkineticsMapXML->pressurelog_reactions(pos_classic_plog_reaction).Beta()) ;
	}
	
}

// EPLR
void OpenSMOKEDirectApplicInterface::ChangelnA_EPLR(int r_index, double coefficient_new, std::string bath_gas)
{
	int pos_EPLR          = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();
	int pos_BathGas_EPLR  = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).get_bath_gas_position(bath_gas);

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_A(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR));
		
		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size(); i++)
		{
			std::cout<<"After change lnA_EPLR["<< r_index <<"]["<< i <<"] = "<< ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR)[i] << " " <<std::endl;
		}

	} else
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_A(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR));
	}
}

void OpenSMOKEDirectApplicInterface::Change_ER_EPLR(int r_index, double coefficient_new, std::string bath_gas)
{		
	int pos_EPLR          = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();
	int pos_BathGas_EPLR  = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).get_bath_gas_position(bath_gas);

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_E_over_R(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).E_over_R(pos_BathGas_EPLR));
		
		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size(); i++)
		{
			std::cout<<"After change E_over_R_EPLR["<< r_index <<"]["<< i <<"] = "<< ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR)[i] << " " <<std::endl;
		}
		
	} else
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_E_over_R(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).E_over_R(pos_BathGas_EPLR));
	}
} 

void OpenSMOKEDirectApplicInterface::Change_Beta_EPLR(int r_index, double coefficient_new, std::string bath_gas)
{
	int pos_EPLR          = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelog_reactions().begin();
	int pos_BathGas_EPLR  = ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).get_bath_gas_position(bath_gas);

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_Beta(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).Beta(pos_BathGas_EPLR));
		
		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR).size(); i++)
		{
			std::cout<<"After change Beta_EPLR["<< r_index <<"]["<< i <<"] = "<< ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).lnA(pos_BathGas_EPLR)[i] << " " <<std::endl;
		}
		
	} else
	{
		ObjectInput2.kineticsMapXML->extendedpressurelog_reactions(pos_EPLR).modify_Beta(pos_BathGas_EPLR, coefficient_new, ObjectInput2.nominalkineticsMapXML->extendedpressurelog_reactions(pos_EPLR).Beta(pos_BathGas_EPLR));
	}
}
// EXTENDED PLOGs
void OpenSMOKEDirectApplicInterface::ChangelnA_ExtPLOG(int r_index, double lnA_plog_NEW, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species)
{
	int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();

	std::vector<std::string> species_temp;
	for (int j=0; j < list_of_reactions.size(); j++){
		if (list_of_reactions[j] ==r_index){
			species_temp.push_back(list_of_species[j]);
		}			
	}

	std::vector<int> index_species_temp(species_temp.size(),1);
	for (int j=0; j < species_temp.size(); j++){
			index_species_temp[j] = Extract_SpeciesPos_ExtPLOG(pos_extended_plog_reaction, species_temp[j]);
	}

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_A(index_species_temp, lnA_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0));
		
		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0).size(); i++)
		{
			std::cout<<"After change lnA_ExtPLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0)[i] << " " <<std::endl;
		}
	} else
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_A(index_species_temp, lnA_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).lnA(0));
	}
}

void OpenSMOKEDirectApplicInterface::Change_ER_ExtPLOG(int r_index, double E_over_R_plog_NEW, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species)
{		
	int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();

	std::vector<std::string> species_temp;
	for (int j=0; j < list_of_reactions.size(); j++){
		if (list_of_reactions[j] ==r_index){
			species_temp.push_back(list_of_species[j]);
		}			
	}

	std::vector<int> index_species_temp(species_temp.size(),1);
	for (int j=0; j < species_temp.size(); j++){
			index_species_temp[j] = Extract_SpeciesPos_ExtPLOG(pos_extended_plog_reaction, species_temp[j]);
	}

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_E_over_R(index_species_temp, E_over_R_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0));
		
		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0).size(); i++)
		{
			std::cout<<"After change E_over_R_ExtPLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0)[i] << " " <<std::endl;
		}
	} else
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_E_over_R(index_species_temp, E_over_R_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).E_over_R(0));
	}
} 

void OpenSMOKEDirectApplicInterface::Change_Beta_ExtPLOG(int r_index, double Beta_plog_NEW, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species)
{
	int pos_extended_plog_reaction = std::find(ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin(),ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().end(),r_index)-ObjectInput2.kineticsMapXML->indices_of_extendedpressurelogopt_reactions().begin();

	std::vector<std::string> species_temp;
	for (int j=0; j < list_of_reactions.size(); j++){
		if (list_of_reactions[j] ==r_index){
			species_temp.push_back(list_of_species[j]);
		}			
	}

	std::vector<int> index_species_temp(species_temp.size(),1);
	for (int j=0; j < species_temp.size(); j++){
			index_species_temp[j] = Extract_SpeciesPos_ExtPLOG(pos_extended_plog_reaction, species_temp[j]);
	}

	if (Debug_ChangeParam ==true)
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_Beta(index_species_temp, Beta_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0));

		for (int i=0; i < ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0).size(); i++)
		{
			std::cout<<"After change Beta_ExtPLOG["<< r_index <<"]["<< i <<"] = "<<ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0)[i] << " " <<std::endl;
		}
	} else
	{
		ObjectInput2.kineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).modify_Beta(index_species_temp, Beta_plog_NEW, ObjectInput2.nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).Beta(0));		
	}
}
// EXTENDED PLOGs - UTILITIES
int OpenSMOKEDirectApplicInterface::Extract_SpeciesPos_ExtPLOG(int extPLOG_pos, std::string species){
	
	int pos;
	for (int k=0; k < ObjectInput2.kineticsMapXML->extendedplogopt_reactions(extPLOG_pos).species().size(); k++ )
	{			
		std::string temp_string = ObjectInput2.kineticsMapXML->extendedplogopt_reactions(extPLOG_pos).species()[k];
		if(!ObjectInput2.kineticsMapXML->extendedplogopt_reactions(extPLOG_pos).species()[k].compare(0,temp_string.length(),species))
		{	
			pos = k;
			break;     
		}
	}
	return pos;
}

void OpenSMOKEDirectApplicInterface::BootStrapping_exp_data(std::vector<std::vector<std::vector<std::vector<double>>>> Exp_data)
{
		// ******************************************************** //
		//															//
		//					PUT BS variations equal to				//
		//					ordinates in Exp_data					//
		//															//
		// ******************************************************** //
		
		bootstrapExp.resize(Exp_data.size());

    	for (int a=0; a<Exp_data.size(); ++a) {

			

			bootstrapExp[a].resize(Exp_data[a].size()); 
			for (int b=0; b<Exp_data[a].size(); ++b) {
				
				
				bootstrapExp[a][b].resize(ObjectInput2.numberOfBootstrapVariations); 
				for (int c=0; c<ObjectInput2.numberOfBootstrapVariations; ++c) {
					
					bootstrapExp[a][b][c].resize(Exp_data[a][b][1].size()); 
					for (int d=0; d <Exp_data[a][b][1].size(); ++d) {

						bootstrapExp[a][b][c][d]=Exp_data[a][b][1][d];
					}
				}
			}		
		}
		
		// ******************************************************** //
		//															//
		//				  Replace BS variations with				//
		//					samples from Gaussians					//
		//															//
		// ******************************************************** //
		//std::cout<< "Number of sigma in distribution "<< ObjectInput2.SigmaExpDistribution <<std::endl; 
		
		for (int i=0; i<bootstrapExp.size(); ++i) {
		
			for (int c=0; c<bootstrapExp[i].size(); ++c){

				for (int a=0; a<bootstrapExp[i][c][0].size(); ++a) {

					
					// create an object of class std::default_random_engine, which generates pseudo-random numbers
					std::default_random_engine generator;
					// initializes the seed of the random_engine_generator
					generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
					// initializes a normal distribution for this experimental point, with mean Exp_data[i][c][1][a] and standard deviation standard_deviations[i][c][a]
					// if (dist=="normal"){
					// 	std::normal_distribution<double> distribution(Exp_data[i][c][1][a],ObjectInput2.standard_deviations[i][c][a]);
					// } else {
					// 	std::normal_distribution<double> distribution(Exp_data[i][c][1][a],ObjectInput2.standard_deviations[i][c][a]);
					// }
					std::normal_distribution<double> distribution(Exp_data[i][c][1][a],ObjectInput2.standard_deviations[i][c][a]);	

					// loop over the number of Bootstrap variations
					for (int b=1; b<ObjectInput2.numberOfBootstrapVariations; ++b) 
					{
						// generate a random number from the distribution
						double number = distribution(generator);
						// if negative replace with 0
						if (number < 0 && possibleNegativeOrdinates == false)
						{
						//std::cout<< "it enters the if" <<std::endl;
							number = 0;
						}
						// if logScale is false, returns the number
						// otherwise, return the log of the number to bootstrap.
						if (logScale == false)
						{
							// if the data is 0, then always put 0, otherwise, replace the nominal value with the sampled number
							if (Exp_data[i][c][1][a]==0)
							{
								bootstrapExp[i][c][b][a] = 0;
							} else
							{
								bootstrapExp[i][c][b][a] = number;

							}
						}else
						{
							// IN PRINCIPLES IT WILL NEVER GO HERE
							bootstrapExp[i][c][b][a] = log(number);
						}	
						
					}
				}
			}
		}
		// ******************************************************** //
		//															//
		//				  Print out the BS object 					//
		//															//
		// ******************************************************** //
		// print out bootstrap
		if(ObjectInput2.print_bootstrap==true)
		{
			for (int k = 0; k < bootstrapExp.size(); k++){
				for (int l = 0; l < bootstrapExp[k].size(); l++){
					for (int m = 0; m < bootstrapExp[k][l].size(); m++){
						for (int t = 0; t < bootstrapExp[k][l][m].size(); t++) {
						std::cout<< "in dataset	"<< k << " for the t. sp. " << l <<" at row " << m << "in column " << t <<" the value is 	"<< bootstrapExp[k][l][m][t] <<std::endl;
						}
					}
				}
			}
		}
}

// Function for calculating the objective function
double OpenSMOKEDirectApplicInterface::Calculate_ObjFunc(std::vector<std::vector<std::vector<double>>> Sim_values)
{
		Dakota::Real Obj_Func_val=0;
		std::vector<std::vector<std::vector<double>>> Diff_obj;

		if (ObjectInput2.Objective_Function == "CurveMatching") {
			
			
			std::vector<std::vector<std::vector<double>>> CM_values;

			indexes.resize(Sim_values.size());
			CM_values.resize(Sim_values.size());
			//std::cout<<"It is up to computing the average score"<<std::endl;
			double tStart_CKI = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			//std::cout <<"Print Indexes? " <<ObjectInput2.print_indexes<< " Print Splines? "<< ObjectInput2.print_splines<<std::endl;
			for (int i=0; i < Sim_values.size(); ++i)
			{	
				CM_values[i].resize(Sim_values[i].size());
				indexes[i].resize(Sim_values[i].size());
				

				for (int j=0; j < Sim_values[i].size(); ++j) {

					CM_values[i][j].resize(ObjectInput2.numberOfBootstrapVariations);
					indexes[i][j].resize(ObjectInput2.numberOfBootstrapVariations);
					
					for (int a=0; a< ObjectInput2.numberOfBootstrapVariations; ++a) {
						CM_values[i][j][a] = 0;
						
						CM_values[i][j][a] = indexes[i][j][a].solve(false, false, splinesExp[i][j][a], Sim_values[i][j], ObjectInput2.Exp_data[i][j][0], i, ObjectInput2.print_indexes, ObjectInput2.print_splines, ObjectInput2.QoI[i]);
						
					}
					
				}
			}
			//std::cout<<"It has finished the scores computation"<<std::endl;
			// AVERAGE OF INDEXES - THAN (1 - AV) - GIVE BACK TO DAKOTA
			
			
			double final_index = 0;
			for (int i =0; i<indexes.size(); ++i)
			{
				
				double i_th_index = 0;

				for (int j=0; j<indexes[i].size(); ++j) {
					
					double j_th_index = 0;

					for (int a=0; a<ObjectInput2.numberOfBootstrapVariations; ++a)
					{
						j_th_index  = j_th_index +CM_values[i][j][a]/ObjectInput2.numberOfBootstrapVariations;	
					}

					i_th_index = i_th_index  + j_th_index/indexes[i].size();
					
				}

				std::cout<<"The CM index of the " << i<<"-th dataset is 	"<<i_th_index/ObjectInput2.numberOfBootstrapVariations<< std::endl;	
				final_index = final_index+i_th_index/indexes.size();
			}
			

			std::cout<<"The final CM index is "<<final_index<< std::endl;
			double tEnd_CKI = OpenSMOKE::OpenSMOKEGetCpuTime();

			std::cout << "Time to compute the Curve Matching using "<<ObjectInput2.numberOfBootstrapVariations<<" bootstrap variations: 	" << tEnd_CKI - tStart_CKI << std::endl;			
			Obj_Func_val = 1 - final_index;
			
		} else if (ObjectInput2.Objective_Function == "L1-norm")
		{
			// first dimension will have the dimensionality of the Datasets
			Diff_obj.resize(Sim_values.size());
			// loop over this dimension
			for (int i = 0; i < Sim_values.size(); i++)
			{
				Diff_obj[i].resize(Sim_values[i].size());
				
				if (ObjectInput2.QoI[i]=="IDT")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
                		    Diff_obj[i][z][m] = std::abs(std::log(ObjectInput2.Exp_data[i][z][1][m])-std::log(Sim_values[i][z][m]))/std::abs(std::log(ObjectInput2.standard_deviations[i][z][m]));
						}
					}
				} else if (ObjectInput2.QoI[i]=="Species" || ObjectInput2.QoI[i]=="LFS" || ObjectInput2.QoI[i]=="m_SP" || ObjectInput2.QoI[i]=="BSF")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
                	        Diff_obj[i][z][m] = std::abs(ObjectInput2.Exp_data[i][z][1][m]-Sim_values[i][z][m])/ObjectInput2.standard_deviations[i][z][m];
						}
					}
				} else if (ObjectInput2.type_of_reactor[i]=="KTP")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
 			                Diff_obj[i][z][m] = std::abs(ObjectInput2.Exp_data[i][z][1][m]-Sim_values[i][z][m])/std::abs(ObjectInput2.Exp_data[i][z][1][m]);
						}
					}
				}
			}
			for (int i=0; i < Diff_obj.size(); i++)
			{
				for(int z=0; z< Diff_obj[i].size(); z++)
				{
					Dakota::Real Obj_Func_val_temp;
					Obj_Func_val_temp=0;
					for (int m=0; m < Diff_obj[i][z].size(); m++)
					{		
						Obj_Func_val_temp+=Diff_obj[i][z][m];
					}

					if (isnan(Obj_Func_val_temp)){
						Obj_Func_val= std::pow(10,9);
						return Obj_Func_val;
						
					} else {
						Obj_Func_val+=(Obj_Func_val_temp/Diff_obj[i][z].size());
					}

					std::cout<<"Objective function value for target "<<ObjectInput2.what_2_calc[i][z]<< " in data set "<<i+1<<":	 	"<< Obj_Func_val_temp/Diff_obj[i][z].size() <<std::endl;
				}
			}
		
		} else if (ObjectInput2.Objective_Function == "L2-norm")
		{
			// first dimension will have the dimensionality of the Datasets
			Diff_obj.resize(Sim_values.size());
			// loop over this dimension
			for (int i = 0; i < Sim_values.size(); i++)
			{
				Diff_obj[i].resize(Sim_values[i].size());

				if (ObjectInput2.QoI[i]=="IDT")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
			                Diff_obj[i][z][m] = std::pow(std::log(ObjectInput2.Exp_data[i][z][1][m])-std::log(Sim_values[i][z][m]),2)/std::pow(std::log(ObjectInput2.standard_deviations[i][z][m]),2);
						}
					}
				} else if (ObjectInput2.QoI[i]=="Species" || ObjectInput2.QoI[i]=="LFS" || ObjectInput2.QoI[i]=="m_SP" || ObjectInput2.QoI[i]=="BSF")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
 			                Diff_obj[i][z][m] = std::pow(ObjectInput2.Exp_data[i][z][1][m]-Sim_values[i][z][m],2)/std::pow(ObjectInput2.standard_deviations[i][z][m],2);
						}
					}
				} else if (ObjectInput2.type_of_reactor[i]=="KTP")
				{
					for (int z=0; z < Sim_values[i].size(); z++)
					{
						Diff_obj[i][z].resize(Sim_values[i][z].size());
						for (int m = 0; m < Sim_values[i][z].size(); m++)
						{	
 			                Diff_obj[i][z][m] = std::pow(ObjectInput2.Exp_data[i][z][1][m]-Sim_values[i][z][m],2)/std::pow(ObjectInput2.Exp_data[i][z][1][m],2);
						}
					}
				}
			}
			for (int i=0; i < Diff_obj.size(); i++)
			{
				for(int z=0; z< Diff_obj[i].size(); z++)
				{
					Dakota::Real Obj_Func_val_temp;

					Obj_Func_val_temp=0;
					for (int m=0; m < Diff_obj[i][z].size(); m++)
					{		
						Obj_Func_val_temp+=Diff_obj[i][z][m];
					}

					if (isnan(Obj_Func_val_temp)){
						Obj_Func_val= std::pow(10,9);
						return Obj_Func_val;
						
					} else {
						Obj_Func_val+=(Obj_Func_val_temp/Diff_obj[i][z].size());
					}

					

					std::cout<<"Objective function value for target "<<ObjectInput2.what_2_calc[i][z]<< " in data set "<<i+1<<":	 	"<< Obj_Func_val_temp/Diff_obj[i][z].size() <<std::endl;
				}
			}
		
		}
		else 
		{
			OpenSMOKE::FatalErrorMessage("Only modified L1-norm and L2-norm, and Curve Matching Index currently available");		
		
		}
		
		return Obj_Func_val;

}

} // namespace SIM
