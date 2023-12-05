namespace SIM {

	int SerialDakotaInterface::derived_map_ac(const Dakota::String& ac_name)
	{
        // Brutto coglione se è l'interfaccia seriale perchè c'è una call a MPI DEBUG?
		#ifdef MPI_DEBUG
			Cout << "analysis server " << analysisServerId << " invoking " << ac_name
				<< " within SIM::OpenSMOKEDirectApplicInterface." << std::endl;
		#endif // MPI_DEBUG

		if (multiProcAnalysisFlag) {
			Cerr << "Error: plugin serial direct fn does not support multiprocessor "
				<< "analyses." << std::endl;
            Dakota::abort_handler(-1);
        }
  
        int fail_code = 0;
        if (ac_name == "opensmoke_plugin") {
            fail_code = simulations_interface(xC, directFnASV[0], fnVals[0]);
        }
        else {
            Cerr << ac_name << " is not available as an analysis within "
                << "OptiSMOKE::SerialDakotaInterface." << std::endl;
            Dakota::abort_handler(Dakota::INTERFACE_ERROR);
        }

        // Failure capturing
        if (fail_code) {
            std::string err_msg("Error evaluating plugin analysis_driver ");
            err_msg += ac_name;
            throw Dakota::FunctionEvalFailure(err_msg);
        }

        return 0;
    }

    void SerialDakotaInterface::wait_local_evaluations(Dakota::PRPQueue& prp_queue)
    {
        if (multiProcAnalysisFlag) {
            Cerr << "Error: plugin serial direct fn does not support multiprocessor "
                << "analyses." << std::endl;
            Dakota::abort_handler(-1);
        }
    
        for (Dakota::PRPQueueIter prp_iter = prp_queue.begin(); prp_iter != prp_queue.end(); prp_iter++) 
        {
            // For each job in the processing queue, evaluate the response
            int fn_eval_id = prp_iter->eval_id();
            const Dakota::Variables& vars = prp_iter->variables();
            const Dakota::ActiveSet& set  = prp_iter->active_set();
            Dakota::Response         resp = prp_iter->response(); // shared rep
            if (outputLevel > Dakota::SILENT_OUTPUT)
                Cout << "SerialDakotaInterface: evaluating function evaluation "
                    << fn_eval_id << " in batch mode." << std::endl;
            short asv = set.request_vector()[0];
            Dakota::Real& fn_val = resp.function_value_view(0);
            simulations_interface(vars.continuous_variables(), asv, fn_val);

            // indicate completion of job to ApplicationInterface schedulers
            completionSet.insert(fn_eval_id);
        }
    }

    int SerialDakotaInterface::simulations_interface(const Dakota::RealVector& c_vars, short asv, Dakota::Real& fn_val)
    {
        eval_nr++;
        std::vector<double> b(c_vars.length());
        for(unsigned int j = 0; j < c_vars.length(); j++)
            b[j] = c_vars[j];

        sim_iface_->SubstituteKineticParameters(b);

		if(data_.optimization_setup().penalty_function())
			violated_uncertainty = sim_iface_->CheckKineticConstasts();
        
		if(violated_uncertainty){
			if (data_.optimization_setup().objective_function_type() == "CurveMatching")
				fn_val = 1;
			else
				fn_val = 10000000;
		}
        else{
            sim_iface_->run();
            fn_val = sim_iface_->ComputeObjectiveFunction();
        
            if(eval_nr == 1)
                prev_fn_val = fn_val;
        
            if(prev_fn_val > fn_val) {
                prev_fn_val = fn_val;
                opti_kinetics_->WriteOptimizedMechanism();
                std::cout << " * Wrote optimized mechanism" << std::endl;
            }
        }
        
        return 0;
    }
} // namespace SIM
