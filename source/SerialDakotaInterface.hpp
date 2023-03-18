namespace SIM {

	int SerialDakotaInterface::derived_map_ac(const Dakota::String& ac_name)
	{
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
        unsigned int count = 0;
        OptiSMOKE::SimulationsInterface sim_iface(data_);
        sim_iface.Setup();
        sim_iface.run();
        return 0;
    }

} // namespace OptiSMOKE
