namespace OptiSMOKE
{
  OptimizedKinetics::OptimizedKinetics(
      const OptiSMOKE::InputManager& data,
      const OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
      OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML
  ) : data_(data), thermodynamicsMapXML_(thermodynamicsMapXML), kineticsMapXML_(kineticsMapXML)
  {
    isChemkinNameSet_ = false;
    MemoryAllocation();
  }

  OptimizedKinetics::~OptimizedKinetics(){}

  void OptimizedKinetics::MemoryAllocation()
  {
    NS_ = thermodynamicsMapXML_->NumberOfSpecies();
    NR_ = kineticsMapXML_->NumberOfReactions();

    species_list_.resize(NS_);
    iSpeciesList_.resize(NS_);

    for (int i = 0; i < NS_; i++)
    {
      iSpeciesList_[i] = true;
      species_list_[i] = thermodynamicsMapXML_->NamesOfSpecies()[i];
    }
  }

  void OptimizedKinetics::WriteOptimizedMechanism() 
  {
    if (isChemkinNameSet_ == false)
      OptiSMOKE::FatalErrorMessage("Error: chemkin path was not set for reduced kinetics");

    if(!fs::exists(chemkin_path_.parent_path()))
      fs::create_directories(chemkin_path_.parent_path());

    // Log file (Setup)
    boost::filesystem::path file_ascii_log_ = chemkin_path_.parent_path() / "log";

    std::ofstream flog;
    flog.open(std::string(file_ascii_log_.string()).c_str(), std::ios::out);
    OpenSMOKE::CheckIfFileIsOpen(flog, file_ascii_log_.string());
    flog.setf(std::ios::scientific);
    std::cout.setstate(std::ios_base::failbit); // Disable video output

    //Preprocessing kinetic file
    ThermoReader_CHEMKIN* thermoreader;

    // Reading thermodynamic database
    thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
    thermoreader->ReadFromFile(data_.kinetics_data().chemkin_thermodynamics().string());

    //Preprocessing Chemical Kinetics
    preprocessor_kinetics_ = new PreProcessorKinetics_CHEMKIN(flog);
    preprocessor_kinetics_->ReadFromASCIIFile(data_.kinetics_data().chemkin_kinetics().string());

    // Preprocessing the thermodynamics
    PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport;
    preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, *preprocessor_kinetics_, flog);
    CheckForFatalError(preprocessor_species_without_transport->Setup());

    // Read kinetics from file
    CheckForFatalError(preprocessor_kinetics_->ReadKineticsFromASCIIFile(preprocessor_species_without_transport->AtomicTable()));
    std::cout.clear(); // Re-enable video output

    delete thermoreader;
    thermoreader = NULL;

    delete preprocessor_species_without_transport;
    preprocessor_species_without_transport = NULL;

    fChemKin_.open(chemkin_path_.c_str(), std::ios::out);

    //Writing elements
    fChemKin_ << "ELEMENTS" << std::endl;
    for (int i = 0; i < thermodynamicsMapXML_->elements().size(); i++)
      fChemKin_ << thermodynamicsMapXML_->elements()[i] << std::endl;
    fChemKin_ << "END" << std::endl << std::endl;

    //Writing species list
    fChemKin_ << "SPECIES" << std::endl;
    for (int i = 0; i < species_list_.size(); i++) {
      fChemKin_ << std::setw(25) << std::left << species_list_[i];
      if ((i + 1) % 5 == 0 && i + 1 != species_list_.size()) {
        fChemKin_ << std::endl;
      }
    }
    fChemKin_ << std::endl << "END" << std::endl << std::endl;

    fChemKin_ << std::endl << "REACTIONS" << std::endl << std::endl;
    std::vector<unsigned int> indices_of_classic_plog = kineticsMapXML_->IndicesOfPLOGReactions();
    std::vector<unsigned int> indices_of_falloff_reactions = kineticsMapXML_->IndicesOfFalloffReactions();

    for (unsigned int k = 0; k < NR_; k++)
    {
      // initialize the stringstream for the reaction data, a string for the specific reaction
      std::stringstream reaction_data;
      std::string reaction_string;

      // Get the reaction string and erase white spaces in the string
      preprocessor_kinetics_->reactions()[k].GetReactionString(thermodynamicsMapXML_->NamesOfSpecies(), reaction_string);
      boost::erase_all(reaction_string, " ");

      // Set precision to the stringstream
      reaction_data.precision(6);
      // if the reaction tag is not one of these special formats
      if(	preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_LINDEMANN_FALLOFF && 
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_LINDEMANN_CABR && 
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_TROE_FALLOFF && 
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_TROE_CABR && 
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_SRI_FALLOFF && 
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_SRI_CABR &&
          preprocessor_kinetics_->reactions()[k].Tag() != PhysicalConstants::REACTION_EXTENDEDFALLOFF )
      {
        // Print normal reaction! or even dummy values for PLOG!
        reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML_->A(k) / preprocessor_kinetics_->reactions()[k].A_conversion();
        reaction_data.precision(6);
        reaction_data.width(12);
        reaction_data << std::fixed << std::right << kineticsMapXML_->Beta(k);
        reaction_data.precision(6);
        reaction_data << std::setw(17) << std::right << kineticsMapXML_->E_over_R(k) * PhysicalConstants::R_cal_mol << std::endl;
        // If the reaction is duplicate
        if(preprocessor_kinetics_->reactions()[k].IsDuplicate() == true)
          reaction_data << " DUPLICATE" << std::endl;

        if(preprocessor_kinetics_->reactions()[k].IsExplicitlyReversible() == true)
        {
          reaction_data << " REV /  ";
          reaction_data.precision(4);
          reaction_data << std::scientific << preprocessor_kinetics_->reactions()[k].A_reversible() / preprocessor_kinetics_->reactions()[k].Arev_conversion() << "  ";
          reaction_data.precision(3);
          reaction_data <<std::fixed <<  preprocessor_kinetics_->reactions()[k].Beta_reversible() << "  ";
          reaction_data.precision(2);
          reaction_data << preprocessor_kinetics_->reactions()[k].E_over_R_reversible() * PhysicalConstants::R_cal_mol;
          reaction_data << "  /" << std::endl;
        }
        // If it is PLOG
        if(preprocessor_kinetics_->reactions()[k].IsPressureLog() == true)
        {
          int pos_classic_plog_reaction = std::find(indices_of_classic_plog.begin(),indices_of_classic_plog.end(),k+1)-indices_of_classic_plog.begin();
          reaction_data.unsetf(std::ios_base::floatfield);
          reaction_data.precision(6);

          for(unsigned int l = 0; l < preprocessor_kinetics_->reactions()[k].plog_coefficients().size() - 2; l++)
          {
            if(l % 4 == 0)
              reaction_data << " PLOG /  ";

            reaction_data << std::showpoint << std::setw(16) << std::scientific << std::left << kineticsMapXML_->pressurelog_reactions(pos_classic_plog_reaction).GetAdjustedCoefficients()[l];

            if((l+1) % 4 == 0 || l == preprocessor_kinetics_->reactions()[k].plog_coefficients().size() - 3)
            {
              reaction_data << "/" << std::endl;
            }
          }                    
        }

        if(preprocessor_kinetics_->reactions()[k].IsJanevLanger() == true)
        {
          reaction_data.unsetf(std::ios_base::floatfield);
          reaction_data.precision(6);
          for(unsigned int l = 0; l < preprocessor_kinetics_->reactions()[k].janev_langer_coefficients().size(); l++)
          {
            if(l % 5 == 0)
              reaction_data << " JAN /  ";

            reaction_data << std::showpoint << preprocessor_kinetics_->reactions()[k].janev_langer_coefficients()[l]<< " ";
            if((l+1) % 5 == 0 || l == preprocessor_kinetics_->reactions()[k].janev_langer_coefficients().size() - 1)
              reaction_data << "/" << std::endl;
          }
        }

        if(preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_THIRDBODY)
        {
          std::vector<double> third_body_efficiencies_ = preprocessor_kinetics_->reactions()[k].third_body_efficiencies();
          for(unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
          {
            unsigned int third_body_index = preprocessor_kinetics_->reactions()[k].third_body_indices()[j];
            reaction_data << species_list_[third_body_index] << "/ ";
            reaction_data.precision(2);
            reaction_data << std::showpoint << std::fixed << std::left << kineticsMapXML_->ThirdBody(k, third_body_index) << "/ ";
          }    
          if(preprocessor_kinetics_->reactions()[k].third_body_efficiencies().size() != 0)
            reaction_data << std::endl;
        }

        if(preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_CHEBYSHEV)
        {
          reaction_data.unsetf(std::ios_base::floatfield);          

          if( preprocessor_kinetics_->reactions()[k].chebyshev_temperature_limits()[0] != 300 || preprocessor_kinetics_->reactions()[k].chebyshev_temperature_limits()[1] != 2500 )
          {
            reaction_data << " TCHEB/ ";
            reaction_data.precision(1);
            reaction_data << std::showpoint <<std::fixed << preprocessor_kinetics_->reactions()[k].chebyshev_temperature_limits()[0] << " " << preprocessor_kinetics_->reactions()[k].chebyshev_temperature_limits()[1];
            reaction_data << " /" << std::endl;
          }

          if(preprocessor_kinetics_->reactions()[k].chebyshev_pressure_limits()[0] != 0.001 ||  preprocessor_kinetics_->reactions()[k].chebyshev_pressure_limits()[1] != 100)
          {
            reaction_data << " PCHEB/ ";
            reaction_data.precision(4);
            reaction_data << std::showpoint <<std::fixed << std::left << preprocessor_kinetics_->reactions()[k].chebyshev_pressure_limits()[0] << " " << preprocessor_kinetics_->reactions()[k].chebyshev_pressure_limits()[1];
            reaction_data << " /" << std::endl;
          }

          unsigned int chebyshev_size = boost::lexical_cast<unsigned int>(preprocessor_kinetics_->reactions()[k].chebyshev_coefficients()[0]) * boost::lexical_cast<unsigned int>(preprocessor_kinetics_->reactions()[k].chebyshev_coefficients()[1]);

          reaction_data.unsetf(std::ios_base::floatfield);
          reaction_data.precision(6);
          for(unsigned int l=0; l < chebyshev_size+2; l++)
          {                    
            if(l%6 == 0)
              reaction_data << " CHEB/ ";

            if(l < 2)
              reaction_data << std::noshowpoint << preprocessor_kinetics_->reactions()[k].chebyshev_coefficients()[l] << " ";
            else
              reaction_data << std::showpoint << preprocessor_kinetics_->reactions()[k].chebyshev_coefficients()[l] << " ";

            if((l+1)%6 == 0)
              reaction_data << " /" << std::endl;

            if(l == chebyshev_size + 1 && (l+1)%6 != 0 )
              reaction_data << " /" << std::endl;
          }
        }
      }
      else if (preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
      {
        // High-pressure kinetic parameters
        reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML_->A_falloff_inf(k) / preprocessor_kinetics_->reactions()[k].A_inf_conversion();
        reaction_data.precision(3);
        reaction_data.width(9);
        reaction_data << std::fixed << std::right << kineticsMapXML_->Beta_falloff_inf(k);
        reaction_data.precision(2);
        reaction_data << std::setw(13) << std::right << kineticsMapXML_->E_over_R_falloff_inf(k) * PhysicalConstants::R_cal_mol << std::endl;
        reaction_data.unsetf(std::ios_base::floatfield);
        // Low-pressure parameters
        OpenSMOKE::ExtendedFallOff extendedFallOff;
        extendedFallOff.Setup(preprocessor_kinetics_->reactions()[k].extendedfalloff_coefficients());
        extendedFallOff.WriteCHEMKINOnASCIIFile(reaction_data);

        // Add third body efficiencies
        bool iThirdBody_ = false;
        std::vector<double> third_body_efficiencies_ = preprocessor_kinetics_->reactions()[k].third_body_efficiencies();
        for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
        {
          int third_body_index = preprocessor_kinetics_->reactions()[k].third_body_indices()[j];
          reaction_data << species_list_[third_body_index] << "/ ";
          reaction_data.precision(2);
          reaction_data << std::fixed << std::showpoint << kineticsMapXML_->ThirdBody(k, third_body_index) << "/  ";
          iThirdBody_ = true;
        }

        if (iThirdBody_ == true)
          reaction_data << std::endl;
      }
      else
      {
        int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),k+1)-indices_of_falloff_reactions.begin();
        reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << kineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction) / preprocessor_kinetics_->reactions()[k].A_inf_conversion();
        reaction_data.precision(6);
        reaction_data.width(12);
        reaction_data <<std::fixed << std::right << kineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);reaction_data.precision(4);
        reaction_data.precision(4);
        reaction_data << std::setw(15) << std::right << kineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction) * PhysicalConstants::R_cal_mol << std::endl;

        if(preprocessor_kinetics_->reactions()[k].IsDuplicate() == true)
          reaction_data << " DUPLICATE" << std::endl;

        if(preprocessor_kinetics_->reactions()[k].IsExplicitlyReversible() == true)
        {
          reaction_data << " REV /  ";
          reaction_data.precision(4);
          reaction_data << std::scientific << preprocessor_kinetics_->reactions()[k].A_reversible() / preprocessor_kinetics_->reactions()[k].Arev_conversion()<< "  ";
          reaction_data.precision(3);
          reaction_data <<std::fixed << preprocessor_kinetics_->reactions()[k].Beta_reversible() << "  ";
          reaction_data.precision(2);
          reaction_data << preprocessor_kinetics_->reactions()[k].E_over_R_reversible() * PhysicalConstants::R_cal_mol;
          reaction_data << "  /" << std::endl;
        }            

        if(	preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_FALLOFF      ||
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ||
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_FALLOFF)
        {
          reaction_data << " LOW/";
        }
        else if(preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_CABR ||
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ||
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_CABR)
        {
          reaction_data << " HIGH/";
        }

        // print out the parameter values with the desired formatting and precision // For LOW or HIGH
        reaction_data.width(15);
        reaction_data.precision(6);
        reaction_data << std::right << std::scientific << kineticsMapXML_->A(k) / preprocessor_kinetics_->reactions()[k].A_conversion();
        reaction_data.precision(6);
        reaction_data.width(14);
        reaction_data <<std::fixed << std::right << kineticsMapXML_->Beta(k);
        reaction_data.precision(4);
        reaction_data << std::setw(16) << std::right << kineticsMapXML_->E_over_R(k) * PhysicalConstants::R_cal_mol << " /" << std::endl;

        // prints out the TROE parameters for the blending function F
        if( preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_FALLOFF || 
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_TROE_CABR)
        {
          reaction_data << "TROE/";
          reaction_data.width(11);
          reaction_data.unsetf(std::ios_base::floatfield);
          for(unsigned int j = 0; j < preprocessor_kinetics_->reactions()[k].troe().size(); j++)
          {
            reaction_data.precision(4);
            reaction_data << std::showpoint << "   " << preprocessor_kinetics_->reactions()[k].troe()[j];
          }
          reaction_data << "/" << std::endl;                        
        }

        if(	preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ||
            preprocessor_kinetics_->reactions()[k].Tag() == PhysicalConstants::REACTION_SRI_CABR)
        {
          reaction_data << "SRI/ ";
          for(unsigned int j = 0; j < preprocessor_kinetics_->reactions()[k].sri().size(); j++)
          {
            reaction_data.precision(4);
            reaction_data << std::showpoint << "  " << preprocessor_kinetics_->reactions()[k].sri()[j];
          }
          reaction_data << "/" << std::endl;
        }

        bool iThirdBody_ = false;
        std::vector<double> third_body_efficiencies_ = preprocessor_kinetics_->reactions()[k].third_body_efficiencies();
        for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
        {
          int third_body_index = preprocessor_kinetics_->reactions()[k].third_body_indices()[j];
          reaction_data << species_list_[third_body_index] << "/ ";
          reaction_data.precision(5);
          reaction_data << std::fixed << std::showpoint << kineticsMapXML_->ThirdBody(k, third_body_index) << "/  ";
          iThirdBody_ = true;
        }
        if(iThirdBody_ == true)
          reaction_data << std::endl;  
      }
      if(preprocessor_kinetics_->reactions()[k].IsFit1())
      {
        reaction_data.unsetf(std::ios_base::floatfield);
        reaction_data.precision(6);
        for(unsigned int l = 0; l < preprocessor_kinetics_->reactions()[k].fit1_coefficients().size(); l++)
        {
          reaction_data << " FIT1 /  ";
          reaction_data << std::showpoint << preprocessor_kinetics_->reactions()[k].fit1_coefficients()[l]<< " ";
        }
        reaction_data << "/" << std::endl;
      }
      if(preprocessor_kinetics_->reactions()[k].IsLandauTeller())
      {
        reaction_data.unsetf(std::ios_base::floatfield);
        reaction_data.precision(3);
        for(unsigned int l = 0; l < preprocessor_kinetics_->reactions()[k].landau_teller_coefficients().size(); l++)
        {
          reaction_data << " LT /  ";
          reaction_data << std::showpoint << preprocessor_kinetics_->reactions()[k].landau_teller_coefficients()[l]<< " ";
        }
        reaction_data << "/" << std::endl;
      }
      if(preprocessor_kinetics_->reactions()[k].IsFORD())
      {
        reaction_data.unsetf(std::ios_base::floatfield);
        reaction_data.precision(4);
        std::vector<double> reactant_lambda_ = preprocessor_kinetics_->reactions()[k].reactant_lambda();
        std::vector<unsigned int> reactant_lambda_indices_ = preprocessor_kinetics_->reactions()[k].reactant_lambda_indices();

        for(unsigned int l = 0; l < reactant_lambda_.size(); l++)
        {
          reaction_data << " FORD /  ";
          int index = reactant_lambda_indices_[l];
          reaction_data << species_list_[index] << "  ";
          reaction_data << std::showpoint <<std::fixed << reactant_lambda_[l];
          reaction_data << "/" << std::endl;
        }
      }
      fChemKin_ << reaction_data.str();
      fChemKin_ << std::endl;
    }
    // close the reactions part
    fChemKin_ << "END" << std::endl;
    fChemKin_ << std::endl;
    fChemKin_.close();

    delete preprocessor_kinetics_;
    preprocessor_kinetics_ = NULL;
  }

  void OptimizedKinetics::SetChemkinName(const fs::path& path) {
    chemkin_path_ = path;
    isChemkinNameSet_ = true;
  }
} // namspace OptiSMOKE
