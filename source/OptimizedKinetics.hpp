namespace OptiSMOKE
{
    OptimizedKinetics::OptimizedKinetics(
        const OptiSMOKE::InputManager& data,
        const OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
        const OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML
    ) : data_(data),
        thermodynamicsMapXML_(thermodynamicsMapXML),
        kineticsMapXML_(kineticsMapXML)
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

    void OptimizedKinetics::WriteOptimizedMechanism() {
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
        
        //Writing reactions list
        fChemKin_ << "REACTIONS" << std::endl;
        for (int i = 0; i < NR_; i++) {

			// TODO
		}
        fChemKin_ << "END";
        fChemKin_.close();
        
        delete preprocessor_kinetics_;
        preprocessor_kinetics_ = NULL;
    }

    void OptimizedKinetics::SetChemkinName(const fs::path& path) {
        chemkin_path_ = path;
        isChemkinNameSet_ = true;
    }
} // namspace OptiSMOKE