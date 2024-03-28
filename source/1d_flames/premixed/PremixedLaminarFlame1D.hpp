namespace OptiSMOKE {
void PremixedLaminarFlame1D::Setup(const std::string input_file_name_,
                                   OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
                                   OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML,
                                   OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML) {

  std::cout << "Ciao sono in setup 1" << std::endl;
  // Pointers
  thermodynamicsMapXML_ = thermodynamicsMapXML;
  kineticsMapXML_       = kineticsMapXML;
  transportMapXML_      = transportMapXML;

  std::cout << "Ciao sono in setup 2" << std::endl;

  // Defines the grammar rules
  OptiSMOKE::Grammar_PremixedLaminarFlame grammar_laminarflame;

  // Define the dictionaries
  const std::string main_dictionary_name_ = "PremixedLaminarFlame1D";
  OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
  dictionaries.ReadDictionariesFromFile(input_file_name_);
  dictionaries(main_dictionary_name_).SetGrammar(grammar_laminarflame);

  // Species bundling
  double species_bundling = 0.;
  if (dictionaries(main_dictionary_name_).CheckOption("@SpeciesBundling") == true) {
    dictionaries(main_dictionary_name_).ReadDouble("@SpeciesBundling", species_bundling);
  }

  // Read kinetic modifier
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@KineticsModifier") == true) {
      std::string name_of_kinetics_modifier_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@KineticsModifier", name_of_kinetics_modifier_subdictionary);

      // OpenSMOKE::KineticsModifier modifier;
      // modifier.SetupFromDictionary(dictionaries(name_of_kinetics_modifier_subdictionary));
      // modifier.Setup(*thermodynamicsMapXML, *kineticsMapXML);
    }
  }

  double inlet_velocity  = 0.;
  double inlet_mass_flux = 0.;

  // Read inlet conditions
  {
    std::vector<std::string> list_of_strings;
    if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
      dictionaries(main_dictionary_name_).ReadOption("@InletStream", list_of_strings);

    // If multiple inlet streams are specified
    if (list_of_strings.size() != 1) {
      inlet_T.resize(list_of_strings.size());
      P_Pa.resize(list_of_strings.size());
      inlet_omega.resize(list_of_strings.size());
      equivalence_ratios.resize(list_of_strings.size());
      for (unsigned int i = 0; i < list_of_strings.size(); i++){
        GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML, inlet_T[i], P_Pa[i],
                                   inlet_omega[i]);
      }
    }
    // If a single inlet stream is defined
    else {
      GetGasStatusFromDictionary(dictionaries(list_of_strings[0]), *thermodynamicsMapXML, inlet_T, P_Pa, inlet_omega,
                                 equivalence_ratios);
    }
  }

  // Read outlet conditions
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@OutletStream") == true) {
      double P_Pa_outlet;
      std::string name_of_gas_status_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@OutletStream", name_of_gas_status_subdictionary);

      GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, outlet_T,
                                 P_Pa_outlet, outlet_omega);

      if (P_Pa_outlet != P_Pa[0]) {
        OpenSMOKE::FatalErrorMessage("The pressure of outlet stream does not match with the inlet stream");
      }
    } else {  // Equilibrium calculation
      std::cout << "Calculation of thermodynamic equilibrium at the outlet section..." << std::endl;
      // Construction of the main object
      OpenSMOKE::ThermodynamicEquilibrium thermodynamic_equilibrium(*thermodynamicsMapXML);
      thermodynamic_equilibrium.SetVerbosityLevel(1);

      // Initialize the thermodynamic equilibrium (pressure and composition)
      const double TFirstGuess = 2000.;
      const double TMinimum    = 1800.;
      thermodynamic_equilibrium.SetTemperaturePressureAndMassFractions(inlet_T[0], P_Pa[0], inlet_omega[0]);

      // Calculate the equilibrium
      double mw;
      OpenSMOKE::OpenSMOKEVectorDouble inlet_x(thermodynamicsMapXML->NumberOfSpecies());
      thermodynamicsMapXML->SetPressure(P_Pa[0]);
      thermodynamicsMapXML->SetTemperature(inlet_T[0]);
      thermodynamicsMapXML->MoleFractions_From_MassFractions(inlet_x.GetHandle(), mw, inlet_omega[0].GetHandle());

      const double H        = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(inlet_x.GetHandle());  // [J/kmol]
      const double H_over_R = H / PhysicalConstants::R_J_kmol;                                               // [K]
      const int flag = thermodynamic_equilibrium.CalculateEquilibriumFixedEnthalpyAndPressure(H_over_R, TFirstGuess);

      outlet_T = thermodynamic_equilibrium.EquilibriumTemperature();
      OpenSMOKE::ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &outlet_omega, true);
      thermodynamic_equilibrium.EquilibriumMassFractions(&outlet_omega);
      std::cout << " * Equilibrium temperature [K]: " << outlet_T << std::endl;
      outlet_T = std::max(TMinimum, outlet_T);
      std::cout << " * Outlet temperature [K]:      " << outlet_T << " (after correction)" << std::endl;
      std::cout << std::endl;
    }
  }

  // Read inlet velocity
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@InletVelocity") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@InletVelocity", value, units);
      if (units == "m/s") inlet_velocity = value;
      else if (units == "cm/s") inlet_velocity = value / 100.;
      else if (units == "mm/s") inlet_velocity = value / 1000.;
      else if (units == "km/h") inlet_velocity = value * 10. / 36.;
      else OpenSMOKE::FatalErrorMessage("Unknown velocity units");
    }
  }

  // Read inlet mass flux
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@InletMassFlux") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@InletMassFlux", value, units);
      if (units == "kg/m2/s") inlet_mass_flux = value;
      else if (units == "kg/m2/min") inlet_mass_flux = value / 60.;
      else if (units == "g/cm2/s") inlet_mass_flux = value * 10.;
      else if (units == "g/cm2/min") inlet_mass_flux = value * 10. / 60.;
      else OpenSMOKE::FatalErrorMessage("Unknown mass flux units");
    }
  }

  std::cout << "Ciao sono in setup 3" << std::endl;
  Eigen::VectorXd w;
  // Adaptive grid
  // OpenSMOKE::Grid1D* grid;
  {
    std::string name_of_adaptive_grid_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@Grid") == true)
      dictionaries(main_dictionary_name_).ReadDictionary("@Grid", name_of_adaptive_grid_subdictionary);

    grid = new OpenSMOKE::Grid1D(dictionaries(name_of_adaptive_grid_subdictionary), w);
  }

  std::cout << "Ciao sono in setup 4" << std::endl;
  flame_premixed =
      new OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *grid);
  std::cout << "Ciao sono in setup 5" << std::endl;

  // Output folder
  if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true) {
    boost::filesystem::path output_folder;
    dictionaries(main_dictionary_name_).ReadPath("@Output", output_folder);
    flame_premixed->SetOutputFolder(output_folder);
  }

  // Solver type
  {
    std::string solver_type;
    if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true) {
      dictionaries(main_dictionary_name_).ReadString("@Type", solver_type);
      flame_premixed->SetSolverType(solver_type);
    }
  }

  // Soret effect
  {
    bool soret = true;
    if (dictionaries(main_dictionary_name_).CheckOption("@Soret") == true)
      dictionaries(main_dictionary_name_).ReadBool("@Soret", soret);
    flame_premixed->SetSoret(soret);
  }

  // Use Kuo's correlation for thermal diffusion coefficients instead of Chapman-Cowling
  // theory
  {
    bool flag = false;
    if (dictionaries(main_dictionary_name_).CheckOption("@SoretKuoCorrelation") == true)
      dictionaries(main_dictionary_name_).ReadBool("@SoretKuoCorrelation", flag);
    flame_premixed->SetSoretKuoCorrelation(flag);
  }

  // Frozen mass diffusivities
  {
    bool frozen_mass_diffusivities = false;
    if (dictionaries(main_dictionary_name_).CheckOption("@FrozenMassDiffusivities") == true)
      dictionaries(main_dictionary_name_).ReadBool("@FrozenMassDiffusivities", frozen_mass_diffusivities);
    flame_premixed->SetFrozenMassDiffusivities(frozen_mass_diffusivities);
  }

  // Radiative heat transfer
  bool radiative_heat_transfer = false;
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@Radiation") == true)
      dictionaries(main_dictionary_name_).ReadBool("@Radiation", radiative_heat_transfer);
    flame_premixed->SetRadiativeHeatTransfer(radiative_heat_transfer);
  }

  // Read environment temperature
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
      if (units == "K") value *= 1.;
      else if (units == "C") value += 273.15;
      else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

      flame_premixed->SetEnvironmentTemperature(value);
    }
  }

  // Simplified fluxes on the boundaries
  {
    bool simplified_boundary_fluxes = false;
    if (dictionaries(main_dictionary_name_).CheckOption("@SimplifiedBoundaryFluxes") == true)
      dictionaries(main_dictionary_name_).ReadBool("@SimplifiedBoundaryFluxes", simplified_boundary_fluxes);
    flame_premixed->SetSimplifiedBoundaryFluxes(simplified_boundary_fluxes);
  }

  // Testing option
  {
    std::string testing_option;
    if (dictionaries(main_dictionary_name_).CheckOption("@TestingOption") == true) {
      dictionaries(main_dictionary_name_).ReadString("@TestingOption", testing_option);
      flame_premixed->SetTestingOption(testing_option);
    }
  }

  // Polimi soot
  // OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true) {
      std::string name_of_polimisoot_analyzer_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot", name_of_polimisoot_analyzer_subdictionary);
      polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML,
                                                       dictionaries(name_of_polimisoot_analyzer_subdictionary));
      // polimi_soot->ClassesFromXMLFile(path_kinetics_output / "kinetics.xml");
      //
      // // Kinetic modifier (if requested)
      // {
      //   std::vector<unsigned int> indices;
      //   std::vector<double> coefficients;
      //   if (polimi_soot->ClassesCorrectionCoefficients(indices, coefficients) == true) {
      //     OpenSMOKE::KineticsModifier modifier;
      //     modifier.SetupFromIndices(indices, coefficients);
      //     modifier.Setup(*thermodynamicsMapXML, *kineticsMapXML);
      //   }
      // }
      //
      // if (polimi_soot->number_sections() != 0) flame_premixed->SetPolimiSoot(polimi_soot);
    }
  }

  // Flammability limits
  {
    flammability_limits = new OpenSMOKE::FlammabilityLimits(*thermodynamicsMapXML, flame_premixed->output_folder());

    // if (dictionaries(main_dictionary_name_).CheckOption("@FlammabilityLimits") ==
    // true) {
    //   std::string name_of_options_subdictionary;
    //   dictionaries(main_dictionary_name_)
    //       .ReadDictionary("@FlammabilityLimits", name_of_options_subdictionary);
    //   flammability_limits->SetupFromDictionary(
    //       dictionaries(name_of_options_subdictionary));
    //   flame_premixed->SetFlammabilityLimits(flammability_limits);
    //
    //   // Read inlet conditions
    //   {
    //     std::vector<std::string> list_of_strings;
    //     if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
    //       dictionaries(main_dictionary_name_)
    //           .ReadOption("@InletStream", list_of_strings);
    //
    //     if (list_of_strings.size() == 1) {
    //       std::vector<std::string> fuel_names;
    //       std::vector<double> moles_fuel;
    //       std::vector<std::string> oxidizer_names;
    //       std::vector<double> moles_oxidizer;
    //       GetGasStatusFromDictionary(dictionaries(list_of_strings[0]),
    //                                  *thermodynamicsMapXML, fuel_names, moles_fuel,
    //                                  oxidizer_names, moles_oxidizer);
    //
    //       flammability_limits->SetFlame(fuel_names, moles_fuel, oxidizer_names,
    //                                     moles_oxidizer);
    //     } else
    //       OpenSMOKE::FatalErrorMessage(
    //           "If the @FlammabilityLimits option is turned on, no multiple
    //           @InletStream " "option is allowed");
    //   }
    // }
  }

  // On the fly PostProcessing
  // OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
  {
    on_the_fly_post_processing =
        new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, flame_premixed->output_folder());

    if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true) {
      std::string name_of_options_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
      on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
      flame_premixed->SetOnTheFlyPostProcessing(on_the_fly_post_processing);
    }
  }

  // Fixed temperature profile
  {
    std::string name_of_fixed_temperature_profile_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@FixedTemperatureProfile") == true) {
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@FixedTemperatureProfile", name_of_fixed_temperature_profile_subdictionary);

      OpenSMOKE::OpenSMOKEVectorDouble x, y;
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_fixed_temperature_profile_subdictionary), x, y, x_variable,
                                 y_variable);

      if (x_variable != "length")
        OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be defined versus spacee");
      if (y_variable != "temperature")
        OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be define the temperature profile");

      flame_premixed->SetFixedTemperatureProfile(x, y);
    }
  }

  // Fixed specific (i.e. per unit area) mass flow rate profile
  {
    std::string name_of_fixed_specific_mass_flow_rate_profile_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@FixedSpecificMassFlowRateProfile") == true) {
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@FixedSpecificMassFlowRateProfile",
                          name_of_fixed_specific_mass_flow_rate_profile_subdictionary);

      OpenSMOKE::OpenSMOKEVectorDouble x, y;
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_fixed_specific_mass_flow_rate_profile_subdictionary), x, y,
                                 x_variable, y_variable);

      if (x_variable != "length")
        OpenSMOKE::FatalErrorMessage("The @FixedSpecificMassFlowRateProfile must be defined versus spacee");
      if (y_variable != "specific-mass-flow-rate")
        OpenSMOKE::FatalErrorMessage(
            "The @FixedSpecificMassFlowRateProfile must be define the specific (i.e. "
            "per unit area) mass flow rate profile");

      flame_premixed->SetFixedSpecificMassFlowRateProfile(x, y);
    }
  }

  // Read outlet temperature
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@FixedOutletTemperature") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@FixedOutletTemperature", value, units);
      if (units == "K") value = value;
      else if (units == "C") value = value + 273.15;
      else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

      flame_premixed->SetFixedOutletTemperature(value);
    }
  }

  // Internal diameter
  double internal_diameter = 0.01;
  if (dictionaries(main_dictionary_name_).CheckOption("@InternalDiameter") == true) {
    std::string units_internal_diameter;
    dictionaries(main_dictionary_name_).ReadMeasure("@InternalDiameter", internal_diameter, units_internal_diameter);
    if (units_internal_diameter == "m") internal_diameter = internal_diameter;
    else if (units_internal_diameter == "cm") internal_diameter = internal_diameter * 0.01;
    else if (units_internal_diameter == "mm") internal_diameter = internal_diameter * 0.001;
    else OpenSMOKE::FatalErrorMessage("Unknown length units");
  }

  // Taylor-Aris correction for mass and heat diffusion
  if (dictionaries(main_dictionary_name_).CheckOption("@TaylorArisCorrection") == true) {
    std::string value;
    OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::TaylorArisCorrection_Type taylorArisCorrection;
    dictionaries(main_dictionary_name_).ReadString("@TaylorArisCorrection", value);
    if (value == "none")
      taylorArisCorrection = OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::TaylorArisCorrection_Type::NONE;
    else if (value == "mass-heat")
      taylorArisCorrection = OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::TaylorArisCorrection_Type::MASS_HEAT;
    else if (value == "mass")
      taylorArisCorrection = OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::TaylorArisCorrection_Type::MASS;
    else if (value == "heat")
      taylorArisCorrection = OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::TaylorArisCorrection_Type::HEAT;
    else
      OpenSMOKE::FatalErrorMessage(
          "Unknown @TaylorArisCorrection value. Available options: none | mass-heat | "
          "mass | heat");

    flame_premixed->SetTaylorArisCorrection(taylorArisCorrection);
  }

  // Wall heat exchange
  if (dictionaries(main_dictionary_name_).CheckOption("@WallHeatExchangeCoefficient") == true ||
      dictionaries(main_dictionary_name_).CheckOption("@WallHeatNusseltNumber") == true) {
    double wall_heat_exchange_coefficient = 0.;
    double wall_heat_nusselt_number       = 0.;

    if (dictionaries(main_dictionary_name_).CheckOption("@WallHeatExchangeCoefficient") == true) {
      std::string units_wall_heat_exchange_coefficient;
      dictionaries(main_dictionary_name_)
          .ReadMeasure("@WallHeatExchangeCoefficient", wall_heat_exchange_coefficient,
                       units_wall_heat_exchange_coefficient);
      if (units_wall_heat_exchange_coefficient == "W/m2/K")
        wall_heat_exchange_coefficient = wall_heat_exchange_coefficient;
      else if (units_wall_heat_exchange_coefficient == "W/m2/C")
        wall_heat_exchange_coefficient = wall_heat_exchange_coefficient;
      else if (units_wall_heat_exchange_coefficient == "kcal/m2/K")
        wall_heat_exchange_coefficient = wall_heat_exchange_coefficient * 4186.8;
      else if (units_wall_heat_exchange_coefficient == "kcal/m2/C")
        wall_heat_exchange_coefficient = wall_heat_exchange_coefficient * 4186.8;

      else OpenSMOKE::FatalErrorMessage("Unknown heat exchange coefficient");
    } else if (dictionaries(main_dictionary_name_).CheckOption("@WallHeatNusseltNumber") == true) {
      dictionaries(main_dictionary_name_).ReadDouble("@WallHeatNusseltNumber", wall_heat_nusselt_number);
    }

    OpenSMOKE::OpenSMOKEVectorDouble wall_temperature_profile_x;
    OpenSMOKE::OpenSMOKEVectorDouble wall_temperature_profile_y;
    {
      std::string name_of_wall_temperature_profile_subdictionary;
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@WallTemperatureProfile", name_of_wall_temperature_profile_subdictionary);
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_wall_temperature_profile_subdictionary),
                                 wall_temperature_profile_x, wall_temperature_profile_y, x_variable, y_variable);

      if (x_variable != "length")
        OpenSMOKE::FatalErrorMessage("The @WallTemperatureProfile must be defined versus space");
      if (y_variable != "temperature")
        OpenSMOKE::FatalErrorMessage("The @WallTemperatureProfile must be define the temperature profile");
    }

    flame_premixed->SetWallHeatExchange(wall_heat_exchange_coefficient, wall_heat_nusselt_number, internal_diameter,
                                        wall_temperature_profile_x, wall_temperature_profile_y);
  }

  // Lewis numbers
  if (dictionaries(main_dictionary_name_).CheckOption("@LewisNumbers") == true) {
    std::string name_of_lewis_numbers_subdictionary;
    dictionaries(main_dictionary_name_).ReadDictionary("@LewisNumbers", name_of_lewis_numbers_subdictionary);

    std::vector<double> lewis_numbers;
    OpenSMOKE::GetLewisNumbersFromDictionary(dictionaries(name_of_lewis_numbers_subdictionary), thermodynamicsMapXML[0],
                                             lewis_numbers);
    flame_premixed->SetLewisNumbers(lewis_numbers);
  }

  // Use Stefan-Maxwell
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@StefanMaxwell") == true) {
      std::string name_of_stefanmaxwell_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@StefanMaxwell", name_of_stefanmaxwell_subdictionary);
      transportMapXML->SetupStefanMaxwellFromDictionary(dictionaries(name_of_stefanmaxwell_subdictionary),
                                                        thermodynamicsMapXML->NamesOfSpecies());
      flame_premixed->SetStefanMaxwell(true);
    }
  }

  // Use multicomponent transport based on MuTLib
  // OpenSMOKE::MulticomponentTransportLibrary* mutlib;
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@MulticomponentTransport") == true) {
      std::string name_of_mutlib_subdictionary;
      dictionaries(main_dictionary_name_).ReadDictionary("@MulticomponentTransport", name_of_mutlib_subdictionary);
      mutlib = new OpenSMOKE::MulticomponentTransportLibrary();
      mutlib->Initialize(*transportMapXML, *thermodynamicsMapXML);
      mutlib->SetupFromDictionary(dictionaries(name_of_mutlib_subdictionary));
      flame_premixed->SetMulticomponentTransport(mutlib);
    }
  }

  // Initialize from backup
  bool use_userdefined_grid_for_backup = false;
  if (dictionaries(main_dictionary_name_).CheckOption("@Backup") == true) {
    if (dictionaries(main_dictionary_name_).CheckOption("@DontUseBackupGrid") == true) {
      dictionaries(main_dictionary_name_).ReadBool("@DontUseBackupGrid", use_userdefined_grid_for_backup);
    }

    if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED) {
      std::cout << "Siamo in backup" << std::endl;
      fs::path path_backup;
      dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

      // Set inlet and outlet values and first guess velocity according to user values
      flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
      flame_premixed->SetOutlet(outlet_T, outlet_omega);
      if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
      else flame_premixed->SetInletMassFlux(inlet_mass_flux);

      // Setup the solution, accordingly to backup file
      flame_premixed->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
    } else if (flame_premixed->solver_type() ==
               OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED) {
      fs::path path_backup;
      dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

      // Set inlet and outlet values and first guess velocity according to user values
      flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
      flame_premixed->SetOutlet(outlet_T, outlet_omega);
      if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
      else flame_premixed->SetInletMassFlux(inlet_mass_flux);

      // Setup the solution, accordingly to backup file
      flame_premixed->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
    }
  } else {
    if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED) {
      flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
      flame_premixed->SetOutlet(outlet_T, outlet_omega);
      if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
      else flame_premixed->SetInletMassFlux(inlet_mass_flux);
      flame_premixed->SetupForFlameSpeed(w);
    } else if (flame_premixed->solver_type() ==
               OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED) {
      flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
      flame_premixed->SetOutlet(outlet_T, outlet_omega);
      if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
      else flame_premixed->SetInletMassFlux(inlet_mass_flux);
      flame_premixed->SetupForBurnerStabilized(w);
    }
  }

  // Use NLS Solver
  {
    bool use_nls_solver = true;
    if (dictionaries(main_dictionary_name_).CheckOption("@UseNlsSolver") == true)
      dictionaries(main_dictionary_name_).ReadBool("@UseNlsSolver", use_nls_solver);
    flame_premixed->SetUseNlsSolver(use_nls_solver);
  }

  // Use DAE Solver
  {
    bool use_dae_solver = true;
    if (dictionaries(main_dictionary_name_).CheckOption("@UseDaeSolver") == true)
      dictionaries(main_dictionary_name_).ReadBool("@UseDaeSolver", use_dae_solver);
    flame_premixed->SetUseDaeSolver(use_dae_solver);
  }

  // Sensitivity Options
  // OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
  if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true) {
    sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
    std::string name_of_sensitivity_options_subdictionary;
    dictionaries(main_dictionary_name_)
        .ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
    sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
    flame_premixed->EnableSensitivityAnalysis(*sensitivity_options);
  }

  // Dae Options
  // DaeSMOKE::DaeSolver_Parameters* dae_parameters;
  dae_parameters = new DaeSMOKE::DaeSolver_Parameters();
  if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true) {
    std::string name_of_subdictionary;
    dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
    dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
  }

  // Nls Options
  // NlsSMOKE::NonLinearSolver_Parameters* nls_parameters;
  nls_parameters = new NlsSMOKE::NonLinearSolver_Parameters();
  if (dictionaries(main_dictionary_name_).CheckOption("@NlsParameters") == true) {
    std::string name_of_subdictionary;
    dictionaries(main_dictionary_name_).ReadDictionary("@NlsParameters", name_of_subdictionary);
    nls_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
  }

  // Pseudo Transient Options
  // NlsSMOKE::FalseTransientSolver_Parameters* false_transient_parameters;
  false_transient_parameters = new NlsSMOKE::FalseTransientSolver_Parameters();
  if (dictionaries(main_dictionary_name_).CheckOption("@FalseTransientParameters") == true) {
    std::string name_of_subdictionary;
    dictionaries(main_dictionary_name_).ReadDictionary("@FalseTransientParameters", name_of_subdictionary);
    false_transient_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
  }

  // Write RHS contributions
  {
    bool write_rhs = false;
    if (dictionaries(main_dictionary_name_).CheckOption("@WriteRHS") == true) {
      dictionaries(main_dictionary_name_).ReadBool("@WriteRHS", write_rhs);
    }
    flame_premixed->SetWriteRHS(write_rhs);
  }

  if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasTemperature") == true) {
    std::string value;
    dictionaries(main_dictionary_name_).ReadString("@DerivativeGasTemperature", value);
    if (value == "upwind") flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
    else if (value == "backward") flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
    else if (value == "forward") flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
    else if (value == "centered") flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
    else OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas temperature");
  }

  if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasMassFractions") == true) {
    std::string value;
    dictionaries(main_dictionary_name_).ReadString("@DerivativeGasMassFractions", value);
    if (value == "upwind") flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
    else if (value == "backward") flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
    else if (value == "forward") flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
    else if (value == "centered") flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
    else OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas mass fractions");
  }

  // Hybrid Method of Moments
  // OpenSMOKE::HMOM* hmom;
  if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true) {
    // hmom = new OpenSMOKE::HMOM();
    // std::string name_of_subdictionary;
    // dictionaries(main_dictionary_name_).ReadDictionary("@HMOM",
    // name_of_subdictionary);
    // hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));
    //
    // flame_premixed->SolveHMOMFromExistingSolution(
    //     *hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);
    //
    // OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_PremixedLaminarFlame1D",
    //                           "Alberto Cuoci (alberto.cuoci@polimi.it)");
    //
    // return OPENSMOKE_SUCCESSFULL_EXIT;
  }

  if (flammability_limits->is_active() == true) {
    // boost::filesystem::path output_folder_root = flame_premixed->output_folder();
    //
    // std::ofstream fOutput((output_folder_root / "FlameSpeeds.out").string().c_str(),
    //                       std::ios::out);
    // fOutput.setf(std::ios::scientific);
    // {
    //   fOutput << std::left;
    //   fOutput << std::setw(8) << "Case(1)";
    //   fOutput << std::setw(20) << "Speed[cm/s](2)";
    //   fOutput << std::setw(20) << "Eq.Ratio(3)";
    //   fOutput << std::setw(20) << "Pressure[atm](4)";
    //   fOutput << std::setw(20) << "TempInlet[K](5)";
    //   fOutput << std::setw(20) << "TempMax[K](6)";
    //   fOutput << std::setw(20) << "Fuel_x(7)";
    //   fOutput << std::setw(20) << "Ox_x(8)";
    //   fOutput << std::setw(20) << "Fuel_w(9)";
    //   fOutput << std::setw(20) << "Ox_w(10)";
    //
    //   {
    //     std::vector<double> X(thermodynamicsMapXML->NumberOfSpecies());
    //     std::vector<double> Y(thermodynamicsMapXML->NumberOfSpecies());
    //     flammability_limits->CompositionFromEquivalenceRatio(1., X, Y);
    //
    //     unsigned int count = 11;
    //     for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies(); j++) {
    //       if (X[j] > 1e-16) {
    //         std::string title = thermodynamicsMapXML->NamesOfSpecies()[j] + "_x";
    //         OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, title, count);
    //       }
    //     }
    //     for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies(); j++) {
    //       if (Y[j] > 1e-16) {
    //         std::string title = thermodynamicsMapXML->NamesOfSpecies()[j] + "_w";
    //         OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, title, count);
    //       }
    //     }
    //   }
    //
    //   fOutput << std::endl;
    // }
    //
    // std::ofstream fFlammability(
    //     (output_folder_root / "FlammabilityLimits.out").string().c_str(),
    //     std::ios::out);
    // fFlammability.setf(std::ios::scientific);
    // {
    //   fFlammability << std::left;
    //   fFlammability << std::setw(8) << "Dummy(1)";
    //   fFlammability << std::setw(20) << "Speed[cm/s](2)";
    //   fFlammability << std::setw(20) << "Eq.Ratio(3)";
    //   fFlammability << std::setw(20) << "Pressure[atm](4)";
    //   fFlammability << std::setw(20) << "TempInlet[K](5)";
    //   fFlammability << std::setw(20) << "TempMax[K](6)";
    //   fFlammability << std::setw(20) << "Fuel_x(7)";
    //   fFlammability << std::setw(20) << "Ox_x(8)";
    //   fFlammability << std::setw(20) << "Fuel_w(9)";
    //   fFlammability << std::setw(20) << "Ox_w(10)";
    //
    //   {
    //     std::vector<double> X(thermodynamicsMapXML->NumberOfSpecies());
    //     std::vector<double> Y(thermodynamicsMapXML->NumberOfSpecies());
    //     flammability_limits->CompositionFromEquivalenceRatio(1., X, Y);
    //
    //     unsigned int count = 11;
    //     for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies(); j++) {
    //       if (X[j] > 1e-16) {
    //         std::string title = thermodynamicsMapXML->NamesOfSpecies()[j] + "_x";
    //         OpenSMOKE::PrintTagOnASCIILabel(20, fFlammability, title, count);
    //       }
    //     }
    //     for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies(); j++) {
    //       if (Y[j] > 1e-16) {
    //         std::string title = thermodynamicsMapXML->NamesOfSpecies()[j] + "_w";
    //         OpenSMOKE::PrintTagOnASCIILabel(20, fFlammability, title, count);
    //       }
    //     }
    //   }
    //
    //   fFlammability << std::endl;
    // }
  }
}

void PremixedLaminarFlame1D::Solve() {
  if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED) {
    if (inlet_omega.size() == 1) {  // Solve only for a single flame
      // time_t timerStart;
      // time_t timerEnd;
      //
      // time(&timerStart);
      std::cout << "Eccoci" << std::endl;
      fs::path output_folder_root = "/dev/null";
      flame_premixed->SetOutputFolder(output_folder_root);
      flame_premixed->SolveFlameSpeedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
      LFS_ = flame_premixed->flame_speed()*100.;
      // time(&timerEnd);
      //
      // std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" <<
      // std::endl;
    } else {  // Solve for several flames
      fs::path output_folder_root = "/dev/null";
      for (unsigned int i = 0; i < inlet_omega.size(); i++) {
        flame_premixed->SetOutputFolder(output_folder_root);
        flame_premixed->ChangeInletConditions(inlet_T[i], P_Pa[i], inlet_omega[i]);
        flame_premixed->SolveFlameSpeedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
      }
    }
  } else if (flame_premixed->solver_type() ==
             OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED) {
    flame_premixed->SolveBurnerStabilizedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
  }
}
}  // namespace OptiSMOKE
