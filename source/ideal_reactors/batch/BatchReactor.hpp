namespace OptiSMOKE {

void BatchReactor::Setup(const std::string input_file_name,
                         OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML,
                         OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML) {
  // Pointers
  thermodynamicsMapXML_ = thermodynamicsMapXML;
  kineticsMapXML_       = kineticsMapXML;

  // Defines the grammar rules
  OptiSMOKE::Grammar_BatchReactor grammar_BatchReactor;

  // Define the dictionaries
  std::string main_dictionary_name_ = "BatchReactor";
  OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
  dictionaries.ReadDictionariesFromFile(input_file_name);
  dictionaries(main_dictionary_name_).SetGrammar(grammar_BatchReactor);

  // Read thermodynamics and kinetics maps
  double T, P_Pa;
  OpenSMOKE::OpenSMOKEVectorDouble omega;
  tEnd_         = 0.;
  tStart_       = 0.;  // default 0
  double volume = 1.;  // default value [1 m3]

  // Read initial conditions
  {
    std::string name_of_gas_status_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@InitialStatus") == true)
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@InitialStatus", name_of_gas_status_subdictionary);

    GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary),
                               *thermodynamicsMapXML, T, P_Pa, omega);
  }

  // Read end time
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@EndTime") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@EndTime", value, units);
      if (units == "s") tEnd_ = value;
      else if (units == "ms") tEnd_ = value / 1000.;
      else if (units == "min") tEnd_ = value * 60.;
      else if (units == "h") tEnd_ = value * 3600.;
      else OpenSMOKE::FatalErrorMessage("Unknown time units");
    }
  }

  // Read start time
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@StartTime") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@StartTime", value, units);
      if (units == "s") tStart_ = value;
      else if (units == "ms") tStart_ = value / 1000.;
      else if (units == "min") tStart_ = value * 60.;
      else if (units == "h") tStart_ = value * 3600.;
      else OpenSMOKE::FatalErrorMessage("Unknown time units");
    }
  }

  // Read volume
  {
    if (dictionaries(main_dictionary_name_).CheckOption("@Volume") == true) {
      double value;
      std::string units;

      dictionaries(main_dictionary_name_).ReadMeasure("@Volume", value, units);
      if (units == "m3") volume = value;
      else if (units == "dm3") volume = value / 1.e3;
      else if (units == "cm3") volume = value / 1.e6;
      else if (units == "mm3") volume = value / 1.e9;
      else if (units == "l") volume = value / 1.e3;
      else OpenSMOKE::FatalErrorMessage("Unknown volume units");
    }
  }

  // Read exchange area
  double exchange_area = 0.;
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@ExchangeArea") == true) {
      dictionaries(main_dictionary_name_).ReadMeasure("@ExchangeArea", value, units);
      if (units == "m2") exchange_area = value;
      else if (units == "dm2") exchange_area = value / 1.e2;
      else if (units == "cm2") exchange_area = value / 1.e4;
      else if (units == "mm2") exchange_area = value / 1.e6;
      else OpenSMOKE::FatalErrorMessage("Unknown area units");
    }
  }

  // Read global thermal exchange coefficient
  double global_thermal_exchange_coefficient = 0.;
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_)
            .CheckOption("@GlobalThermalExchangeCoefficient") == true) {
      dictionaries(main_dictionary_name_)
          .ReadMeasure("@GlobalThermalExchangeCoefficient", value, units);
      if (units == "W/m2/K") global_thermal_exchange_coefficient = value;
      else if (units == "W/m2/C") global_thermal_exchange_coefficient = value;
      else if (units == "kcal/m2/K")
        global_thermal_exchange_coefficient = value * 4186.8;
      else if (units == "kcal/m2/C")
        global_thermal_exchange_coefficient = value * 4186.8;
      else
        OpenSMOKE::FatalErrorMessage(
            "Unknown global thermal exchange coefficient units");
    }
  }

  // Environment temperature
  double T_environment = T;
  {
    double value;
    std::string units;
    if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") ==
        true) {
      dictionaries(main_dictionary_name_)
          .ReadMeasure("@EnvironmentTemperature", value, units);
      if (units == "K") T_environment = value;
      else if (units == "C") T_environment = value + 273.15;
      else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
    }
  }

  // Type
  {
    std::string value;
    if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true) {
      dictionaries(main_dictionary_name_).ReadString("@Type", value);
      if (value == "Isothermal-ConstantVolume")
        type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV;
      else if (value == "Isothermal-ConstantPressure")
        type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP;
      else if (value == "NonIsothermal-ConstantVolume")
        type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV;
      else if (value == "NonIsothermal-ConstantPressure")
        type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP;
      else if (value == "NonIsothermal-UserDefinedVolume")
        type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME;
      else OpenSMOKE::FatalErrorMessage("Unknown batch reactor type: " + value);
    }
  }

  // Options
  {
    batch_options = new OpenSMOKE::BatchReactor_Options();
    if (dictionaries(main_dictionary_name_).CheckOption("@Options") == true) {
      std::string name_of_options_subdictionary;
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@Options", name_of_options_subdictionary);
      batch_options->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
    }
    fs::path tmp = "/dev/null";
    batch_options->SetOutputPath(tmp);
    batch_options->SetVerboseVideo(false);
    batch_options->SetVerboseOutput(false);
    batch_options->SetSensitivityAnalysis(false);
  }

  // ODE Parameters
  {
    ode_parameters = new OpenSMOKE::ODE_Parameters();
    if (dictionaries(main_dictionary_name_).CheckOption("@OdeParameters") == true) {
      std::string name_of_ode_parameters_subdictionary;
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@OdeParameters", name_of_ode_parameters_subdictionary);
      ode_parameters->SetupFromDictionary(
          dictionaries(name_of_ode_parameters_subdictionary));
    }
  }

  // Sensitivity Options
  // To check non sensitivty analysis is allowed
  {
    sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
    if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true) {
      //	std::string name_of_sensitivity_options_subdictionary;
      //	dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis",
      // name_of_sensitivity_options_subdictionary);
      //
      //	batch_options->SetSensitivityAnalysis(true);
      //	sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
    }
  }

  // On the fly ROPA
  {
    onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML, *kineticsMapXML);
    if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyROPA") == true) {
      // std::string name_of_options_subdictionary;
      // dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyROPA",
      // name_of_options_subdictionary);
      //  No ROPA (disabled)
      //  onTheFlyROPA->SetupFromDictionary(dictionaries(name_of_options_subdictionary),
      //  path_kinetics_output);
    }
  }

  // On the fly CEMA
  {
    onTheFlyCEMA = new OpenSMOKE::OnTheFlyCEMA(*thermodynamicsMapXML, *kineticsMapXML,
                                               batch_options->output_path());
    if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyCEMA") == true) {
      // std::string name_of_options_subdictionary;
      // dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyCEMA",
      // name_of_options_subdictionary);
      // onTheFlyCEMA->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
    }
  }

  // On the fly PostProcessing
  {
    on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(
        *thermodynamicsMapXML, *kineticsMapXML, batch_options->output_path());

    if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") ==
        true) {
      // std::string name_of_options_subdictionary;
      // dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing",
      // name_of_options_subdictionary);
      // on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
    }
  }

  // Ignition Delay Times
  {
    idt = new OpenSMOKE::IgnitionDelayTimes_Analyzer();
    if (dictionaries(main_dictionary_name_).CheckOption("@IgnitionDelayTimes") == true) {
      if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV ||
          type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
        OpenSMOKE::FatalErrorMessage(
            "The @IgnitionDelayTimes can be used only for NonIsothermal reactors");

      std::string name_of_idt_subdictionary;
      dictionaries(main_dictionary_name_)
          .ReadDictionary("@IgnitionDelayTimes", name_of_idt_subdictionary);
      idt->SetupFromDictionary(dictionaries(name_of_idt_subdictionary),
                               *thermodynamicsMapXML);
    }
  }

  {
    volume_profile_      = false;
    temperature_profile_ = false;
    // Read pressure coefficient
    if (dictionaries(main_dictionary_name_).CheckOption("@PressureCoefficient") ==
        true) {
      batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();
      batchreactor_volumeprofile->SetType(
          OpenSMOKE::BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_Type::
              BatchReactor_VolumeProfile_PressureCoefficient);

      if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
        OpenSMOKE::FatalErrorMessage(
            "The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume "
            "reactors");

      double value;
      std::string units;

      dictionaries(main_dictionary_name_)
          .ReadMeasure("@PressureCoefficient", value, units);
      if (units == "Pa/s") batchreactor_volumeprofile->SetPressureCoefficient(value);
      else if (units == "bar/s")
        batchreactor_volumeprofile->SetPressureCoefficient(value * 1.e5);
      else if (units == "atm/s")
        batchreactor_volumeprofile->SetPressureCoefficient(value * 101325.);
      else if (units == "Pa/ms")
        batchreactor_volumeprofile->SetPressureCoefficient(value * 1000.);
      else if (units == "bar/ms")
        batchreactor_volumeprofile->SetPressureCoefficient(value * 1.e5 * 1000.);
      else if (units == "atm/ms")
        batchreactor_volumeprofile->SetPressureCoefficient(value * 101325. * 1000.);
      else
        OpenSMOKE::FatalErrorMessage(
            "Unknown pressure coefficient units. Available: Pa/s || Pa/ms || bar/s || "
            "bar/ms || atm/s || atm/ms");
      volume_profile_ = true;
    }

    // Read volume profile
    if (dictionaries(main_dictionary_name_).CheckOption("@VolumeProfile") == true) {
      batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

      std::string name_of_profile_subdictionary;

      if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
        OpenSMOKE::FatalErrorMessage(
            "The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume "
            "reactors");

      dictionaries(main_dictionary_name_)
          .ReadDictionary("@VolumeProfile", name_of_profile_subdictionary);

      OpenSMOKE::OpenSMOKEVectorDouble x, y;
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y,
                                 x_variable, y_variable);

      if (x_variable != "time")
        OpenSMOKE::FatalErrorMessage(
            "The @VolumeProfile must be defined versus the time");
      if (y_variable != "volume")
        OpenSMOKE::FatalErrorMessage(
            "The @VolumeProfile must be defined versus the volume");

      if (y[1] != volume)
        OpenSMOKE::FatalErrorMessage(
            "The @VolumeProfile and the @Volume options must be consistent");

      batchreactor_volumeprofile->SetType(
          OpenSMOKE::BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_Type::
              BatchReactor_VolumeProfile_VolumeHistory);
      batchreactor_volumeprofile->SetProfile(1, x, y);
      volume_profile_ = true;
    }

    // Read pressure profile
    if (dictionaries(main_dictionary_name_).CheckOption("@PressureProfile") == true) {
      batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

      std::string name_of_profile_subdictionary;

      if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
        OpenSMOKE::FatalErrorMessage(
            "The @PressureProfile can be used only for NonIsothermal-UserDefinedVolume "
            "reactors");

      dictionaries(main_dictionary_name_)
          .ReadDictionary("@PressureProfile", name_of_profile_subdictionary);

      OpenSMOKE::OpenSMOKEVectorDouble x, y;
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y,
                                 x_variable, y_variable);

      if (x_variable != "time")
        OpenSMOKE::FatalErrorMessage(
            "The @PressureProfile must be defined versus the time");
      if (y_variable != "pressure")
        OpenSMOKE::FatalErrorMessage(
            "The @PressureProfile must be defined versus the volume");

      if (y[1] != P_Pa)
        OpenSMOKE::FatalErrorMessage(
            "The @PressureProfile and the initial pressure of mixture must be "
            "consistent");

      batchreactor_volumeprofile->SetType(
          OpenSMOKE::BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_Type::
              BatchReactor_VolumeProfile_FromPressureProfile);
      batchreactor_volumeprofile->SetProfile(1, x, y);
      volume_profile_ = true;
    }

    // Read temperature profile
    if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true) {
      std::string name_of_profile_subdictionary;

      if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
        OpenSMOKE::FatalErrorMessage(
            "The @TemperatureProfile can be used only for "
            "NonIsothermal-UserDefinedVolume reactors");

      dictionaries(main_dictionary_name_)
          .ReadDictionary("@TemperatureProfile", name_of_profile_subdictionary);

      OpenSMOKE::OpenSMOKEVectorDouble x, y;
      std::string x_variable, y_variable;
      GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y,
                                 x_variable, y_variable);

      if (x_variable != "time")
        OpenSMOKE::FatalErrorMessage(
            "The @TemperatureProfile must be defined versus the time");
      if (y_variable != "temperature")
        OpenSMOKE::FatalErrorMessage(
            "The @TemperatureProfile must be defined versus the temperature");

      if (y[1] != T)
        OpenSMOKE::FatalErrorMessage(
            "The @TemperatureProfile and the initial temperature of mixture must be "
            "consistent");

      // Assign profiles
      if (dictionaries(main_dictionary_name_).CheckOption("@PressureProfile") == false) {
        batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();
        OpenSMOKE::OpenSMOKEVectorDouble xp(2);
        xp[1] = tStart_;
        xp[2] = tEnd_;
        OpenSMOKE::OpenSMOKEVectorDouble yp(2);
        yp[1] = P_Pa;
        yp[2] = P_Pa;
        batchreactor_volumeprofile->SetProfile(1, xp, yp);
      }

      batchreactor_volumeprofile->SetType(
          OpenSMOKE::BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_Type::
              BatchReactor_VolumeProfile_FromPressureAndTemperatureProfiles);
      batchreactor_volumeprofile->SetProfile(2, x, y);
      temperature_profile_ = true;
    }
  }

  // Polimi soot
  // No polimi soot at the moment
  polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);
  {
    std::string name_of_polimisoot_analyzer_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true) {
      //	dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot",
      // name_of_polimisoot_analyzer_subdictionary);
      //	polimi_soot->SetupFromDictionary(dictionaries(name_of_polimisoot_analyzer_subdictionary));
    }
  }

  // Solve the ODE system: NonIsothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV) {
    batch_nonisothermal_constantv_ =
        new OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume(
            *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
            *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt,
            *polimi_soot, volume, T, P_Pa, omega, global_thermal_exchange_coefficient,
            exchange_area, T_environment);
  }

  // Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME) {
    batch_nonisothermal_userdefinedv_ =
        new OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume(
            *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
            *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt,
            *polimi_soot, volume, T, P_Pa, omega, tStart_,
            global_thermal_exchange_coefficient, exchange_area, T_environment);

    batch_nonisothermal_userdefinedv_->SetVolumeProfile(*batchreactor_volumeprofile);
  }

  // Solve the ODE system: Isothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV) {
    batch_isothermal_constantv_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantVolume(
        *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
        *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot,
        volume, T, P_Pa, omega);
  }

  // Solve the ODE system: NonIsothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP) {
    batch_nonisothermal_constantp_ =
        new OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure(
            *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
            *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt,
            *polimi_soot, volume, T, P_Pa, omega, global_thermal_exchange_coefficient,
            exchange_area, T_environment);
  }

  // Solve the ODE system: Isothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP) {
    batch_isothermal_constantp_ =
        new OpenSMOKE::BatchReactor_Isothermal_ConstantPressure(
            *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
            *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt,
            *polimi_soot, volume, T, P_Pa, omega);
  }
}

void BatchReactor::Solve() {
  // std::cout.setstate(std::ios_base::failbit); // Disable video output
  // Solve the ODE system: NonIsothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
    batch_nonisothermal_constantv_->Solve(tStart_, tEnd_);

  // Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
    batch_nonisothermal_userdefinedv_->Solve(tStart_, tEnd_);

  // Solve the ODE system: Isothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
    batch_isothermal_constantv_->Solve(tStart_, tEnd_);

  // Solve the ODE system: NonIsothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
    batch_nonisothermal_constantp_->Solve(tStart_, tEnd_);

  // Solve the ODE system: Isothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
    batch_isothermal_constantp_->Solve(tStart_, tEnd_);
  // std::cout.clear(); // Re-enable video-output
}

double BatchReactor::GetIgnitionDelayTime(std::string criterion) {
  double tau_ign_temp;
  std::string delimiter = "-";

  // --- Think about rewrting this piece
  int pos = 0;
  std::vector<std::string> token;
  std::string tmp = criterion;  // This because I don't want to overwrite criterion
  while ((pos = tmp.find(delimiter)) != std::string::npos) {
    token.push_back(tmp.substr(0, pos));
    tmp.erase(0, pos + delimiter.length());
  }
  token.push_back(tmp);

  if (token[1] == "increase") {
    if (token[0] == "Temperature")
      // Temperature-increase
      tau_ign_temp = idt->temperature_increase_tau();
    else if (token[0] == "Pressure")
      // Pressure-increase
      tau_ign_temp = idt->pressure_increase_tau();
    else
      OptiSMOKE::FatalErrorMessage(
          criterion + " unknown for the evaluation of Ignition Delay Time!");
  } else if (token[1] == "max" && token.size() == 2) {
    if (token[0] == "Temperature")
      // Temperature-max
      tau_ign_temp = idt->temperature_max_tau();
    else if (token[0] == "Pressure")
      // Pressure-max
      tau_ign_temp = idt->pressure_max_tau();
    else {
      // <Species>-max
      std::vector<unsigned int> species_index_temp;
      species_index_temp = idt->species_index();
      int species_index =
          std::distance(species_index_temp.begin(),
                        std::find(species_index_temp.begin(), species_index_temp.end(),
                                  thermodynamicsMapXML_->IndexOfSpecies(token[0]))) -
          1;
      tau_ign_temp = idt->species_max_tau()[species_index];
    }
  } else if (token[1] == "max" && token.size() == 3) {
    if (token[0] == "Temperature")
      // Temperature-max-slope
      tau_ign_temp = idt->temperature_slope_tau();
    else if (token[0] == "Pressure")
      // Pressure-max-slope
      tau_ign_temp = idt->pressure_slope_tau();
    else if (token[2] == "intercept") {
      // <Species>-max-intercept
      std::vector<unsigned int> species_index_temp;
      species_index_temp = idt->species_intercept_max_index();

      if (std::find(species_index_temp.begin(), species_index_temp.end(),
                    thermodynamicsMapXML_->IndexOfSpecies(token[0]) - 1) ==
          species_index_temp.end())
        OptiSMOKE::FatalErrorMessage(
            "Species: " + token[0] +
            " not specified in the OpenSMOKE++ input file properly!");

      int species_index =
          std::distance(species_index_temp.begin(),
                        std::find(species_index_temp.begin(), species_index_temp.end(),
                                  thermodynamicsMapXML_->IndexOfSpecies(token[0]) - 1));
      tau_ign_temp = idt->species_intercept_max_tau()[species_index];
    } else {
      // <Species>-max-slope
      std::vector<unsigned int> species_index_temp;
      species_index_temp = idt->species_index();
      int species_index =
          std::distance(species_index_temp.begin(),
                        std::find(species_index_temp.begin(), species_index_temp.end(),
                                  thermodynamicsMapXML_->IndexOfSpecies(token[0]))) -
          1;
      tau_ign_temp = idt->species_slope_tau()[species_index];
    }
  } else if (token[1] == "min" && token.size() == 3) {
    // <Species>-min-intercept
    std::vector<unsigned int> species_index_temp;
    species_index_temp = idt->species_intercept_min_index();

    if (std::find(species_index_temp.begin(), species_index_temp.end(),
                  thermodynamicsMapXML_->IndexOfSpecies(token[0]) - 1) ==
        species_index_temp.end())
      OptiSMOKE::FatalErrorMessage(
          "Species: " + token[0] +
          " not specified in the OpenSMOKE++ input file properly!");

    int species_index =
        std::distance(species_index_temp.begin(),
                      std::find(species_index_temp.begin(), species_index_temp.end(),
                                thermodynamicsMapXML_->IndexOfSpecies(token[0]) - 1));
    tau_ign_temp = idt->species_intercept_min_tau()[species_index];
  } else
    OptiSMOKE::FatalErrorMessage(
        criterion + " not yet implemented for the evaluation of Ignition Delay Time!");

  CleanMemory();
  return tau_ign_temp;
}

void BatchReactor::CleanMemory() {
  delete batch_options;
  batch_options = NULL;

  delete ode_parameters;
  ode_parameters = NULL;

  delete sensitivity_options;
  sensitivity_options = NULL;

  delete onTheFlyROPA;
  onTheFlyROPA = NULL;

  delete onTheFlyCEMA;
  onTheFlyCEMA = NULL;

  delete on_the_fly_post_processing;
  on_the_fly_post_processing = NULL;

  delete polimi_soot;
  polimi_soot = NULL;

  if (volume_profile_ || temperature_profile_) {
    delete batchreactor_volumeprofile;
    batchreactor_volumeprofile = NULL;
  }

  delete idt;
  idt = NULL;

  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV) {
    delete batch_nonisothermal_constantv_;
    batch_nonisothermal_constantv_ = NULL;
  }

  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME) {
    delete batch_nonisothermal_userdefinedv_;
    batch_nonisothermal_userdefinedv_ = NULL;
  }

  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV) {
    delete batch_isothermal_constantv_;
    batch_isothermal_constantv_ = NULL;
  }

  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP) {
    delete batch_nonisothermal_constantp_;
    batch_nonisothermal_constantp_ = NULL;
  }

  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP) {
    delete batch_isothermal_constantp_;
    batch_isothermal_constantp_ = NULL;
  }
}

}  // namespace OptiSMOKE
