namespace OptiSMOKE {

void ErrorMessage(const std::string functionName, const std::string errorMessage) {
  std::cout << "Function:     " << functionName << std::endl;
  std::cout << "Fatal error:  " << errorMessage << std::endl;
  std::cout << "Press enter to exit..." << std::endl;
  getchar();
  exit(OPTISMOKE_FATAL_ERROR_EXIT);
}

int FatalErrorMessage(const std::string errorMessage) {
  std::cout << "Fatal error:  " << errorMessage << std::endl;
  std::cout << "Press enter to exit..." << std::endl;
  getchar();
  exit(OPTISMOKE_FATAL_ERROR_EXIT);
  return OPTISMOKE_FATAL_ERROR_EXIT;
}

void OptiSMOKE_logo(const std::string author_name) {
  std::string current_time = __TIME__;
  std::string current_date = __DATE__;
  std::string author_complete = "Authors: " + author_name;
  std::string compilation_time = "Compilation date: " + current_date + " at " + current_time;
  std::string version = "Version: ";
  version += __OPTISMOKE_VERSION__;

  std::cout << "-----------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "          ____        __  _ _____ __  _______  __ __ ______                     " << std::endl;
  std::cout << "         / __ \\____  / /_(_) ___//  |/  / __ \\/ //_// ____/___  ____            " << std::endl;
  std::cout << "        / / / / __ \\/ __/ /\\__ \\/ /|_/ / / / / ,<  / __/ / __ \\/ __ \\           " << std::endl;
  std::cout << "       / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           " << std::endl;
  std::cout << "       \\____/ .___/\\__/_//____/_/  /_/\\____/_/ |_/_____/ .___/ .___/            " << std::endl;
  std::cout << "           /_/                                        /_/   /_/                 " << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "       [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>             " << std::endl;
  std::cout << "           Department of Chemistry, Materials and Chemical Engineering          " << std::endl;
  std::cout << "           Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano      " << std::endl;
  std::cout << "                                                                                " << std::endl;
  std::cout << "       [2] BRITE Research Group <https://brite-research.be>                     " << std::endl;
  std::cout << "           Brussels Institute for Thermal-fluid systems and clean Energy        " << std::endl;
  std::cout << "           Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                 " << std::endl;
  std::cout << std::endl;
  std::cout << "       " << version << std::endl;
  std::cout << "       " << author_complete << std::endl;
  std::cout << "       " << compilation_time << std::endl;
  std::cout << std::endl;
  std::cout << "-----------------------------------------------------------------------------" << std::endl;
}
}  // namespace OptiSMOKE
