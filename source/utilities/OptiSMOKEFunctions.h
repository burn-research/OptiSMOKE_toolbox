#pragma once

#define __OPTISMOKE_VERSION__ "3.0.0-beta"
#define OPTISMOKE_FATAL_ERROR_EXIT -1
#define OPTISMOKE_SUCCESSFULL_EXIT 0

namespace OptiSMOKE {

void ErrorMessage(const std::string functionName, const std::string errorMessage);

int FatalErrorMessage(const std::string errorMessage);

void OptiSMOKE_logo(const std::string application_name, const std::string author_name);
}  // namespace OptiSMOKE

#include "OptiSMOKEFunctions.hpp"
