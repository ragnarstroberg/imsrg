#include "version.hh"

#include <string>

#ifndef BUILDVERSION
#define BUILDVERSION "unversioned"
#endif

std::string BuildVersion() {
    return BUILDVERSION;
}
