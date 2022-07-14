#include "version.hh"

#include <string>

#ifndef BUILDVERSION
#define BUILDVERSION "unversioned"
#endif

namespace version {
std::string BuildVersion() {
    return BUILDVERSION;
}
}
