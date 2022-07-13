#include "version.hh"

#ifndef BUILDVERSION
#define BUILDVERSION "unversioned"
#endif

std::string BuildVersion() {
    return BUILDVERSION;
}
