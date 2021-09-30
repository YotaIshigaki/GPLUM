#ifndef AMBERFILE_DEBUG_HH
#define AMBERFILE_DEBUG_HH

#include <stdint.h>
#include <string>
#include <vector>

#include "ReadAmber.h"

namespace amberfile {
namespace debug {

void DumpSectionIndex(const amberfile::readamb::SectionIndex& index);
void DumpVector(const char* flag, const std::vector<std::string>& vec);
void DumpVector(const char* flag, const std::vector<int32_t>& vec);
void DumpVector(const char* flag, const std::vector<int64_t>& vec);
void DumpVector(const char* flag, const std::vector<double>& vec);

}  // namespace debug
}  // namespace amberfile

#endif  // AMBERFILE_DEBUG_HH
