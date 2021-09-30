#include "debug.h"

#include <cstdio>

namespace amberfile {
namespace debug {

void DumpSectionIndex(const amberfile::readamb::SectionIndex& index) {
#ifndef NDEBUG
  for (amberfile::readamb::SectionIndex::const_iterator idx = index.begin();
       idx != index.end(); ++idx) {
    printf("<%s, %d?%d.%d>\n", idx->first.c_str(), idx->second.number,
           idx->second.length, idx->second.precision);
  }
#endif  // NDEBUG
}

void DumpVector(const char* flag, const std::vector<std::string>& vec) {
#ifndef NDEBUG
  if (vec.size() > 0) {
    printf("%s\t[%s, ..., %s]\n", flag, vec[0].c_str(),
           vec[vec.size() - 1].c_str());
  } else {
    printf("%s\t[]\n", flag);
  }
#endif  // NDEBUG
}
void DumpVector(const char* flag, const std::vector<int32_t>& vec) {
#ifndef NDEBUG
  if (vec.size() > 0) {
    printf("%s\t[%d, ..., %d]\n", flag, vec[0], vec[vec.size() - 1]);
  } else {
    printf("%s\t[]\n", flag);
  }
#endif  // NDEBUG
}
void DumpVector(const char* flag, const std::vector<int64_t>& vec) {
#ifndef NDEBUG
  if (vec.size() > 0) {
    printf("%s\t[%lld, ..., %lld]\n", flag, vec[0], vec[vec.size() - 1]);
  } else {
    printf("%s\t[]\n", flag);
  }
#endif  // NDEBUG
}
void DumpVector(const char* flag, const std::vector<double>& vec) {
#ifndef NDEBUG
  if (vec.size() > 0) {
    printf("%s\t[%g, ..., %g]\n", flag, vec[0], vec[vec.size() - 1]);
  } else {
    printf("%s\t[]\n", flag);
  }
#endif  // NDEBUG
}

}  // namespace debug
}  // namespace amberfile
