#ifndef AMBERFILE_READAMBER_H
#define AMBERFILE_READAMBER_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace amberfile {
namespace readamb {

struct VersionInfo {
  float version;
  char  date[9];
  char  time[9];
};

struct SectionInfo {
  long  offset;
  int   number;
  int   length;
  int   precision;
};

typedef std::map<std::string, SectionInfo> SectionIndex;

extern const SectionInfo info_a;
extern const SectionInfo info_e;
extern const SectionInfo info_f;
extern const SectionInfo info_i;

void BuildIndex(SectionIndex* index, VersionInfo* version, FILE* fp);

void ReadTitle(const VersionInfo& version, const SectionIndex& index,
               std::string* title, FILE* fp);

template <typename T>
void ReadData(const VersionInfo& version, const SectionIndex& index,
              const char* flag, int max_data, const char type,
              std::vector<T>* data, FILE* fp);

}  // namespace readamb
}  // namespace amberfile

#endif  // AMBERFILE_READAMBER_H
