#include "ReadAmber.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <sstream>
#include <string>

#include "getline.h"

namespace amberfile {
namespace readamb {

const SectionInfo info_a = { -1, 20,  4, 0 };
const SectionInfo info_e = { -1,  5, 16, 8 };
const SectionInfo info_f = { -1,  6, 12, 7 };
const SectionInfo info_i = { -1, 12,  6, 0 };

void BuildIndex(SectionIndex* index, VersionInfo* version, FILE* fp) {
  using namespace std;

  int line_num = 0;
  long offset = 0;
  string flag;
  SectionInfo section_info;
  enum ReadState { kStart, kFlag, kFormat } state = kStart;

  version->version = 0000.000;
  strncpy(version->date, "00/00/00", 9), version->date[8] = '\0';
  strncpy(version->time, "00:00:00", 9), version->time[8] = '\0';

  char* line = NULL;
  size_t len = 0;
  for (ssize_t read; (read = getline(&line, &len, fp)) != -1; offset += read) {
    ++line_num;
    switch (state) {
      case kStart:
        if (strncmp(line, "%VERSION ", 9) == 0) {
          assert(version->version <= 0.0);
#ifndef NDEBUG
          int match =
#endif  // NDEBUG
              sscanf(line,"%%VERSION VERSION_STAMP = V%f DATE = %s %s",
                     &version->version, version->date, version->time);
          assert(match > 0);
        } else if (strncmp(line, "%FLAG ", 6) == 0) {
          char *buf = new char[len];
#ifndef NDEBUG
          int match =
#endif  // NDEBUG
              sscanf(line + 6, "%s", buf);
          assert(match == 1);
          flag = buf;
          state = kFlag;
          delete []buf;
        }
        break;
      case kFlag:
        if (strncmp(line, "%FORMAT", 7) == 0) {
          section_info.offset = section_info.number = 0;
          section_info.length = section_info.precision = 0;
#ifndef NDEBUG
          int match =
#endif  // NDEBUG
              sscanf(line + 7, "(%d%*c%d.%d)", &section_info.number,
                     &section_info.length, &section_info.precision);
          assert(match >= 2);
          state = kFormat;
        } else if (read > 0 && *line == '%') {
        } else {
          assert(false);
        }
        break;
      case kFormat:
        if (read > 0 && *line == '%') {
          assert(false);
        } else {
          section_info.offset = offset;
          assert(index->find(flag) == index->end());
          (*index)[flag] = section_info;
          state = kStart;
        }
        break;
      default:
        assert(false);
    }
  }
  if (line) free(line);
}

void ReadTitle(const VersionInfo& version, const SectionIndex& index,
               std::string* title, FILE* fp) {
  using namespace std;
  if (version.version <= 0.0) {
    rewind(fp);
  } else {
    SectionIndex::const_iterator idx = index.find(string("TITLE"));
    if (idx == index.end()) return;
    const SectionInfo& info = idx->second;
    fseek(fp, info.offset, SEEK_SET);
  }
  char* line = NULL;
  size_t len = 0;
  ssize_t read = getline(&line, &len, fp);
  if (read != -1) {
    if (line[read - 1] == '\n') line[read - 1] = '\0';
    *title = line;
  }
  if (line) free(line);
}

}  // namespace readamb
}  // namespace amberfile
