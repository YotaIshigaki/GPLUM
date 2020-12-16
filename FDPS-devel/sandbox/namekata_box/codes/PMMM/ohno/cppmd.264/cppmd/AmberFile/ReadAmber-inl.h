#ifndef AMBERFILE_READAMBER_INL_H
#define AMBERFILE_READAMBER_INL_H

#include "ReadAmber.h"

#include <sstream>
#include "getline.h"

namespace amberfile {
namespace readamb {

template <typename T>
  void ReadData(const VersionInfo& version, const SectionIndex& index,
                const char* flag, int max_data, char type,
                std::vector<T>* data, FILE* fp) {
    using namespace std;
    data->clear();

    if (max_data > 0) {
      data->reserve(max_data);

      const SectionInfo* info = &info_a;
      if (version.version <= 0.0) {
        switch (type) {
          case 'A': info = &info_a; break;
          case 'E': info = &info_e; break;
          case 'F': info = &info_f; break;
          case 'I': info = &info_i; break;
          default: assert(false);
        }
      } else {
        SectionIndex::const_iterator idx = index.find(string(flag));
        if (idx == index.end()) return;
        info = &idx->second;
        fseek(fp, idx->second.offset, SEEK_SET);
      }

      int count = 0;
      char *buf = new char[info->length + 1];
      char* line = NULL;
      size_t len = 0;
      for (ssize_t read;
           count < max_data && (read = getline(&line, &len, fp)) != -1; ) {
        if (*line == '%') break;
        if (line[read - 1] == '\n') --read;
        for (char* p = line; p + info->length <= line + read;
             p += info->length) {
          strncpy(buf, p, info->length);
          buf[info->length] = '\0';
          T datum;
          stringstream ss; ss << buf; ss >> datum;
          data->push_back(datum);
          if (++count >= max_data) break;
        }
      }
      if (line) free(line);
      delete []buf;
    }
    assert(data->size() <= static_cast<size_t>(max_data));
  }

}  // namespace readamb
}  // namespace amberfile

#endif  // AMBERFILE_READAMBER_INL_H
