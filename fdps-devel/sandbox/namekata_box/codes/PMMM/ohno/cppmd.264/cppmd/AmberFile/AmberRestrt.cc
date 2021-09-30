#include "AmberRestrt.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include "ReadAmber-inl.h"

namespace {

void ReadNatomTime(const amberfile::readamb::VersionInfo& /*version*/,
                   const amberfile::readamb::SectionIndex& /*index*/,
                   amberfile::AmberRestrt* rst, FILE* fp) {
  char* line = NULL;
  size_t len = 0;
  ssize_t read = getline(&line, &len, fp);
  if (read != -1) {
#ifdef CPPMD_ENABLE_LARGEMODEL
# if SIZEOF_LONG == 8
    sscanf(line, "%ld %lf", &rst->natom, &rst->current_time);
# else  // SIZEOF_LONG == 8
    sscanf(line, "%ld %llf", &rst->natom, &rst->current_time);
# endif  // SIZEOF_LONG == 8
#else  // CPPMD_ENABLE_LARGEMODEL
    sscanf(line, "%d %lf", &rst->natom, &rst->current_time);
#endif  // CPPMD_ENABLE_LARGEMODEL
  }
  if (line) free(line);
}

}  // namespace

namespace amberfile {

int ReadFromFile(const std::string& file_name,
                 AmberRestrt* rst,
                 int irest) {
  FILE* fp = fopen(file_name.c_str(), "r");
  if (fp == NULL) return 1;
  int ret = ReadFromFile(fp, rst, irest);
  fclose(fp);
  return ret;
}

int ReadFromFile(FILE* fp, AmberRestrt* rst, int irest) {
  using namespace std;

  readamb::VersionInfo version;
  readamb::SectionIndex index;

  readamb::BuildIndex(&index, &version, fp);

  rst->version = version.version;
  rst->date = version.date;
  rst->time = version.time;
  rst->irest = irest;

  readamb::ReadTitle(version, index, &rst->title, fp);
  ReadNatomTime(version, index, rst, fp);

# define READ(f,n,v,t) do { \
  readamb::ReadData(version, index, (f), (n), (t), &rst->v, fp); } while (0)

  READ("COORDINATES", rst->natom * 3, coordinates, 'F');
  READ("VELOCITIES", rst->natom * 3, velocities, 'F');
  READ("BOX_DIMENSIONS", 6, box, 'F');
  if (irest == 0) {
    if (rst->velocities.size() >= 6 &&
        rst->velocities.size()
        < static_cast<std::vector<double>::size_type>(rst->natom * 3)) {
      rst->box.resize(6);
      for (int i = 0; i < 6; ++i) {
        rst->box[i] = rst->velocities[i];
      }
    }
    rst->velocities.clear();
    rst->velocities.resize(rst->natom * 3);
  }
  if (fabs(rst->box[3]) < numeric_limits<double>::epsilon()
      && fabs(rst->box[4]) < numeric_limits<double>::epsilon()
      && fabs(rst->box[5]) < numeric_limits<double>::epsilon()) {
    rst->box[3] = rst->box[4] = rst->box[5] = 90.0;
  }

#undef READ

  // returns 1 when error: rst->coordinates.size() != rst->velocities.size()
  return rst->coordinates.size() != rst->velocities.size();
}

int WriteToFile(const AmberRestrt& rst, const std::string& file_name, int ntb) {
  if (file_name.empty()) return 0;
  FILE* fp = fopen(file_name.c_str(), "w");
  if (fp == NULL) return 1;
  int ret = WriteToFile(rst, fp, ntb);
  fclose(fp);
  return ret;
}

int WriteToFile(const AmberRestrt& rst, FILE* fp, int ntb) {
#ifdef CPPMD_ENABLE_LARGEMODEL
# if SIZEOF_LONG == 8
  fprintf(fp, "%s\n% ld% 15.7E",
          rst.title.c_str(), rst.natom, rst.current_time);
# else  // SIZEOF_LONG == 8
  fprintf(fp, "%s\n% lld% 15.7E",
          rst.title.c_str(), rst.natom, rst.current_time);
# endif  // SIZEOF_LONG == 8
#else  // CPPMD_ENABLE_LARGEMODEL
  fprintf(fp, "%s\n% 5d% 15.7E",
          rst.title.c_str(), rst.natom, rst.current_time);
#endif  // CPPMD_ENABLE_LARGEMODEL
  for (int i = 0; i < rst.natom * 3; ++i)
    fprintf(fp, "%s%12.7f", ((i % 6) ? "" : "\n"), rst.coordinates[i]);
  for (int i = 0; i < rst.natom * 3; ++i)
    fprintf(fp, "%s%12.7f", ((i % 6) ? "" : "\n"), rst.velocities[i]);
  if (ntb > 0) {
    for (int i = 0; i < 6; ++i)
      fprintf(fp, "%s%12.7f", ((i % 6) ? "" : "\n"), rst.box[i]);
  }
  fputc('\n', fp);
  return 0;
}

FILE* PrepareTitledFile(const AmberRestrt& rst, const std::string& file_name) {
  FILE* fp = fopen(file_name.c_str(), "w");
  if (fp) fprintf(fp, "%s", rst.title.c_str());
  return fp;
}

void CloseTitledFile(FILE* fp) {
  if (fp) {
    fputc('\n', fp);
    fclose(fp), fp = 0;
  }
}

int AppendCoordinatesToFile(const AmberRestrt& rst, FILE* fp, int ntb,
                            int iwrap) {
  if (fp == NULL) return 0;
  int i = 0;
  for (; i < rst.natom * 3; ++i) {
    double c = rst.coordinates[i];
    if (iwrap) {
      c = fmod(c, rst.box[i % 3]);
      if (c < 0.0) c += rst.box[i % 3];
    }
    fprintf(fp, "%s%8.3f", ((i % 10) ? "" : "\n"), c);
  }
  if (ntb > 0) {
    fputc('\n', fp);
    for (int i = 0; i < 3; ++i) fprintf(fp, "%8.3f", rst.box[i]);
    //for (int i = 3; i < 6; ++i) fprintf(fp, "%8.3f", rst.box[i]);
  }
  return 0;
}

int AppendVelocitiesToFile(const AmberRestrt& rst, FILE* fp, int ntb) {
  if (fp == NULL) return 0;
  for (int i = 0; i < rst.natom * 3; ++i)
    fprintf(fp, "%s%8.3f", ((i % 10) ? "" : "\n"), rst.velocities[i]);
  if (ntb > 0) {
    fputc('\n', fp);
    for (int i = 0; i < 3; ++i) fprintf(fp, "%8.3f", rst.box[i]);
    //for (int i = 3; i < 6; ++i) fprintf(fp, "%8.3f", rst.box[i]);
  }
  return 0;
}

}  // namespace amberfile
