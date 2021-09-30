#ifndef ERRORPOS_H
#define ERRORPOS_H

#include <string>
#include <sstream>

#ifndef NDEBUG
#ifdef __FCC_VERSION
#define errorPos(x) makeErrorPos(x,__FILE__,__LINE__,"FCC __FUNCTION__ not work")
#else
#define errorPos(x) makeErrorPos(x,__FILE__,__LINE__,__FUNCTION__)
#endif
#else
#define errorPos(x) (x)
#endif

namespace {
const std::string makeErrorPos(const std::string s, char const *file, int line, char const* func)
{
  std::ostringstream os;
  os << s << ": " << file << ":" << line << ":" << func;
  return os.str();
}
}

#endif
