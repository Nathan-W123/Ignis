// SPDX-License-Identifier: MIT
#include "ignis/core/Version.hpp"

#include <sstream>

namespace ignis {

std::string versionBanner() {
  std::ostringstream os;
  os << "Ignis " << kVersion << " (" << kGitDescribe << ")\n"
     << "  build type : " << (kBuildType[0] ? kBuildType : "unspecified") << "\n"
     << "  compiler   : " << kCompilerId << " " << kCompilerVersion << "\n"
     << "  data dir   : " << kDefaultDataDir;
  return os.str();
}

}  // namespace ignis
