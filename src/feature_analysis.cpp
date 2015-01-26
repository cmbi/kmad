#include "feature_analysis.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <fstream>

namespace fs = boost::filesystem;
feature_analysis::CodesMap feature_analysis::parse_mapfile(
    std::string filename) {
  feature_analysis::CodesMap c;
  fs::path p(filename);
  if (!fs::exists(p)) {
    throw std::invalid_argument("File not found: " + filename);
  }
  std::ifstream mapfile(filename.c_str());
  std::string line;
  while (std::getline(mapfile, line)) {
      std::vector<std::string> result;
      boost::split(result, line, boost::is_any_of("\t "));
      if (result.size() != 2) {
        throw std::runtime_error("Invalid feature map format: " + line);
      }
      c[result[0]] = result[1];
  }
  mapfile.close();
  return c;
}
