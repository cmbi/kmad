#include <iostream>
#include <map>

namespace feature_analysis {
  typedef std::map<std::string, std::string> CodesMap;

  CodesMap parse_mapfile(std::string filename);
}
