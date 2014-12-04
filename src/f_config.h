#include "types.h"
#include <iostream>
#include <libconfig.h++>

namespace f_config {
  struct FeaturePositions {
    int seq;
    std::vector<int> positions;
  };
  struct FeatureSettings {
    std::string tag;
    int add_score;
    int subtract_score;
    FeatureNamesList add_features;
    FeatureNamesList add_tags;
    FeatureNamesList add_exceptions;
    FeatureNamesList subtract_features;
    FeatureNamesList subtract_tags;
    FeatureNamesList subtract_exceptions;
    std::vector<FeaturePositions> positions;
  };
  typedef std::map<std::string, FeatureSettings> UsrFeatureMap;
  class ConfParser {
    public:
      static UsrFeatureMap parse_conf_file(const char* filename);
    private:
      static UsrFeatureMap process_config(const libconfig::Config& cnfg);
  };

}
