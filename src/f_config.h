#ifndef F_CONFIG_H
#define F_CONFIG_H
#include "types.h"

#include <iostream>
#include <libconfig.h++>


namespace f_config
{
  struct FeaturePositions
  {
    int seq;
    std::vector<int> positions;
  };

  struct FeatureSettings
  {
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

  typedef std::map<std::string, FeatureSettings> FeatureSettingsMap;

  class ConfParser
  {
    public:
      static FeatureSettingsMap parse_conf_file(const std::string& filename);
    private:
      static FeatureSettingsMap process_config(const libconfig::Config& cnfg);
  };
}

#endif /* F_CONFIG_H */
