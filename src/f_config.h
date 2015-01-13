#ifndef F_CONFIG_H
#define F_CONFIG_H
#include "types.h"

#include <iostream>
#include <libconfig.h++>


namespace f_config
{
  struct FeaturePositions
  {
    int seq_no;
    std::vector<int> positions;
  };

  struct RawFeatureSettings {
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

  struct FeatureSettings {
    std::string tag;
    int add_score;
    int subtract_score;
    FeatureNamesList add_features;
    FeatureNamesList subtract_features;
    std::vector<FeaturePositions> positions;
  };

  typedef std::map<std::string, RawFeatureSettings> RawFeatureSettingsMap;
  typedef std::map<std::string, FeatureSettings> FeatureSettingsMap;

  class ConfParser
  {
    public:
      static FeatureSettingsMap parse_conf_file(
          const std::string& filename);
    private:
      static RawFeatureSettingsMap process_config(
          const libconfig::Config& cnfg);
      static FeatureSettingsMap process_settings(
          const RawFeatureSettingsMap raw_feat_set);
  };
}

#endif /* F_CONFIG_H */
