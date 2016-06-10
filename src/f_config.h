#ifndef F_CONFIG_H
#define F_CONFIG_H

#include <libconfig.h++>
#include <unordered_map>
#include <vector>

typedef std::vector<std::string> FeatureNamesList;
namespace f_config
{
  struct FeaturePositions
  {
    int seq_no;
    std::vector<int> positions;
  };
  ///
  /// Struct to hold the feature settings(for a single feature) in the same
  /// structure as in the
  /// config file (get parsed to FeatureSettings)
  ///
  struct RawFeatureSettings {
    std::string tag;
    int add_score;
    int subtract_score;
    std::string pattern;
    FeatureNamesList add_features;
    FeatureNamesList add_tags;
    FeatureNamesList add_exceptions;
    FeatureNamesList subtract_features;
    FeatureNamesList subtract_tags;
    FeatureNamesList subtract_exceptions;
    std::vector<FeaturePositions> positions;
  };
  ///
  /// Processed feature settings (from RawFeatureSettings)
  /// *_features, *_tags, and *_exceptions elements from raw... are now
  /// combined in one *_features
  ///
  struct FeatureSettings {
    int add_score;
    int subtract_score;
    std::string pattern;
    FeatureNamesList add_features;
    FeatureNamesList subtract_features;
    std::vector<FeaturePositions> positions;
  };

  typedef std::unordered_map<std::string, RawFeatureSettings> RawFeatureSettingsMap;
  typedef std::unordered_map<std::string, FeatureSettings> FeatureSettingsMap;

  class ConfParser
  {
    public:
      ///
      /// take a configuration file name and returns processed map of feature
      /// settings - calls process_config(...) and process_settings(...)
      ///
      static FeatureSettingsMap parse_conf_file(const std::string& filename);
    private:
      static RawFeatureSettingsMap process_config(
          const libconfig::Config& cnfg);
      ///
      /// converts a map of RawFeatureSettings structs to FeatureSettings
      /// by combining the different add_* and subtract_* elements
      ///
      static FeatureSettingsMap process_settings(
          const RawFeatureSettingsMap& raw_feat_set);
  };
}

#endif /* F_CONFIG_H */
