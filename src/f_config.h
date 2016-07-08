#ifndef F_CONFIG_H
#define F_CONFIG_H

#include <libconfig.h++>
#include <unordered_map>
#include <vector>

#include "fasta.h"


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

  fasta::FastaData get_conf_data(const fasta::FastaData& fasta_data,
                          const FeatureSettingsMap& f_set, bool gapped);

  /// \brief add features based on the provided regular expression @pattern
  void assign_feature_by_pattern(fasta::SequenceList& sequences,
                                 const std::string& pattern,
                                 const std::string& feat_name);

  /// \brief find position in the alignment based on position in the sequence
  ///  e.g. for @sequence '--A' and @position 0 will return 2
  int find_real_pos(const std::string& sequence, int position);


  /// creates a list of all features present in @sequences + default ptm
  ///         features
  FeatureNamesList make_feature_list(const fasta::SequenceList& sequences);

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
