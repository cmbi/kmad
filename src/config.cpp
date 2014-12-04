#include "config.h"
#include "types.h"

#include <iostream>
#include <libconfig.h++>

namespace lcg = libconfig;

config::UsrFeatureMap config::parse_conf_file(const char* filename) {
  lcg::Config cnfg;
  try
  {
    cnfg.readFile(filename);
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return process_config(cnfg);
}


UsrFeatureMap config::process_config(const Config& cnfg) {
  config::UsrFeatureMap feat_config;
  try
  {
    const Setting &root  = cnfg.getRoot();
    const Setting &features = root["feature_settings"]["usr_features"];
    int count = features.getLength();
    for (unsigned i = 0; i < count; ++i) {
      FeatureSettings feat_set;
      const Setting& feature = features[i];

      std::string name, tag;
      double add_score, subtract_score;
      FeatureNamesList add_features, add_tags, add_exceptions;
      FeatureNamesList subtract_features, subtract_tags, subtract_exceptions;
      config::FeaturePositions positions; 

      bool found_add = feature.lookupValue("add_score", 
                                           feat_set.add_score)
                       && (feature.lookupValue("add_features", 
                                               feat_set.add_features)
                           || feature.lookupValue("add_tags", 
                                                  feat_set.add_tags));

      bool found_subtract = feature.lookupValue("subtract_score", 
                                                 feat_set.subtract_score)
                            && (feature.lookupValue("subtract_features", 
                                                    feat_set.subtract_features)
                               || feature.lookupValue("subtract_tags", 
                                                      feat_set.subtract_tags));

      if (!(feature.lookupValue("name", name) 
            && feature.lookupValue("positions.seq", feat_set.positions.seq)
            && feature.lookupValue("positions.pos", feat_set.positions.pos)
            && (found_add || found_subtract))) {
        continue;
      }
      feat_config[name] = feat_set;
    }
  }
  catch(const SettingNotFoundException &nfex)
  {
    // Ignore.
  }
}




