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
      if (!(result.size() == 3 && result[0].substr(0, 2) == "m_")
          && !(result.size() == 2 && result[0].substr(0, 2) != "m_")){
        throw std::runtime_error("Invalid feature map format: " + line);
      }
      c[result[0]] = {result[1]};
      if (result[0].substr(0, 2) == "m_") {
        c[result[0]].push_back(result[2]);
      }
  }
  mapfile.close();
  return c;
}

feature_analysis::ConsensusSequence feature_analysis::analyze_alignment(
    feature_analysis::CodesMap codes_map,
    std::vector<fasta::SequenceList> alignment,
    double conservation_cutoff) {
  feature_analysis::ConsensusSequence cs;
  return cs;


}

void feature_analysis::write_consensus_to_file(
    feature_analysis::ConsensusSequence cons_seq,
    std::string out_cons_filename) {

}
