#include "sequences.h"
#include "residue.h"
#include "features_profile.h"
#include "profile.h"
#include "scoring_matrix.h"
#include "substitution_matrix.h"
#include "vec_util.h"

#include <iostream>
#include <tuple>
#include <string>
#include <vector>
#include <ctime>


Sequences::Sequences(CodonSeqWithNamesList &s) {
  for (unsigned int i = 0; i < s.size(); i++) {
    m_sequence_names.push_back(s[i][0][0]);
    ResidueSequence new_seq;
    for (unsigned int j = 0; j < s[i][1].size(); j++) {
      Residue newRes(s[i][1][j]);
      new_seq.push_back(newRes);
    }
    m_sequences_aa.push_back(new_seq);
  }
  m_seq_nr = s.size(); 
  m_first_sequence_size = s[0][1].size();
}


Sequences::Sequences() {
}


StringSequences Sequences::PerformMSAfirstRound(Profile& output_profile, 
                                                FeaturesProfile& output_features_profile, 
                                                double gap_open_pen, 
                                                double end_pen, 
                                                double gap_ext_pen, 
                                                int codon_length, 
                                                IdentitiesList& identities) {
  output_profile = Profile(substitution_matrix::ConvertToProfileFormat(m_sequences_aa[0])); 
  //working alignment - without lowercase around cut out residues
  SequenceList alignment_without_lowercase;  
  //lowercase before and after cut out residues -- final result 
  SequenceList alignment_with_lowercase;    
  alignment_without_lowercase.push_back(m_sequences_aa[0]);
  alignment_with_lowercase.push_back(m_sequences_aa[0]);
  // identity of the 1st one to itself
  //to build the first profile based only on the first seqeunce
  identities.push_back(1); 
  output_features_profile.ExpandListOfFeatures(m_sequences_aa);
  //create features profile based on the 1st seq
  output_features_profile.CreateProfile(alignment_without_lowercase, codon_length);   
  add_feature_indexes(output_features_profile);
  //pairwise alignment without lowercase characters
  ResidueSequence al_without_lower; 
  //pairwise alignment with lowercase characters where chars were removed
  ResidueSequence al_with_lower; 
  for (auto &seqI: m_sequences_aa) {
    AlignPairwise(al_without_lower, al_with_lower, seqI, output_profile, 
                  output_features_profile, gap_open_pen, end_pen, gap_ext_pen, 
                  codon_length);
    double identity = CalcIdentity(al_without_lower);
    identities.push_back(identity);
  }
  return vec_util::Flatten(alignment_with_lowercase);    
}


void Sequences::PerformMSAnextRound(StringSequences& prev_alignment, 
                                    Profile& output_profile,
                                    FeaturesProfile& output_features_profile, 
                                    double gap_open_pen, 
                                    double end_pen, 
                                    double gap_ext_pen,
                                    double identity_cutoff,
                                    int codon_length, 
                                    IdentitiesList& identities, 
                                    int& prev_alignments) {
  int next_alignments = CountAlignments(identity_cutoff, identities);
  if (next_alignments > prev_alignments) {
    //working alignment - without lowercase around cut out residues
    //would make latter aligning more complicated
    SequenceList alignment_without_lowercase;  
    //lowercase before and after cut out residues -- final result 
    SequenceList alignment_with_lowercase;    
    alignment_without_lowercase.push_back(m_sequences_aa[0]);
    alignment_with_lowercase.push_back(m_sequences_aa[0]);
    // tmp pairwise alignment (and so is al_with_lower)
    ResidueSequence al_without_lower; 
    ResidueSequence al_with_lower;
    for (int i = 1; i < m_seq_nr; i++) {
      if (identities[i] > identity_cutoff) {
        // NW alignment of the ith seq against the profile
        AlignPairwise(al_without_lower, al_with_lower, m_sequences_aa[i], output_profile,
                      output_features_profile, gap_open_pen, end_pen,
                      gap_ext_pen, codon_length); 
        alignment_without_lowercase.push_back(al_without_lower);
        alignment_with_lowercase.push_back(al_with_lower);
      }
    }
    //create features profile based on the 1st seq
    output_features_profile.CreateProfile(alignment_without_lowercase, 
                                          codon_length);
    output_profile.ProcessProfile(alignment_without_lowercase);
    prev_alignment = vec_util::Flatten(alignment_with_lowercase);
    //update number of performed alignments
    prev_alignments = next_alignments; 
  }
}


double Sequences::CalcIdentity(const ResidueSequence& aligned_sequence) {
  double identical_residues = 0;
  for (unsigned int i = 0; i < aligned_sequence.size(); i++) {
    if (aligned_sequence[i].get_aa() == m_sequences_aa[0][i].get_aa()) {
      identical_residues++;
    }
  }
  return identical_residues/double(m_first_sequence_size);
}


void Sequences::RemoveGaps(ResidueSequence& alignment_with_lowercase, 
                           ResidueSequence& alignment_without_lowercase, 
                           SequenceList& alignment) {
  ResidueSequence s1 = alignment[0];
  ResidueSequence s2 = alignment[1];
  ResidueSequence new_s2;
  ResidueSequence new_s2_lower;
  char gap = '-';
  bool lower_flag = false;
  for (unsigned int i = 0; i < alignment[0].size(); i++) {
    char s1char = s1[i].get_aa();
    if (s1char == gap) {
      if (new_s2_lower.size() > 0) {
        //change previous character to lowercase
        new_s2_lower[new_s2_lower.size()-1].change_to_lowercase(); 
      }
      // flag to true so that the next character is also lowercase
      lower_flag = true; 
    } else {
      if (lower_flag) {   //lowercase char
        Residue new_residue = s2[i];
        new_residue.change_to_lowercase();
        //add lowercase char to the alignment with lowercases
        new_s2_lower.push_back(new_residue); 
        //add uppercase alignment to the alignment without lowercases
        new_s2.push_back(s2[i]);        
        lower_flag = false;
      } else {     
        //uppercase char
        // adds the same uppercase char to both alignments (with lowercases and 
        // without lowercases)
        new_s2_lower.push_back(s2[i]); 
        new_s2.push_back(s2[i]);
      }
    }
  }
  alignment_with_lowercase = new_s2_lower;
  alignment_without_lowercase = new_s2;
}


void Sequences::AlignPairwise(ResidueSequence& al_without_lower, 
                              ResidueSequence& al_with_lower, 
                              ResidueSequence& seq2, 
                              Profile& prf, 
                              FeaturesProfile& feat_prf, 
                              double gap_open_pen, double end_pen, 
                              double gap_ext_pen, 
                              int codon_length) {
  int profile_length = prf.get_matrix()[0].size();
  SequenceList alignment;
  ScoringMatrix scores(profile_length, seq2.size(), gap_open_pen, 
                       end_pen, gap_ext_pen);
  scores.CalculateScores(seq2, prf, feat_prf, 
                         codon_length);
  scores.PerformNWAlignment(&alignment, seq2, prf, feat_prf, 
                            codon_length);

  RemoveGaps(al_with_lower, al_without_lower, alignment); 
}


int Sequences::CountAlignments(double identity_cutoff, 
                               IdentitiesList& identities) {
  int count = 0;
  for (auto &item: identities) {
    if (item > identity_cutoff) {
      count++;
    }
  }
  return count;
}


void Sequences::add_usr_features(RuleTuplesList& feature_rules) {
  for (auto &rule: feature_rules) {
    std::string feat_name = std::string("USR_")
                            + std::get<0>(rule)
                            + std::string("_")
                            + std::get<1>(rule);
    int sequence_no = std::get<2>(rule);
    int start = std::get<3>(rule);
    int end = std::get<4>(rule)+1;
    signed int seq_length = m_sequences_aa[sequence_no].size();
    if (sequence_no < (signed)m_sequences_aa.size()) {
      for (int j = start; j < end && j < seq_length; j++) {
        m_sequences_aa[sequence_no][j].add_feature(feat_name);
      }
    }
  }
}


SeqNames Sequences::get_names() {
  return m_sequence_names;
}


void Sequences::add_feature_indexes(FeaturesProfile& fprf) {
  std::string nothing = "AA";
  for (auto &seq: m_sequences_aa) {
    for (auto &res: seq) {
        FeatureNamesList features = res.get_features();
        FeaturesList indexes;
        for (auto &feat: features) {
          if (feat != nothing) {
            indexes.push_back(fprf.FindFeaturesIndex(feat));
          }
        }
        res.set_feat_indexes(indexes);
    }
  }
}
