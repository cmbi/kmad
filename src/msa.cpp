#include "msa.h"

#include "features_profile.h"
#include "profile.h"
#include "scoring_matrix.h"
#include "substitution_matrix.h"
#include "vec_util.h"

#include <iostream>
#include <vector>
#include <tuple>


std::vector<std::string> msa::run_msa(fasta::FastaData fasta_data,
                             double gap_open_pen, double gap_ext_pen,
                             double end_pen, int domain_score, int motif_score,
                             int phosph_score, int codon_length) {

      Profile prf;
      FeaturesProfile fprf(domain_score, phosph_score, motif_score,
                           fasta_data.probabilities);


      // first round of the alignment - all vs 1st
      // TODO: first round is used to set identities. Return value is ignored.
      // Can this be improved?
      std::vector<double> identities;
      msa::PerformMSAfirstRound(fasta_data, prf, fprf, gap_open_pen, end_pen,
          gap_ext_pen, codon_length, identities);


      std::vector<std::string> alignment;
      //int prev_alignments = 0;
      //for (int i = 8; i >= 0; i--) {
        //double cutoff = double(i) / 10;
        //sequences.PerformMSAnextRound(alignment, prf, fprf, gap_open_pen,
                                      //end_pen, gap_ext_pen,
                                      //cutoff, codon_length,
                                      //identities, prev_alignments);
        ////prev_alignments - number of alignments performed in the previous
        ////rounds - to omit this round if the number of aligned sequences is the
        ////same as in the previous round
      //}
      //prev_alignments = 0;  // to align (again) all sequences to the profile
      //sequences.PerformMSAnextRound(alignment, prf, fprf,
                                    //gap_open_pen, end_pen, gap_ext_pen, 0,
                                    //codon_length, identities,
                                    //prev_alignments);
      return alignment;
}


std::vector<std::string> msa::PerformMSAfirstRound(
    fasta::FastaData fasta_data, Profile& output_profile,
    FeaturesProfile& output_features_profile, double gap_open_pen,
    double end_pen, double gap_ext_pen, int codon_length,
    IdentitiesList& identities)
{
  output_profile = Profile(substitution_matrix::ConvertToProfileFormat(
        fasta_data.sequences[0].residues));

  std::vector<fasta::Sequence> alignment_without_lowercase; //working alignment - without lowercase around cut out residues
  std::vector<fasta::Sequence> alignment_with_lowercase; //lowercase before and after cut out residues -- final result
  alignment_without_lowercase.push_back(fasta_data.sequences[0]);
  alignment_with_lowercase.push_back(fasta_data.sequences[0]);

  // identity of the 1st one to itself
  // to build the first profile based only on the first seqeunce
  identities.push_back(1);
  output_features_profile.ExpandListOfFeatures(fasta_data.sequences);

  //create features profile based on the 1st seq
  output_features_profile.CreateProfile(alignment_without_lowercase,
      codon_length);
  add_feature_indexes(output_features_profile);
  fasta::Sequence al_without_lower; //pairwise alignment without lowercase characters
  fasta::Sequence al_with_lower; //pairwise alignment with lowercase characters where chars were removed

  // TODO: Should this be a const loop?
  for (auto& sequence: fasta_data.sequences) {
    msa::AlignPairwise(
        al_without_lower, al_with_lower, sequence, output_profile,
        output_features_profile, gap_open_pen, end_pen, gap_ext_pen,
        codon_length);

    double identity = msa::CalcIdentity(al_without_lower);
    identities.push_back(identity);
  }
  return vec_util::Flatten(alignment_with_lowercase);
}


// TODO: Implement
double msa::CalcIdentity(const fasta::Sequence& aligned_sequence) {
  //double identical_residues = 0;
  //for (unsigned i = 0; i < aligned_sequence.size(); i++) {
    //if (aligned_sequence[i].get_aa() == m_sequences_aa[0][i].get_aa()) {
      //identical_residues++;
    //}
  //}
  //return identical_residues / double(m_first_sequence_size);
  return 0.0;
}


// TODO: Implement
// TODO: Return the result as normal
void msa::RemoveGaps(fasta::Sequence& alignment_with_lowercase,
                     fasta::Sequence& alignment_without_lowercase,
                     std::vector<fasta::Sequence>& alignment) {
  //fasta::Sequence s1 = alignment[0];
  //fasta::Sequence s2 = alignment[1];
  //fasta::Sequence new_s2;
  //fasta::Sequence new_s2_lower;
  //char gap = '-';
  //bool lower_flag = false;
  //for (unsigned int i = 0; i < alignment[0].size(); i++) {
    //char s1char = s1[i].get_aa();
    //if (s1char == gap) {
      //if (new_s2_lower.size() > 0) {
        ////change previous character to lowercase
        //new_s2_lower[new_s2_lower.size()-1].change_to_lowercase();
      //}
      //// flag to true so that the next character is also lowercase
      //lower_flag = true;
    //} else {
      //if (lower_flag) {   //lowercase char
        //Residue new_residue = s2[i];
        //new_residue.change_to_lowercase();
        ////add lowercase char to the alignment with lowercases
        //new_s2_lower.push_back(new_residue);
        ////add uppercase alignment to the alignment without lowercases
        //new_s2.push_back(s2[i]);
        //lower_flag = false;
      //} else {
        ////uppercase char
        //// adds the same uppercase char to both alignments (with lowercases and
        //// without lowercases)
        //new_s2_lower.push_back(s2[i]);
        //new_s2.push_back(s2[i]);
      //}
    //}
  //}
  //alignment_with_lowercase = new_s2_lower;
  //alignment_without_lowercase = new_s2;
}


// TODO: Implement
// TODO: Return a value as normal
void msa::AlignPairwise(fasta::Sequence& al_without_lower,
                        fasta::Sequence& al_with_lower,
                        fasta::Sequence& seq2,
                        Profile& prf,
                        FeaturesProfile& feat_prf,
                        double gap_open_pen, double end_pen,
                        double gap_ext_pen,
                        int codon_length) {
  //int profile_length = prf.get_matrix()[0].size();
  //std::vector<fasta::Sequence> alignment;
  //ScoringMatrix scores(profile_length, seq2.size(), gap_open_pen,
                       //end_pen, gap_ext_pen);
  //scores.CalculateScores(seq2, prf, feat_prf,
                         //codon_length);
  //scores.PerformNWAlignment(alignment, seq2, prf, feat_prf,
                            //codon_length);

  //msa::RemoveGaps(al_with_lower, al_without_lower, alignment);
}


// TODO: Implement
void msa::PerformMSAnextRound(
    std::vector<std::string>& prev_alignment,
    Profile& output_profile,
    FeaturesProfile& output_features_profile,
    double gap_open_pen,
    double end_pen,
    double gap_ext_pen,
    double identity_cutoff,
    int codon_length,
    IdentitiesList& identities,
    int& prev_alignments)
{
  //int next_alignments = CountAlignments(identity_cutoff, identities);
  //if (next_alignments > prev_alignments) {
    ////working alignment - without lowercase around cut out residues
    ////would make latter aligning more complicated
    //std::vector<fasta::Sequence> alignment_without_lowercase;
    ////lowercase before and after cut out residues -- final result
    //std::vector<fasta::Sequence> alignment_with_lowercase;
    //alignment_without_lowercase.push_back(m_sequences_aa[0]);
    //alignment_with_lowercase.push_back(m_sequences_aa[0]);
    //// tmp pairwise alignment (and so is al_with_lower)
    //fasta::Sequence al_without_lower;
    //fasta::Sequence al_with_lower;
    //for (int i = 1; i < m_seq_nr; i++) {
      //if (identities[i] > identity_cutoff) {
        //// NW alignment of the ith seq against the profile
        //msa::AlignPairwise(al_without_lower, al_with_lower, m_sequences_aa[i],
                      //output_profile, output_features_profile, gap_open_pen,
                      //end_pen, gap_ext_pen, codon_length);
        //alignment_without_lowercase.push_back(al_without_lower);
        //alignment_with_lowercase.push_back(al_with_lower);
      //}
    //}
    ////create features profile based on the 1st seq
    //output_features_profile.CreateProfile(alignment_without_lowercase,
                                          //codon_length);
    //output_profile.ProcessProfile(alignment_without_lowercase);
    //prev_alignment = vec_util::Flatten(alignment_with_lowercase);
    ////update number of performed alignments
    //prev_alignments = next_alignments;
  //}
}


// TODO: Implement. Can probably be moved somewhere else as well.
void msa::add_feature_indexes(FeaturesProfile& fprf) {
  //std::string nothing = "AA";
  //for (auto &seq : m_sequences_aa) {
    //for (auto &res : seq) {
        //FeatureNamesList features = res.get_features();
        //FeaturesList indexes;
        //for (auto &feat : features) {
          //if (feat != nothing) {
            //indexes.push_back(fprf.FindFeaturesIndex(feat));
          //}
        //}
        //res.set_feat_indexes(indexes);
    //}
  //}
}


int msa::CountAlignments(double identity_cutoff,
                         IdentitiesList& identities) {
  int count = 0;
  for (auto& item: identities) {
    if (item > identity_cutoff) {
      count++;
    }
  }
  return count;
}
