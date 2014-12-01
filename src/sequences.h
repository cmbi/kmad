#ifndef SEQUENCES_H
#define SEQUENCES_H


#include "types.h"
#include <iostream>
#include <vector>


class Residue;
class Profile;
class FeaturesProfile;
class Sequences{
public:
  ///
  /// constructor; creates a list of residues out of a list of sequences with
  /// names
  ///
  Sequences(CodonSeqWithNamesList& s);
  ///
  /// constructor; creates an empty object Sequences
  /// 
  Sequences();
  ///
  /// returns a list of sequence names (fasta headers)
  ///
  SeqNames get_names();
  ///
  /// performs the first round of alignments, 
  /// all vs query seq (first calculates profile 
  /// based only on the query seq, then 
  /// aligns all sequences and calculates 
  /// identity of each sequence to the query seq.)
  ///
	StringSequences PerformMSAfirstRound(Profile& output_profile, 
                                        FeaturesProfile& output_features_profile, 
                                        double penalty, double endPenalty, 
                                        double extensionPenalty, 
                                        int codon_length, 
                                        IdentitiesList& identities);
  ///
  /// performs next round of MSA (good for all rounds except for the first one
  /// - you need a profile)
  ///
	void PerformMSAnextRound(StringSequences& prevAlignment,
                           Profile& output_profile, 
                           FeaturesProfile& output_features_profile, 
                           double penalty, double endPenalty, 
                           double extensionPenalty, 
                           double identityCutoff, 
                           int codon_length, IdentitiesList& identities, 
                           int& prev_alignments);
  ///
  /// adds features from the tuple 'feature_rules'(usr defined) to relevant 
  /// residues (also specified in 'feature_rules')
  ///
  void add_usr_features(RuleTuplesList& feature_rules);
private:
  ///
  /// takes pairwise alignment, removes 
  /// characters from the 2nd sequence that match gaps from 1st seq and returns 
  /// vector<string> of 2 elements, where the 1st one is 2nd sequence with cut 
  /// out chars and 2nd one
  /// is 2nd sequence with cut out chars and lowercase chars 
  /// before and after that
  ///
	void RemoveGaps(ResidueSequence& alignment_with_lowercase, 
                  ResidueSequence& alignment_without_lowercase, 
                  SequenceList& alignment);
  ///
  /// takes a sequence and profiles, returns an
  /// alignment of the two, with gaps cut out
  ///
	void AlignPairwise(ResidueSequence& al_without_lower, 
                     ResidueSequence& al_with_lower, 
                     ResidueSequence& seq2,
                     Profile& prf, FeaturesProfile& featPrf,
                     double penalty, double endPenalty, double extensionPenalty,
                     int codon_length);
  ///
  /// calculates identity with the query sequence 
  /// @param alignedSequence sequence aligned to the profile with the gaps cut
  /// out (its length is equal to the profile's length) 
  ///
	double CalcIdentity(const ResidueSequence& alignedSequence);
  ///
  /// count alignments that will be performed in this round
  /// 
	int CountAlignments(double identity_cutoff, IdentitiesList& identities);
  void add_feature_indexes(FeaturesProfile& fprf);
  ///
  ///m_seq_nr - number of sequences
  ///
	int m_seq_nr;
	int m_first_sequence_size;
	SequenceList m_sequences_aa;
	SeqNames m_sequence_names;
};

#endif /* SEQUENCES_H */
