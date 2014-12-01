#include "residue.h"
#include <iostream>
#include <vector>
#include <tuple>

class Residue;


typedef std::vector<std::string> StringSequences;
typedef std::vector<std::string> SeqNames;
typedef std::vector<std::string> CodonSeq;
typedef std::vector<CodonSeq> CodonSeqList;
typedef std::vector<CodonSeq> CodonSeqWithName;
typedef std::vector<CodonSeqList> CodonSeqWithNamesList;
typedef std::vector<std::string> IDsList;
typedef std::vector<double> ProbsList;
typedef std::vector<int> FeaturesList;
typedef std::vector<std::string> FeatureNamesList;
typedef std::vector<Residue> ResidueSequence;
typedef std::vector<ResidueSequence> SequenceList;
typedef std::tuple<std::string, std::string, int, int, int, double, double,
                   double, double, std::string, std::string> RulesTuple;
typedef std::vector<RulesTuple> RuleTuplesList;
typedef	std::tuple<std::string, int, double, double, double,
                   double, std::vector<int>,
                   std::vector<int> > ProcessedRules;
typedef std::vector<ProcessedRules> PrcRulesList;
typedef std::vector<double> IdentitiesList;
typedef std::vector<double> ProfileMatrixRow;
typedef std::vector<double> ProfileMatrixColumn;
typedef std::vector<ProfileMatrixRow> ProfileMatrix;
typedef std::vector<std::vector<double>> Matrix2D;
typedef std::vector<std::vector<int>> SbstMatrix;
typedef std::vector<int> SbstMatColumn;
typedef std::vector<SbstMatColumn> SbstMatrixColumns;
typedef std::vector<int> IndexList;
typedef std::tuple<std::string, std::string, int, int, int> DefaultRuleTuple;
typedef std::vector<DefaultRuleTuple> DefaultRulesList;
typedef std::vector<char> AlphabetVec;
