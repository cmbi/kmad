#include "Residue.h"
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
                   double, double, std::string, std::string> rulesTuple;
typedef std::vector<rulesTuple> rulesTuplesList;
typedef	std::tuple<std::string, int, double, double, double,
                   double, std::vector<int>,
                   std::vector<int> > processedRules;
typedef std::vector<processedRules> prcRulesList;
typedef std::vector<double> identitiesList;
typedef std::vector<double> profileMatrixRow;
typedef std::vector<double> profileMatrixColumn;
typedef std::vector<profileMatrixRow> profile_matrix;
typedef std::vector<std::vector<double>> matrix2d;
typedef std::vector<std::vector<int>> sbst_matrix;
typedef std::vector<std::vector<int>> sbst_matrix_columns;
typedef std::vector<int> indexList;
typedef std::vector<int> sbstMatColumn;
typedef std::tuple<std::string, std::string, int, int, int> defaultRuleTuple;
typedef std::vector<defaultRuleTuple> defaultRulesList;
typedef std::vector<char> alphabetVec;
