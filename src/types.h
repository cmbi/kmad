#include "Residue.h"
#include <iostream>
#include <vector>
#include <tuple>

typedef std::vector<std::string> string_sequences;
typedef std::vector<std::string> ids_list;
typedef std::vector<double> probs_list;
typedef std::vector<int> featuresList;
typedef std::vector<std::string> featureNamesList;
typedef std::vector<Residue> sequence;
typedef std::vector<sequence> sequenceList;
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
typedef std::vector<int> indexList;
