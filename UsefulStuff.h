#include <iostream>
#include <vector>
class UsefulStuff{
	public:
		//text processing
		int convertStringToInt(std::string);
		std::vector< std::vector<std::string> > processFASTA(std::string);
		std::vector< std::vector< std::vector<std::string> > > processFASTA(std::string, int);
		void writeAlignmentToFile(std::vector<std::string>,std::vector< std::vector<std::string> >, std::string);
		//max values
		int maxValue(int[],int);
		int maxValue(int,int,int); 		
		std::vector< std::vector<int> > nMaxValues(std::vector<int>, int);
		int getMaxDoubleValuesIndex(std::vector<double>);
		//operations on matrices
		void transposeVec(std::vector< std::vector<int> >&);
		void transposeVec(std::vector< std::vector<double> >&);
		//std::vector<double> divideVectorByAScalar(std::vector<double>, int);
		void divideVectorByAScalar(std::vector<double>&, int);
		void multiplyVectorByAScalar(std::vector<double>&, double);
		std::vector<double> addUp(std::vector< std::vector<double> >);
		std::vector<double> convertIntVectorToDoubleVector(std::vector<int>);
		//
		int findAminoAcidsNo(char);
		void addSequenceToMultipleAlignment(std::vector<std::string>&,std::vector<std::string>);
		int countTrueValuesInVector(const std::vector<bool>&);
		void printDoubleVector(const std::vector<double>&);
};
