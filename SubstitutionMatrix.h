#include <iostream>
#include <vector>
class SubstitutionMatrix{
public:
	SubstitutionMatrix(); //constructor
	void printSbstMatrix();
	std::vector< std::vector<double> > convertToProfileFormat(std::string);
	//getters
	char getLetter(int);
	int getElement(char,char);
	int getElement(int,int);
	std::vector<int> getColumn(int);
private:
	//functions
	//variables
/*
	std::vector<char> alphabet;
	std::vector< std::vector<int> > simScores;*/
	int testInt;
};
