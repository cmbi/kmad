#include "txtProc.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
//function convertStringToInt
int txtProc::convertStringToInt(std::string s){
	int convertedInt;
	std::istringstream iss(s);
	iss >> convertedInt;	
	return convertedInt;
}
/* function processFASTA
reads 'filename' file, returns an array of sequences and their names - [[sequenceName,sequence],...] */
std::vector< std::vector<std::string> > txtProc::processFASTA(std::string filename){
	std::vector< std::vector<std::string> > result;
	std::string line;
	std::string fastaSymbol = ">";
	std::vector<std::string> newEntry;
	std::ifstream fastafile (filename.c_str());
	int seqNo = -1;
	std::string newSeq = "";
	if (fastafile.is_open()){
		while(fastafile){
			getline(fastafile,line);
			std::string firstChar = line.substr(0,1);
			if (firstChar == fastaSymbol){
				seqNo++;
				result.push_back(newEntry);
				result.at(seqNo).push_back(line);
				result.at(seqNo).push_back(newSeq);
			}
			else{
				result.at(seqNo).at(1)=result.at(seqNo).at(1).append(line);	
			}

		}
		fastafile.close();
	}
	else{
		std::cout << "Where is the file????";
	}
	return result;
}
//function writeAlignmentToFile
void txtProc::writeAlignmentToFile(std::vector<std::string> sequences,std::vector< std::vector<std::string> > sequencesWithNames, std::string filename){
	std::stringstream sstr;
	sstr << filename << "_al";
	std::ofstream outputFile(sstr.str().c_str(),std::ios::out);
	for (int i = 0; i < sequences.size() ;i++){
		outputFile << sequencesWithNames.at(i).at(0)<< "\n" << sequences.at(i) << "\n";
	}
}
