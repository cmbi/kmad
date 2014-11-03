#include "txtProc.h"
#include "vecUtil.h"
#include "Sequences.h"
#include "FeaturesProfile.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <tuple>
#include <algorithm>
#include <iterator>
namespace {
	static const std::vector<char> acceptable_characters = { 'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4','5','6','7','8','9'};
}
int txtProc::convertStringToInt(std::string s){
	int convertedInt;
	std::istringstream iss(s);
	iss >> convertedInt;
	return convertedInt;
}
double txtProc::convertStringToDouble(std::string s){
	double convertedDouble = atof(s.c_str());;
	return convertedDouble;
}
/*
//alternative version of split
std::vector<std::string> &txtProc::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
*/
//splits a string by the delimiter, returns a vector of strings
std::vector<std::string> txtProc::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    //split(s, delim, elems);
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
//function processFASTA - reads fasta file with encoded sequence
//writes sequences + seqNames to vector<vector<string>>; motifs' ids to vector<string> ids and motifs' probabilities to vector<double> probs
std::vector< std::vector< std::vector<std::string> > > txtProc::processFASTA(std::string filename,int codonLength, std::vector<std::string>* ids, std::vector<double>* probs){
	std::vector< std::vector< std::vector<std::string> > > resultSequences;
	std::string fastaSymbol = ">";
	std::ifstream fastafile (filename.c_str());
	std::vector<std::string> newName;
	std::vector<std::string> newSequence;
	std::vector< std::vector<std::string> > newEntry;
	bool sequences = true;
	int seqNo = -1;
	std::string newSeq = "";
	std::string line;
	while(!safeGetline(fastafile, line).eof()){
			std::string firstChar = line.substr(0,1);
			if (line != std::string("## PROBABILITIES")){
				if (sequences && firstChar == fastaSymbol){
					seqNo++;
					resultSequences.push_back(newEntry);
					newName.push_back(line);
					resultSequences[seqNo].push_back(newName);
					newName.clear();
					resultSequences[seqNo].push_back(newSequence);
				}
				else if (sequences){
					for (int i = 0; i < line.size();i++){
						if (i % codonLength == 0){
							std::string newResidue = "";
							//j for goes through all codon postions of this residue
							for (int j = i;j < i + codonLength; j++){
									if (acceptableChar(line[j])){
										newResidue += line[j];
									}
									else{
										std::cout << "I found a weird character (" << line[j] << "), so you probbaly want to go through your files and double check them"  << std::endl;
										std::exit(0);
									}
							}
							resultSequences[seqNo][1].push_back(newResidue);
						}
					}
				}
				//else means we're already in the motifs probs section
				else{
					std::istringstream iss(line);
					std::vector<std::string> motif{std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>{}};
					if (motif.size() == 2){
						ids->push_back(motif[0]);
						probs->push_back(convertStringToDouble(motif[1]));
					}
				}
			}
			else{
				sequences = false;
			}
		}
	fastafile.close();
	return resultSequences;
}
//function writeAlignmentToFile
void txtProc::writeAlignmentToFile(std::vector<std::string> sequences,std::vector< std::vector< std::vector<std::string> > > sequencesWithNames, std::string filename){
	std::stringstream sstr;
	sstr << filename << "_al";
	std::ofstream outputFile(sstr.str().c_str(),std::ios::out);
	for (int i = 0; i < sequences.size() ;i++){
		outputFile << sequencesWithNames[i][0][0]<< "\n" << sequences[i] << "\n";
	}
}
//write alignment to file as a nonencoded fasta
void txtProc::writeAlignmentWithoutCodeToFile(std::vector<std::string> sequences,std::vector< std::vector<std::vector<std::string> > > sequencesWithNames, std::string filename){
	std::stringstream sstr;
	sstr << filename << "_al";
	std::ofstream outputFile(sstr.str().c_str(),std::ios::out);
	for (int i = 0; i < sequences.size() ;i++){
		outputFile << sequencesWithNames[i][0][0]<< "\n";
		std::string sequence="";
		for (int j = 0; j < sequences[i].size(); j+=4){
			sequence += sequences[i][j];
		}
		outputFile << sequence << std::endl;
	}
}
//write to file (2d vector as a 1d vector)
void txtProc::writeVector(std::vector<std::vector<double>> vec, std::string filename){
	std::stringstream sstr;
	sstr << filename;
	std::ofstream outputFile(sstr.str().c_str(),std::ios::out);
	for (int i = 0; i < vec.size() ;i++){
		for (int j = 0; j < vec[0].size(); j++){
			outputFile << vec[i][j] << std::endl;
		}
	}
}
std::string txtProc::charToString(char mychar){
	return std::string(1,mychar);
}
std::string txtProc::charToString(char mychar1, char mychar2){
	std::string newstring = std::string(1,mychar1);
	newstring.push_back(mychar2);
	return newstring;
}
std::istream& txtProc::safeGetline(std::istream& is, std::string& t)
{
	t.clear();
	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();
        for(;;) {
	        int c = sb->sbumpc();
		switch (c) {
			case '\n':
				return is;
			case '\r':
				if(sb->sgetc() == '\n')
				sb->sbumpc();
				return is;
			case EOF:
				// Also handle the case when the last line has no line ending
				if(t.empty())
					is.setstate(std::ios::eofbit);
				return is;
			default:
				t += (char)c;
		}
	}
}
//check if the character is supported
bool txtProc::acceptableChar(char my_char){
	bool result = false;
	for (int i = 0; i < acceptable_characters.size(); i++){
		if (acceptable_characters[i]==my_char){
			result = true;
			break;
		}
	}
	return result;
}
void txtProc::process_conf_file(std::string filename, FeaturesProfile& feat_profile, Sequences& sequences_aa){
	std::ifstream conf_file(filename.c_str());
	std::string line;
	std::vector<std::tuple<std::string, std::string, int, int, int, double, double, double, double, std::string, std::string> > usr_feature_rules;
	std::vector<std::tuple<std::string, std::string, int, int, int> > feature_rules;
    bool features = true;
	while(!safeGetline(conf_file, line).eof()){
       // if (line[0] != '#' && features){

       // }
		else if (line[0] != '#'){
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			std::vector<std::string> tmp_vector = split(line,';');
			usr_feature_rules.push_back(std::make_tuple(tmp_vector[0], tmp_vector[1], std::stoi(tmp_vector[2]), std::stoi(tmp_vector[3]), std::stoi(tmp_vector[4]),std::stod(tmp_vector[5]),std::stod(tmp_vector[6]),std::stod(tmp_vector[7]),std::stod(tmp_vector[8]),tmp_vector[9],tmp_vector[10]));
		}
	}
	feat_profile.add_USR_features(feature_rules);
	sequences_aa.add_features(feature_rules);
	feat_profile.setRules(feature_rules);
}
// converts the string form the conf file to vector of positions of features to be scored
std::vector<int> txtProc::unfold(std::string conf_string, std::vector<std::string> listOfFeatures){
	std::vector<std::string> tmp_vector = split(conf_string,',');
	std::vector<int> out_vector;
	for (int i = 0; i < tmp_vector.size(); i++){
		if (split(tmp_vector[i],'_').size() > 1){						// this is a single feature entry, e.g. 'PF_A'
			out_vector.push_back(vecUtil::findIndex(std::string("USR_")+tmp_vector[i], listOfFeatures));
		}
		else if (split(tmp_vector[i],'[').size() == 1){						// this is an entry with only the tag specified (without any exceptions)
			for (int j = 0; j < listOfFeatures.size(); j++){
				std::vector<std::string> singlefeat = split(listOfFeatures[j],'_');
				if (singlefeat.size() > 1 && singlefeat[1] == tmp_vector[i]){out_vector.push_back(j);}
			}
		}
		else{//TAG with exceptions
			std::vector<std::string> tagfeat = split(tmp_vector[i],'[');
			std::string tag = tagfeat[0];
			std::vector<std::string> exceptions = split(split(tagfeat[1],']')[0],'.');
			for (int j = 0; j < listOfFeatures.size(); j++){
				std::vector<std::string> singlefeat = split(listOfFeatures[j],'_');
				if (singlefeat.size() > 1 && singlefeat[1] == tag && !vecUtil::contains(exceptions, singlefeat[2])){
					out_vector.push_back(j);
				}
			}
		}
	}
	return out_vector;
}
