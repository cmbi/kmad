//UsefulStuff class implementation
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "UsefulStuff.h"
using namespace std;
// function maxValue - finds maximum value from an array of ints
int UsefulStuff::maxValue(int arrayOfInts[], int arraySize){
	int m = arrayOfInts[0];
	for (int i = 0; i < arraySize; i++){
		if (arrayOfInts[i]>m){
			m=arrayOfInts[i];
		}
	}
	return m;
}
//function maxValueOf3 - finds maximum value from 3 ints
int UsefulStuff::maxValue(int int1,int int2,int int3){
	int m = int1;
	if (int2 > int1 && int2 > int3){
		m = int2;
	}
	else if (int3 > int1){
		m = int3;
	}
	return m;
}
/* function maxValuesFromVector
finds n max values from vector<int>  vectorA */
vector< vector<int> > UsefulStuff::nMaxValues(vector<int> vectorA, int n){
	int i = 0;
	vector< vector<int> > result;
	vector<int> newEntry;
	if (vectorA.size() < n){
		n = vectorA.size();
	}
	while (i < n){
		newEntry.push_back(vectorA.at(i));
		newEntry.push_back(i);
		result.push_back(newEntry);
		newEntry.clear();
		i++;
	}
	for (int i = n; i < vectorA.size(); i++){
		for (int j = 0; j < n; j++){
			if (vectorA.at(i) > result.at(j).at(0) ){
				int min = 100000;
				int minIndex = -1;
				for (int k = 0; k < n; k++){
					if (min > result.at(k).at(0)){
						min = result.at(k).at(0);
						minIndex = k;	
					}	
				}
				result.at(minIndex).at(0) = vectorA.at(i);
				result.at(minIndex).at(1) = i;
				break;
			}
		}
	}
	return result;
}
/* function processFASTA
reads 'filename' file, returns an array of sequences and their names - [[sequenceName,sequence],...] */
vector< vector<string> > UsefulStuff::processFASTA(string filename){
	vector< vector<string> > result;
	string line;
	string fastaSymbol = ">";
	vector<string> newEntry;
	ifstream fastafile (filename.c_str());
	int seqNo = -1;
	string newSeq = "";
	if (fastafile.is_open()){
		while(fastafile){
			getline(fastafile,line);
			string firstChar = line.substr(0,1);
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
		cout << "Where is the file????";
	}
	return result;
}
//function processFASTA - reads fasta file with encoded sequence
vector< vector< vector<string> > > UsefulStuff::processFASTA(string filename,int codonLength){
	vector< vector< vector<string> > > result;
	string line;
	string fastaSymbol = ">";
	ifstream fastafile (filename.c_str());
	vector<string> newName;
	vector<string> newSequence;
	vector< vector<string> > newEntry;
	int seqNo = -1;
	string newSeq = "";
	if (fastafile.is_open()){
		while(fastafile){
			getline(fastafile,line);
			string firstChar = line.substr(0,1);
			if (firstChar == fastaSymbol){
				seqNo++;
				result.push_back(newEntry);
				newName.push_back(line);
				result.at(seqNo).push_back(newName);
				newName.clear();
				result.at(seqNo).push_back(newSequence);
			}
			else{
				for (int i = 0; i < line.size();i++){
					if (i % codonLength == 0){
						string newResidue = "";
						for (int j = i;j < i + codonLength; j++){
							newResidue += line[j];
						}
						result.at(seqNo).at(1).push_back(newResidue);
					}
				}
			}

		}
		fastafile.close();
	}
	else{
		cout << "Where is the file????";
	}
	return result;
}
/* function combinePairwiseAlignments
arguments: 
	vector<string> al1 - multiple alignment
	vector<string> al2 - pairwise alignment
returns:
	vector<string>  - multiple alignment with one added sequence
*/
void UsefulStuff::addSequenceToMultipleAlignment( vector<string>& al1,  vector<string> al2){    ///al1 - multiple alignment, al2 - pairwise alignment
	int i = 0;
	int j = 0;
	vector<string> result;
	string seq="";
	for (int i =0; i<al1.size()+1;i++){
		result.push_back(seq);
	}
	char gap = '-';
	while (i < al1.at(0).size() || j < al2.at(0).size()){
		char al20j = al2.at(0).at(j);
		char al10i = al1.at(0).at(i);
		if (al10i != gap && al20j != gap){
			for (int k = 0; k< al1.size();k++){
				result.at(k)+=al1.at(k).at(i);
			}
			result.at(al1.size())+=al2.at(1).at(j);		
			i++;
			j++;
		}
		else if (al10i == gap){
			for (int k = 0; k< al1.size();k++){
				result.at(k)+=al1.at(k).at(i);
			}
			result.at(al1.size())+="-";
			i++;
		}
		else{
			for (int k = 0; k< al1.size();k++){
				result.at(k)+="-";
			}
			result.at(al1.size())+=al2.at(1).at(j);
			j++;
		}
	}
	al1 = result;
}
//function convertStringToInt
int UsefulStuff::convertStringToInt(string s){
	int convertedInt;
	istringstream iss(s);
	iss >> convertedInt;	
	return convertedInt;
}
//function transposeVec
void UsefulStuff::transposeVec(vector< vector<int> >& vec){
	vector< vector<int> > newVec;
	vector<int> newRow;
	for (int i = 0; i < vec.at(0).size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec.at(j).at(i));	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
	//return newVec;
}
//function transposeVec  - transposes vector< vector<double> > and returns vector< vector<double> >
void UsefulStuff::transposeVec(vector< vector<double> >& vec){
	vector< vector<double> > newVec;
	vector<double> newRow;
	for (int i = 0; i < vec.at(0).size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec.at(j).at(i));	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}
//function findAminoAcidsNo - finds index of the given char aa amino acid
int UsefulStuff::findAminoAcidsNo(char aa){
	char alphabetArr[] = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V' };
	int aAcidint = -1;
	for (int i = 0; i < 20; i++){
		if (aa == alphabetArr[i]){
			aAcidint = i;
			break;
		}
	}
	return aAcidint;
}
//function getMaxDoubleValuesIndex - returns index of the maximum value from vector<double> someVector
int UsefulStuff::getMaxDoubleValuesIndex(vector<double> someVector){
	int max = -100000;
	int maxIndex;
	for (int i = 0; i < someVector.size(); i++){
		if (someVector.at(i) > max){
			max = someVector.at(i);
			maxIndex = i;
		}
	}
	return maxIndex;
}
//function divideVectorByAScalar
void UsefulStuff::divideVectorByAScalar(vector<double>& vec, int scalar){
	vector<double> result;
	for (int i = 0; i < vec.size(); i++){
		result.push_back(vec[i]/scalar);	
	}
	vec = result;
}
void UsefulStuff::multiplyVectorByAScalar(vector<double>& vec, double scalar){
	vector<double> result;
	for (int i = 0; i < vec.size(); i++){
		result.push_back(vec[i]*scalar);
	}
	vec = result;
}
//function writeAlignmentToFile
void UsefulStuff::writeAlignmentToFile(vector<string> sequences,vector< vector<string> > sequencesWithNames, string filename){
	stringstream sstr;
	sstr << filename << "_al";
	ofstream outputFile(sstr.str().c_str(),ios::out);
	for (int i = 0; i < sequences.size() ;i++){
		outputFile << sequencesWithNames.at(i).at(0)<< endl << sequences.at(i) << endl;
	}
}
//convertIntVectorToDoubleVector
vector<double> UsefulStuff::convertIntVectorToDoubleVector(vector<int> vec){
	vector<double> result;
	for (int i = 0; i < vec.size();i++){
		result.push_back(double(vec.at(i)));
	}
	return result;
}
//function addUp - takes 2D matrix, adds up elements from each column, returns a 1D vector
vector<double> UsefulStuff::addUp(vector< vector<double> > vec){
	vector<double> newVec;
	for (int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for (int j = 0; j < vec.size(); j++){
			sum += vec.at(j).at(i);
		}
		newVec.push_back(sum);
	}
	return newVec;
}
//function countTrueValuesInVector
int UsefulStuff::countTrueValuesInVector(const vector<bool>& vec){
	int result = 0;
	for (int i = 0; i < vec.size(); i++){
		if (vec[i]){
			result++;	
		}
	}
	return result;
}
//function printDoubleVector
void UsefulStuff::printDoubleVector(const vector<double>& vec){
	for (int i = 0; i < vec.size(); i++){
		cout << vec[i];
	}
	cout << "\n";
}

