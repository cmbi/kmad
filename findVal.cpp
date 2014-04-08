#include "findVal.h"
#include <iostream>
#include <vector>
#include <string>
//function maxValueOf3 - finds maximum value from 3 doubles
double findVal::maxValueDoubles(double double1,double double2,double double3){
	double m = double1;
	if (double2 > double1 && double2 > double3){
		m = double2;
	}
	else if (double3 > double1){
		m = double3;
	}
	return m;
}
/* function maxValuesFromVector
finds n max values from vector<int>  vectorA */
std::vector< std::vector<int> > findVal::nMaxValues(std::vector<int> vectorA, int n){
	int i = 0;
	std::vector< std::vector<int> > result;
	std::vector<int> newEntry;
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
//function getMaxDoubleValuesIndex - returns index of the maximum value from vector<double> someVector
int findVal::getMaxDoubleValuesIndex(std::vector<double> someVector){
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
