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
std::vector< std::vector<int> > findVal::nMaxValues(std::vector<int>& vectorA, int n){
	int i = 0;
	std::vector< std::vector<int> > result;
	std::vector<int> newEntry;
	if ((signed)vectorA.size() < n){
		n = vectorA.size();
	}
	while (i < n){
		newEntry.push_back(vectorA[i]);
		newEntry.push_back(i);
		result.push_back(newEntry);
		newEntry.clear();
		i++;
	}
	for (unsigned int i = n; i < vectorA.size(); i++){
		for (int j = 0; j < n; j++){
			if (vectorA[i] > result[j][0] ){
				int min = 100000;
				int minIndex = -1;
				for (int k = 0; k < n; k++){
					if (min > result[k][0]){
						min = result[k][0];
						minIndex = k;	
					}	
				}
				result[minIndex][0] = vectorA[i];
				result[minIndex][1] = i;
				break;
			}
		}
	}
	return result;
}
//function getMaxDoubleValuesIndex - returns index of the maximum value from vector<double> someVector
int findVal::getMaxDoubleValuesIndex(std::vector<double>& someVector){
	int max = -100000;
	int maxIndex = -1;
	for (unsigned int i = 0; i < someVector.size(); i++){
		if (someVector[i] > max){
			max = someVector[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}
