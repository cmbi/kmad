#include <iostream>
#include <vector>
using namespace std;

int main(){
	vector<double> hello;
	hello.push_back(0);
	hello.push_back(0);
	hello.push_back(3);
	hello[-1]++;
	for (int i = 0; i < hello.size();i++){
		cout << hello.at(i)<< "\n";
		cout << hello[i]<< "\n";
	}
	cout << hello[-1] << endl;
	
}
