#include <iostream>
#include <vector>
namespace{
	int a = 8;
}
int main(){
	std::vector<std::string> myvector ={"aasia"};
	std::reverse(myvector.begin(),myvector.end());
	for (int i = 0; i < myvector.size(); i++) std::cout << myvector.at(i) << "\n";
	a++;
	std::cout << a << "\n"; 
	char newchar = myvector[0][2];
	char newchar2 = 's';
	if (newchar==newchar2){
		std::cout << "passed" << std::endl;
	}
	else std::cout << myvector[0][2] << std::endl;

}
