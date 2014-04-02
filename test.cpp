#include <iostream>

void AddOne(int &y)
{
	    y++;
}
 
int main()
{
	    int x = 5;
	    std::cout << "x = " << x << std::endl;
	AddOne(x);
	std::cout << "x = " << x << std::endl;
			 
			    return 0;
}
