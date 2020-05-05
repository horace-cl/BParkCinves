
/*
Anything that begins with a # is a pre-processing statement
They happen juste before the actual compilation

In this case include find a file (IOSTREAM) take all the contents of this file and paste them in 
the current file

Are typically called header files
*/


#include<iostream>

int main(){

	// << Overloaded Operator 
	std::cout << "Hello wordl!" << std::endl;
	std::cin.get();
}