#include "CSystemOfLinearEquations.h"
#include <iostream>

int main(void) {
	int rang;
	std::cout << "Enter the rang of the matrix: ";
	std::cin >> rang;
	CSystemOfLinearEquations mySystem(rang);
	std::cin >> mySystem;
	std::cout << "\nYour system of equations: \n";
	std::cout << mySystem;

	if (mySystem.checkIfConvergent()) {
		std::cout << "The system is convergent\n";
	}
	else {
		std::cout << "The system is not convergent!";
		return -1;
	}

	double* mySolutions, precision;
	int iterations;
	
	std::cout << "\nEnter precision to solve the system: ";
	std::cin >> precision;
	std::cout << '\n';

	mySolutions = mySystem.solveUsingYakobi(precision, iterations);
	
	std::cout << "\nSolutions by Yakobi: ";
	for (int i = 0; i < rang; i++) {
		std::cout << "x" << i + 1 << " = " << mySolutions[i] << ", ";
	}
	std::cout << "Iterations: " << iterations << "\n\n";

	mySolutions = mySystem.solveUsingZeidel(0.001, iterations);

	std::cout << "\nSolutions by Zeidel: ";
	for (int i = 0; i < rang; i++) {
		std::cout << "x" << i + 1 << " = " << mySolutions[i] << ", ";
	}
	std::cout << "Iterations: " << iterations;
	std::cout << "\n\n";

	return 0;
}