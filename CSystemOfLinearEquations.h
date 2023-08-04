#pragma once

#include <iostream>

class CSystemOfLinearEquations {
	int countOfUnknown;
	double** matrixOfCoefficients;
	double* vectorOfSolutions;

public:
	CSystemOfLinearEquations(int countOfUnknown);
	friend std::istream& operator>> (std::istream& in, CSystemOfLinearEquations object);
	friend std::ostream& operator<< (std::ostream& out, CSystemOfLinearEquations object);
	double* solveUsingYakobi(double precision, int& countOfIterations);
	double* solveUsingZeidel(double precision, int& countOfIterations);
	bool checkIfConvergent();
};

