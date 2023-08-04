#include "CSystemOfLinearEquations.h"
#include <iomanip>
#include <vector>


double findMaxInArray(double* array, int count) {
	double max = array[0];
	for (int i = 0; i < count; i++) {
		if (max > array[i]) {
			max = array[i];
		}
	}
	return max;
}

double sumOf(double** matrixOfCoefficients, int countOfUnknown, int row, double* tempSolutions) {
	double sum = 0;
	for (int i = 0; i < countOfUnknown; i++) {
		if (i != row) {
			sum += (matrixOfCoefficients[row][i] / matrixOfCoefficients[row][row]) * tempSolutions[i];
		}
	}
	return sum;
}

double findingDelta(double* tempSolutions, double* vectorOfSolutions, int countOfUnknown) {
	double sum = 0;
	for (int i = 0; i < countOfUnknown; i++) {
		sum += (vectorOfSolutions[i] - tempSolutions[i]) * (vectorOfSolutions[i] - tempSolutions[i]);
	}
	return sqrt(sum);
}

CSystemOfLinearEquations::CSystemOfLinearEquations(int countOfUnknown) {
	this->countOfUnknown = countOfUnknown;
	matrixOfCoefficients = new double* [countOfUnknown];
	for (int i = 0; i < countOfUnknown; i++) {
		matrixOfCoefficients[i] = new double[countOfUnknown + 1];
	}
	vectorOfSolutions = new double[countOfUnknown];
	for (int i = 0; i < countOfUnknown; i++) {
		vectorOfSolutions[i] = 0;
	}
}

double* CSystemOfLinearEquations::solveUsingYakobi(double precision, int &countOfIterations) {
	countOfIterations = 0;

	double delta = precision + 1;
	double* tempSolutions = new double[countOfUnknown];
	//starting precisions are free members vector divided by i i element
	for (int i = 0; i < countOfUnknown; i++) {
		tempSolutions[i] = matrixOfCoefficients[i][countOfUnknown] / matrixOfCoefficients[i][i];
	}

	//printting first precisions
	for (int i = 0; i < countOfUnknown; i++) {
		std::cout << std::setw(7) << std::setprecision(4) << tempSolutions[i] << " ";
	}
	std::cout << std::endl;

	while (precision < delta) {
		for (int i = 0; i < countOfUnknown; i++) { //For each x: x1, x2, x3...
			vectorOfSolutions[i] = matrixOfCoefficients[i][countOfUnknown] / matrixOfCoefficients[i][i] - sumOf(matrixOfCoefficients, countOfUnknown, i, tempSolutions);    /*b/a - sum (aij/aii)*x */
		}
		delta = findingDelta(tempSolutions, vectorOfSolutions, countOfUnknown);

		//making tempsolutions the real solutions
		for (int i = 0; i < countOfUnknown; i++) {
			tempSolutions[i] = vectorOfSolutions[i];
			std::cout << std::setw(7) << std::setprecision(4) << tempSolutions[i] << " ";
		}
		countOfIterations++;

		std::cout << std::setw(7) << std::setprecision(4) << delta << "\n";
	}

	delete[] tempSolutions;

	return vectorOfSolutions;
}

double* CSystemOfLinearEquations::solveUsingZeidel(double precision, int& countOfIterations) {
	countOfIterations = 0;
	double delta = precision + 1;

	//We need a temp vector, which we will be using to find delta
	double* vectorOfBefore = new double[countOfUnknown];

	//starting precisions are free members vector / i i element
	for (int i = 0; i < countOfUnknown; i++) {
		vectorOfSolutions[i] = matrixOfCoefficients[i][countOfUnknown] / matrixOfCoefficients[i][i];
		vectorOfBefore[i] = vectorOfSolutions[i];
	}

	//printing starting precisions
	for (int i = 0; i < countOfUnknown; i++) {
		std::cout << std::setw(7) << std::setprecision(4) << vectorOfSolutions[i] << " ";
	}
	std::cout << std::endl;

	while (precision < delta) {
		for (int i = 0; i < countOfUnknown; i++) { //For each x: x1, x2, x3...
			vectorOfSolutions[i] = matrixOfCoefficients[i][countOfUnknown] / matrixOfCoefficients[i][i] - sumOf(matrixOfCoefficients, countOfUnknown, i, vectorOfSolutions);    /*b/a - sum (aij/aii)*x */
		}
		delta = findingDelta(vectorOfBefore, vectorOfSolutions, countOfUnknown);

		//setting back the vectorOfBefore
		for (int i = 0; i < countOfUnknown; i++) {
			vectorOfBefore[i] = vectorOfSolutions[i];
			std::cout << std::setw(7) << std::setprecision(4) << vectorOfBefore[i] << " ";
		}
		countOfIterations++;
		std::cout << std::setw(7) << std::setprecision(4) << delta << "\n";
	}

	delete[] vectorOfBefore;

	return vectorOfSolutions;
}

bool CSystemOfLinearEquations::checkIfConvergent() {
	std::vector<std::vector<double>> matrixAlpha(countOfUnknown, std::vector<double> (countOfUnknown, 0));

	//finding Alpha matrix
	for (int i = 0; i < countOfUnknown; i++) { //Rows
		for (int j = 0; j < countOfUnknown; j++) {  //Columns
			if (i != j) {
				matrixAlpha[i][j] = matrixOfCoefficients[i][j] / matrixOfCoefficients[i][i];
			}
		}
	}

	//printing Alpha matrix
	std::cout << "\nAlpha matrix:\n";
	for (int i = 0; i < countOfUnknown; i++) {
		for (int j = 0; j < countOfUnknown; j++) {
			std::cout << std::setw(7) << std::setprecision(4) << matrixAlpha[i][j] << " ";
		}
		std::cout << std::endl;
	}

	//finding and printing beta vector
	std::cout << "\nBeta-vector:\n";
	for (int i = 0; i < countOfUnknown; i++) {
		std::cout << matrixOfCoefficients[i][countOfUnknown] / matrixOfCoefficients[i][i] << ' ';
	}

	double* firstNormOfMatrix = new double[countOfUnknown], * secondNormOfMatrix = new double[countOfUnknown];

	//nulling matrixes
	for (int i = 0; i < countOfUnknown; i++) {
		firstNormOfMatrix[i] = 0;
		secondNormOfMatrix[i] = 0;
	}

	double thirdNormOfMatrix = 0;
	for (int i = 0; i < countOfUnknown; i++) {
		for (int j = 0; j < countOfUnknown; j++) {
			firstNormOfMatrix[i] += fabs(matrixAlpha[i][j]);
			secondNormOfMatrix[i] += fabs(matrixAlpha[j][i]);
			thirdNormOfMatrix += matrixAlpha[i][j] * matrixAlpha[i][j];
		}
	}

	thirdNormOfMatrix = sqrt(thirdNormOfMatrix);

	std::cout << "\n\nFirst norm = " << findMaxInArray(firstNormOfMatrix, countOfUnknown) << "\nSecond norm = " << findMaxInArray(secondNormOfMatrix, countOfUnknown) << "\nThird norm = " << thirdNormOfMatrix << "\n";

	if (findMaxInArray(firstNormOfMatrix, countOfUnknown) < 1 || findMaxInArray(secondNormOfMatrix, countOfUnknown) < 1 || thirdNormOfMatrix < 1) {
		delete[] firstNormOfMatrix; delete[] secondNormOfMatrix;
		return true;
	}
	else {
		delete[] firstNormOfMatrix; delete[] secondNormOfMatrix;
		return false;
	}
}

std::istream& operator>>(std::istream& in, CSystemOfLinearEquations object) {
	for (int i = 0; i < object.countOfUnknown; i++) {
		std::cout << "Enter coefficients of " << i + 1 << " line: ";
		for (int j = 0; j < object.countOfUnknown + 1; j++) {
			in >> object.matrixOfCoefficients[i][j];
		}
	}
	return in;
}

std::ostream& operator<<(std::ostream& out, CSystemOfLinearEquations object) {
	out << std::setprecision(4);
	for (int i = 0; i < object.countOfUnknown; i++) {
		for (int j = 0; j < object.countOfUnknown + 1; j++) {
			if (j == 0) {
				out << std::setw(8) << object.matrixOfCoefficients[i][j] << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << "x" << j + 1;
			}
			else if (j > 0 && j != object.countOfUnknown) {
				out << std::showpos << std::setw(8) << object.matrixOfCoefficients[i][j] << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << "x" << j + 1;
			}
			else if (j = object.countOfUnknown) {
				out << "=" << object.matrixOfCoefficients[i][j];
			}
		}
		out << std::endl;
	}
	return out;
}