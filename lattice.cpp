#include "lattice.h"
#include <iostream>
#include <time.h>
#include <math.h>

lattice::lattice(int n, double givenBeta, int givenJ) {
	this->averageEnergy = 0;
	this->averageMagnetisation = 0;
	this->averageSpecificHeat = 0;
	this->susceptibility = 0;
	this->dimension = n;
	this->beta = givenBeta;
	this->J = givenJ;
	this->latticeGrid = new bool*[n];

	for (int i = 0; i < n; i++) {
		this->latticeGrid[i] = new bool[n];
	}

	this->randomInitialise(); // default initialize to random

	srand(time(NULL)); // random seed
};

lattice::lattice(const lattice &copy, double givenBeta) {
	this->averageEnergy = 0;
	this->averageMagnetisation = 0;
	this->averageSpecificHeat = 0;
	this->susceptibility = 0;
	this->dimension = copy.dimension;
	this->beta = givenBeta;
	this->J = copy.J;

	this->latticeGrid = new bool*[this->dimension];

	for (int i = 0; i < this->dimension; i++) {
		this->latticeGrid[i] = new bool[this->dimension];
		for (int j = 0; j < this->dimension; j++) {
			this->latticeGrid[i][j] = copy.latticeGrid[i][j];
		}
	}
}

lattice::~lattice() {
	for (int i = 0; i < this->dimension; i++) {
		delete[] this->latticeGrid[i];
	}

	delete[] this->latticeGrid;
}

void lattice::randomInitialise() {
	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			int myRandom = rand() % 2;
			if (myRandom == 0) {
				latticeGrid[i][j] = true;
			} else {
				latticeGrid[i][j] = false;
			}
		}
	}
}

void lattice::upInitialise() {
	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			latticeGrid[i][j] = true;
		}
	}
}

void lattice::downInitialise() {
	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			latticeGrid[i][j] = false;
		}
	}
}

void lattice::printLattice() {
	std::cout << "Beta = " << this->beta << std::endl;
	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			if (latticeGrid[i][j] == true) {
				//std::cout << "↑";
				std::cout << "*";
			} else {
				//std::cout << "↓";
				std::cout << ".";
			}
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}

int getSpin(bool state) { // gives spin value from boolean
	if (state) {
		return 1;
	}

	return -1;
}

int lattice::getMicroEnergy(int x, int y, bool flip) { // calculates energy of microstate
	int x_plus, x_minus, y_plus, y_minus;

	if ((x + 1) == this->dimension) { // at edge
		x_plus = 0;
	} else {
		x_plus = x + 1;
	}

	if ((x - 1) < 0) { // at edge
		x_minus = this->dimension - 1;
	} else {
		x_minus = x - 1;
	}

	if ((y + 1) == this->dimension) { // at edge
		y_plus = 0;
	} else {
		y_plus = y + 1;
	}

	if ((y - 1) < 0) { // at edge
		y_minus = this->dimension - 1;
	} else {
		y_minus = y - 1;
	}

	int s;
	if (flip) { // case of energy for flipped state
		s = getSpin(!(latticeGrid[x][y]));
	} else {
		s = getSpin(latticeGrid[x][y]);
	}

	int s1 = getSpin(latticeGrid[x_plus][y]);
	int s2 = getSpin(latticeGrid[x_minus][y]);
	int s3 = getSpin(latticeGrid[x][y_plus]);
	int s4 = getSpin(latticeGrid[x][y_minus]);

	return -1*(this->J)*((s*s1) + (s*s2) + (s*s3) + (s*s4));
}

void lattice::reachEquilibrium() {
	double energyTotal = 0;
	double energyTotal_squared = 0;
	double averageTotal = 0;
	double spinTotal = 0;
	double spinTotal_squared = 0;

	for (int runs = 0; runs < 100000; runs++) {
		int x = rand() % this->dimension; // random x-coordinate
		int y = rand() % this->dimension; // random y-coordinate

		int E_i = this->getMicroEnergy(x, y, false);
		int E_j = this->getMicroEnergy(x, y, true); // true to flip

		if (E_j < E_i) {
			latticeGrid[x][y] = !(latticeGrid[x][y]); // flip spin
		} else {
			double r = (double) rand() / (RAND_MAX); // random between 0 < r < 1
			double probability = exp(-1*(this->beta)*(E_j - E_i));

			if (r < probability) {
				latticeGrid[x][y] = !(latticeGrid[x][y]); // flip spin
			}
			// std::cout << "Get accepted by the rrrr...\n";
		}

		if (runs > 60000) { // get averages after equilibrium
			averageTotal++;

			int currentEnergy = this->getTotalEnergy();
			energyTotal += currentEnergy;
			energyTotal_squared += currentEnergy*currentEnergy;

			double currentSpinTotal = this->getTotalSpin();
			spinTotal += currentSpinTotal;
			spinTotal_squared += currentSpinTotal*currentSpinTotal;
		}
	}

	this->averageEnergy = energyTotal/averageTotal;

	double averageEnergy_squared = energyTotal_squared/averageTotal;
	this->averageSpecificHeat = (J*beta)*(J*beta)*((averageEnergy_squared - (averageEnergy*averageEnergy))/(dimension*dimension));

	double averageSpin = spinTotal/averageTotal;
	double averageSpin_squared = spinTotal_squared/averageTotal;
	this->susceptibility = (J*beta)*(averageSpin_squared - (averageSpin*averageSpin))/(dimension*dimension);
	this->averageMagnetisation = averageSpin/(dimension*dimension);
}

int lattice::getTotalEnergy() { // returns total energy of current system state
	int total = 0;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			total += getMicroEnergy(i, j, false);
		}
	}

	return total;
}

double lattice::getTotalSpin() { // returns sum off all spins for current system state
	double total = 0;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			total += getSpin(latticeGrid[i][j]);
		}
	}

	return total;
}

double lattice::getAverageEnergy() {
	return this->averageEnergy;
}

double lattice::getAverageMagnetisation() {
	return this->averageMagnetisation;
}

double lattice::getAverageSpecificHeat() {
	return this->averageSpecificHeat;
}

double lattice::getSusceptibility() {
	return this->susceptibility;
}