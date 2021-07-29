#include "lattice.h"
#include <iostream>

int main () {
	// Changeable parameters
	int N = 20; // side-length of lattice
	int J = 1; // -1 for paramagnetic

	int index = 0;
	double energy[16];
	double magnetisation[16];
	double specificheat[16];
	double susceptibility[16];
	// beta range: 0.05 < B < 0.8
	lattice* prevGrid;
	lattice* currGrid;
	for (double beta = 0.05; beta < 0.85; beta += 0.05) {
		if (beta == 0.05) {
			currGrid = new lattice(N, beta, J);
			prevGrid = currGrid;
		} else {
			currGrid = new lattice(*prevGrid, beta);
			delete prevGrid;
			prevGrid = currGrid;
		}

		currGrid->reachEquilibrium();

		currGrid->printLattice();

		energy[index] = currGrid->getAverageEnergy();
		magnetisation[index] = currGrid->getAverageMagnetisation();
		specificheat[index] = currGrid->getAverageSpecificHeat();
		susceptibility[index] = currGrid->getSusceptibility();

		index++;
	}

	delete currGrid;

	std::cout << "Average system energy:\n";
	for (int i = 0; i < 16; i++) {
		std::cout << energy[i];
		if (i != 15) {
			std::cout << ",";
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Magnetisation of system:\n";
	for (int i = 0; i < 16; i++) {
		std::cout << magnetisation[i];
		if (i != 15) {
			std::cout << ",";
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Specific Heat of system:\n";
	for (int i = 0; i < 16; i++) {
		std::cout << specificheat[i];
		if (i != 15) {
			std::cout << ",";
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Susceptibility of system:\n";
	for (int i = 0; i < 16; i++) {
		std::cout << susceptibility[i];
		if (i != 15) {
			std::cout << ",";
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;

	return 0;
}