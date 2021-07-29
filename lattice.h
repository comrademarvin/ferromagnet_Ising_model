#pragma once

class lattice {
	private:
		bool** latticeGrid; // 2D dynamic array - true for up and false for down
		int dimension;
		double beta;
		int J;
		double averageEnergy;
		double averageMagnetisation;
		double averageSpecificHeat;
		double susceptibility;

		int getMicroEnergy(int x, int y, bool flip);
		int getTotalEnergy();
		double getTotalSpin();
	public:
		lattice(int n, double givenBeta, int givenJ);
		lattice(const lattice &copy, double givenBeta);
		~lattice();
		void randomInitialise();
		void upInitialise();
		void downInitialise();
		void printLattice();
		void reachEquilibrium();
		double getAverageEnergy();
		double getAverageMagnetisation();
		double getAverageSpecificHeat();
		double getSusceptibility();
};