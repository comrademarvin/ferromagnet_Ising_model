lattice: lattice.o main.o
	g++ -o lattice main.o lattice.o

lattice.o: lattice.cpp main.cpp
	g++ -c main.cpp lattice.cpp

run:
	./lattice

clean:
	rm *.o lattice