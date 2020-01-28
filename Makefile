all: minSurf

DCO_INC=dco_cpp/include
DCO_LIB=dco_cpp/lib

minSurf:
	g++ -std=c++14 -Ofast -march=native -I${DCO_INC} -L${DCO_LIB} minSurf.cpp solver.cpp cartesianGrid.cpp -ldcoc -o minSurf.out

hard: minSurf
	export NAG_KUSARI_FILE=nag_key.txt && time ./minSurf.out 0

hand: minSurf
	export NAG_KUSARI_FILE=nag_key.txt && time ./minSurf.out 1

dco: minSurf
	export NAG_KUSARI_FILE=nag_key.txt && time ./minSurf.out 2

clean:
	rm -f minSurf.out z0 z1 z2
