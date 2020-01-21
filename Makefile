all: minSurf

DCO_INC=dco_cpp/include
DCO_LIB=dco_cpp/lib

minSurf:
	g++ -std=c++14 -O3 -I${DCO_INC} -L${DCO_LIB} minSurf.cpp solver.cpp cartesianGrid.cpp -ldcoc -o minSurf.out

minSurfHand:
	g++ -std=c++14 -O3 -I${DCO_INC} -L${DCO_LIB} minSurf.cpp solver.cpp cartesianGrid.cpp -ldcoc -o minSurfHand.out

run: minSurf
	export NAG_KUSARI_FILE=nag_key.txt && ./minSurf.out | tee z

hand: minSurfHand
	export NAG_KUSARI_FILE=nag_key.txt && ./minSurfHand.out | tee z1

clean:
	rm -f minSurf.out minSurfHand.out z z1
