CC = g++
CFLAG = -std=c++11 -W -Wall -Wextra -g
TESTGEN = cxxtestgen
TESTS = test.h
FLAGTEST = --error-printer
LIB = -O1 -larmadillo -lm
OBJ = ../src/poly.o ../src/basis.o ../src/sol.o test.o


test: $(OBJ) 
	$(CC) $(CFLAG) $^ -o $@ $(LIB) 

%.o: %.cpp %.h
	$(CC) $(CFLAG) $< -c -o $@ $(LIB)

test.cpp: 
	$(TESTGEN) $(FLAGTEST) $(TESTS) -o $@ 


clean:
	rm  test.cpp
	rm  test
	rm test.o
