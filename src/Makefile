CC = g++

FLAGS = -larmadillo

all: mysrg
	
OBJ = ModelSpace.o ReadWrite.o Operator.o HartreeFock.o

mysrg: main.cc $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(FLAGS) 

%.o: %.cc
	$(CC) -c $^ -o $@ $(FLAGS) 

clean:
	rm -f *.o
