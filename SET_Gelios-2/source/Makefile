CC = g++
CFLAGS = -c -O3 -std=c++17 -fopenmp
LDFLAGS = -fopenmp
SRC = Gelios_2D.cpp Gran.cpp Setka.cpp Cell.cpp Point.cpp Rail.cpp Solvers.cpp Help.cpp sensor.cpp sensor2.cpp MKmethod.cpp Dist_func.cpp Interpol_Setka.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = 1set_progect  # или проект_сет, проект_цуй-шен, проект_эол

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXEC)

distclean: clean
	rm -f *~ *.out *.dat *.log

rebuild: clean all

.PHONY: all clean distclean rebuild