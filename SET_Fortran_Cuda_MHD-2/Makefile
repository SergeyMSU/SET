# Компилятор и флаги
FC = nvfortran
FCFLAGS = -O3 -fopenmp -gpu=rdc -cuda
LDFLAGS = -O3 -fopenmp -gpu=rdc -cuda

# Исходные файлы
SOURCES = Storage_module.f90 Solvers.f90 cuf_kernel.cuf cuf_solvers.cuf FCMHD.f90
OBJECTS = $(SOURCES:.f90=.o)
OBJECTS := $(OBJECTS:.cuf=.o)

# Цели
all: my_program

my_program: $(OBJECTS)
	$(FC) -o $@ $^ $(LDFLAGS)

# Правило для .f90 файлов
%.o: %.f90
	$(FC) -c $(FCFLAGS) $< -o $@

# Правило для .cuf файлов  
%.o: %.cuf
	$(FC) -c $(FCFLAGS) $< -o $@

clean:
	rm -f *.o *.mod my_program

.PHONY: all clean