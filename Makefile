FC = gfortran
FFLAGS = -g -O0 -fbacktrace -fbounds-check -fcheck=all -Wall
SRC = fullmc_funcs.f90 fullmc.f90
OBJ = $(SRC:.f90=.o)
TARGET = exe_fullmc

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(OBJ) -o $@ $(FFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# ensure main is rebuilt after the module
fullmc.o: fullmc_funcs.o

clean:
	rm -f $(OBJ) $(TARGET) *.mod