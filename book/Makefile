FC = gfortran
FFLAGS = -g -O0 -fbacktrace -fbounds-check -fcheck=all -Wall
SRC = fullmc.f90
OBJ = $(SRC:.f90=.o)
TARGET = exe_fullmc

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(OBJ) -o $@ $(FFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# run: $(TARGET)
# 	./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)