# Compiler/Linker settings
FC = gfortran
# FLAGS =  -O3 -g -fno-automatic  -fbounds-check  -ffpe trap=invalid,zero,overflow  -ggdb3
FLFLAGS = ../VTKFortran/static/libvtkfortran.a
FCFLAGS = -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5

IMOD_LIB = -I../VTKFortran/static/mod
IOBJ_LIB = -I../VTKFortran/static/obj

# project directories
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin


# project sources
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES))

MOD_FILES=$(wildcard $(SRC_DIR)/mod*.f90)
MOD_OBJS=$(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MOD_FILES))


PROGRAM = riemann
PRG_OBJ = $(OBJ_DIR)/$(PROGRAM).o


# make without parameters will make first target found.
all : $(BIN_DIR)/$(PROGRAM)

# Compiler steps for all objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -J$(OBJ_DIR) $(IMOD_LIB) $(IOBJ_LIB) -c -o $@ $<

#$@ -c 
# Linker
$(BIN_DIR)/$(PROGRAM) : $(OBJ_FILES)
	$(FC)  -o $@ $^ -I$(OBJ_DIR) $(IMOD_LIB) $(IOBJ_LIB) $(FLFLAGS)

# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "SRCS = $(SRC_FILES)"
	@echo "OBJS = $(OBJ_FILES)"
	@echo "MODS = $(MOD_FILES)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(BIN_DIR)/$(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"


.PHONY: all clean debug



clean:
	rm -rf $(BIN_DIR)/$(PROGRAM) $(OBJ_FILES) $(patsubst %.o,%.mod,$(MOD_OBJS))


# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Module dependencies (manual)
$(OBJ_DIR)/mod_thermodynamics.o: $(OBJ_DIR)/mod_constants.o
$(OBJ_DIR)/mod_write_vtk.o: $(OBJ_DIR)/mod_thermodynamics.o $(OBJ_DIR)/mod_constants.o
