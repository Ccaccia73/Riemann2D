# Compiler/Linker settings
FC = gfortran
# FLAGS =  -O3 -g -fno-automatic  -fbounds-check  -ffpe trap=invalid,zero,overflow  -ggdb3
FLFLAGS = -g
FCFLAGS = -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5


# project directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin


# project sources
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES))



PROGRAM = riemann
PRG_OBJ = $(OBJ_DIR)/$(PROGRAM).o


# make without parameters will make first target found.
all : $(BIN_DIR)/$(PROGRAM)

# Compiler steps for all objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

#$@ -c 
# Linker
$(BIN_DIR)/$(PROGRAM) : $(OBJ_FILES)
	$(FC) $(FLFLAGS) -o $@ $^

# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "SRCS = $(SRC_FILES)"
	@echo "OBJS = $(OBJ_FILES)"
#	@echo "MODS = $(MODS)"
#	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = bin/$(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"


.PHONY: all clean debug



clean:
	rm -rf $(BIN_DIR)/$(PROGRAM) $(OBJ_FILES)


# Dependencies

# Main program depends on all modules
#$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
#mod_blocks.o mod_allocations.o : mod_shared.o