# Makefile for TB-tbG Fortran project
# Compiler and flags
FC = mpiifx
FFLAGS = -O2 -qmkl -module build -Ibuild -xHost
LDFLAGS = -qmkl

# Target executable
TARGET = tbsolver

# Source directory
SRCDIR = src

# Build directory
BUILDDIR = build

# Fortran source files (all .f90 files in src/)
SRCS = $(wildcard $(SRCDIR)/*.f90)

# Object files (placed in build/)
OBJS = $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(SRCS))

# Module files (placed in build/)
MODS = $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.mod, $(filter-out $(SRCDIR)/main.f90, $(SRCS)))

# Determine module dependencies (simplified)
# We'll define a fixed order to ensure modules are compiled before their dependents
# This list should be ordered from least dependent to most dependent
MODULE_ORDER = constants interface types utils parser tbmodel solver eels ioutils mpi_solver

# Convert module names to source files
ORDERED_SRCS = $(foreach mod,$(MODULE_ORDER), $(SRCDIR)/$(mod).f90)

# Add main.f90 at the end
ORDERED_SRCS += $(SRCDIR)/main.f90

# Ordered object files
ORDERED_OBJS = $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(ORDERED_SRCS))

# Default target
all: $(BUILDDIR) $(TARGET)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Link target executable
$(TARGET): $(ORDERED_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compilation rule for .f90 files
# We use a pattern rule that places .o and .mod in build/
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Special handling for module files: ensure build directory exists
# The .mod files are generated during compilation of .o files
# This rule ensures the build directory exists before compilation
$(ORDERED_OBJS): | $(BUILDDIR)

# Clean up
clean:
	rm -rf $(BUILDDIR) $(TARGET)

# Phony targets
.PHONY: all clean

# Dependencies (simplified - in practice should be generated automatically)
# Rule to generate .mod files (they are created when compiling .o files)
$(BUILDDIR)/%.mod: $(SRCDIR)/%.f90
	@touch $@

# constants.f90 is independent
$(BUILDDIR)/constants.o: $(SRCDIR)/constants.f90

# interface.f90 is independent
$(BUILDDIR)/interface.o: $(SRCDIR)/interface.f90

# types.f90 depends on constants.mod
$(BUILDDIR)/types.o: $(SRCDIR)/types.f90 $(BUILDDIR)/constants.mod

# utils.f90 depends on constants.mod and interface.mod
$(BUILDDIR)/utils.o: $(SRCDIR)/utils.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/interface.mod

# For other modules, we'll assume they depend on constants and types at least
# In a real project, you might want to generate dependencies automatically
$(BUILDDIR)/parser.o: $(SRCDIR)/parser.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/types.mod $(BUILDDIR)/utils.mod
$(BUILDDIR)/tbmodel.o: $(SRCDIR)/tbmodel.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/types.mod $(BUILDDIR)/utils.mod
$(BUILDDIR)/eels.o: $(SRCDIR)/eels.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/interface.mod $(BUILDDIR)/types.mod $(BUILDDIR)/utils.mod $(BUILDDIR)/tbmodel.mod $(BUILDDIR)/solver.mod
$(BUILDDIR)/ioutils.o: $(SRCDIR)/ioutils.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/types.mod $(BUILDDIR)/utils.mod
$(BUILDDIR)/solver.o: $(SRCDIR)/solver.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/interface.mod $(BUILDDIR)/utils.mod $(BUILDDIR)/tbmodel.mod
$(BUILDDIR)/mpi_solver.o: $(SRCDIR)/mpi_solver.f90 $(BUILDDIR)/constants.mod $(BUILDDIR)/types.mod $(BUILDDIR)/utils.mod $(BUILDDIR)/tbmodel.mod $(BUILDDIR)/solver.mod $(BUILDDIR)/eels.mod $(BUILDDIR)/parser.mod $(BUILDDIR)/ioutils.mod

# main.f90 depends on almost everything
$(BUILDDIR)/main.o: $(SRCDIR)/main.f90 \
	$(BUILDDIR)/constants.mod \
	$(BUILDDIR)/types.mod \
	$(BUILDDIR)/utils.mod \
	$(BUILDDIR)/parser.mod \
	$(BUILDDIR)/tbmodel.mod \
	$(BUILDDIR)/solver.mod \
	$(BUILDDIR)/mpi_solver.mod \
	$(BUILDDIR)/eels.mod \
	$(BUILDDIR)/ioutils.mod

# Note: .mod files are generated when compiling the corresponding .o files
# The dependency on .mod files ensures correct compilation order
