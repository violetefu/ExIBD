# Set installation path and src/util paths
SF = src
EXEDIR = bin
REF = reference/FDRref.txt
Beagle = beagle_v3.3.2/beagle.jar

# Flags
CFLAGS = -Wall -O3 -fopenmp
PFLAGS = -F

# Libaries
LIB = -fopenmp

# Set the compiler
CC = g++
PY = pyinstaller

#list all programs to be installed by make
CPROGS = ExIBD ExIBD_Refined
PPROGS = ExIBD_Candidate ExIBD_Filtered

all: introduce $(CPROGS) $(PPROGS) install
	@echo "Finished compiling"
	@echo ""
	
Conly: introduce ExIBD_Refined install_Conly
	@echo "Finished compiling"
	@echo ""

introduce:
	@echo "Building ExIBD program files"

install: 
	@echo "Installing ExIBD binaries in $(EXEDIR)"
	@mkdir -p $(EXEDIR)
	@cp $(REF) $(EXEDIR)
	@cp $(Beagle) $(EXEDIR)
	@mv $(CPROGS) $(EXEDIR)
	@mv ./dist/* $(EXEDIR)
	@echo "Finished installation."
	@echo ""

install_Conly: 
	@echo "Installing ExIBD C binaries in $(EXEDIR)"
	@mkdir -p $(EXEDIR)
	@cp $(REF) $(EXEDIR)
	@cp $(Beagle) $(EXEDIR)
	@cp $(SF)/*.py $(EXEDIR)
	@mv ExIBD_Refined $(EXEDIR)
	@echo "Finished installation."
	@echo ""

clean: 
	@echo "Removing object files"
	@rm -f *.o
	@rm -f $(SF)/*.o
	@rm -r ./dist
	@rm -r ./build
	@rm -f *.spec
	@echo "Finished cleanup."
	@echo ""

clean_Conly: 
	@echo "Removing object files"
	@rm -f *.o
	@rm -f $(SF)/*.o
	@echo "Finished cleanup."
	@echo ""

ExIBD : ExIBD.o
	$(CC) $(CFLAGS) $(LIB) ExIBD.o -o ExIBD

ExIBD.o : $(SF)/ExIBD.cpp
	$(CC) $(CFLAGS) -c $<
	
ExIBD_Candidate : 
	$(PY) $(PFLAGS) $(SF)/ExIBD_Candidate.py

ExIBD_Refined : ExIBD_Refined.o
	$(CC) $(CFLAGS) $(LIB) ExIBD_Refined.o -o ExIBD_Refined

ExIBD_Refined.o : $(SF)/ExIBD_Refined.cpp
	$(CC) $(CFLAGS) -c $<
	
ExIBD_Filtered : 
	$(PY) $(PFLAGS) $(SF)/ExIBD_Filtered.py

