
# paths
highs_path = $(HOME)/Documents/HiGHS
metis_path = $(HOME)/Documents/METIS
local_path = $(HOME)/local
openblas_path = /opt/homebrew/Cellar/openblas/0.3.26/

cpp_sources = \
	Analyse.cpp \
	Auxiliary.cpp \
	Factorise.cpp \
	Numeric.cpp \
	Symbolic.cpp \
	main.cpp

c_sources = \
	hsl_wrapper.c \
	DenseFact.c 

# binary file name
binary_name = fact

# object files directory
objdir = obj

# compilers
CPP = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++
CC = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang

# compiler flags
CPPFLAGS = -std=c++11 -O3 -g3 -Wno-deprecated #-fsanitize=address
CFLAGS = -O3 -g3 #-fsanitize=address

# includes and libraries
includes = -I$(highs_path)/build -I$(highs_path)/src/ -I$(metis_path)/include -I$(local_path)/include
libs_path = -L$(highs_path)/build/lib -L$(metis_path)/build/libmetis -L$(local_path)/lib
libs = -lhighs -lmetis -lGKlib -llapack -lblas -lhsl_ma86 -lhsl_ma87 -lhsl_ma97 -lhsl_ma57 -lhsl_mc68 -lfakemetis

# mess to link openmp on mac
OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp


# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(objdir)/%.o)
c_objects = $(c_sources:%.c=$(objdir)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(objdir)/%.d)


# link
$(binary_name): $(cpp_objects) $(c_objects)
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(OPENMP_FLAGS) $(libs_path) $(libs) $^ -o $@

# manage dependencies
-include $(dep)

# compile cpp
$(cpp_objects): $(objdir)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

# compile c
$(c_objects): $(objdir)/%.o: %.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) -MMD -c $(CFLAGS) $(includes) $< -o $@


.PHONY : clean
clean: 
	rm $(objdir)/*.o
	rm $(objdir)/*.d
	rm $(binary_name)
