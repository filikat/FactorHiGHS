
# paths
highs_path = $(HOME)/Documents/HiGHS
metis_path = $(HOME)/Documents/METIS
local_path = $(HOME)/local
blas_path  = /Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers


cpp_sources = \
		Analyse.cpp \
		Auxiliary.cpp \
		Factorise.cpp \
		Numeric.cpp \
		Symbolic.cpp \
		FormatHandler.cpp \
		FullFormatHandler.cpp \
		HybridPackedFormatHandler.cpp \
		HybridHybridFormatHandler.cpp \
		PackedPackedFormatHandler.cpp \
		SolveHandler.cpp \
		FullSolveHandler.cpp \
		PackedSolveHandler.cpp \
		HybridSolveHandler.cpp \
		DataCollector.cpp \
		DenseFactKernel.cpp \
		DenseFactFull.cpp \
		DenseFactHybrid.cpp \
		CallAndTimeBlas.cpp

# binary file name
binary_name = fact

# object files directory
objdir = obj

# compilers
CPP = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++

# compiler flags
CPPFLAGS = -std=c++11 -O3 -g3 -Wno-deprecated #-fsanitize=address

# includes and libraries
includes = -I$(highs_path)/build -I$(highs_path)/src/ -I$(metis_path)/include -I$(local_path)/include -I$(blas_path)
libs_path = -L$(highs_path)/build/lib -L$(metis_path)/build/libmetis -L$(local_path)/lib
libs = -lhighs -lmetis -lGKlib -lblas

# mess to link openmp on mac
# OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp


# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(objdir)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(objdir)/%.d)


# link
$(binary_name): $(cpp_objects)
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(OPENMP_FLAGS) $(libs_path) $(libs) $^ -o $@

# manage dependencies
-include $(dep)

# compile cpp
$(cpp_objects): $(objdir)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

.PHONY : clean
clean: 
	rm $(objdir)/*.o
	rm $(objdir)/*.d
	rm $(binary_name)
