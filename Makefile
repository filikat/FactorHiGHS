
# paths
highs_path = $(HOME)/Documents/HiGHS
metis_path = $(HOME)/Documents/METIS
local_path = $(HOME)/local

cpp_sources = \
	Analyze.cpp \
	Aux_analyze.cpp \
	test.cpp

# binary file name
binary_name = fact

# object files directory
objdir = obj

# compilers
CPP = clang++

# compiler flags
CPPFLAGS = -std=c++11 -g

# includes and libraries
includes = -I$(highs_path)/build -I$(highs_path)/src/ -I$(metis_path)/include -I$(local_path)/include
libs_path = -L$(highs_path)/build/lib -L$(metis_path)/build/libmetis -L$(local_path)/lib
libs = -lhighs -lmetis -lGKlib

# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(objdir)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(objdir)/%.d)


# link
$(binary_name): $(cpp_objects)
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(libs_path) $(libs) $^ -o $@

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
	rm $(binary_name)
	rm $(objdir)/*.d
