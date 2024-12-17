# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/nematsfolder/chem279/homework3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nematsfolder/chem279/homework3/build

# Utility rule file for run.

# Include any custom commands dependencies for this target.
include CMakeFiles/run.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/run.dir/progress.make

CMakeFiles/run: /Users/nematsfolder/chem279/homework3/bin/MoleculeSetup
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Running MoleculeSetup with C2H2.txt, C2H4.txt, and H2.txt and generating output for each"
	/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup /Users/nematsfolder/chem279/homework3/sample_input/C2H2.txt > /Users/nematsfolder/chem279/homework3/calculated_output/C2H2_output.txt
	/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup /Users/nematsfolder/chem279/homework3/sample_input/C2H4.txt > /Users/nematsfolder/chem279/homework3/calculated_output/C2H4_output.txt
	/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup /Users/nematsfolder/chem279/homework3/sample_input/H2.txt > /Users/nematsfolder/chem279/homework3/calculated_output/H2_output.txt

run: CMakeFiles/run
run: CMakeFiles/run.dir/build.make
.PHONY : run

# Rule to build all files generated by this target.
CMakeFiles/run.dir/build: run
.PHONY : CMakeFiles/run.dir/build

CMakeFiles/run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run.dir/clean

CMakeFiles/run.dir/depend:
	cd /Users/nematsfolder/chem279/homework3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nematsfolder/chem279/homework3 /Users/nematsfolder/chem279/homework3 /Users/nematsfolder/chem279/homework3/build /Users/nematsfolder/chem279/homework3/build /Users/nematsfolder/chem279/homework3/build/CMakeFiles/run.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/run.dir/depend

