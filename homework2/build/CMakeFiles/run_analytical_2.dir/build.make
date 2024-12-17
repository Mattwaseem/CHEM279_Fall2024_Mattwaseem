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
CMAKE_SOURCE_DIR = /Users/nematsfolder/chem279/homework2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nematsfolder/chem279/homework2/build

# Utility rule file for run_analytical_2.

# Include any custom commands dependencies for this target.
include CMakeFiles/run_analytical_2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/run_analytical_2.dir/progress.make

CMakeFiles/run_analytical_2: /Users/nematsfolder/chem279/homework2/bin/homework2_exec
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/Users/nematsfolder/chem279/homework2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Running homework2_exec with analytical input file 2.txt"
	/Users/nematsfolder/chem279/homework2/bin/homework2_exec /Users/nematsfolder/chem279/homework2/sample_input/analytical/2.txt > /Users/nematsfolder/chem279/homework2/calculated_output_p2/2-output.txt

run_analytical_2: CMakeFiles/run_analytical_2
run_analytical_2: CMakeFiles/run_analytical_2.dir/build.make
.PHONY : run_analytical_2

# Rule to build all files generated by this target.
CMakeFiles/run_analytical_2.dir/build: run_analytical_2
.PHONY : CMakeFiles/run_analytical_2.dir/build

CMakeFiles/run_analytical_2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run_analytical_2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run_analytical_2.dir/clean

CMakeFiles/run_analytical_2.dir/depend:
	cd /Users/nematsfolder/chem279/homework2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nematsfolder/chem279/homework2 /Users/nematsfolder/chem279/homework2 /Users/nematsfolder/chem279/homework2/build /Users/nematsfolder/chem279/homework2/build /Users/nematsfolder/chem279/homework2/build/CMakeFiles/run_analytical_2.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/run_analytical_2.dir/depend

