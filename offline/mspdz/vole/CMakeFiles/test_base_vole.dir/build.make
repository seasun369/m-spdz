# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ztx/SPDZ/mspdz

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ztx/SPDZ/mspdz

# Include any dependencies generated for this target.
include vole/CMakeFiles/test_base_vole.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include vole/CMakeFiles/test_base_vole.dir/compiler_depend.make

# Include the progress variables for this target.
include vole/CMakeFiles/test_base_vole.dir/progress.make

# Include the compile flags for this target's objects.
include vole/CMakeFiles/test_base_vole.dir/flags.make

vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o: vole/CMakeFiles/test_base_vole.dir/flags.make
vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o: vole/base_vole.cpp
vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o: vole/CMakeFiles/test_base_vole.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztx/SPDZ/mspdz/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o"
	cd /home/ztx/SPDZ/mspdz/vole && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o -MF CMakeFiles/test_base_vole.dir/base_vole.cpp.o.d -o CMakeFiles/test_base_vole.dir/base_vole.cpp.o -c /home/ztx/SPDZ/mspdz/vole/base_vole.cpp

vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_base_vole.dir/base_vole.cpp.i"
	cd /home/ztx/SPDZ/mspdz/vole && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ztx/SPDZ/mspdz/vole/base_vole.cpp > CMakeFiles/test_base_vole.dir/base_vole.cpp.i

vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_base_vole.dir/base_vole.cpp.s"
	cd /home/ztx/SPDZ/mspdz/vole && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ztx/SPDZ/mspdz/vole/base_vole.cpp -o CMakeFiles/test_base_vole.dir/base_vole.cpp.s

# Object files for target test_base_vole
test_base_vole_OBJECTS = \
"CMakeFiles/test_base_vole.dir/base_vole.cpp.o"

# External object files for target test_base_vole
test_base_vole_EXTERNAL_OBJECTS =

bin/test_base_vole: vole/CMakeFiles/test_base_vole.dir/base_vole.cpp.o
bin/test_base_vole: vole/CMakeFiles/test_base_vole.dir/build.make
bin/test_base_vole: /usr/local/lib/libemp-tool.so
bin/test_base_vole: /usr/lib/x86_64-linux-gnu/libssl.so
bin/test_base_vole: /usr/lib/x86_64-linux-gnu/libcrypto.so
bin/test_base_vole: vole/CMakeFiles/test_base_vole.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ztx/SPDZ/mspdz/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/test_base_vole"
	cd /home/ztx/SPDZ/mspdz/vole && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_base_vole.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
vole/CMakeFiles/test_base_vole.dir/build: bin/test_base_vole
.PHONY : vole/CMakeFiles/test_base_vole.dir/build

vole/CMakeFiles/test_base_vole.dir/clean:
	cd /home/ztx/SPDZ/mspdz/vole && $(CMAKE_COMMAND) -P CMakeFiles/test_base_vole.dir/cmake_clean.cmake
.PHONY : vole/CMakeFiles/test_base_vole.dir/clean

vole/CMakeFiles/test_base_vole.dir/depend:
	cd /home/ztx/SPDZ/mspdz && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ztx/SPDZ/mspdz /home/ztx/SPDZ/mspdz/vole /home/ztx/SPDZ/mspdz /home/ztx/SPDZ/mspdz/vole /home/ztx/SPDZ/mspdz/vole/CMakeFiles/test_base_vole.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : vole/CMakeFiles/test_base_vole.dir/depend

