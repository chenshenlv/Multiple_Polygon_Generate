# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex

# Include any dependencies generated for this target.
include CMakeFiles/multi_polygon.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/multi_polygon.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/multi_polygon.dir/flags.make

CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o: CMakeFiles/multi_polygon.dir/flags.make
CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o: multi_polygon.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o -c /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/multi_polygon.cpp

CMakeFiles/multi_polygon.dir/multi_polygon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multi_polygon.dir/multi_polygon.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/multi_polygon.cpp > CMakeFiles/multi_polygon.dir/multi_polygon.cpp.i

CMakeFiles/multi_polygon.dir/multi_polygon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multi_polygon.dir/multi_polygon.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/multi_polygon.cpp -o CMakeFiles/multi_polygon.dir/multi_polygon.cpp.s

CMakeFiles/multi_polygon.dir/polygon.cpp.o: CMakeFiles/multi_polygon.dir/flags.make
CMakeFiles/multi_polygon.dir/polygon.cpp.o: polygon.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/multi_polygon.dir/polygon.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/multi_polygon.dir/polygon.cpp.o -c /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/polygon.cpp

CMakeFiles/multi_polygon.dir/polygon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multi_polygon.dir/polygon.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/polygon.cpp > CMakeFiles/multi_polygon.dir/polygon.cpp.i

CMakeFiles/multi_polygon.dir/polygon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multi_polygon.dir/polygon.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/polygon.cpp -o CMakeFiles/multi_polygon.dir/polygon.cpp.s

CMakeFiles/multi_polygon.dir/mesh.cpp.o: CMakeFiles/multi_polygon.dir/flags.make
CMakeFiles/multi_polygon.dir/mesh.cpp.o: mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/multi_polygon.dir/mesh.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/multi_polygon.dir/mesh.cpp.o -c /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/mesh.cpp

CMakeFiles/multi_polygon.dir/mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multi_polygon.dir/mesh.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/mesh.cpp > CMakeFiles/multi_polygon.dir/mesh.cpp.i

CMakeFiles/multi_polygon.dir/mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multi_polygon.dir/mesh.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/mesh.cpp -o CMakeFiles/multi_polygon.dir/mesh.cpp.s

# Object files for target multi_polygon
multi_polygon_OBJECTS = \
"CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o" \
"CMakeFiles/multi_polygon.dir/polygon.cpp.o" \
"CMakeFiles/multi_polygon.dir/mesh.cpp.o"

# External object files for target multi_polygon
multi_polygon_EXTERNAL_OBJECTS =

multi_polygon: CMakeFiles/multi_polygon.dir/multi_polygon.cpp.o
multi_polygon: CMakeFiles/multi_polygon.dir/polygon.cpp.o
multi_polygon: CMakeFiles/multi_polygon.dir/mesh.cpp.o
multi_polygon: CMakeFiles/multi_polygon.dir/build.make
multi_polygon: /usr/local/lib/libCGAL.so.13.0.3
multi_polygon: /usr/lib/x86_64-linux-gnu/libmpfr.so
multi_polygon: /usr/lib/x86_64-linux-gnu/libgmp.so
multi_polygon: CMakeFiles/multi_polygon.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable multi_polygon"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/multi_polygon.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/multi_polygon.dir/build: multi_polygon

.PHONY : CMakeFiles/multi_polygon.dir/build

CMakeFiles/multi_polygon.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/multi_polygon.dir/cmake_clean.cmake
.PHONY : CMakeFiles/multi_polygon.dir/clean

CMakeFiles/multi_polygon.dir/depend:
	cd /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex /media/jason/UData/home/soundpadlab/jason/pipeline_nonconvex/CMakeFiles/multi_polygon.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/multi_polygon.dir/depend
