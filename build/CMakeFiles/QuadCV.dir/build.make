# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.1

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build

# Include any dependencies generated for this target.
include CMakeFiles/QuadCV.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/QuadCV.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/QuadCV.dir/flags.make

CMakeFiles/QuadCV.dir/main.cpp.o: CMakeFiles/QuadCV.dir/flags.make
CMakeFiles/QuadCV.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/QuadCV.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/QuadCV.dir/main.cpp.o -c /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/main.cpp

CMakeFiles/QuadCV.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/QuadCV.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/main.cpp > CMakeFiles/QuadCV.dir/main.cpp.i

CMakeFiles/QuadCV.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/QuadCV.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/main.cpp -o CMakeFiles/QuadCV.dir/main.cpp.s

CMakeFiles/QuadCV.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/QuadCV.dir/main.cpp.o.requires

CMakeFiles/QuadCV.dir/main.cpp.o.provides: CMakeFiles/QuadCV.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/QuadCV.dir/build.make CMakeFiles/QuadCV.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/QuadCV.dir/main.cpp.o.provides

CMakeFiles/QuadCV.dir/main.cpp.o.provides.build: CMakeFiles/QuadCV.dir/main.cpp.o

CMakeFiles/QuadCV.dir/calib.cpp.o: CMakeFiles/QuadCV.dir/flags.make
CMakeFiles/QuadCV.dir/calib.cpp.o: ../calib.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/QuadCV.dir/calib.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/QuadCV.dir/calib.cpp.o -c /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/calib.cpp

CMakeFiles/QuadCV.dir/calib.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/QuadCV.dir/calib.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/calib.cpp > CMakeFiles/QuadCV.dir/calib.cpp.i

CMakeFiles/QuadCV.dir/calib.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/QuadCV.dir/calib.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/calib.cpp -o CMakeFiles/QuadCV.dir/calib.cpp.s

CMakeFiles/QuadCV.dir/calib.cpp.o.requires:
.PHONY : CMakeFiles/QuadCV.dir/calib.cpp.o.requires

CMakeFiles/QuadCV.dir/calib.cpp.o.provides: CMakeFiles/QuadCV.dir/calib.cpp.o.requires
	$(MAKE) -f CMakeFiles/QuadCV.dir/build.make CMakeFiles/QuadCV.dir/calib.cpp.o.provides.build
.PHONY : CMakeFiles/QuadCV.dir/calib.cpp.o.provides

CMakeFiles/QuadCV.dir/calib.cpp.o.provides.build: CMakeFiles/QuadCV.dir/calib.cpp.o

# Object files for target QuadCV
QuadCV_OBJECTS = \
"CMakeFiles/QuadCV.dir/main.cpp.o" \
"CMakeFiles/QuadCV.dir/calib.cpp.o"

# External object files for target QuadCV
QuadCV_EXTERNAL_OBJECTS =

QuadCV: CMakeFiles/QuadCV.dir/main.cpp.o
QuadCV: CMakeFiles/QuadCV.dir/calib.cpp.o
QuadCV: CMakeFiles/QuadCV.dir/build.make
QuadCV: /usr/lib/libopencv_videostab.so.2.4.10
QuadCV: /usr/lib/libopencv_ts.a
QuadCV: /usr/lib/libopencv_superres.so.2.4.10
QuadCV: /usr/lib/libopencv_stitching.so.2.4.10
QuadCV: /usr/lib/libopencv_contrib.so.2.4.10
QuadCV: /lib64/libGLU.so
QuadCV: /lib64/libGL.so
QuadCV: /lib64/libSM.so
QuadCV: /lib64/libICE.so
QuadCV: /lib64/libX11.so
QuadCV: /lib64/libXext.so
QuadCV: /usr/lib/libopencv_nonfree.so.2.4.10
QuadCV: /usr/lib/libopencv_ocl.so.2.4.10
QuadCV: /usr/lib/libopencv_gpu.so.2.4.10
QuadCV: /usr/lib/libopencv_photo.so.2.4.10
QuadCV: /usr/lib/libopencv_objdetect.so.2.4.10
QuadCV: /usr/lib/libopencv_legacy.so.2.4.10
QuadCV: /usr/lib/libopencv_video.so.2.4.10
QuadCV: /usr/lib/libopencv_ml.so.2.4.10
QuadCV: /usr/lib/libopencv_calib3d.so.2.4.10
QuadCV: /usr/lib/libopencv_features2d.so.2.4.10
QuadCV: /usr/lib/libopencv_highgui.so.2.4.10
QuadCV: /usr/lib/libopencv_imgproc.so.2.4.10
QuadCV: /usr/lib/libopencv_flann.so.2.4.10
QuadCV: /usr/lib/libopencv_core.so.2.4.10
QuadCV: CMakeFiles/QuadCV.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable QuadCV"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/QuadCV.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/QuadCV.dir/build: QuadCV
.PHONY : CMakeFiles/QuadCV.dir/build

CMakeFiles/QuadCV.dir/requires: CMakeFiles/QuadCV.dir/main.cpp.o.requires
CMakeFiles/QuadCV.dir/requires: CMakeFiles/QuadCV.dir/calib.cpp.o.requires
.PHONY : CMakeFiles/QuadCV.dir/requires

CMakeFiles/QuadCV.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/QuadCV.dir/cmake_clean.cmake
.PHONY : CMakeFiles/QuadCV.dir/clean

CMakeFiles/QuadCV.dir/depend:
	cd /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1 /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1 /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build /home/andi/Studium/Master/1.Semester/PDS/comvis-10-1/build/CMakeFiles/QuadCV.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/QuadCV.dir/depend
