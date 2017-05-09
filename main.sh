#!/usr/bin/env bash
# By using the above line (shebang) we can launch this script by
# "./main.sh". In case if you don't have permission to do it, run: 
# "chmod 744 ./main.sh". But you can always run it by ". main.sh"


# This shell script is designed to compile/build multimodule fortran
# projects. It uses a script "run-fortran" written in Python, so 
# if you want to run it, you will have to install Python3.5 or higher 
# and "click" package for it. If you don't have permission to install
# "click" automatically (by pip) you can download the source at 
# https://pypi.python.org/pypi/click and install it manually by running:
# "python3.6 setup.py install --user".
# P.S: If you don't want to depend on Python and its own dependencies, 
# you can try to make executable from run-fortran.py by "pyinstaller",
# but be sure to do it on a system with the lowest GLIBC possible,
# as building the executable on a system with too high GLIBC will 
# result in a bad portability to older systems.


# Defining variables for paths of the project:
FOLDER=.  # Current folder
SRC=./modules  # Folder for the most common modules like "math", "I/O"
TMP=./temp
OUTPUTS=./outputs  # Data for plotting
CODE=./code  # Modules relevant only to this particular project


# As .mod files are compiler-dependent, everytime we build the program 
# we should delete them and compile new ones:
rm $TMP/*.mod  # Deleting old .mod files in temp folder
rm $FOLDER/*.exe  # Deleting old executable


# "run-fortran" creates an ordered list of fortran files to be compiled:
files_to_be_compiled=$(python3.6 run-fortran.py run --path $FOLDER)


# Determining compiler options of gfortran. Some other useful ones:
#   -Wall - generate warnings about many common sources of bugs
#   -Wextra - warnings about even more potential problems
#   -Wconversion - warnings about implicit conversions
#   -fmax-errors=n - max number of errors on output, 0 - no limit
#   -pedantic - warn about not fully supported f95 features
#   -std=f95 - same as pedantic but gives errors instead of warnings
#   -g or -g3 - generates debugging info to use it by GDB
#   -fbacktrace - tells what function was called last during the crash
#   -fbounds-check - check if array indexes are within the bounds
#   -ffpe-trap=zero,overflow,underflow - kill program if divide by 0 etc
#   -Olevel - code optimization, level=0(no opt.),1,2,3(max opt.)
# For the full list see: 
#   https://gcc.gnu.org/onlinedocs/gfortran/Invoking-GNU-Fortran.html
Fortran_compiler_options=$(echo "-Wall"\
                                "-Wextra"\
                                "-Wconversion"\
                                "-fbounds-check")


# Compiling all files:
compiled_object_files=""  # Here we will put paths to .o files (ordered)
for file in $files_to_be_compiled; do
    filename=${file##*/}  # Leaving only the name of the file
	filename=${filename%.*}  # Deleting .f90 from the name
    filename="$filename.o"  # Adding ".o" - extension of compiled files
	# Idir - specifies where to look for includes
    # Jdir - specifies where to put .mod file
    # -c - compile the following file
    # -o - where to put compiled file and how to name it
    gfortran -I$TMP -J$TMP -c $file -o $TMP/$filename $Fortran_compiler_options
    compiled_object_files="$compiled_object_files $TMP/$filename"
done


# Building executable:
gfortran -o main.exe $compiled_object_files

 
# Test execution:
time $FOLDER/main.exe -s


# Deleting unnecessary files:
rm $TMP/*.o  # Deleting compiled object files, they always appear
# Sometimes linter can autogenerate garbage files. Deleting them:
for file in $FOLDER/* $SRC/* $CODE/*; do
	filename=${file##*/}  # Getting file name
	extension="${filename##*.}"  # Getting extension of the file
	if [ "$extension" = "mod" ]
	then
		rm $file
	fi
	if [ "$extension" = "o" ]
	then
		rm $file
	fi
done
