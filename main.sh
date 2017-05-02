FOLDER=.
SRC=./modules
TMP=./temp
OUTPUTS=./outputs
CODE=./code

# Deleting old mod files in temp folder
rm $TMP/*.mod

# NOTE: next part should be uncommented when run-fortran is ready
# Creating sorted list of path/to/modules.f90 to be compiled 
# by run-fortran program
# modules=$(./run-fortran/run-fortran run -p .)
# echo "$modules"
# extension=".o" # Of compiled file
# for file in $modules; do
# 	filename=${file%.*}  # Getting 'path/to/module' string
# 	filename="$filename$extension"  # Adding '.o'
# 	# echo "Compiling $filename"
# 	# gfortran -I$TMP -J$TMP -c $file -o $filename -Wall -Wextra -Wconversion -fbounds-check
# done

# NOTE: This part should be deleted when run-fortran is ready
# Compilation of module files:
modules=""  # Here we put paths to compiled module files
extension=".o"  # of compiled module file 
for file in $SRC/*.f90 $CODE/*.f90; do
	filename=${file##*/}  # Getting 'myModule.f90' string
	filename=${filename%.*}  # Getting 'myModule' string
	filename="$filename$extension"  # Get 'myModule.o' string
# 	# Compile module file. Jdir - specifies where to put .mod file
	gfortran -J$TMP -c $file -o $TMP/$filename -Wall -Wextra -Wconversion -fbounds-check
	modules="$modules$TMP/$filename "
done

# Compile main file. Idir - specifies where to look for includes
gfortran -I$TMP -c main.f90 -o $TMP/main.o -Wall -Wextra -Wconversion -fbounds-check

# Collecting all compiled files together
gfortran -o main.exe $TMP/main.o $modules 
 
# Test execution
time $FOLDER/main.exe -s -l

rm $TMP/*.o # Deleting compiled module files

# Sometimes .mod and .o files can be generated automatically in the same 
# folder with modules. In order to keep workplace clean we delete them.
# Loop is needed in order to avoid errors in case of absence of files
for file in $FOLDER/* $SRC/* $CODE/*; do
	filename=${file##*/}  # Getting 'myModule.f90' string
	extension="${filename##*.}"  # Getting 'myModule' string
	if [ "$extension" = "mod" ]
	then
		rm $file
	fi
	if [ "$extension" = "o" ]
	then
		rm $file
	fi
done