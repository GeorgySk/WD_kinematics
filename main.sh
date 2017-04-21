FOLDER=.
SRC=./modules
TMP=./temp
OUTPUTS=./outputs
CODE=./code

# NOTE: For some reason if I have my old .mod files from different
# compiler, I will get errors saying that they are from different 
# version. So this next command doesn't work? I should fix it, if not - 
# I have to clean temp folder by my own.
# Deleting old mod files in temp folder
rm $TMP/*.mod

# TODO: figure out what to do with logic here
# Compilation of module files:
#	1st goes the SRC folder
#   2nd goes the folder with dependent modules
#   3rd goes module parts of program
modules=""  # Here we put paths to compiled module files
extension=".o"  # of compiled module file 
for file in $SRC/*.f90 $SRC/dependants/*.f90 $CODE/*.f90 $CODE/dependants/*.f90; do
	filename=${file##*/}  # Getting 'myModule.f90' string
	filename=${filename%.*}  # Getting 'myModule' string
	filename="$filename$extension"  # Get 'myModule.o' string
	# Compile module file. Jdir - specifies where to put .mod file
	gfortran -J$TMP -c $file -o $TMP/$filename -Wall -Wextra -Wconversion -fbounds-check
	modules="$modules$TMP/$filename "
done

# Compile main file. Idir - specifies where to look for includes
gfortran -I$TMP -c main.f90 -o $TMP/main.o -Wall -Wextra -Wconversion -fbounds-check

# Collecting all compiled files together
gfortran -o main.exe $TMP/main.o $modules 
 
# Test execution
time $FOLDER/main.exe -o -l

rm $TMP/*.o # Deleting compiled module files

# Sometimes .mod files can be generated automatically in the same 
# folder with modules. In order to keep workplace clean we delete them.
# Loop is needed in order to avoid errors in case of absence of files
for file in $SRC/* $SRC/dependants/* $CODE/* $CODE/dependants/*; do
	filename=${file##*/}  # Getting 'myModule.f90' string
	extension="${filename##*.}"  # Getting 'myModule' string
	if [ "$extension" = "mod" ]
	then
		rm $file
	fi
done