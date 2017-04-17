SRC=./modules
TMP=./temp
OUTPUTS=./outputs
CODE=./code

# Deleting old mod files in temp folder
rm $TMP/*.mod

# Compilation of module files:
#	1st goes the SRC folder
#   2nd goes the folder with dependent modules
#   3rd goes module parts of program
modules=""  # Here we put paths to compiled module files
extension=".o"  # of compiled module file 
for file in $SRC/*.f90 $SRC/dependants/*.f90 $CODE/*.f90; do
	filename=${file##*/}  # Getting 'myModule.f90' string
	filename=${filename%.*}  # Getting 'myModule' string
	filename="$filename$extension"  # Get 'myModule.o' string
	# Compile module file. Jdir - specifies where to put .mod file
	gfortran -J$TMP -c $file -o $TMP/$filename  
	modules="$modules$TMP/$filename "
done

# Compile main file. Idir - specifies where to look for includes
gfortran -I$TMP -c main.f90 -o $TMP/main.o

# Collecting all compiled files together
gfortran -o $TMP/main.exe $TMP/main.o $modules 

ln -s $OUTPUTS/vw.dat fort.11
ln -s $OUTPUTS/uw.dat fort.12
ln -s $OUTPUTS/uv.dat fort.13
 
time $TMP/main.exe

rm fort.* # Deleting fortran link files
rm $TMP/main.exe
rm $TMP/*.o # Deleting compiled module files

# Sometimes .mod files can be generated automatically in the same 
# folder with modules. In order to keep workplace clean we delete them.
for file in $SRC/*.mod $SRC/dependants/*.mod $CODE/*.mod
# rm $SRC/*.mod
# rm $SRC/dependants/*.mod
# rm $CODE/*.mod
	rm file
done