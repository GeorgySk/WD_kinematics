kinematics.sublime-build - If you want to build a program right from Sublime 
	Text 3, place this file in folder with packages (Preferences -> Browse 
	packages). Change shell_cmd argument with whatever command you use to build
	the program in command line (in my case it was shell script execution)

SublimeLinter.sublime-settings - If you want SublimeLinter to see dependencies
	located in separate folders, then you should add the folder with them to
	gfortran calling flags. (see "gfortranmodern" -> "args") 