 STEPS USED TO CREATE PROJECT

1) Download and unzip CImg library

2) Add CImg.h file to project folder

3) (LINUX)
   // Be sure the following line is commented out "system("pause");"
   // This is only for windows. Also make sure that the following line
   // is uncommented "std::cin.ignore.get();" this is for pausing executable
   // so that when we run it it doesn't immediately close out. Simply run the
   // following make command or type out the commands in the makefile.

   make
   

   main: main.cpp
	
	   g++ -o capHW1 main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
	
	   chmod +x capHW1
	
	   ./capHW1

   (WINDOWS)
   
   // uncomment the line of code at the very bottom of the program
   // "system("pause");" Also be sure to comment out "std::cin.ignore.get();"
   // as this is only for linux.
   
   cd eclipse-workspace
   cd capHW1
   g++ -o capHW1.exe main.cpp -O2 -lgdi32

