main: main.cpp
	g++ -o capHW1 main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
	chmod +x capHW1
	./capHW1
