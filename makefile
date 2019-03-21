graph.o: graph.cpp
	c++ -o graph graph.cpp myLib.cc `root-config --cflags --glibs`
	c++ -o risol risol.cpp myLib.cc `root-config --cflags --glibs`

