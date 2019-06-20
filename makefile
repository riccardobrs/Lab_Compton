graph.o: graph.cpp
	c++ -o graph graph.cpp myLib.cc `root-config --cflags --glibs`
	c++ -o fit fit.cpp myLib.cc `root-config --cflags --glibs`
	c++ -o compton compton.cpp myLib.cc `root-config --cflags --glibs`
	c++ -o eff eff.cpp

