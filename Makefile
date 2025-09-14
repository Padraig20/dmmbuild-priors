CXXFLAGS = -Wall -O3 -g

all: dummer-build dummer-build-lib

dummer-build: dummer-build.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer-build.cc

dummer-build-lib: dummer-build.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -fPIC -shared dummer-build.cc -o dummer-build.so

clean:
	rm -f *.o
