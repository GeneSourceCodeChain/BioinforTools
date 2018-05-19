CXXFLAGS=-I. -std=c++14 -O2
LIBS=
OBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: alignment

alignment: Alignment.o main.o
	$(CXX) $^ $(LIBS) -o ${@}

clean:
	$(RM) $(OBJS) alignment

