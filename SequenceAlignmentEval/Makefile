CXXFLAGS=-I. -std=c++14 -O2
LIBS=-lboost_filesystem -lboost_system -lboost_program_options -lboost_regex
OBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: evaluation

evaluation: main.o Evaluation.o
	$(CXX) $^ $(LIBS) -o ${@}

clean:
	$(RM) evaluation $(OBJS)
