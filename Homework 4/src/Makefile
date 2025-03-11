CXX      := g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -I../include
EXT_OBJS := $(wildcard ../obj/*.o)
EXEC     := ../bin/hw4
PARALLEL_EXEC := ../bin/hw4_parallel
SRCS     := ExampleFunction.cpp GlobalPlacer.cpp main.cpp
PARALLEL_SRCS := ExampleFunction.cpp GlobalPlacer.cpp main_parallel.cpp
OBJS     := $(SRCS:.cpp=.o)
PARALLEL_OBJS := $(PARALLEL_SRCS:.cpp=.o)

all: $(EXEC) $(PARALLEL_EXEC)

$(EXEC): $(OBJS) $(EXT_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(PARALLEL_EXEC): $(PARALLEL_OBJS) $(EXT_OBJS)
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

main_parallel.o: main_parallel.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

ifeq (test, $(firstword $(MAKECMDGOALS)))
  TESTCASE := $(word 2, $(MAKECMDGOALS))
  $(eval $(TESTCASE):;@:)
endif

test: $(EXEC)
	@echo test on ../testcase/$(TESTCASE)/$(TESTCASE).aux
	./$(EXEC) ../testcase/$(TESTCASE)/$(TESTCASE).aux ../output/$(TESTCASE).gp.pl
	./../verifier/verify ../testcase/$(TESTCASE)/$(TESTCASE).aux ../output/$(TESTCASE).gp.pl

EXAMPLE := ../bin/example

example: $(EXAMPLE)

$(EXAMPLE): Example.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(EXT_OBJS)

clean:
	rm -f $(EXEC) $(PARALLEL_EXEC) $(OBJS) $(PARALLEL_OBJS) $(EXAMPLE)

.PHONY: all test example clean
