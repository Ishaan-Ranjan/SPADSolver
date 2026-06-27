CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude -Itests
LDFLAGS  =

SRCS = src/GaussianElimination.cpp src/interpolation.cpp
OBJS = $(SRCS:.cpp=.o)
LIB  = libspadsolver.a

TEST_INTERP_SRC = tests/test_interpolation.cpp
TEST_INTERP_BIN = test_interpolation

TEST_GE_SRC = tests/test_gaussian_elimination.cpp
TEST_GE_BIN = test_gaussian_elimination

.PHONY: all clean test test-interpolation test-gaussian-elimination

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_INTERP_BIN): $(TEST_INTERP_SRC) src/interpolation.o
	$(CXX) $(CXXFLAGS) $(TEST_INTERP_SRC) src/interpolation.o -o $@

$(TEST_GE_BIN): $(TEST_GE_SRC) src/GaussianElimination.o
	$(CXX) $(CXXFLAGS) $(TEST_GE_SRC) src/GaussianElimination.o -o $@

test: test-interpolation test-gaussian-elimination

test-interpolation: $(TEST_INTERP_BIN)
	./$(TEST_INTERP_BIN)

test-gaussian-elimination: $(TEST_GE_BIN)
	./$(TEST_GE_BIN)

clean:
	rm -f $(OBJS) $(LIB) $(TEST_INTERP_BIN) $(TEST_GE_BIN)
