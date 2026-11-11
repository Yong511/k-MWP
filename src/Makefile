CXX := g++
CXXFLAGS := -std=c++17 -O3 -march=native -pipe -Wall -Wextra -Wpedantic -DNDEBUG
LDFLAGS :=

TARGET := MWPlex
SRCS := main.cpp Graph.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
