CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
TARGET = a
SRC = a.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)