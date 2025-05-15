CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
SRCS = a.cpp b.cpp
TARGETS = $(basename $(SRCS))

$(TARGETS): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

check-%:
	python3 check.py ./$* test_data/

clean:
	rm -f $(TARGETS)