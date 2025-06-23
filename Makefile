CXX = g++
DEBUG = 0
ifeq ($(DEBUG), 1)
    CXXFLAGS = -std=c++11 -Wall -Wextra -fsanitize=address -g -O0 -DDEBUG
else
    CXXFLAGS = -std=c++11 -Wall -Wextra -O2
endif


SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
TARGETS = $(basename $(SRCS))
TARGETS_NAMES = $(notdir $(TARGETS))

default:
	@echo "Use make check-<target>"

$(TARGETS_NAMES): %: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $(SRC_DIR)/$@ $<

check-%: %
	python3 scripts/check.py $(SRC_DIR)/$* tests/

clean:
	rm -f $(TARGETS)