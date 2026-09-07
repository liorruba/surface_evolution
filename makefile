# Build type: "release" (default, optimized) or "debug" (no optimization, debug symbols,
# AddressSanitizer + UndefinedBehaviorSanitizer). Each type has its own object directory and
# binary, so `make` and `make debug` can coexist without a `make clean` in between.
BUILD_TYPE ?= release

CXX      := g++
WARNINGS := -pedantic-errors -Wall -Wextra -Werror -std=c++17
BUILD    := ./build
APP_DIR  := $(BUILD)/apps
INCLUDE  := -Iinclude/
LDFLAGS  := -lm

ifeq ($(BUILD_TYPE),debug)
CXXFLAGS := $(WARNINGS) -g -O0 -DDEBUG -fsanitize=address,undefined -fno-omit-frame-pointer
OBJ_DIR  := $(BUILD)/obj/debug
TARGET   := regolit_main_debug.run
else
CXXFLAGS := $(WARNINGS) -O2
OBJ_DIR  := $(BUILD)/obj/release
TARGET   := regolit_main.run
endif

SRC 	 := $(wildcard src/*.cpp)
OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
HEADERS  := $(wildcard include/*.hpp)

all: main

main: build $(APP_DIR)/$(TARGET)
$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LDFLAGS)

.PHONY: all main build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug:
	$(MAKE) BUILD_TYPE=debug main

release:
	$(MAKE) BUILD_TYPE=release main

clean:
	-@rm -rvf $(BUILD)/obj
	-@rm -rvf $(APP_DIR)/*
	-@rm -f "log/log.txt"
	-@rm -rf "output/"*
