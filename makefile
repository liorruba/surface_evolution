EXE = regolit_main
SRC_DIR = src
OBJ_DIR = obj
UNAME := $(shell uname)

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CPPFLAGS += -std=c99
CFLAGS += -Wall
LDLIBSLINUX += -lm -lrt
LDLIBSMAC += -lm

all: $(EXE)

$(EXE): $(OBJ)
ifeq ($(UNAME),Linux)
	$(CC) $(LDFLAGS) $^ $(LDLIBSLINUX) -o $@
endif

ifeq ($(UNAME),Darwin)
	$(CC) $(LDFLAGS) $^ $(LDLIBSMAC) -o $@
endif
	# Create dirs:
	mkdir log
	mkdir output

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJ) regolit_main
	rm -r log
	rm -r output

.PHONY: all clean
