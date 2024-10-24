CC = gcc
CFLAGS += -Wall -O3 -g -fsanitize=address
TARGET = gnwa
CWD = $(shell pwd)

# Directories
SRC_DIR = $(CWD)/src
BIN_DIR = $(CWD)/bin
OBJ_DIR = $(CWD)/obj
INC_DIR = $(CWD)/include
LIB_DIR = $(CWD)/lib

# Source and Object files
SRCS = $(SRC_DIR)/gnwa.c $(SRC_DIR)/main.c
OBJS = $(OBJ_DIR)/gnwa.o $(OBJ_DIR)/main.o

# Executable
EXE = $(BIN_DIR)/$(TARGET)

.PHONY:all clean tests lib

all:$(LIB_DIR)/libgnwa.a $(EXE)

lib:$(LIB_DIR)/libgnwa.a

$(EXE): $(OBJS)
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(OBJS) -o $@ $(CFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Static library creation
$(LIB_DIR)/libgnwa.a: $(OBJ_DIR)/gnwa.o
	@mkdir -p $(@D)
	ar rcs $@ $^

tests: all
	$(shell $(EXE))

clean:
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)