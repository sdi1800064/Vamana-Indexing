# Makefile for compiling the fvecs project

# Compiler
CC = gcc

# Directories
SRC_DIR = src
OBJ_DIR = obj
HEADER_DIR = headers

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.c)

# Object files
OBJ_FILES = $(SRC_FILES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Executable name
EXEC = my_program

# Compiler flags
CFLAGS = -I$(HEADER_DIR) -Wall -Wextra

# Default rule
all: $(EXEC)

# Rule to link the object files into an executable
$(EXEC): $(OBJ_FILES)
	$(CC) -o $@ $^

# Rule to compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR) # Create obj directory if it doesn't exist
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule to remove object files and executable
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC) output.txt

# Phony targets
.PHONY: all clean
