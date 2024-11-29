# Makefile for compiling the fvecs project

# Compiler
CC = gcc

# Directories
SRC_DIR = src
OBJ_DIR = obj
HEADER_DIR = headers

# List all .c source files in the source directory
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)

# List of files to exclude
EXCLUDE_FILES := $(SRC_DIR)/groundtruth.c

# Filter out the excluded files
INCLUDE_FILES := $(filter-out $(EXCLUDE_FILES), $(SRC_FILES))

# Generate object files list from the filtered source files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(INCLUDE_FILES))

# Executable name
EXEC = my_program

# Compiler flags (for normal build)
CFLAGS = -I$(HEADER_DIR) -Wall -Wextra

# Compiler flags for debug build (-g for debug info, -O0 to disable optimization)
DEBUG_FLAGS = -I$(HEADER_DIR) -Wall -Wextra -g -O0

# Linker flags (include -lm to link the math library)
LDFLAGS = -lm

# Default rule (normal build)
all: $(EXEC)

# Rule for debug build
debug: CFLAGS = $(DEBUG_FLAGS)
debug: $(EXEC)

# Rule to link the object files into an executable
$(EXEC): $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

# Rule to compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR) # Create obj directory if it doesn't exist
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to compile and link the test executable
test: test.c $(SRC_DIR)/dataset.c
	$(CC) -o test test.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(CFLAGS) $(LDFLAGS)

groundtruth: $(SRC_DIR)/groundtruth.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c
	$(CC) -o groundtruth $(SRC_DIR)/groundtruth.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(CFLAGS) $(LDFLAGS)

# Clean rule to remove object files and executables
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC) test output.txt

# Phony targets
.PHONY: all debug clean test
