# Makefile for compiling the fvecs project

# Compiler
CC = gcc

# Directories
SRC_DIR = src
OBJ_DIR = obj
HEADER_DIR = headers
TEST_DIR = tests/unit
MOCK_DIR = tests/mocks

# List all .c source files in the source directory
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)

# Test files
TEST_FILES := $(wildcard $(TEST_DIR)/*.c)
TEST_OBJ_FILES := $(TEST_FILES:$(TEST_DIR)/%.c=$(OBJ_DIR)/%.o)

# Executables
TEST_TARGET = run_tests
EXEC = exclude

# # List of files to exclude
EXCLUDE_FILES := $(SRC_DIR)/groundtruth.c $(SRC_DIR)/stitchedVamana.c $(SRC_DIR)/recallStitchedVamana.c $(SRC_DIR)/recallFilteredVamana.c

# # Filter out the excluded files
SRC_FILES := $(filter-out $(EXCLUDE_FILES), $(SRC_FILES))

OBJ_FILES := $(SRC_FILES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Executable name

# Compiler flags (for normal build)
CFLAGS = -mavx2 -O3 -I$(HEADER_DIR) -Wall -Wextra -pthread

# Compiler flags for debug build (-g for debug info, -O0 to disable optimization)
DEBUG_FLAGS = -I$(HEADER_DIR) -Wall -Wextra -g -O0 -pthread

# Linker flags (include -lm to link the math library)
LDFLAGS = -lcheck -lsubunit -lm -pthread

# Default rule (normal build)
all: $(EXEC)


#Rule to link the object files into an executable
$(EXEC): $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

# Rule to compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR) # Create obj directory if it doesn't exist
	$(CC) $(CFLAGS) -c $< -o $@

recalls: $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/recallStitchedVamana.c 
	$(CC) -o recallSticthedVamana $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/recallStitchedVamana.c $(CFLAGS)

recallf: $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/recallFilteredVamana.c
	$(CC) -o recallFilteredVamana $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/recallFilteredVamana.c $(CFLAGS)

stitched: $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/stitchedVamana.c
	$(CC) -o stitchedVamana $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/stitchedVamana.c $(CFLAGS)

#create the  filtered vamana
createFiltered: $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/CreateFilteredVamana.c
	$(CC) -mavx2 -O3 -o createFilteredVamana $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(SRC_DIR)/CreateFilteredVamana.c $(CFLAGS)

# Rule to compile and link the test executable
smth: smth.c $(SRC_DIR)/dataset.c
	$(CC) -o smth smth.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(CFLAGS) $(LDFLAGS)

# Rule to compile and link the test executable
test: $(TEST_OBJ_FILES) $(OBJ_FILES)
	@mkdir -p $(OBJ_DIR) # Ensure obj directory exists before compiling test files
	$(CC) $(CFLAGS) $(TEST_OBJ_FILES) $(OBJ_FILES) $(LDFLAGS) -o $(TEST_TARGET)
	./$(TEST_TARGET)

# Compile test files into object files
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(OBJ_DIR) # Ensure obj directory exists for test files
	$(CC) $(CFLAGS) -c $< -o $@

groundtruth: $(SRC_DIR)/groundtruth.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c
	$(CC) -o groundtruth $(SRC_DIR)/groundtruth.c $(SRC_DIR)/dataset.c $(SRC_DIR)/graph.c $(CFLAGS) $(LDFLAGS)

# Clean rule to remove object files and executables
clean:
	rm -rf $(OBJ_DIR)/*.o recallSticthedVamana recallFilteredVamana smth groundtruth output.txt stitchedVamana run_tests

# Phony targets
.PHONY: all debug clean test
