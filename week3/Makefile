CC = gcc

CFLAGS = -Wall -Wextra -std=c99 -O2 -I./src  # Added -O2 for optimization
LDFLAGS = -lm  # Link against math library

SRC_DIR = ./src
BUILD_DIR = ./build
RUNS_DIR = ./runs

SRC_FILES = $(SRC_DIR)/gravity.c
OBJ_FILES = $(BUILD_DIR)/gravity.o

EXEC = $(BUILD_DIR)/gravity_sim

all: $(EXEC)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJ_FILES)
	$(CC) $(OBJ_FILES) -o $(EXEC) $(LDFLAGS)

# Run the simulation
run: $(EXEC)
	./$(EXEC)

clean:
	rm -rf $(BUILD_DIR)/*.o $(EXEC)
	rm -rf $(RUNS_DIR)/*

.PHONY: all run clean
