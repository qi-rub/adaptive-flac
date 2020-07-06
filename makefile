PROGRAM := encoder

SRC_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(patsubst src/%.cpp,build/%.o,$(SRC_FILES))

COMPILER_FLAGS = -std=c++17 -Wall -Ofast -MD \
		  	-I"src"

main: $(OBJ_FILES)
	g++ $^ -o $(PROGRAM)

build/%.o: src/%.cpp
	@echo 'Building $(patsubst $src/%,%,$<).'
	@mkdir -p ${@D}
	g++ $< -o $@ -c $(COMPILER_FLAGS)

debug: $(OBJ_FILES)
	g++ -g $^ -o $(PROGRAM) -fsanitize=address

profiler: $(OBJ_FILES)
	g++ -g $^ -o $(PROGRAM) -pg -O3
	./$(PROGRAM)
	gprof $(PROGRAM) > $(PROGRAM)_profile

run:
	make main
	./$(PROGRAM)

clean:
	rm -rf build
	rm -f $(PROGRAM)

