CC = gcc -g
C++FLAGS = -std=c++14
GCOVFLAGS = -fprofile-arcs -ftest-coverage

OS := $(shell uname)

ifeq ($(OS), Darwin)
    FLAGSYSTEM = -lgtest -lgcov -lstdc++ -lm
else
    FLAGSYSTEM = -lgtest -lgcov -pthread -lstdc++ -lrt -lm
endif

all: s21_matrix_oop.a test

s21_matrix_oop.a:
	$(CC) $(C++FLAGS) -c s21_matrix_oop.cc
	ar rcs s21_matrix_oop.a *.o
	ranlib s21_matrix_oop.a

test:  
	$(CC) $(C++FLAGS) s21_matrix_oop.cc s21_test_matrix.cc -o s21_test $(FLAGSYSTEM)
	./s21_test

rebuild: clean all

clean: 
	rm -rf *.o s21_matrix_test s21_matrix_oop.a s21_test
cf:
	clang-format --style=Google -n $(shell find . -name '*.c') $(shell find . -name '*.h') $(shell find . -name '*.cpp')
	clang-format --style=Google -i $(shell find . -name '*.c') $(shell find . -name '*.h') $(shell find . -name '*.cpp')
