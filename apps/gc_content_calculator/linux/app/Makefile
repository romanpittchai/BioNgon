CC = gcc
OBJ =  obj/main.o obj/DB.o obj/interface_functions.o obj/nucleotide_functions.o obj/sqlite3.o obj/util.o obj/about.o
FLAGS_OBJ = -c -Wall -Wextra `pkg-config --libs --cflags gtk+-3.0` -ldl -lpthread -o
FLAGS_BIN = `pkg-config --libs --cflags gtk+-3.0` -ldl -lpthread
GC = gc_content_calculator
SH = run.sh
DIR_O = obj
DIR_B = bin

all: $(GC)

$(GC): $(OBJ)
	mkdir -p $(GC)
	$(CC) -o $(DIR_B)/$(GC) $^ $(FLAGS_BIN)
	cp bin/$(GC) $(GC) && cp -r interface $(GC)
	chmod +x $(GC)/$(GC)
	@echo "./$(GC)" >> $(GC)/$(SH)
	chmod +x $(GC)/$(SH)
	clear
	@echo "Build successful! Check for the gc_content_calculator folder. It should contain an executable file and a folder with the interface. Have a nice job!"
	
clean:
	rm -rf $(DIR_O) $(DIR_B)
	clear
	@echo "Cleanup complete!"

$(DIR_O)/main.o: main.c interface_functions.h nucleotide_functions.h sqlite3.h structures.h util.h DB.h about.h
	mkdir -p bin && mkdir -p obj
	$(CC) $< $(FLAGS_OBJ) $@

$(DIR_O)/DB.o: DB.c structures.h util.h sqlite3.h
	$(CC) $< $(FLAGS_OBJ) $@

$(DIR_O)/interface_functions.o: interface_functions.c structures.h util.h 
	$(CC) $< $(FLAGS_OBJ) $@

$(DIR_O)/nucleotide_functions.o: nucleotide_functions.c structures.h util.h
	$(CC) $< $(FLAGS_OBJ) $@

$(DIR_O)/sqlite3.o: sqlite3.c
	$(CC) -c $< -Wall -Wextra -o $@

$(DIR_O)/util.o: util.c structures.h sqlite3.h
	$(CC) $< $(FLAGS_OBJ) $@

$(DIR_O)/about.o: about.c structures.h
	$(CC) $< $(FLAGS_OBJ) $@
	
