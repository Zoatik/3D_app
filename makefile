main: main.o M_matrices.o Color_picker.o
	g++ obj/main.o obj/M_matrices.o obj/Color_picker.o -o main

main.o: main.cpp M_matrices.h Color_picker.h
	g++ -c main.cpp -o obj/main.o

M_matrices.o: M_matrices.cpp M_matrices.h
	g++ -c M_matrices.cpp -o obj/M_matrices.o

Color_picker.o: Color_picker.cpp Color_picker.h
	g++ -c Color_picker.cpp -o obj/Color_picker.o

clean:
	del obj\*.o