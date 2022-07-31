all: botbin

botbin: .obj/main.o .obj/Graph.o .obj/Utility.o .obj/HashSet.o .obj/BitMap.o
	g++ .obj/main.o .obj/Graph.o .obj/Utility.o .obj/HashSet.o  .obj/BitMap.o Config.h -o botbin -w
	rm .obj/*.o

.obj/main.o: main.cpp
	g++ -c -O3 -o .obj/main.o main.cpp -w

.obj/Graph.o: Graph.cpp
	g++ -c -O3 -o .obj/Graph.o Graph.cpp -w

.obj/Utility.o: Utility.cpp
	g++ -c -O3 -o .obj/Utility.o Utility.cpp -w

.obj/HashSet.o: HashSet.cpp
	g++ -c -O3 -o .obj/HashSet.o HashSet.cpp -w

.obj/BitMap.o: BitMap.cpp
	g++ -c -O3 -o .obj/BitMap.o BitMap.cpp -w

clean:
	rm -rf *o .obj/
	mkdir .obj