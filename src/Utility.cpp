#include "Utility.h"

FILE* open_file(const char* file_name, const char* mode) {
	FILE* f = fopen(file_name, mode);
	if (f == NULL) {
		printf("Can not open file: %s\n", file_name);
		exit(1);
	}

	return f;
}

Utility::Utility() {

}

Utility::~Utility() {

}

void Utility::format_graph(string dest, string src) {
	format_graph(dest, src, 100, false);
}

void Utility::format_graph(string dest, string src, int percent, bool isSampleEdge) {
	printf("Format file: %s\n", src.c_str());
	percent = percent / 10 - 1;

	int memEdges = 100000000;
	dir_ = dest;

	FILE* fp = open_file(src.c_str(), "r");


	Edge* edges = new Edge[memEdges];


	//vertex_map_ = new int[MAX_ID];
	//memset(vertex_map_, -1, sizeof(int) * MAX_ID);

	// save every line read from txtFile
	char line[100];

	int u, v;

	unsigned long es = 0;

	int size = 0, num = 0, tmpFile = 0;
	//FILE* f_edges = open_file((dir_ + "dynamic_edges").c_str(), "wb");
	//fwrite(&DYNAMIC_EDGE_NUM, sizeof(int), 1, f_edges);
	//fclose(f_edges);
	printf("Separating: \n");
	while (fgets(line, 100, fp)) {

		if (line[0] < '0' || line[0] > '9')
			continue;

		sscanf(line, "%d\t%d", &u, &v);
		//printf("%d\t%d\n", u, v);
		//cout << u << " " << v << endl;
		//exit(0);
		//printf("%d\t%d\n", u, v);
		//while (true);
		if (u == v) continue;
		++es;
		if (isSampleEdge) {
			if (es % 10 > percent) continue;
		}
		else {
			if (u % 10 > percent || v % 10 > percent) continue;
		}

		u = get_vertex_id(u, num);
		v = get_vertex_id(v, num);


		edges[size].a = u;
		edges[size].b = v;
		++size;

		edges[size].a = v;
		edges[size].b = u;
		++size;
		//if (_ytb_len - size < 2*DYNAMIC_EDGE_NUM) {
		//	fwrite(&u, sizeof(int), 1, f_edges);
		//	fwrite(&v, sizeof(int), 1, f_edges);

		//	//fwrite(edges_list, sizeof(int), DYNAMIC_EDGE_NUM * 2, f_edges);

		//}
		if (size >= memEdges) {
			printf("load %d edges\n", size);
			save_tmp_edges(edges, size, tmpFile);
			size = 0;
			++tmpFile;
		}

	}
	printf("load %d edges\n", size);
	save_tmp_edges(edges, size, tmpFile);
	//fclose(f_edges);
	int* new2old = new int[num];
	//int cur = 0;
	//int i = 0;
	//while (cur < num) {
	//	if (vertex_map_[i] == -1) {
	//		++i;
	//		continue;
	//	}
	//	new2old[vertex_map_[i]] = i;
	//	++cur;
	//	++i;
	//}

	// write match file
	// printf("write match file\n");
	string mts = dest + "match.st";
	FILE* mch = fopen(mts.c_str(), "wb");
	fwrite(new2old, sizeof(int), num, mch);
	fclose(mch);
	
	delete[] new2old;
	//delete[] vertex_map_;
	delete[] edges;

	fclose(fp);

	merge(tmpFile + 1);
}

bool Utility::edge_compare(const Edge& e1, const Edge& e2) {
	if (e1.a < e2.a) {
		return true;
	}
	if (e1.a > e2.a) {
		return false;
	}
	return e1.b < e2.b;
}

// sort edges and save
void Utility::save_tmp_edges(Edge* edges, int size, int tmpFile) {
	printf("sort edges...\n");
	sort(edges, edges + size, edge_compare);

	string destDir = dir_ + "sort_edge_tmp";

	if (access(destDir.c_str(), 0) == -1) {
		// use for LINUX
#ifdef _LINUX_
		int flag = mkdir(destDir.c_str(), 0777);
#else
		int flag = mkdir(destDir.c_str());
#endif
		if (flag == 0) printf("Directory \"%s\" is created\n", destDir.c_str());
		else {
			printf("Create directory failed\n");
			exit(1);
		}
	}


	char fileName[200];
	sprintf(fileName, "%ssort_edge_tmp/edges_tmp_%d", dir_.c_str(), tmpFile);
	printf("Creating tmp_file_%d: %s\n", tmpFile, fileName);

	FILE* fo = fopen(fileName, "wb");
	for (int i = 0; i < size; ++i) {
		fwrite(edges + i, sizeof(Edge), 1, fo);
		// printf("edge[%d,%d]\n",edges[i].a,edges[i].b );
	}
	printf("------\n\n");
	fclose(fo);

}

// merge_sort edges from all edges_tmp files
void Utility::merge(int size) {
	FILE** frl = new FILE * [size];
	Edge* es = new Edge[size];
	string idxPath = dir_ + "graph.idx";
	string datPath = dir_ + "graph.dat";
	FILE* fIdx = fopen(idxPath.c_str(), "wb");
	FILE* fDat = fopen(datPath.c_str(), "wb");

	long pos;

	for (int i = 0; i < size; ++i) {
		char fileName[200];
		sprintf(fileName, "%ssort_edge_tmp/edges_tmp_%d", dir_.c_str(), i);

		frl[i] = fopen(fileName, "rb");
		// get first edge of all edges_tmp files
		fread(&es[i], sizeof(Edge), 1, frl[i]);
	}

	unsigned int edgeNum = 0;

	int minIndex, previousA = -1, previousB = -1;

	int i = 0;

	int maxDegree = 0;

	int degree = -1;
	// printf("start merge\n");
	int x = 0;
	int f = 0;
	while (is_merge_finished(es, size)) {
		minIndex = get_min_edge(es, size);
		int u = es[minIndex].a;
		int v = es[minIndex].b;

		if (u != previousA) {
			// u != previousA demonstrates that all previousA's neighbors have been writen
			if (u % 100000 == 0) {
				printf("[%d]\n", u);
			}

			if (degree != -1) {

				// write the vertex degree in .idx file
				fwrite(&degree, sizeof(int), 1, fIdx);
				edgeNum += degree;
				maxDegree = degree > maxDegree ? degree : maxDegree;
				// printf("max degree: %d\n",maxDegree );
			}

			// write the vertex beginning position in .dat file to .idx file
			pos = ftell(fDat);
			fwrite(&pos, sizeof(long), 1, fIdx);

			degree = 1;

			fwrite(&v, sizeof(int), 1, fDat);


		}
		else if (v != previousB) {

			fwrite(&v, sizeof(int), 1, fDat);
			++degree;

		}

		// if u==previousA & v==previousB, ignore edge(u,v) cause it is same as previous one.

		previousA = u;
		previousB = v;

		// replace es[minIndex] by picking up the first edge from file edges_tmp_minIndex
		if (!fread(&es[minIndex], sizeof(Edge), 1, frl[minIndex])) {
			es[minIndex].a = MAX_ID;
		}


	}

	fwrite(&degree, sizeof(int), 1, fIdx);
	edgeNum += degree;
	edgeNum /= 2;
	maxDegree = degree > maxDegree ? degree : maxDegree;

	// write the vertex num and max degree
	int vertexNum = previousA + 1;
	printf("vertex num: %d\n", vertexNum);
	printf("edge num: %d\n", edgeNum);

	string infoPath = dir_ + "graph.info";
	FILE* fInfo = fopen(infoPath.c_str(), "wb");
	fwrite(&vertexNum, sizeof(int), 1, fInfo);
	fwrite(&edgeNum, sizeof(unsigned int), 1, fInfo);
	fwrite(&maxDegree, sizeof(int), 1, fInfo);
	fclose(fInfo);

	for (int i = 0; i < size; ++i) {
		fclose(frl[i]);
	}

	fclose(fIdx);
	fclose(fDat);

	delete[] frl;
	delete[] es;
}

bool Utility::is_merge_finished(Edge* es, int size) {

	for (int i = 0; i < size; ++i) {
		if (es[i].a != MAX_ID) {
			return true;
		}
	}

	return false;
}

// get minimum edge from edge list
int Utility::get_min_edge(Edge* es, int size) {
	int min = 0;
	for (int i = 1; i < size; ++i) {
		if (es[i].a < es[min].a) {
			min = i;
		}
		else if (es[i].a > es[min].a) {
			continue;
		}
		else if (es[i].b < es[min].b) {
			min = i;
		}
	}

	return min;
}

// get final vertexID from map array
int Utility::get_vertex_id(int u, int& num) {
	if (vertex_map_.find(u) == vertex_map_.end()) {
		vertex_map_[u] = num++;
		return num;
	}
	return vertex_map_[u];
}
