#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>

#define _LINUX_
#define _DEBUG_

// #define _CN_HASH_SET_
// #define _CN_SORT_MERGE_
#define _CN_BITMAP_

const int THER = 1;
const int THER1 = 1000000;

const int DYNAMIC_EDGE_NUM = 10000;

#ifdef _LINUX_
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>  
#include <sys/stat.h>
#include <sys/resource.h>
#include <unordered_map>
#else
#include <io.h>
#include <direct.h>
#endif

using namespace std;

FILE* open_file(const char* file_name, const char* mode);


const int MAX_ID = 2000000000;
const int HASH_SCALE = 2;

struct Edge {
	int a;
	int b;
};

class Utility {

	//int* vertex_map_;
	unordered_map<int, int> vertex_map_;
	string dir_;

	static bool edge_compare(const Edge& e1, const Edge& e2);
	bool is_merge_finished(Edge* es, int size);

	void save_tmp_edges(Edge* edges, int size, int tmpFile);
	void merge(int size);

	int get_min_edge(Edge* es, int size);
	int get_vertex_id(int u, int& num);

public:
	Utility();
	~Utility();
	// create binary graph files from .txt file
	void format_graph(string dest, string src);
	void format_graph(string dest, string src, int percent, bool isSampleEdge);
};


#endif
