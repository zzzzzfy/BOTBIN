#ifndef _GRAPH_H_
#define _GRAPH_H_
//Thanks to the GS-index author for some code

#include "Utility.h"
#include "HashSet.h"
//#include "HashMap.h"
#include "BitMap.h"
#include <queue>
#include <vector>
#include <stack>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>

class Graph {


	int tot_core = 0;
	double bucket_width_ = 0.01;
	vector<vector<int>> BOTBIN_cnt;
	int len_eps_;
	vector<float> eps_idx;

	int iter_update_core = 0;
	int* sort_buf_;
	double pf = 0.001;
	bool dice_sim = false;
	int* pa_;
	int* cluster_belong;
	int* rank_;
	int* non_core_vis;
	float q_eps_;
	int q_miu_;
	string dir_; 
	int* deg_;
	long long* tot_deg_;
	int* edges_; 
	unsigned int* reverse_; 
	static unsigned int* offset_; 
	int* common_nbr_;
	int n_;	
	unsigned int m_; 
	int max_degree_;
	static vector<int> rnk;
	bool* tiny_deg_bm;
	int* BOTBIN_core_tag; // for query
	vector<vector<int> > sketch;

	BitMap* dynamic_bm_;

	vector<set<pair<int, int>, greater<pair<int, int>>>> BOTBIN_bucket_;
	vector<set<pair<float, int>, greater<pair<float, int>>>> BOTBIN_nbr_mp_;



	struct ed {
		ed() {}
		ed(int _cn, int _kth) {
			edcn = _cn;
			edkth = _kth;
		}
		int edcn, edkth;
	};
	long long qk_vis = 0;
	unordered_map<int, ed>* BOTBIN_edge_; // store edge

	static int core_comp_deg_;
	static unsigned int nbr_comp_offset_;

	vector<pair<int, int>> result_non_core_;
	vector<vector<pair<float, int>>> simialrity_mp_;

	BitMap* core_bm_;

	int* BOTBIN_core_bm_;
	//HashSet** cn_hs_;
	BitMap* cn_bm_;

	bool BOTBIN_check_nbr_mp_miu_eps(int u, int miu, float eps);
	void BOTBIN_quick_query(const float& eps, int& miu, string& result_file, bool non_core, const int& time_stamp);
	void BOTBIN_query(const float& eps, int& miu, string& result_file, const int& time_stamp, int non_core);
	void calc_k(double rho);
	void link_reverse_edges();

	bool static cmp_rnk(const int& u, const int& v);
	void sketch_common_nbr(const int& u, const int& v);
	bool sketch_search(int &u, int &pos);

	void sketch_insert(vector<int>& s, const int& w);
	void BOTBIN_update_nbr_order(int u, int root, float new_ss, float old__s);
	int find_eps_idx(float& sim);
	void BOTBIN_update_core_order(int u, float new_ss, float old_ss);
	bool BOTBIN_add_sub(int u, int v);
	bool BOTBIN_add_sub_naive(int u, int v);
	void BOTBIN_add_sub_old(int u, int v);

	void BOTBIN_update_all_core_order_del(int u, float sim, bool opt);
	void BOTBIN_update_all_core_order_add(int u, float sim, bool opt);

	bool BOTBIN_del_sub(int u, int v);
	bool BOTBIN_del_sub_naive(int u, int v);
	void BOTBIN_add(int u, int v);
	void BOTBIN_add_naive(int u, int v);

	void BOTBIN_del_naive(int a, int b);
	void BOTBIN_del(int a, int b);

public:
	int k_ = 6000;
	double rho_ = 0;
	Graph(string dir, double rho);
	~Graph();
	void sample_edges();
	void init();
	void BOTBIN_construct();
	void read_BOTBIN(int delta);
	void BOTBIN_query_old(const float& eps, int& miu, string& result_file, const int& time_stamp);

	void BOTBIN_load();

	bool is_BOTBIN_exist();

	void load();

	void print_nbr_mp(int u, int k);

	void print_sketch(int u);
	unsigned int binary_search_edge(const int u, const int v);

	void dynamic_BOTBIN(bool is_query, string result_file, vector<int>& eps_list, vector<int>& mu_list, 
		bool isbuild, int delta, bool is_naive);

	void dfs_clr(int x, int* group);
	void cluster_ari(string st, vector<int>& eps_list, vector<int>& mu_list);
};

#endif