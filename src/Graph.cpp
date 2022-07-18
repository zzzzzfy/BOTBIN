#include "Graph.h"
#include<iostream>
#include<cstdio>
#include <cstring>
#include <string.h>
#include <vector>
#include <numeric>
#include <random>
#include <chrono>
#include <ratio>
#include "assert.h"
#include<algorithm>
#include <fstream>
using namespace std::chrono;

const int INF = 1000000000;
int Graph::core_comp_deg_ = 0;
unsigned Graph::nbr_comp_offset_ = 0;
unsigned int* Graph::offset_ = NULL;
vector<int> Graph::rnk = vector<int>();

Graph::Graph(string dir, double rho ) {

	edges_ = NULL;   
	reverse_ = NULL;
	offset_ = NULL;
	common_nbr_ = NULL;
	cn_bm_ = NULL;
	sort_buf_ = NULL;
	dir_ = dir;
	n_ = 0;
	m_ = 0;

	// idx0_sd_ = NULL;
	// idx0_ed_ = NULL;
	pa_ = NULL;
	rank_ = NULL;
	rho_ = rho;
	calc_k(rho);
	if (dice_sim)
		cerr << "Dice sim:";
	cout << dir << " " << k_<<" "<<rho << "\n";


}


Graph::~Graph() {

	if (edges_ != NULL) delete[] edges_;
	if (reverse_ != NULL) delete[] reverse_;
	if (offset_ != NULL) {
		delete[] offset_;
		offset_ = NULL;
	}
}

void Graph::load() {
	// initialize number of vertices and edges
	printf("Begin to load graph\n");
	FILE* f_info = open_file((dir_ + "graph.info").c_str(), "rb");
	fread(&n_, sizeof(int), 1, f_info);
	fread(&m_, sizeof(unsigned int), 1, f_info);
	fread(&max_degree_, sizeof(int), 1, f_info);
	m_ *= 2;
	fclose(f_info);
	cout << "N:" << n_ << " M:" << m_ << "\n";
	edges_ = new int[m_];

	deg_ = new int[m_];

	offset_ = new unsigned int[n_ + 1];
	reverse_ = new unsigned int[m_];
	common_nbr_ = new int[m_];

	// initialize edges array
	FILE* f_dat = open_file((dir_ + "graph.dat").c_str(), "rb");
	fread(edges_, sizeof(int), m_, f_dat);
	fclose(f_dat);
	for (int i = 0; i < m_; i++) {
		if (edges_[i] > n_)
			cerr << i << "\n";
	}
	// initialize vertex degree and neighbors
	FILE* f_idx = open_file((dir_ + "graph.idx").c_str(), "rb");
	long pos = 0;
	unsigned int cur = 0;
	int degree;
	long long tot_deg = 0;
	int max_degree = 0;
	for (int i = 0; i < n_; ++i) {
		offset_[i] = cur;
		fread(&pos, sizeof(long), 1, f_idx);
		fread(&degree, sizeof(int), 1, f_idx);

		deg_[i] = degree;
		cur += degree;
		max_degree = max(max_degree, degree);
		tot_deg += (long long)(degree) * degree;
	}

	cout << "average_degree_through_edge:" << " " << double(tot_deg)/m_ <<" max_degree:"<<max_degree<< "\n";
	offset_[n_] = m_;
	fclose(f_idx);
	link_reverse_edges();
	printf("load graph finished\n");
}

void Graph::init() {

#if defined(_CN_HASH_SET_) || defined(_CN_SORT_MERGE_)
#ifdef _CN_HASH_SET_
	// _CN_HASH_SET_
	build_hash_set();
#endif
#else
	// _CN_BITMAP_
	cn_bm_ = new BitMap(n_);
#endif
}

unsigned int Graph::binary_search_edge(const int u, const int v) {
	unsigned int left = offset_[u];
	unsigned int right = offset_[u + 1] - 1;
	unsigned int mid;
	while (left <= right) {
		mid = left + (right - left) / 2;
		if (edges_[mid] < v) left = mid + 1;
		else if (edges_[mid] > v) right = mid - 1;
		else return mid;
	}

#ifdef _DEBUG_
	printf("search edge [%d, %d] failed\n", u, v);
#endif

	return mid;
}

void Graph::link_reverse_edges() {

	for (int v = 0; v < n_; ++v) {
		for (unsigned int i = offset_[v]; i < offset_[v + 1]; ++i) {
			int u = edges_[i];
			if (offset_[v + 1] - offset_[v] < offset_[u + 1] - offset_[u]) continue;

			unsigned int reverse_edge_id = binary_search_edge(u, v);

			reverse_[i] = reverse_edge_id;
			reverse_[reverse_edge_id] = i;
		}
	}

}



bool Graph::is_BOTBIN_exist() {
	return access((dir_ + to_string(k_) + "_rnk").c_str(), 0) == 0;
}


bool Graph::cmp_rnk(const int& u, const int& v) {
	return rnk[u] < rnk[v];
}

void Graph::sketch_common_nbr(const int& u, const int& v) {
	vector<int> uv;
	set_intersection(sketch[u].begin(), sketch[u].end(),
		sketch[v].begin(), sketch[v].end(), back_inserter(uv));
}

bool Graph::sketch_search(int& u, int& pos) {
	int l = 0, r = sketch[u].size() - 1;
	while (l <= r) {
		int mid = (l + r) >> 1;
		if (sketch[u][mid] == pos) return true;
		if (sketch[u][mid] < pos) {
			l = mid + 1;
		}
		else r = mid - 1;
	}
	return false;
}

int find_pos_deg(long long* a, int len, long long value) {
	int l = 0, r = len - 1;
	int res = len - 1;
	while (l <= r) {
		int mid = (l + r) >> 1;
		if (a[mid] >= value) {
			res = mid;
			r = mid - 1;
		}
		else l = mid + 1;
	}
	return res;
}

int tot_500, tot_1000, tot_2000, tot_4000, tot_5k, tot_1w5;

void Graph::sample_edges() {
	load();
	mt19937 rng(random_device{}());
	int* edges_list = new int[DYNAMIC_EDGE_NUM * 2];
	int* id = new int[m_];

	for (int u = 0; u < n_;u++) {
		for (unsigned int i = offset_[u]; i < offset_[u + 1]; ++i) {
			id[i] = u;
		}
	}
	
	for (int i = 0; i < DYNAMIC_EDGE_NUM; i++) {

		unsigned int m_id = rng() % m_;
		while (id[m_id] == -1) m_id = rng() % m_;
		int u = id[m_id];
		int v = edges_[m_id];
		edges_list[2 * i] = u;
		edges_list[2 * i + 1] = v;


		id[m_id] = -1; 
		id[reverse_[m_id]] = -1;
		int u_deg = offset_[u + 1] - offset_[u];
		int v_deg = offset_[v + 1] - offset_[v];
		if (u_deg >= 500 || v_deg >= 500) tot_500++;
		if (u_deg >= 1000 || v_deg >= 1000) tot_1000++;
		if (u_deg >= 2000 || v_deg >= 2000) tot_2000++;
		if (u_deg >= 4000 || v_deg >= 4000) tot_4000++;
		if (u_deg >= 5000 || v_deg >= 5000) tot_5k++;
		if (u_deg >= 15000 || v_deg >= 15000) tot_1w5++;
	}

	cout << "tot" << tot_500 << " " << double(tot_500) / DYNAMIC_EDGE_NUM << "\n";
	cout << "tot" << tot_1000 << " " << double(tot_1000) / DYNAMIC_EDGE_NUM << "\n";
	cout << "tot" << tot_2000 << " " << double(tot_2000) / DYNAMIC_EDGE_NUM << "\n";
	cout << "tot" << tot_4000 << " " << double(tot_4000) / DYNAMIC_EDGE_NUM << "\n";
	cout << "tot" << tot_5k << " " << double(tot_5k) / DYNAMIC_EDGE_NUM << "\n";
	cout << "tot" << tot_1w5 << " " << double(tot_1w5) / DYNAMIC_EDGE_NUM << "\n";

	FILE* f_edges = open_file((dir_ + "dynamic_edges").c_str(), "wb");
	fwrite(&DYNAMIC_EDGE_NUM, sizeof(int), 1, f_edges);
	fwrite(edges_list, sizeof(int), DYNAMIC_EDGE_NUM * 2, f_edges);
	fclose(f_edges);
}

int intersect(vector<int>& a, vector<int>& b, vector<int>& c) {
	vector<int> result;
	vector<int> ab;
	set_intersection(a.begin(), a.end(),
		b.begin(), b.end(),
		std::back_inserter(ab));
	set_intersection(ab.begin(), ab.end(),
		c.begin(), c.end(),
		std::back_inserter(result));
	return result.size();
}

int insert_find_kthelement_ab(const vector<int>& a, const vector<int>& b, const int& k, const int& old_k) {
	if (k > a.size() + b.size())
		return INF;
	int cnt = 0, res = 0;
	auto cur1 = a.begin();
	auto cur2 = b.begin();
	while (cnt != k) {
		if (*cur1 < *cur2) {
			cnt++;
			res = *cur1;
			cur1++;

		}
		else if (*cur1 == *cur2) {
			cnt++;
			res = *cur1;
			cur1++;
			cur2++;
		}
		else {
			cnt++;
			res = *cur2;
			cur2++;
		}
		if (cur1 == a.end() || cur2 == b.end())
			break;
	}
	//cout << cnt << "\n";

	if (cnt < k) {
		if (cur1 == a.end()) {
			for (; cnt != k && cur2 != b.end(); cnt++) {
				res = *cur2;
				cur2++;
			}
		}
		else {
			for (; cnt != k && cur1 != a.end(); cnt++) {
				res = *cur1;
				cur1++;
			}
		}
	}
	if (cnt < k)
		return INF;
	return res;
}

void find_kthelementcn_ab(const vector<int>& a, const vector<int>& b, const int& k, int &tot, int& res) {
	//if (k > a.size() + b.size())
	//	return INF;
	res = 0;
	tot = 0;
	int cnt = 0;
	auto cur1 = a.begin();
	auto cur2 = b.begin();
	while (cnt != k) {
		if (*cur1 < *cur2) {
			cnt++;
			res = *cur1;
			cur1++;

		}
		else if (*cur1 == *cur2) {
			cnt++;
			res = *cur1;
			tot++;
			cur1++;
			cur2++;
		}
		else {
			cnt++;
			res = *cur2;
			cur2++;
		}
		if (cur1 == a.end() || cur2 == b.end())
			break;
	}
	//cout << cnt << "\n";

	if (cnt < k) {
		if (cur1 == a.end()) {
			for (; cnt != k && cur2 != b.end(); cnt++) {
				res = *cur2;
				cur2++;
			}
		}
		else {
			for (; cnt != k && cur1 != a.end(); cnt++) {
				res = *cur1;
				cur1++;
			}
		}
	}
	if (cnt < k)
		res = INF;
	//return make_pair(cnt, res);
}

int find_kthelement_ab(const vector<int>& a, const vector<int>& b, const int& k) {
	if (k > a.size() + b.size()) 
		return INF;
	int cnt = 0, res = 0;
	auto cur1 = a.begin();
	auto cur2 = b.begin();
	while (cnt != k) {
		if (*cur1 < *cur2) {
			cnt++;
			res = *cur1;
			cur1++;

		}
		else if (*cur1 == *cur2) {
			cnt++;
			res = *cur1;
			cur1++;
			cur2++;
		}
		else {
			cnt++;
			res = *cur2;
			cur2++;
		}
		if (cur1 == a.end() || cur2 == b.end())
			break;
	}
	//cout << cnt << "\n";

	if (cnt < k) {
		if (cur1 == a.end()) {
			for (; cnt != k && cur2 != b.end(); cnt++) {
				res = *cur2;
				cur2++;
			}
		}
		else {
			for (; cnt != k && cur1 != a.end(); cnt++) {
				res = *cur1;
				cur1++;
			}
		}
	}
	if (cnt < k)
		return INF;
	return res;
}

void Graph::BOTBIN_construct() {
	// rewrite sort by vector

	printf("begin to construct BOTBIN\n");
	FILE* f_info = open_file((dir_ + "graph.info").c_str(), "rb");
	fread(&n_, sizeof(int), 1, f_info);
	fread(&m_, sizeof(unsigned int), 1, f_info);
	fread(&max_degree_, sizeof(int), 1, f_info);
	m_ *= 2;
	fclose(f_info);
	cout << n_ << " " << m_ << "\n";

	int* tmp_edges_ = new int[m_];

	tiny_deg_bm = new bool[n_];
	memset(tiny_deg_bm, 0, sizeof(bool) * n_);

	//nbr_order_.resize(n_);
	BOTBIN_edge_ = new unordered_map<int, ed>[n_];
	//BOTBIN_nbr_v_ = new vector<pair<int, ed>>[n_];
	simialrity_mp_.resize(n_);
	//BOTBIN_nbr_mp_.resize(n_, set<pair<float, int>, greater<pair<float, int>>>());


	FILE* f_dat = open_file((dir_ + "graph.dat").c_str(), "rb");
	fread(tmp_edges_, sizeof(int), m_, f_dat);
	fclose(f_dat);

	FILE* f_idx = open_file((dir_ + "graph.idx").c_str(), "rb");
	long pos = 0;
	unsigned int cur = 0, pre_cur = 0;
	int degree;
	for (int i = 0; i < n_; ++i) {
		fread(&pos, sizeof(long), 1, f_idx);
		fread(&degree, sizeof(int), 1, f_idx);
		pre_cur = cur;
		cur += degree;
		for (unsigned int j = pre_cur; j < cur; j++) {
			BOTBIN_edge_[i].insert(make_pair(tmp_edges_[j], ed(-1, -1)));
		}

	}
	fclose(f_idx);
	//printf("%d\n", cur);
	cout << cur << "\n";

	delete[] tmp_edges_;

	mt19937 rng(random_device{}());
	rnk.resize(n_);
	iota(rnk.begin(), rnk.end(), 0);
	shuffle(rnk.begin(), rnk.end(), rng);

	int* rnk_id = new int[n_];
	for (int i = 0; i < n_; i++) {
		rnk_id[i] = i;
	}
	sort(rnk_id, rnk_id + n_, cmp_rnk);
	sketch.resize(n_);

	for (int i = 0; i < n_; i++) {

		int u = rnk_id[i];
		if (sketch[u].size() < k_)
			sketch[u].push_back(i);

		for (auto j = BOTBIN_edge_[u].begin(); j != BOTBIN_edge_[u].end(); j++) {
			int v = j->first;
			if (sketch[v].size() >= k_) {
				continue;
			}
			sketch[v].push_back(i);
		}
	}
	int tot_deg = 0, deg_50 = 0;
	for (int u = 0; u < n_; u++) {
		for (int i = 0; i < sketch[u].size(); i++) {
			tiny_deg_bm[sketch[u][i]] = true;
			//tiny_deg_bm[sketch[u][i]] = u;
		}

		int u_deg = sketch[u].size();
		if (u % 1000000 == 0) {
			printf("Process Vertex id:%d deg:%d\n", u, BOTBIN_edge_[u].size());

			struct rusage rUsage;
			getrusage(RUSAGE_SELF, &rUsage);
			long ms = rUsage.ru_maxrss;
			float gms = ms / 1024 / 1024;
			printf("Memory usage = %ldKB, %ldMB, %fGB\n", ms, ms / 1024, gms);
		}
		//assert(u_deg>sketch[u].size());
		float j_sim = 0;

		for (auto j = BOTBIN_edge_[u].begin(); j != BOTBIN_edge_[u].end(); j++) {
			int v = j->first, cn = 0;
			if (v < u) {
				continue;
			}


			int cur_k = find_kthelement_ab(sketch[u], sketch[v], k_);
			int v_deg = sketch[v].size();

			if (cur_k == INF) {
				for (int i = 0; i < v_deg; i++) {
					cn += tiny_deg_bm[sketch[v][i]];
				}
				int uv_deg = u_deg + v_deg - cn;
				j->second.edcn = cn;
				j->second.edkth = INF;
				j_sim = float(cn) / uv_deg;
			}
			else {

				for (int i = 0; i < v_deg; i++) {
					int val = sketch[v][i];
					//cout << i <<" "<<sketch[v][i] << "\n";
					if (val <= cur_k)
						cn += tiny_deg_bm[sketch[v][i]];
					else break;
				}
				j->second.edcn = cn;
				j->second.edkth = cur_k;
				j_sim = float(cn) / k_;
			}
			auto it_v = BOTBIN_edge_[v].find(u);
			it_v->second.edcn = j->second.edcn;
			it_v->second.edkth = j->second.edkth;
			if (dice_sim)
				j_sim = 2 - (2 / (1 + j_sim));
			simialrity_mp_[u].push_back(make_pair(j_sim, v));
			simialrity_mp_[v].push_back(make_pair(j_sim, u));
			//BOTBIN_nbr_mp_[u].insert(make_pair(j_sim, v));
		}
		for (int i = 0; i < sketch[u].size(); i++) {
			tiny_deg_bm[sketch[u][i]] = false;
		}
	}


	ofstream tmpf;
	cout << "store_rnk\n";
	freopen((dir_ + to_string(k_) + "_rnk").c_str(), "w", stdout);
	for (int i = 0; i < n_; i++)
		cout << rnk[i] << " ";
	fclose(stdout);

	FILE* f_BOTBIN = open_file((dir_ + to_string(k_) + "_nbr_mp").c_str(), "wb");
	//freopen((dir_ + "_nbr_mp").c_str(), "w", stdout);
	for (int i = 0; i < n_; i++) {
		int a = BOTBIN_edge_[i].size();
		fwrite(&a, sizeof(int), 1, f_BOTBIN);
		//cout << BOTBIN_nbr_[i].size() << "\n";
		float _s;
		int _v1, _v2, _v3;
		for (auto j = BOTBIN_edge_[i].begin(); j != BOTBIN_edge_[i].end(); j++) {
			_v1 = j->first; _v2 = j->second.edcn; _v3 = j->second.edkth;
			fwrite(&_v1, sizeof(int), 1, f_BOTBIN);
			fwrite(&_v2, sizeof(int), 1, f_BOTBIN);
			fwrite(&_v3, sizeof(int), 1, f_BOTBIN);
		}
		sort(simialrity_mp_[i].begin(), simialrity_mp_[i].end(), greater<pair<float, int>>());
		for (auto j = simialrity_mp_[i].begin(); j != simialrity_mp_[i].end(); j++) {
			_s = j->first; _v1 = j->second;
			fwrite(&_s, sizeof(float), 1, f_BOTBIN);
			fwrite(&_v1, sizeof(int), 1, f_BOTBIN);
		}
	}
	fclose(f_BOTBIN);


	// calculate construction time
	delete[] rnk_id;
}


void Graph::read_BOTBIN(int delta) {


	printf("begin to read BOTBIN\n");
	FILE* f_info = open_file((dir_ + "graph.info").c_str(), "rb");
	fread(&n_, sizeof(int), 1, f_info);
	fread(&m_, sizeof(unsigned int), 1, f_info);
	fread(&max_degree_, sizeof(int), 1, f_info);
	m_ *= 2;
	fclose(f_info);
	cout << n_ << " " << m_ << "\n";
	struct rusage rUsage;
	cout << "delta:" << delta << "\n";
	BOTBIN_edge_ = new unordered_map<int, ed>[n_];
	BOTBIN_core_bm_ = new int[n_];
	BOTBIN_core_tag = new int[n_];
	for (int i = 0; i < n_; i++) BOTBIN_core_bm_[i] = BOTBIN_core_tag[i] = 0;

	BOTBIN_nbr_mp_.resize(n_, set<pair<float, int>, greater<pair<float, int>>>());

	bucket_width_ = 1.0/delta;

	len_eps_ = delta;
	cout << "len_delta:" << len_eps_ << "\n";
	BOTBIN_bucket_.resize(len_eps_, set<pair<int, int>, greater<pair<int, int>>>());
	eps_idx.resize(len_eps_);

	double st = 1 - bucket_width_;
	for (int i = 0; i < len_eps_; i++, st -= bucket_width_){
		eps_idx[i] = 1 - ((double)(i+1))/(delta);
	}
	eps_idx[len_eps_ - 1] = 0;

	BOTBIN_cnt.resize(n_);

	ifstream rnkf;
	rnkf.open((dir_ + to_string(k_) + "_rnk").c_str(), ios::app);
	for (int i = 0; i < n_; i++) {
		int u;
		rnkf >> u;
		rnk.push_back(u);
	}
	rnkf.close();
	FILE* f_nbr = open_file((dir_ + to_string(k_)+"_nbr_mp").c_str(), "rb");
	int size_;
	float sim_;
	for (int i = 0; i < n_; i++) {
		//cin >> size_;
		fread(&size_, sizeof(int), 1, f_nbr);
		BOTBIN_cnt[i].resize(len_eps_);
		int u1, u2, u3;
		for (int j = 0; j < size_; j++) {
			fread(&u1, sizeof(int), 1, f_nbr);
			fread(&u2, sizeof(int), 1, f_nbr);
			fread(&u3, sizeof(int), 1, f_nbr);
			BOTBIN_edge_[i].insert(make_pair(u1, ed(u2, u3)));
		}
		int cur_eps = 0;
		for (int j = 0; j < size_; j++) {

			fread(&sim_, sizeof(float), 1, f_nbr);
			fread(&u1, sizeof(int), 1, f_nbr);
			//cin >> sim_ >> u1;
			//169222 8024

			BOTBIN_nbr_mp_[i].insert(make_pair(sim_, u1));
			if (j == 0) {
				while (sim_ < eps_idx[cur_eps]) cur_eps++;
			}
			else {
				if (sim_ < eps_idx[cur_eps]) {

					while (sim_ < eps_idx[cur_eps]) {
						BOTBIN_cnt[i][cur_eps] = j;
						BOTBIN_bucket_[cur_eps].insert(make_pair(j, i));
						cur_eps++;
					}
				}
			}
		}
		while (sim_ < eps_idx[cur_eps]) {
			cur_eps++;
		}
		while (cur_eps != len_eps_) {
			BOTBIN_cnt[i][cur_eps] = size_;
			BOTBIN_bucket_[cur_eps].insert(make_pair(size_, i));
			cur_eps++;
		}
		if (i % 1000000 == 0) {
			getrusage(RUSAGE_SELF, &rUsage);
			long ms = rUsage.ru_maxrss;
			float gms = ms / 1024 / 1024;
			cout << i << ":";
			printf("Memory usage = %ldKB, %ldMB, %fGB\n", ms, ms / 1024, gms);
		}
	}


	cout << "begin to set sketch\n";
	int* rnk_id = new int[n_];
	for (int i = 0; i < n_; i++) {
		rnk_id[i] = i;
	}
	sort(rnk_id, rnk_id + n_, cmp_rnk);
	sketch.resize(n_);

	for (int i = 0; i < n_; i++) {

		int u = rnk_id[i];
		if (sketch[u].size() < k_)
			sketch[u].push_back(i);

		for (auto j = BOTBIN_edge_[u].begin(); j != BOTBIN_edge_[u].end(); j++) {
			int v = j->first;
			if (sketch[v].size() >= k_) {
				continue;
			}
			sketch[v].push_back(i);
		}
	}
	delete[] rnk_id;

	getrusage(RUSAGE_SELF, &rUsage);

	long ms = rUsage.ru_maxrss;
	float gms = ms / 1024 / 1024;
	printf("Final Memory usage = %ldKB, %ldMB, %fGB\n", ms, ms / 1024, gms);


	ofstream tmpf;

	tmpf.open("BOTBIN_build.txt", ios::app);
	tmpf << "Memory usage:" << ms/1024 << " " << gms << "\n";
	
	tmpf.close();
}


bool Graph::BOTBIN_check_nbr_mp_miu_eps(int u, int miu, float eps) {
	int deg_cur = 1;
	for (auto j = BOTBIN_nbr_mp_[u].begin(); j != BOTBIN_nbr_mp_[u].end(); j++, deg_cur++) {
		if (deg_cur == miu) {
			return j->first >= eps;
		}
	}
	return false;
}

void Graph::BOTBIN_quick_query(const float& eps, int& miu, string& result_file, bool non_core, const int& time_stamp) {
	freopen(result_file.c_str(), "w", stdout);
	queue<int> q;
	if (non_core) {
		cerr << "BOTBIN_query: output non core to result\n";
	}
	if (cluster_belong == NULL)
		cluster_belong = new int[n_];

	if (non_core_vis == NULL){
		non_core_vis = new int[n_];
		for (int i = 0; i < n_; i++) non_core_vis[i] = 0;
	}
	cerr <<result_file<<" " << "quickquery\n";
	vector<vector<int> > core_vec;
	int tot_cluster = 0;
	result_non_core_.clear(); // non-core
	int q_eps_ = 0;
	for (int i = 0; i < len_eps_; i++) {
		//cerr << eps_idx[i] << " ";
		if (eps < eps_idx[i]) continue;
		else{
		// find the first bucket contain eps_
			q_eps_ = i;
			break;
		}
	}
	cerr << eps << " " << miu << '\n';
	cerr << eps_idx[q_eps_] << "\n";
	vector<int> tmp_core;
	auto core_ed = BOTBIN_bucket_[q_eps_].end();
	for (auto i = BOTBIN_bucket_[q_eps_].begin(); i != core_ed; i++) {
		//tot_core++;
		if (miu > i->first)
			break;
		BOTBIN_core_tag[i->second] = time_stamp;
		tmp_core.push_back(i->second);
	}


	for (int i = 0; i < tmp_core.size(); i++) {
		int u = tmp_core[i];
		if (BOTBIN_core_bm_[u] == time_stamp) continue;

		q.push(u);
		core_vec.push_back(vector<int>());
		core_vec[tot_cluster].push_back(u);

		BOTBIN_core_bm_[u] = time_stamp;
		while (!q.empty()) {
			int q_id = q.front();
			cluster_belong[q_id] = tot_cluster;
			q.pop();
			auto tmp_ed = BOTBIN_nbr_mp_[q_id].end();
			for (auto j = BOTBIN_nbr_mp_[q_id].begin(); j != tmp_ed; j++) {
				int v = j->second;
				if (j->first < eps) break;
				if (BOTBIN_core_bm_[v] == time_stamp) continue;
				//nbr_mp_[v][miu - 1] >= eps
				if (BOTBIN_nbr_mp_[v].size() >= miu && BOTBIN_core_tag[v]==time_stamp) {
					BOTBIN_core_bm_[v] = time_stamp;
					core_vec[tot_cluster].push_back(v);
					q.push(v);
				}
				else {
					result_non_core_.push_back(make_pair(u, v));
				}
			}
		}
		tot_cluster++;
	}
	if (non_core) {
		for (auto i : result_non_core_) {
			int u = i.second;
			if (non_core_vis[u]) continue;
			int min_id = 1e9;
			for (auto j : BOTBIN_nbr_mp_[u]) {
				int v = j.second;
				if (j.first < eps) break;
				if (BOTBIN_core_bm_[v] == time_stamp)
					min_id = min(min_id, v);
			}
			core_vec[cluster_belong[min_id]].push_back(u);
			non_core_vis[u] = true;
		}
	}
	cout << tot_cluster << "\n";
	for (int i = 0; i < tot_cluster; i++) {
		cout << core_vec[i].size();
		for (auto j = core_vec[i].begin(); j != core_vec[i].end(); j++) {
			tot_core++;
			non_core_vis[*j] = false;
			cout << " " << *j;
		}
		cout << "\n";
	}
	cout.flush();
	fclose(stdout);
	return;
}

void Graph::BOTBIN_query(const float &eps, int &miu, string &result_file, const int &time_stamp, int non_core) {
	freopen(result_file.c_str(), "w", stdout);
	queue<int> q;
	if (non_core) {
		cerr << "BOTBIN_query: output non core to result\n";
	}
	if (cluster_belong == NULL)
		cluster_belong = new int[n_];

	if (non_core_vis == NULL)
		non_core_vis = new int[n_];

	if (non_core) {
		for (int i = 0; i < n_; i++) non_core_vis[i] = 0;
	}
	result_non_core_.clear();

	vector<vector<int> > core_vec;
	int tot_cluster = 0;

	int q_eps_ = 0;
	for (int i = 0; i < len_eps_; i++) {
		if (abs(eps - eps_idx[i]) < 1e-4)
			q_eps_ = i;
	}
	cerr << eps << " " << miu << '\n';
	//while (eps < eps_idx[q_eps_]) q_eps_++;
	//cerr << eps_idx[q_eps_] << "\n";
	//q_eps_++;
	cerr << eps_idx[q_eps_] << "\n";

	//q_eps_++;
	cerr << q_eps_ << " " << eps_idx[q_eps_]<<"\n";
	for (auto i = BOTBIN_bucket_[q_eps_].begin(); i != BOTBIN_bucket_[q_eps_].end(); i++) {
		//tot_core++;
		if (BOTBIN_core_bm_[i->second] == time_stamp) continue;
		if (miu > i->first)
			break;
		int u = i->second;

		q.push(u);
		core_vec.push_back(vector<int>());
		core_vec[tot_cluster].push_back(u);

		BOTBIN_core_bm_[u] = time_stamp;
		while (!q.empty()) {
			int q_id = q.front();
			cluster_belong[q_id] = tot_cluster;
			q.pop();
			for (auto j = BOTBIN_nbr_mp_[q_id].begin(); j != BOTBIN_nbr_mp_[q_id].end(); j++) {
				int v = j->second;
				if (j->first < eps) break;
				if (BOTBIN_core_bm_[v] == time_stamp) continue;
				//nbr_mp_[v][miu - 1] >= eps
				if (BOTBIN_nbr_mp_[v].size() >= miu && BOTBIN_check_nbr_mp_miu_eps(v, miu, eps)) {
					BOTBIN_core_bm_[v] = time_stamp;
					core_vec[tot_cluster].push_back(v);
					q.push(v);
				}
				else {
					result_non_core_.push_back(make_pair(u, v));
				}
			}
		}
		tot_cluster++;
	}

	if (non_core) {
		for (auto i : result_non_core_) {
			int u = i.second;
			if (non_core_vis[u]) continue;
			int min_id = 1e9;
			for (auto j : BOTBIN_nbr_mp_[u]) {
				int v = j.second;
				if (j.first < eps) break;
				if(BOTBIN_core_bm_[v] == time_stamp)
					min_id = min(min_id, v);
			}
			core_vec[cluster_belong[min_id]].push_back(u);
			non_core_vis[u] = true;
		}
	}


	cout << tot_cluster << "\n";
	for (int i = 0; i < tot_cluster; i++) {
		cout << core_vec[i].size();
		for (auto j = core_vec[i].begin(); j != core_vec[i].end(); j++) {
			tot_core++;
			cout << " " << *j;
		}
		cout << "\n";
	}
	cout.flush();
	fclose(stdout);

}

void Graph::calc_k(double rho) {
	long long M = 2405026092;
	long long avg_deg = 1000; 
	// M and avg_deg should be set in different graph
	k_ = 0;

	long long fz = (long long)M * avg_deg;

	int tmp_k = (1 / (2 * rho * rho)) * log(2  / (pf/M/avg_deg));
	k_ = tmp_k;

}



void Graph::BOTBIN_load() {
	if (is_BOTBIN_exist()) {
		return;
	}
	else {
		//load();
		BOTBIN_construct();
	}
}



int generate_vertex_deg_dist(uniform_int_distribution <long long>& dis, mt19937_64& eng, long long* tot_deg, int n_) {
	long long ran = dis(eng);
	//cout << ran << "\n";
	int lo = 0, hi = n_ - 1;
	while (lo < hi) {
		int mid = (lo + hi) >> 1;
		if (tot_deg[mid] < ran) {
			lo = mid + 1;
		}
		else {
			hi = mid;
		}
	}
	return lo;
}



void Graph::dfs_clr(int x, int* group) {
	queue<int> q;
	q.push(x);
	while (!q.empty()) {
		x = q.front();
		q.pop();
		for (unsigned int i = offset_[x]; i < offset_[x + 1]; ++i) {
			int v = edges_[i];
			//cerr << x << " " << v << "\n";
			if (group[v] == 0) {
				group[v] = group[x];
				q.push(v);
			}
		}
	}
}

void Graph::cluster_ari(string st, vector<int>& eps_list, vector<int>& mu_list) {
	// calculate the ari between two results file _1_xx _2_xx
	load();
	int* group = new int[n_];
	int componet_num = 0;
	for (int i = 0; i < n_; i++) group[i] = 0;
	for (int i = 0; i < n_; i++) {
		if (group[i] == 0) {
			componet_num++;
			group[i] = componet_num;
			dfs_clr(i, group);
		}
	}
	cerr <<"component:" << componet_num << "\n";
	int* mx_cls_sz = new int [componet_num + 5];
	for (int T = 0; T < eps_list.size(); T++) {
		for (int i = 1; i <= componet_num; i++) mx_cls_sz[i] = 0;
		int eps = eps_list[T];
		int mu = mu_list[T];
		cerr << eps << " " << mu << "\n";
		string file1 = st + "_result/_1_" + to_string(eps) + "_" + to_string(mu);
		ifstream tmpf;
		int tot_cluster_1, cluster_size, id, tot_node_1 = 0;
		tmpf.open(file1.c_str(), ios::app);
		tmpf >> tot_cluster_1;
		int mx_cluster = 0;
		map<int, int> cluster_label1;
		map<int, int> cluster_label2;
		//cerr << "cls2:" << mx_cls_sz[2] << "\n";
		for (int i = 0; i < tot_cluster_1; i++) {
			tmpf >> cluster_size;
			tot_node_1 += cluster_size;
			mx_cluster = max(mx_cluster, cluster_size);
			for (int j = 0; j < cluster_size; j++) {
				tmpf >> id;
			}
			if (cluster_size < 0) {
				cerr << i <<" size"<<cluster_size << "~~~~~~\n";
			}
			//if()
			mx_cls_sz[group[id]] = max(mx_cls_sz[group[id]], cluster_size);
		}
		int tot_v = 0;
		for (int i = 1; i <= componet_num; i++) {
			tot_v += mx_cls_sz[i];
			if (mx_cls_sz[i] < 0) {
				cerr << i << " " << mx_cls_sz[i] << "~~~~~~\n";
			}
		}

		cerr << tot_node_1 <<" "<<tot_v<<" "<<n_ <<" "<<sizeof mx_cls_sz << "ratio:" << (double)tot_v / n_ << "\n";
		//string file2 = st + "_result/_2_" + to_string(eps) + "_" + to_string(mu);
	}



}


void Graph::print_nbr_mp(int u, int k) {
	cout << "nbr_mp\n";
	int cnt = 0;
	for (auto it = BOTBIN_nbr_mp_[u].begin(); it != BOTBIN_nbr_mp_[u].end() && cnt!=k; it++, k--) {
		cout <<cnt<<" "<<u<<" "<< it->first << " " << it->second << "\n";
	}
}

void Graph::print_sketch(int u) {
	cout << "sketch" << u << "\n";
	for (int i = 0; i < sketch[u].size(); i++)
		cout << sketch[u][i] << " ";
	cout << '\n';
}



void Graph::dynamic_BOTBIN(bool is_query, string result_file, vector<int>& eps_list, vector<int>& mu_list, 
	bool isbuild, int delta, bool is_naive) {
	if (isbuild || !is_BOTBIN_exist()) {
		BOTBIN_construct();
		return;
	}
	else {
		read_BOTBIN(delta);
	}
	int* edges_list = new int[DYNAMIC_EDGE_NUM * 3];

	if (is_query) {

		cerr << "BOTBIN_query" << dir_ << " " << delta << " " << rho_ << "\n";
		
		for (int T = 0; T < eps_list.size(); T++) {
			int i = eps_list[T], j = mu_list[T];
			string tmp = result_file + "_result/_2_" + to_string(i) + "_" + to_string(j);
			BOTBIN_quick_query(double(i) / 100, j, tmp, false, T+1);
		}

		ofstream tmpf;

		cerr << dir_ << " " << delta << " " << rho_ << " " << "BOTBIN_query:"  << "\n";
		cerr << "tot_core:"<<tot_core << "\n";
		return;
	}


	printf("load edge list from disk\n");
	FILE* f_edges = open_file((dir_ + "dynamic_edges").c_str(), "rb");
	int edge_num;
	fread(&edge_num, sizeof(int), 1, f_edges);
	if (edge_num != DYNAMIC_EDGE_NUM) {
		printf("edge numbers are not same: %d edges in file, DYNAMIC_EDGE_NUM = %d \n", edge_num, DYNAMIC_EDGE_NUM);
		exit(1);
	}
	fread(edges_list, sizeof(int), DYNAMIC_EDGE_NUM * 2, f_edges);
	fclose(f_edges); 
	high_resolution_clock::time_point start_cal = high_resolution_clock::now();

	long long tot_vis = 0;

	for (int i = 0; i < DYNAMIC_EDGE_NUM; ++i) {
		int a = edges_list[2 * i];
		int b = edges_list[2 * i + 1];

		tot_vis += BOTBIN_edge_[a].size() + BOTBIN_edge_[b].size();

		is_naive ? BOTBIN_del_naive(a, b) : BOTBIN_del(a, b);
	}

	high_resolution_clock::time_point end_cal = high_resolution_clock::now();
	duration<double> time_cal_single = duration_cast<duration<double>>(end_cal - start_cal);
	double total_cal = time_cal_single.count();
	cout << "BOTBIN_del:" << total_cal <<" "<<tot_vis<<" "<< qk_vis << "\n";
	ofstream tmpf{};

	tmpf.open("/home/user_name/update_time.txt", ios::app);
	if (is_naive) tmpf << "NAIVE_BOTBIN:";
	tmpf << dir_ << " " << rho_ << " " << delta << "\n";
	tmpf << "BOTBIN_del:" << total_cal << " " << tot_vis << " " << qk_vis << "\n";
	tmpf.close();

	tot_vis = 0;
	qk_vis = 0;
	start_cal = high_resolution_clock::now();
	for (int i = 0; i < DYNAMIC_EDGE_NUM; ++i) {
		int a = edges_list[2 * i];
		int b = edges_list[2 * i + 1];
		tot_vis += BOTBIN_edge_[a].size() + BOTBIN_edge_[b].size();
		//cout <<i<<":"<< a << " " << b << "\n";
		is_naive ? BOTBIN_add_naive(a, b) : BOTBIN_add(a, b);
	}
	end_cal = high_resolution_clock::now();
	time_cal_single = duration_cast<duration<double>>(end_cal - start_cal);
	total_cal = time_cal_single.count();
	cout << "BOTBIN_add:" << total_cal << "\n";
	cout.flush();

	tmpf.open("/home/user_name/update_time.txt", ios::app);
	if (is_naive) tmpf << "NAIVE_BOTBIN:";
	tmpf << dir_ << " " << rho_ << "\n";
	tmpf << "BOTBIN_add:" << total_cal << " " << tot_vis << " " << qk_vis << "\n";
	tmpf.close();

	fclose(stdout);
}


int cal_intersect_ab(vector<int>& a, vector<int>& b) {
	int tot = 0;
	auto i = a.begin(), j = b.begin();
	for (; i != a.end() && j != b.end();) {
		if (*i == *j) {
			tot++;
			i++;
			j++;
		}
		else if (*i > * j) {
			j++;
		}
		else {
			i++;
		}
		//if(*i>)
	}
	return tot;
}

int cal_intersect_ab_k(vector<int>& a, vector<int>& b, const int& k) {
	int tot = 0;
	auto i = a.begin(), j = b.begin();
	for (; i != a.end() && j != b.end();) {
		if (*i > k || *j > k) {
			break;
		}
		if (*i == *j) {
			tot++;
			i++;
			j++;
		}
		else if (*i > * j) {
			j++;
		}
		else {
			i++;
		}

		//if(*i>)
	}
	return tot;
}

void Graph::sketch_insert(vector<int>& s, const int& w) {
	int l = 0, r = s.size() - 1;
	if (s[r] < w) {
		if (r < k_ - 1)
			s.push_back(w);
		return;
	}
	int mid = 0;
	int pos = -1;
	while (l <= r) {
		mid = (l + r) >> 1;
		if (s[mid] < w) {
			pos = mid;
			l = mid + 1;
		}
		else r = mid - 1;
	}
	s.insert(s.begin() + pos + 1, w);
	if (s.size() > k_)
		s.pop_back();
	return;
}

void Graph::BOTBIN_update_nbr_order(int u, int root, float new_ss, float old_ss) {
	auto rit = BOTBIN_nbr_mp_[u].find(make_pair(old_ss, root));
	BOTBIN_nbr_mp_[u].erase(rit);
	BOTBIN_nbr_mp_[u].insert(make_pair(new_ss, root));
}

int Graph::find_eps_idx(float& sim) {
	if (sim == 1)
		return 0;
	else {
		int cur = len_eps_ - floor(sim / bucket_width_) - 1;
		if (sim >= eps_idx[cur]) return cur;
		else if (sim >= eps_idx[cur + 1]) return cur + 1;
		else return cur - 1;
	}
	//else
		//return len_eps_ - floor(sim / eps_) - 1;
}

void Graph::BOTBIN_update_core_order(int u, float new_ss, float old_ss) {

	int old_cur = find_eps_idx(old_ss);
	int new_cur = find_eps_idx(new_ss);
	//if(u== 230260)
	//cout << new_ss << " " << old_ss << " " << old_cur << " " << new_cur << "\n";
	if (old_cur == new_cur) return;
	if (new_ss > old_ss) {
		// new_cur < old_cur
		for (int i = new_cur; i < old_cur; i++) {
			if (BOTBIN_cnt[u][i] == 0) {
				BOTBIN_bucket_[i].insert(make_pair(1, u));
			}
			else {
				auto rit = BOTBIN_bucket_[i].find(make_pair(BOTBIN_cnt[u][i], u));
				BOTBIN_bucket_[i].erase(rit);
				BOTBIN_bucket_[i].insert(make_pair(BOTBIN_cnt[u][i] + 1, u));

			}
			BOTBIN_cnt[u][i]++;
		}
	}
	else {
		for (int i = old_cur; i < new_cur; i++) {
			auto rit = BOTBIN_bucket_[i].find(make_pair(BOTBIN_cnt[u][i], u));
			//cout << u << " " << i << " " << BOTBIN_cnt[u][i] << "\n";
			BOTBIN_bucket_[i].erase(rit);
			BOTBIN_bucket_[i].insert(make_pair(BOTBIN_cnt[u][i] - 1, u));

			BOTBIN_cnt[u][i]--;
		}
	}

}

int find_lower(vector<int>& s1, vector<int>& s2, const int& y) {
	int l = 0, r = s1.size() - 1, mid;
	int res = 0;
	while (l <= r) {
		mid = (l + r) >> 1;
		//cout << mid << " " << s1[mid] << "\n";
		if (s1[mid] < y) {
			res = s1[mid];
			l = mid + 1;
		}
		else r = mid - 1;
	}

	l = 0, r = s2.size() - 1;
	while (l <= r) {
		mid = (l + r) >> 1;
		if (s2[mid] < y) {
			res = max(s2[mid], res);
			l = mid + 1;
		}
		else r = mid - 1;
	}
	return res;
}

bool Graph::BOTBIN_add_sub(int u, int v) {
	if (sketch[u].size() >= k_ ) {
		// must add a element to sketch[u]
		//int pos = find_element_sketch(sketch[u], rnk[v]);
		//if (u == 8024)
		//	cout <<u<<" "<<pos<< "!!\n";
		if (sketch[u][k_-1]<rnk[v]) {
			qk_vis += BOTBIN_edge_[u].size();
			//cout << "quick_del" << " " << u << " " << v << "\n";
			return true;
		}
		else {
			int pre = sketch[u][k_ - 1], x = n_;
			sketch_insert(sketch[u], rnk[v]); // auto pop

			BOTBIN_nbr_mp_[u].clear();
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				int w = i->first;
				if (w == v) continue;
				int cn = i->second.edcn;

				int uw_deg = sketch[u].size() + sketch[w].size();
				float old_sim = float(cn) / (min(uw_deg - cn, k_));
				int y = i->second.edkth;
				int new_y = y;
				if (rnk[v] <= y) {
					if (sketch_search(w, rnk[v]))
						cn += 1;
					else {
						new_y = find_lower(sketch[u], sketch[w], y);
						if ((y == pre || sketch_search(u, y)) && sketch_search(w, y))
							cn -= 1;
						i->second.edkth = new_y;
					}
				}
				float new_sim = float(cn) / (min(uw_deg - cn, k_));
				BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
				auto itw = BOTBIN_edge_[w].find(u);

				i->second.edcn = cn;
				i->second.edkth = new_y;

				itw->second.edcn = cn;
				itw->second.edkth = new_y;
				if (old_sim == new_sim)
					continue;
				BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
				BOTBIN_update_core_order(w, new_sim, old_sim);
			}
		}
	}
	else {
		BOTBIN_nbr_mp_[u].clear();
		//print_sketch(u);
		//cout << rnk[v] << "\n";
		//if (u == 5386) cout << sketch[u].size() << "\n";
		sketch_insert(sketch[u], rnk[v]);
		//int pos = find_element_sketch(sketch[u], rnk[v]);
		//sketch[u].erase(sketch[u].begin() + pos);

		for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {

			int w = i->first;
			if (w == v) continue;


			int cn = i->second.edcn;
			int uw_deg = sketch[u].size() + sketch[w].size();
			float old_sim = float(cn) / (min(uw_deg - cn - 1, k_));
			//float new_sim = old_sim;
			//cout << w <<" "<<BOTBIN_nbr_mp_[w].size()<< " " << old_sim << " " << cn <<" "<<uw_deg<< "\n";
			int y = i->second.edkth;
			int new_y = y;

			if (rnk[v] <= y) {
				if (sketch_search(w, rnk[v]))
					cn += 1;
				else {
					if (y != INF) {
						new_y = find_lower(sketch[u], sketch[w], y);
						if (sketch_search(u, y) && sketch_search(w, y))
							cn -= 1;
					}
					if (y == INF) {
						int new_uw_deg = sketch[u].size() + sketch[w].size() - cn;
						if (new_uw_deg == k_) {
							new_y = min(sketch[u][sketch[u].size() - 1], sketch[w][sketch[w].size() - 1]);
						}
						// can quickly find
						//new_y = find_kthelement_ab(sketch[u], sketch[v], k_);
					}
				}
			}
				

			float new_sim = float(cn) / (min(uw_deg - cn, k_));
			auto itw = BOTBIN_edge_[w].find(u);

			i->second.edcn = cn;
			i->second.edkth = new_y;

			itw->second.edcn = cn;
			itw->second.edkth = new_y;
			BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			if (new_sim == old_sim) continue;
			
			//if()
			BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
			BOTBIN_update_core_order(w, new_sim, old_sim);
		}
	}
	return false;
}

bool Graph::BOTBIN_add_sub_naive(int u, int v) {
	if (sketch[u].size() >= k_) {
		if (sketch[u][k_ - 1] < rnk[v]) {
			qk_vis += BOTBIN_edge_[u].size();
			return true;
		}
		else {
			int pre = sketch[u][k_ - 1], x = n_;
			sketch_insert(sketch[u], rnk[v]); // auto pop

			BOTBIN_nbr_mp_[u].clear();
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				int w = i->first;
				if (w == v) continue;
				int cn = i->second.edcn;

				int uw_deg = sketch[u].size() + sketch[w].size();
				float old_sim = float(cn) / (min(uw_deg - cn, k_));
				int y = i->second.edkth;
				int new_y = y;

				int uw_kth = find_kthelement_ab(sketch[u], sketch[w], k_);
				cn = cal_intersect_ab_k(sketch[u], sketch[w], uw_kth);

				float new_sim = float(cn) / (min(uw_deg - cn, k_));
				BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
				auto itw = BOTBIN_edge_[w].find(u);

				i->second.edcn = cn;
				i->second.edkth = new_y;

				itw->second.edcn = cn;
				itw->second.edkth = new_y;
				if (old_sim == new_sim)
					continue;
				BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
				BOTBIN_update_core_order(w, new_sim, old_sim);
			}
		}
	}
	else {
		BOTBIN_nbr_mp_[u].clear();
		//print_sketch(u);
		//cout << rnk[v] << "\n";
		//if (u == 5386) cout << sketch[u].size() << "\n";
		sketch_insert(sketch[u], rnk[v]);
		//int pos = find_element_sketch(sketch[u], rnk[v]);
		//sketch[u].erase(sketch[u].begin() + pos);

		for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {

			int w = i->first;
			if (w == v) continue;


			int cn = i->second.edcn;
			int uw_deg = sketch[u].size() + sketch[w].size();
			float old_sim = float(cn) / (min(uw_deg - cn - 1, k_));
			//float new_sim = old_sim;
			//cout << w <<" "<<BOTBIN_nbr_mp_[w].size()<< " " << old_sim << " " << cn <<" "<<uw_deg<< "\n";
			int y = i->second.edkth;
			int new_y = y;


			int uw_kth = find_kthelement_ab(sketch[u], sketch[w], k_);
			cn = cal_intersect_ab_k(sketch[u], sketch[w], uw_kth);

			float new_sim = float(cn) / (min(uw_deg - cn, k_));
			auto itw = BOTBIN_edge_[w].find(u);

			i->second.edcn = cn;
			i->second.edkth = new_y;

			itw->second.edcn = cn;
			itw->second.edkth = new_y;
			//if (w == 207538) {
			//	cout << "add_sub:" << w << " " << u << " " << uw_deg << " " << cn
			//		<< " " << sketch[u].size() << " " << sketch[w].size()<<" "<<new_y << "\n";
			//}
			// rnk[v] > y only insert rnk[v] to sketch[u]
			//BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			if (new_sim == old_sim) continue;

			//if()
			BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
			BOTBIN_update_core_order(w, new_sim, old_sim);
			//if(w == 207538)
		}
	}
	return false;
}

void Graph::BOTBIN_add_sub_old(int u, int v) {
	if (sketch[u].size() >= k_) {
		int x = sketch[u][k_ - 1];
		if (x < rnk[v]) {
			return;
		}
		else {
			sketch_insert(sketch[u], rnk[v]);
			BOTBIN_nbr_mp_[u].clear();
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				int w = i->first;
				int cn = i->second.edcn;
				int uw_deg = sketch[u].size() + sketch[w].size();
				float old_sim = float(cn) / (min(uw_deg - cn, k_));
				int y = i->second.edkth;
				if (rnk[v] <= y) {
					if (sketch_search(w, rnk[v])) {
						cn += 1;
						if (x == y) {
							if (sketch_search(w, x)) cn -= 1;
							else {
								i->second.edkth = find_lower(sketch[u], sketch[w], y);
							}
						}
					}
					else {
						if (x == y) {
							if (sketch_search(w, x)) cn -= 1;
							else {
								i->second.edkth = find_lower(sketch[u], sketch[w], y);
							}
						}
						else {
							if (sketch_search(w, y) && sketch_search(u, y)) {
								cn -= 1;
							}
							i->second.edkth = find_lower(sketch[u], sketch[w], y);
							//else 
						}
					}
				}
				float new_sim = float(cn) / (min(uw_deg - cn, k_));
				i->second.edcn = cn;
				// rnk[v] > y only insert rnk[v] to sketch[u]
				BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
				if (new_sim == old_sim) continue;

				auto itw = BOTBIN_edge_[w].find(u);
				itw->second.edcn = cn;
				itw->second.edkth = i->second.edkth;

				BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
				BOTBIN_update_core_order(w, new_sim, old_sim);
			}
		}
	}
	else {
		BOTBIN_nbr_mp_[u].clear();
		sketch_insert(sketch[u], rnk[v]);
		for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
			int w = i->first;
			int cn = i->second.edcn;
			int uw_deg = sketch[u].size() + sketch[w].size();
			float old_sim = float(cn) / (min(uw_deg - cn - 1, k_));
			float new_sim = old_sim;
			int y = i->second.edkth;
			if (sketch_search(w, rnk[v])) {
				if (rnk[v] <= y) {
					cn += 1;
					i->second.edcn++;
					new_sim = float(cn) / min(uw_deg - cn, k_);
				}
			}
			else {
				if (rnk[v] > y) {
					new_sim = float(cn) / min(uw_deg - cn, k_);
				}
				else {
					// INF can't be searched from sketch
					if (sketch_search(u, y) && sketch_search(w, y)) {
						cn -= 1;
						i->second.edcn--;
						new_sim = float(cn) / min(uw_deg - cn, k_);
					}
					if (y == INF) {
						if (uw_deg - cn == k_) {
							i->second.edkth = max(sketch[u][sketch[u].size() - 1], sketch[w][sketch[w].size() - 1]);
						}
					}
					else {
						i->second.edkth = find_lower(sketch[u], sketch[w], y);
					}
				}
			}
			BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			if (new_sim == old_sim) continue;
			auto itw = BOTBIN_edge_[w].find(u);
			itw->second.edcn = cn;
			itw->second.edkth = i->second.edkth;
			BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
			BOTBIN_update_core_order(w, new_sim, old_sim);
		}
	}
}

void Graph::BOTBIN_update_all_core_order_del(int u, float sim, bool opt) {
	if (opt) {
		int cur_eps = find_eps_idx(sim);
		for (; cur_eps != len_eps_; cur_eps++) {
			
			auto rit = BOTBIN_bucket_[cur_eps].find(make_pair(BOTBIN_cnt[u][cur_eps], u));
			BOTBIN_bucket_[cur_eps].erase(rit);
			BOTBIN_cnt[u][cur_eps]--;
			BOTBIN_bucket_[cur_eps].insert(make_pair(BOTBIN_cnt[u][cur_eps], u));
		}
	}
	else {
		for (int i = 0; i < len_eps_; i++) {
			if (BOTBIN_cnt[u][i] == 0) continue;
			auto rit = BOTBIN_bucket_[i].find(make_pair(BOTBIN_cnt[u][i], u));
			BOTBIN_bucket_[i].erase(rit);
			BOTBIN_cnt[u][i] = 0;
		}
		int cur_eps = 0, cnt = 0;
		float sim_;
		if (BOTBIN_nbr_mp_[u].size() == 0)
			return;
		for (auto j = BOTBIN_nbr_mp_[u].begin(); j != BOTBIN_nbr_mp_[u].end();cnt++,j++) {
			sim_ = j->first;
			if (cnt == 0) {
				while (sim_ < eps_idx[cur_eps]) cur_eps++;
			}
			else {
				if (sim_ < eps_idx[cur_eps]) {
					while (sim_ < eps_idx[cur_eps]) {
						BOTBIN_cnt[u][cur_eps] = cnt;
						BOTBIN_bucket_[cur_eps].insert(make_pair(cnt, u));
						cur_eps++;
					}
				}
			}
		}
		while (sim_ < eps_idx[cur_eps]) {
			cur_eps++;
		}
		while (cur_eps != len_eps_) {
			BOTBIN_cnt[u][cur_eps] = BOTBIN_edge_[u].size();
			BOTBIN_bucket_[cur_eps].insert(make_pair(BOTBIN_edge_[u].size(), u));
			cur_eps++;
		}

	}
}


void Graph::BOTBIN_update_all_core_order_add(int u, float sim, bool opt) {

	if (opt) {
		int cur_eps = find_eps_idx(sim);
		for (; cur_eps != len_eps_; cur_eps++) {
			if (BOTBIN_cnt[u][cur_eps]) {
				auto rit = BOTBIN_bucket_[cur_eps].find(make_pair(BOTBIN_cnt[u][cur_eps], u));
				BOTBIN_bucket_[cur_eps].erase(rit);
			}
			BOTBIN_cnt[u][cur_eps]++;
			BOTBIN_bucket_[cur_eps].insert(make_pair(BOTBIN_cnt[u][cur_eps], u));
		}
	}
	else {
		for (int i = 0; i < len_eps_; i++) {
			if (BOTBIN_cnt[u][i] == 0) continue;
			auto rit = BOTBIN_bucket_[i].find(make_pair(BOTBIN_cnt[u][i], u));
			BOTBIN_bucket_[i].erase(rit);
			BOTBIN_cnt[u][i] = 0;
		}
		int cur_eps = 0, cnt = 0;
		float sim_;
		if (BOTBIN_nbr_mp_[u].size() == 0)
			return;
		for (auto j = BOTBIN_nbr_mp_[u].begin(); j != BOTBIN_nbr_mp_[u].end(); cnt++, j++) {
			sim_ = j->first;
			if (cnt == 0) {
				while (sim_ < eps_idx[cur_eps]) cur_eps++;
			}
			else {
				if (sim_ < eps_idx[cur_eps]) {
					while (sim_ < eps_idx[cur_eps]) {
						BOTBIN_cnt[u][cur_eps] = cnt;
						BOTBIN_bucket_[cur_eps].insert(make_pair(cnt, u));
						cur_eps++;
					}
				}
			}
		}
		while (sim_ < eps_idx[cur_eps]) {
			cur_eps++;
		}
		while (cur_eps != len_eps_) {
			BOTBIN_cnt[u][cur_eps] = BOTBIN_edge_[u].size();
			BOTBIN_bucket_[cur_eps].insert(make_pair(BOTBIN_edge_[u].size(), u));
			cur_eps++;
		}

	}
}

int find_element_sketch(vector<int>& v1, int& _s) {
	int l = 0, r = v1.size(), mid;
	while (l <= r) {
		mid = (l + r) >> 1;
		if (v1[mid] == _s)
			return mid;
		else if (v1[mid] < _s)
			l = mid + 1;
		else r = mid - 1;
	}
	return -1;
}

int find_upper(vector<int>& s1, vector<int>& s2, const int& y) {
	int l = 0, r = s1.size() - 1, mid;
	int res = INF;
	while (l <= r) {
		mid = (l + r) >> 1;
		//cout << mid << " " << s1[mid] << "\n";
		if (s1[mid] > y) {
			res = s1[mid];
			r = mid - 1;
		}
		else l = mid + 1;
	}

	l = 0, r = s2.size() - 1;
	while (l <= r) {
		mid = (l + r) >> 1;
		//cout << mid << " " << s2[mid];
		if (s2[mid] > y) {
			//if (res != -1)
			res = min(s2[mid], res);
			//else
			//	res = s2[mid];
			r = mid - 1;
		}
		else
			l = mid + 1;
	}
	return res;
}
//long long qk_vis = 0;
bool Graph::BOTBIN_del_sub(int u, int v) {
	if (sketch[u].size() >= k_ && BOTBIN_edge_[u].size() >= k_) {
		// must add a element to sketch[u]
		int pos = find_element_sketch(sketch[u], rnk[v]);
		//if (u == 8024)
		//	cout <<u<<" "<<pos<< "!!\n";
		if (pos == -1) {
			qk_vis += BOTBIN_edge_[u].size();
			//cout << "quick_del" << " " << u << " " << v << "\n";
			return true;
		}
		else {
			int pre = sketch[u][k_ - 1], x = n_;
			// find the rnk[w] upper than pre
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				if (rnk[i->first] > pre && rnk[i->first] < x)
					x = rnk[i->first];
			}
			if (rnk[u] > pre&& rnk[u] < x)
				x = rnk[u];
			// erase pos & add res
			sketch[u].erase(sketch[u].begin() + pos);
			sketch[u].push_back(x);

			BOTBIN_nbr_mp_[u].clear();
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				int w = i->first;
				if (w == v) continue;
				int cn = i->second.edcn;

				int uw_deg = sketch[u].size() + sketch[w].size();
				float old_sim = float(cn) / (min(uw_deg - cn, k_));
				int y = i->second.edkth;
				int new_y = y;
				if (rnk[v] <= y) {
					if (sketch_search(w, rnk[v])) 
						cn -= 1;
					else {
						new_y = find_upper(sketch[u], sketch[w], y);			
						if (sketch_search(u,new_y) && sketch_search(w, new_y)) 
							cn += 1;
						i->second.edkth = new_y;
					}
				}
				float new_sim = float(cn) / (min(uw_deg - cn, k_));
				BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
				auto itw = BOTBIN_edge_[w].find(u);

				i->second.edcn = cn;
				i->second.edkth = new_y;

				itw->second.edcn = cn;
				itw->second.edkth = new_y;
				if (old_sim == new_sim) 
					continue;
				BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
				BOTBIN_update_core_order(w, new_sim, old_sim);
			}
		}
	}
	else {
		BOTBIN_nbr_mp_[u].clear();
		//print_sketch(u);
		//cout << rnk[v] << "\n";
		//sketch_insert(sketch[u], rnk[v]);
		int pos = find_element_sketch(sketch[u], rnk[v]);
		sketch[u].erase(sketch[u].begin() + pos);

		for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {

			int w = i->first;
			if (w == v) continue;


			int cn = i->second.edcn;
			int uw_deg = sketch[u].size() + sketch[w].size();
			float old_sim = float(cn) / (min(uw_deg - cn + 1, k_));
			//float new_sim = old_sim;
			int y = i->second.edkth;
			int new_y = y;
			if (rnk[v] <= y) {
				if (sketch_search(w, rnk[v]))
					cn -= 1;
				else {
					new_y = find_upper(sketch[u], sketch[w], y);
					//cout << y << " " << new_y << "\n";
					if (new_y != INF) {
						if (sketch_search(u, new_y) && sketch_search(w, new_y))
							cn += 1;
					}
				}
			}

			float new_sim = float(cn) / (min(uw_deg - cn, k_));
			auto itw = BOTBIN_edge_[w].find(u);

			i->second.edcn = cn;
			i->second.edkth = new_y;
			
			itw->second.edcn = cn;
			itw->second.edkth = new_y;

			BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			if (new_sim == old_sim) continue;
			//if (u == 392410)
			////printf("%d %.6f %.6f\n", w, old_sim, new_sim);
			BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
			BOTBIN_update_core_order(w, new_sim, old_sim);
		}
	}
	return false;
}

bool Graph::BOTBIN_del_sub_naive(int u, int v) {
	if (sketch[u].size() >= k_ && BOTBIN_edge_[u].size() >= k_) {
		int pos = find_element_sketch(sketch[u], rnk[v]);
		if (pos == -1) {
			qk_vis += BOTBIN_edge_[u].size();
			return true;
		}
		else {
			int pre = sketch[u][k_ - 1], x = n_;
			// find the rnk[w] upper than pre
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				if (rnk[i->first] > pre && rnk[i->first] < x)
					x = rnk[i->first];
			}
			if (rnk[u] > pre && rnk[u] < x)
				x = rnk[u];
			// erase pos & add res
			sketch[u].erase(sketch[u].begin() + pos);
			sketch[u].push_back(x);

			BOTBIN_nbr_mp_[u].clear();
			for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {
				int w = i->first;
				if (w == v) continue;
				int cn = i->second.edcn;

				int uw_deg = sketch[u].size() + sketch[w].size();
				float old_sim = float(cn) / (min(uw_deg - cn, k_));
				int y = i->second.edkth;
				int new_y = y;

				int uw_kth = find_kthelement_ab(sketch[u], sketch[w], k_);
				cn = cal_intersect_ab_k(sketch[u], sketch[w], uw_kth);

				float new_sim = float(cn) / (min(uw_deg - cn, k_));
				BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
				auto itw = BOTBIN_edge_[w].find(u);

				i->second.edcn = cn;
				i->second.edkth = new_y;

				itw->second.edcn = cn;
				itw->second.edkth = new_y;
				if (old_sim == new_sim)
					continue;
				BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
				BOTBIN_update_core_order(w, new_sim, old_sim);
			}
		}
	}
	else {
		BOTBIN_nbr_mp_[u].clear();

		int pos = find_element_sketch(sketch[u], rnk[v]);
		sketch[u].erase(sketch[u].begin() + pos);

		for (auto i = BOTBIN_edge_[u].begin(); i != BOTBIN_edge_[u].end(); i++) {

			int w = i->first;
			if (w == v) continue;


			int cn = i->second.edcn;
			int uw_deg = sketch[u].size() + sketch[w].size();
			float old_sim = float(cn) / (min(uw_deg - cn + 1, k_));
			//float new_sim = old_sim;
			//cout << w <<" "<<BOTBIN_nbr_mp_[w].size()<< " " << old_sim << " " << cn <<" "<<uw_deg<< "\n";
			int y = i->second.edkth;
			int new_y = y;

			int uw_kth = find_kthelement_ab(sketch[u], sketch[w], k_);
			cn = cal_intersect_ab_k(sketch[u], sketch[w], uw_kth);


			float new_sim = float(cn) / (min(uw_deg - cn, k_));
			auto itw = BOTBIN_edge_[w].find(u);

			i->second.edcn = cn;
			i->second.edkth = new_y;

			itw->second.edcn = cn;
			itw->second.edkth = new_y;

			BOTBIN_nbr_mp_[u].insert(make_pair(new_sim, w));
			if (new_sim == old_sim) continue;

			BOTBIN_update_nbr_order(w, u, new_sim, old_sim);
			BOTBIN_update_core_order(w, new_sim, old_sim);
		}
	}
	return false;
}

void Graph::BOTBIN_del_naive(int u, int v) {
	auto itu = BOTBIN_edge_[u].find(v);

	int uv_kth = find_kthelement_ab(sketch[u], sketch[v], k_);
	int cn = cal_intersect_ab_k(sketch[u], sketch[v], uv_kth);
	
	float uv_sim = float(itu->second.edcn) /
		min(int(sketch[u].size() + sketch[v].size() - itu->second.edcn), k_);

	bool core_u_modify = BOTBIN_del_sub_naive(u, v);
	bool core_v_modify = BOTBIN_del_sub_naive(v, u);

	BOTBIN_edge_[u].erase(itu);
	auto itv = BOTBIN_edge_[v].find(u);
	BOTBIN_edge_[v].erase(itv);
	
	auto it_mpu = BOTBIN_nbr_mp_[u].find(make_pair(uv_sim, v));
	auto it_mpv = BOTBIN_nbr_mp_[v].find(make_pair(uv_sim, u));
	if (it_mpu != BOTBIN_nbr_mp_[u].end())
		BOTBIN_nbr_mp_[u].erase(it_mpu);
	if (it_mpv != BOTBIN_nbr_mp_[v].end())
		BOTBIN_nbr_mp_[v].erase(it_mpv);

	BOTBIN_update_all_core_order_del(u, uv_sim, core_u_modify);
	BOTBIN_update_all_core_order_del(v, uv_sim, core_v_modify);
}



void Graph::BOTBIN_del(int u, int v) {
	auto itu = BOTBIN_edge_[u].find(v);
	//if (u == 3785 || v== 3785) {
	//	auto itv = BOTBIN_nbr_[v].find(u);
	//}
	//cout <<u<<" "<<v<<" "<< sketch[u].size() << " " << sketch[v].size() << " " << itu->second.edcn << "\n";
	float uv_sim = float(itu->second.edcn) /
		min(int(sketch[u].size() + sketch[v].size() - itu->second.edcn), k_);

	bool core_u_modify = BOTBIN_del_sub(u, v);
	bool core_v_modify = BOTBIN_del_sub(v, u);

	BOTBIN_edge_[u].erase(itu);
	auto itv = BOTBIN_edge_[v].find(u);
	BOTBIN_edge_[v].erase(itv);
	//cout << u << " " << v << ":edcn " << itu->second.edcn << " " << itv->second.edcn <<" "<<uv_sim<< "\n";

	auto it_mpu = BOTBIN_nbr_mp_[u].find(make_pair(uv_sim, v));
	auto it_mpv = BOTBIN_nbr_mp_[v].find(make_pair(uv_sim, u));
	if(it_mpu != BOTBIN_nbr_mp_[u].end())
		BOTBIN_nbr_mp_[u].erase(it_mpu);
	if (it_mpv != BOTBIN_nbr_mp_[v].end())
		BOTBIN_nbr_mp_[v].erase(it_mpv);

	BOTBIN_update_all_core_order_del(u, uv_sim, core_u_modify);
	BOTBIN_update_all_core_order_del(v, uv_sim, core_v_modify);
}

void Graph::BOTBIN_add(int u, int v) {
	// remember N[a] = deg[a] + 1
	// nbr_ is different from sketch
	//assert(sketch[u].size() == k_);
	//cout << u << " " << v << "\n";

	bool core_u_modify = BOTBIN_add_sub(u, v);
	bool core_v_modify = BOTBIN_add_sub(v, u);
	// less kth element
	int uv_kth = find_kthelement_ab(sketch[u], sketch[v], k_);
	int cn = cal_intersect_ab_k(sketch[u], sketch[v], uv_kth);
	float uv_sim = float(cn) / min(int(sketch[u].size() + sketch[v].size() - cn), k_);
	BOTBIN_edge_[u].insert(make_pair(v, ed(cn, uv_kth)));
	BOTBIN_edge_[v].insert(make_pair(u, ed(cn, uv_kth)));

	BOTBIN_nbr_mp_[u].insert(make_pair(uv_sim, v));
	BOTBIN_nbr_mp_[v].insert(make_pair(uv_sim, u));

	//nbr_order_[u].push_back(uv_sim);
	//nbr_order_[v].push_back(uv_sim);
	//cout << u << " " << v << " " << uv_sim << "\n";
	BOTBIN_update_all_core_order_add(u, uv_sim, core_u_modify);
	BOTBIN_update_all_core_order_add(v, uv_sim, core_v_modify);

	return;
}


void Graph::BOTBIN_add_naive(int u, int v) {
	bool core_u_modify = BOTBIN_add_sub_naive(u, v);
	bool core_v_modify = BOTBIN_add_sub_naive(v, u);
	// less kth element
	int uv_kth = find_kthelement_ab(sketch[u], sketch[v], k_);
	int cn = cal_intersect_ab_k(sketch[u], sketch[v], uv_kth);
	float uv_sim = float(cn) / min(int(sketch[u].size() + sketch[v].size() - cn), k_);
	BOTBIN_edge_[u].insert(make_pair(v, ed(cn, uv_kth)));
	BOTBIN_edge_[v].insert(make_pair(u, ed(cn, uv_kth)));

	BOTBIN_nbr_mp_[u].insert(make_pair(uv_sim, v));
	BOTBIN_nbr_mp_[v].insert(make_pair(uv_sim, u));

	//nbr_order_[u].push_back(uv_sim);
	//nbr_order_[v].push_back(uv_sim);
	//cout << u << " " << v << " " << uv_sim << "\n";
	BOTBIN_update_all_core_order_add(u, uv_sim, core_u_modify);
	BOTBIN_update_all_core_order_add(v, uv_sim, core_v_modify);

	return;
}


