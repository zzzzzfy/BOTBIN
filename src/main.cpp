#include "Utility.h"
#include "Graph.h"
#include "Config.h"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
void check_abs() {
	map<pair<int, int>, float> mp;
	freopen("/home/user_name/idx1.log", "r", stdin);
	int u, v; float s;
	while (scanf("%d%d%f", &u, &v, &s) != EOF) {
		mp[make_pair(u, v)] = s;
	}
	fclose(stdin);

	freopen("/home/user_name/BOTBIN.log", "r", stdin);
	//int u, v; float s;
	float mx = 0;
	while (scanf("%d%d%f", &u, &v, &s) != EOF) {
		mx = max(abs(s - mp[make_pair(u, v)]), mx);
		cout << u << "," << v << " " << abs(s - mp[make_pair(u, v)]) << "\n";
	}
	cout << mx << endl;
	fclose(stdin);
	exit(0);
}

int intersection(vector<int>& v1,vector<int>& v2) {
	vector<int> v3;
	set_intersection(v1.begin(), v1.end(),
		v2.begin(), v2.end(),
		back_inserter(v3));
	return v3.size();
}

long long calc_2(int a) {
	if (a <= 1) return 0;
	else return ((long long)a * (a - 1)) / 2;
}

double calc_ari(vector<vector<int>> &cluster1, vector<vector<int>> &cluster2, int n) {
	//long long n = 0;
	long long tot1 = 0;//nij 2
	long long tota = 0;
	long long totb = 0;
	map<int, int> cls;
	for (int i = 0; i < cluster1.size(); i++) {
		for (int j = 0; j < cluster1[i].size(); j++)
			cls[cluster1[i][j]] = i;
	}
	map<int, int> cnt;
	for (int i = 0; i < cluster2.size(); i++) {
		cnt.clear();
		for (int j = 0; j < cluster2[i].size(); j++) {
			int ele = cluster2[i][j];
			//if (cls.find(ele) == cls.end()) {
			//	cout << "！！！！" << "\n";
			//}
			if (cnt.find(cls[ele]) != cnt.end()) {
				cnt[cls[ele]]++;
			}
			else cnt[cls[ele]] = 1;
		}
		for (auto it = cnt.begin(); it != cnt.end(); it++) {
			tot1 += calc_2(it->second);
		}
			//cls[cluster2[i][j]] = i;
	}
	
	for (int i = 0; i < cluster1.size(); i++) {
		//n += cluster1[i].size();
		tota += calc_2(cluster1[i].size());

	}
	for (int j = 0; j < cluster2.size(); j++) {
		totb += calc_2(cluster2[j].size());
	}
	//exit(0);
	cout << tot1 << " " << tota << " " << totb << " " << n << "\n";
	if (((tota + totb) * calc_2(n) - 2 * (tota + totb)) == 0) return -1;
	return (tot1 - ((double)tota * totb) / calc_2(n)) / ((tota + totb) * 0.5 -((double)tota * totb) / calc_2(n));
}
vector<vector<int>> cluster1;
vector<vector<int>> cluster2;
vector<int> eps_list;
vector<int> mu_list;
void check_ari(string st, int n_, int rho, int delta) {
	int tot_calc = 0;
	double tot_pre = 0, tot_rec = 0, tot_ari = 0;
	double min_ari = 1;
	pair<int, int> min_ptr;
	map<int, bool> cluster_map;
	long long core_num1 = 0;
	long long core_num2 = 0;
	double ari_3 = 0, ari_5 = 0;
	double avr_ari[100];
	for (int i = 0; i < 100; i++) avr_ari[i] = 0.0;
	//int core_num1, core_num;
	for (int T = 0; T < eps_list.size(); T++) {
		int eps = eps_list[T];
		int mu = mu_list[T];
		string file1 = st + "_result/_1_" + to_string(eps) + "_" + to_string(mu);
		string file2 = st + "_result/_2_" + to_string(eps) + "_" + to_string(mu);
		ifstream tmpf;
		tmpf.open(file1.c_str(), ios::app);
		int tot_cluster_1, cluster_size, id, tot_node_1 = 0, tot_node_2 = 0;
		tmpf >> tot_cluster_1;

		cout << eps << " " << mu <<" "<<file1 << "\n";
		cout << tot_cluster_1 << "\n";
		cluster1.clear();
		cluster2.clear();
		if (tot_cluster_1 <= 1 ) {
			fclose(stdin);
			continue;
		}
		for (int i = 0; i < tot_cluster_1; i++) {
			tmpf >> cluster_size;
			core_num1 += cluster_size;
			tot_node_1 += cluster_size;
			cluster1.push_back(vector<int>());
			for (int j = 0; j < cluster_size; j++) {
				tmpf >> id;
				cluster1[i].push_back(id);
				cluster_map[id] = true;
			}
			sort(cluster1[i].begin(), cluster1[i].end());
			
		}
		tmpf.close();
		tmpf.open(file2.c_str(), ios::app);

		int precision_hit = 0, recall_hit = 0;
		int tot_cluster_2;
		tmpf >> tot_cluster_2;
		vector<int> other_cluster;
		
		//cout << tot_cluster_2 << "\n";
		for (int i = 0; i < tot_cluster_2; i++) {
			tmpf >> cluster_size;
			core_num2 += cluster_size;
			tot_node_2 += cluster_size;
			cluster2.push_back(vector<int>());
			//cout << cluster2[i].size() << "\n";
			for (int j = 0; j < cluster_size; j++) {
				tmpf >> id;
				if (cluster_map.find(id)!=cluster_map.end()) {
					cluster2[i].push_back(id);
					cluster_map[id] = false;
				}
				else {
					// noise can be ignored in different setting 
				}
				//cluster_map[id] = true;
			}

			sort(cluster2[i].begin(), cluster2[i].end());
		}
		//cluster1.push_back(other_cluster);
		other_cluster.clear();
		for (auto it = cluster_map.begin(); it != cluster_map.end(); it++) {
			if (it->second == true) {
				other_cluster.push_back(it->first);
			}
		}
		sort(other_cluster.begin(), other_cluster.end());
		cluster2.push_back(other_cluster);
		//cluster1.push_back(other_cluster);
		//double ari = calc_ari(cluster1, cluster2, tot_node_1 );

		double ari = calc_ari(cluster1, cluster2, tot_node_1);
		//fclose(stdin);
		if (tot_node_1 >= 1 && (ari != -1)) {
			tot_ari += ari;
			tot_calc++;
			if (ari < min_ari) {
				min_ari = ari;
				min_ptr = make_pair(eps, mu);
			}
		}
		if (eps == 3) ari_3 += ari;
		else if (eps == 5) ari_5 += ari;
		avr_ari[eps] += ari;
		cout <<eps<<" "<<mu<<" "<< ari <<" cluster_1:"<<tot_node_1<<" cluster_2:" << tot_node_2 << "\n";
		//delete[] cluster_map;
		cluster_map.clear();
		//exit(0);
	}
	ofstream tmpf;
	tmpf.open("ari_result.txt", ios::app);
	tmpf << st << " rho:" << rho <<" delta:"<<delta << "\n";
	tmpf << "ari:" << tot_ari << " " << tot_calc << " " << tot_ari / tot_calc << "\n";
	tmpf << "min:" << min_ptr.first << " " << min_ptr.second << " " << min_ari << "\n";
	tmpf.close();
	cout << "ari:" << tot_ari << " " << tot_calc << " " << tot_ari / tot_calc << "\n";
	cout << "min:" << min_ptr.first << " " << min_ptr.second << " " << min_ari<<"\n";
	cout << "core_num1:" << core_num1 << " " << "core_num2:" << core_num2 << "\n";
	cout << "ari3:" << ari_3 / 14 << " ari5:" << ari_5 / 14;
	for (int i = 20; i <= 80; i++)
		cout <<i<<":" << avr_ari[i] / 14 << "\n";
	//cout << tot_calc << "pres:" << tot_pre / tot_calc << " recall:" << tot_rec / tot_calc << "\n";
	exit(0);
}

void check_auc(string st,int n_, int rho, int delta) {
	int tot_calc = 0;
	double tot_pre = 0, tot_rec = 0;
	double min_pre = 1, min_rec = 1;
	pair<int, int> min_pre_ptr;
	pair<int, int> min_rec_ptr;
	map<int, bool> cluster_map;
	int core_num1=0, core_num2=0;
	for (int T = 0; T < eps_list.size(); T++) {
	//for (int eps = 1; eps <= 99; eps++) {
		//for (int mu = 2; mu <= 20; mu++) {
		int eps = eps_list[T];
		int mu = mu_list[T];
		string file1 = st + "_result/_1_" + to_string(eps) + "_" + to_string(mu);
		string file2 = st + "_result/_2_" + to_string(eps) + "_" + to_string(mu);
		cout << eps << " " << mu <<" "<<file1<<" "<<file2 << "\n";

		ifstream tmpf;
		tmpf.open(file1.c_str(), ios::app);
		int tot_cluster_1 = 0, cluster_size = 0, id = 0, tot_node_1 = 0;
		tmpf >> tot_cluster_1;
		//bool* cluster_map = new bool[50000000];
		int mx_cluster = 0;
		//memset(cluser_map,false,sizeof clu)
		map<int, int> cluster_label1;
		map<int, int> cluster_label2;
		//int* ari_map = new int[50000000];
		for (int i = 0; i < tot_cluster_1; i++) {
			tmpf >> cluster_size;
			tot_node_1 += cluster_size;
			mx_cluster = max(mx_cluster, cluster_size);
			for (int j = 0; j < cluster_size; j++) {
				tmpf >> id;
				core_num1++;
				cluster_map[id] = true;
			}
		}

		tmpf.close();
		tmpf.open(file2.c_str(), ios::app);
		int precision_hit = 0, recall_hit = 0, tot_node_2 = 0;
		int tot_cluster_2;
		tmpf >> tot_cluster_2;
		for (int i = 0; i < tot_cluster_2; i++) {
			tmpf >> cluster_size;
			tot_node_2 += cluster_size;
			for (int j = 0; j < cluster_size; j++) {
				tmpf >> id; 
				core_num2++;
				if (cluster_map.find(id)!=cluster_map.end()) {
					precision_hit++;
				}
				//cluster_map[id] = true;
			}
		}
		//cout << tot_node_1 << " " << tot_node_2 << "\n";
		if (tot_node_1 >= 1 && tot_node_2 >= 1) {
			double tmp_pre = double(precision_hit) / tot_node_2;
			double tmp_rec = double(precision_hit) / tot_node_1;
			tot_pre += double(precision_hit) / tot_node_2;
			tot_rec += double(precision_hit) / tot_node_1;
			//cout << (double)tot_node_1 / n_<<"\n";
			if (tot_node_1 != tot_node_2)
				cerr << "!!!" << '\n';
			cout << eps << " " << mu << " " <<tot_node_1<<" "<<tot_node_2<<" " << tmp_pre << " " << tmp_rec << " " << (double)mx_cluster / n_ << "\n";
			if (tmp_pre < min_pre) {
				min_pre = tmp_pre;
				min_pre_ptr = make_pair(eps, mu);
			}
			if (tmp_rec < min_rec) {
				min_rec_ptr = make_pair(eps, mu);
				min_rec = tmp_rec;
			}
			tot_calc++;
		}
		cluster_map.clear();
		//delete[] cluster_map;
		
	}	

	ofstream tmpf;
	tmpf.open("prerec_result.txt", ios::app);
	tmpf << st << " rho:" << rho << " delta:" << delta << "\n";
	//tmpf << st << " " << rho << "\n";
	tmpf << tot_calc << "pres:" << tot_pre / tot_calc << " recall:" << tot_rec / tot_calc << "\n";
	tmpf << "recall:" << min_rec_ptr.first << " " << min_rec_ptr.second << " " << min_rec << "\n";
	tmpf << "precision:" << min_pre_ptr.first << " " << min_pre_ptr.second << " " << min_pre << "\n";
	tmpf.close();
	//cout << st << "\n";
	cout <<tot_calc<< "pres:" << tot_pre/tot_calc << " recall:" << tot_rec / tot_calc << "\n";
	cout <<"recall:" << min_rec_ptr.first << " " << min_rec_ptr.second << " " << min_rec << "\n";
	cout <<"precision:" << min_pre_ptr.first << " " << min_pre_ptr.second << " " << min_pre << "\n";
	cout << "core_num1:" << core_num1 << " " << "core_num2:" << core_num2 << "\n";
	//<< " " << min_pre;
	exit(0);
}
int calc_k(double eps, double pf, double M, double avg_deg) {
	return (1 / (2 * eps * eps)) * log(2 / (pf * M * avg_deg));
}


void init_query() {
	int eps, mu;
	freopen("/home/user_name/dataset/query_config", "r", stdin);
	while(cin>>eps>>mu) {
		//cin >> eps >> mu;
		// 
		//cout << eps << " " << mu << "\n";
		eps_list.push_back(eps);
		mu_list.push_back(mu);
	}
	fclose(stdin);
	//exit(0);
}

void create_query(int low_eps, int up_eps, int low_mu, int up_mu) {
	int eps, mu;
	freopen("/home/user_name/dataset/query_config", "w", stdout);
	for (int i = 1; i <= 100; i++) {
		//cin >> eps >> mu;
		eps = rand() % (up_eps - low_eps + 1) + low_eps;
		mu = rand() % (up_mu - low_mu + 1) + low_mu;
		cout << eps << " " << mu << "\n";
		eps_list.push_back(eps);
		mu_list.push_back(mu);
	}
	fclose(stdout);
	exit(0);
}
int cnt = 0;

int main(int argc, char* args[]) {
	init_query();
	Config config;
	int alg;
	int _file = 1;
	int rho = 100;
	int action = 0;
	bool is_build = false;
	bool is_query = false;
	bool is_naive = false;
	int delta = 100; 
	config.action_ = "BOTBIN";
	while (cnt < argc) {
		if (strcmp(args[cnt], "-action") == 0) {
			action = atoi(args[++cnt]);
			if (action == 0) {
				config.action_ = "format";
			}
			if (action == 1) {
				config.action_ = "BOTBIN";
			}
			if (action == 2) {
			   	config.action_ = "check_ari";
			}
			if (action == 3) {
				config.action_ = "check_auc";
			}
			if (action == 4) {
				config.action_ = "sample_edges";
			}
			if (action == 5) {
				config.action_ = "cluster_ari";
			}
		}
		else if (strcmp(args[cnt], "-file") == 0) {
			// Paths that require manual differentiation of filenames
			// Just File ID
			_file = atoi(args[++cnt]);
		}
		else if (strcmp(args[cnt], "-query") == 0) {
			is_query = true;
		}
		else if (strcmp(args[cnt], "-rho") == 0) {
			rho = atoi(args[++cnt]);
		}
		else if (strcmp(args[cnt], "-delta") == 0) {
			delta = atoi(args[++cnt]);
		}
		else if (strcmp(args[cnt], "-eta") == 0) {
			eta = atoi(args[++cnt]);
		}
		else if (strcmp(args[cnt], "-build") == 0) {		
			is_build = true;
		}
		cnt++;
	}

	
	config.set_file(_file);
	cout << config.dest_ << "\n";



	if (config.action_ == "format") {
		// clean the graph data
		Utility ut;
		ut.format_graph(config.dest_, config.src_, 100, false);
		return 0;
	}
	Graph* g = new Graph(config.dest_, double(rho) / 1000);

	if (config.action_ == "check_auc") {
		check_auc(config.dest_, config.n_, rho, delta);
		return 0;
	}
	if (config.action_ == "check_ari") {
		check_ari(config.dest_, config.n_, rho, delta);
		return 0;
	}
	if (config.action_ == "cluster_ari") {
		g->cluster_ari(config.dest_, eps_list, mu_list);
		return 0;
	}
	if (config.action_ == "sample_edges") {
		g->sample_edges();
		return 0;
	}

	if (config.action_ == "BOTBIN") {
		g->dynamic_BOTBIN(is_query, config.dest_, eps_list, mu_list, is_build, delta,
			is_naive);
	}

	return 0;
}