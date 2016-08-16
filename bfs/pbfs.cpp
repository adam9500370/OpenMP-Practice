#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "queue.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "time.h"
#include <omp.h>
#include <math.h>

int N, M;

std::vector<double> ii_arr, jj_arr, norm_arr;
std::vector<bool> ii_bit, jj_bit;
std::vector<int> p_v, p_e;

std::vector<int> edge_list[2]; // 0 for start vertex, 1 for end vertex
std::vector<int> in_deg;
std::vector<int> rand_search_key, search_key;
std::vector<int> merge_queue;
std::vector<Queue_Base*> current_queues, next_queues;
std::vector<bool> visited;
std::vector<std::vector<int> > adjacency_list;
std::vector<int> node_level, node_parent;

std::vector<int> next_queue_size, current_queue_size;

int max_threads = omp_get_max_threads();
bool flag_merge;

void initialize() {
	srand(time(NULL));
	
	in_deg = std::vector<int> (N);
	adjacency_list = std::vector<std::vector<int> > (N);
	//node_level = std::vector<int> (N, -1);
	
	current_queues = std::vector<Queue_Base*> (max_threads);
	for(int i = 0; i < current_queues.size(); i++) {
		current_queues[i] = new Circular_Queue;
		current_queues[i]->clear();
		current_queues[i]->reserve(N + 1);
	}
	next_queues = std::vector<Queue_Base*> (max_threads);
	for(int i = 0; i < next_queues.size(); i++) {
		next_queues[i] = new Circular_Queue;
		next_queues[i]->clear();
		next_queues[i]->reserve(N + 1);
	}
	
	next_queue_size = std::vector<int> (N);
	current_queue_size = std::vector<int> (N);
}

void randperm(int n, std::vector<int> &perm) {
	int j, tmp;
	for(int i = 0; i < n; i++) {
		perm[i] = i;
	}
	
	for(int i = 0; i < n; i++) {
		j = ((double)rand() / (1.0 + RAND_MAX)) * (n-i) + i;
		// swap i, j
		tmp = perm[j];
		perm[j] = perm[i];
		perm[i] = tmp;
	}
}

void rand_arr(int len, std::vector<double> &arr) {
	for(int i = 0; i < len; i++) {
		arr[i] = (double)(rand() % RAND_MAX) / RAND_MAX;
	}
}

void cmp_gt(int len, const std::vector<double> &arr1, const std::vector<double> &arr2, std::vector<bool> &result_arr) { // compare whether arr1 > arr2 or not
	for(int i = 0; i < len; i++) {
		result_arr[i] = (arr1[i] > arr2[i]);
	}
}

void cmp_gt(int len, const std::vector<double> &arr1, double value, std::vector<bool> &result_arr) { // compare whether arr1 > value or not
	for(int i = 0; i < len; i++) {
		result_arr[i] = (arr1[i] > value);
	}
}

void kronecker_generator(int SCALE, int edgefactor) {
	// Set number of vertices.
	N = pow(2, SCALE);
	
	// Set number of edges.
	M = edgefactor * N;
	
	initialize();
	
	// Set initiator probabilities.
	double A = 0.57, B = 0.19, C = 0.19;
	
	// Create index arrays.
	edge_list[0] = std::vector<int> (M);
	edge_list[1] = std::vector<int> (M);
	// Loop over each order of bit.
	double ab, c_norm, a_norm;
	ab = A + B;
	c_norm = C/(1 - (A + B));
	a_norm = A/(A + B);
	
	//printf("N: %d, M: %d, ab: %f, c_norm: %f, a_norm: %f\n", N, M, ab, c_norm, a_norm);
	
	ii_arr = std::vector<double> (M);
	jj_arr = std::vector<double> (M);
	norm_arr = std::vector<double> (M);
	ii_bit = std::vector<bool> (M);
	jj_bit = std::vector<bool> (M);
	for(int ib = 1; ib <= SCALE; ib++) {
		// Compare with probabilities and set bits of indices.
		rand_arr(M, ii_arr);
		cmp_gt(M, ii_arr, ab, ii_bit);
		rand_arr(M, jj_arr);
		for(int k = 0; k < M; k++) {
			norm_arr[k] = c_norm * ii_bit[k] + a_norm * (!ii_bit[k]);
		}
		cmp_gt(M, jj_arr, norm_arr, jj_bit);
		for(int k = 0; k < M; k++) {
			edge_list[0][k] += pow(2, (ib-1)) * ii_bit[k];
			edge_list[1][k] += pow(2, (ib-1)) * jj_bit[k];
		}
	}
	ii_arr.clear();
	jj_arr.clear();
	norm_arr.clear();
	ii_bit.clear();
	jj_bit.clear();
	
	// Permute vertex labels
	p_v = std::vector<int> (N);
	randperm(N, p_v);
	for(int i = 0; i < 2; i++) {
		for(int k = 0; k < M; k++) {
			edge_list[i][k] = p_v[edge_list[i][k]];
		}
	}
	p_v.clear();
	
	// Permute the edge list
	p_e = std::vector<int> (M);
	randperm(M, p_e);
	for(int i = 0; i < 2; i++) {
		for(int k = 0; k < M; k++) {
			edge_list[i][k] = edge_list[i][p_e[k]];
		}
	}
	p_e.clear();
}

void read_edges() {
	freopen("in", "r", stdin);
	
	int start_node, end_node;
	scanf("%d", &N);
	initialize();
	while(scanf("%d%d", &start_node, &end_node) != EOF) {
		edge_list[0].push_back(start_node);
		edge_list[1].push_back(end_node);
	}
}

void read_lists() {
	freopen("in", "r", stdin);
	
	int node_deg, end_node;
	scanf("%d", &N);
	initialize();
	for(int start_node = 0; start_node < N; start_node++) {
		scanf("%d", &node_deg);
		for(int i = 0; i < node_deg; i++) {
			scanf("%d", &end_node);
			edge_list[0].push_back(start_node);
			edge_list[1].push_back(end_node);
		}
	}
}

void read_graph() {
	freopen("in.graph", "r", stdin);
	
	int end_node;
	std::string line;
	scanf("%d%d%c", &N, &M);
	initialize();
	for(int start_node = 0; start_node < N; start_node++) {
		getline(std::cin, line);
		std::istringstream in(line);
		while(in >> end_node) {
			edge_list[0].push_back(start_node);
			edge_list[1].push_back(end_node-1);
		}
	}
}

/*** transform edge list into adjacency list ***/
void construct_graph(bool flag_duplicate) {
	for(int i = 0; i < edge_list[0].size(); i++) {
		if(edge_list[0][i] != edge_list[1][i]) { // no self edge, undirected graph
			adjacency_list[edge_list[0][i]].push_back(edge_list[1][i]);
			in_deg[edge_list[1][i]]++;
			if(flag_duplicate) {
				adjacency_list[edge_list[1][i]].push_back(edge_list[0][i]);
				in_deg[edge_list[0][i]]++;
			}
		}
	}
}

/*** find next_queues elements, each thread maintain their own next queue, directly divide merge queue ***/
void find_next_queues_by_merge_queue(int level) {
	#pragma omp parallel for
	for(int i = 0; i < merge_queue.size(); i++) {
		if(!visited[merge_queue[i]]) {
			for(int j = 0; j < adjacency_list[merge_queue[i]].size(); j++) {
				if(node_parent[adjacency_list[merge_queue[i]][j]] < 0) {
					node_parent[adjacency_list[merge_queue[i]][j]] = merge_queue[i];
					next_queues[ omp_get_thread_num() ]->push(adjacency_list[merge_queue[i]][j]);
				}
			}
		}
		visited[merge_queue[i]] = true;
	}
}

/*** find next_queues elements, each thread maintain their own current and next queue ***/
void find_next_queues_by_current_queues(int level) {
	#pragma omp parallel for
	for(int i = 0; i < max_threads; i++) {
		int queue_size = current_queues[i]->size();
		for(int s = 0; s < queue_size; s++) {
			int node = current_queues[i]->top();
			
			if(!visited[node]) {
				for(int j = 0; j < adjacency_list[node].size(); j++) {
					if(node_parent[adjacency_list[node][j]] < 0) {
						node_parent[adjacency_list[node][j]] = node;
						next_queues[i]->push(adjacency_list[node][j]);
					}
				}
			}
			visited[node] = true;
			
			current_queues[i]->pop();
		}
	}
}

/*** record each next queue size, and determine whether to merge or not ***/
void record_next_queues_size(int level) {
	next_queue_size[level] = 0;
	flag_merge = false;
	//printf("level %d next_queues size: ", level);
	for(int i = 0; i < max_threads; i++) {
		next_queue_size[level] += next_queues[i]->size();
		if(next_queues[i]->size() == 0) {
			flag_merge = true;
		}
		//printf("%d ",next_queues[i]->size());
	}
	//printf(", sum = %d\n", next_queue_size[level]);
	
	int average_size = next_queue_size[level]/max_threads;
	for(int i = 0; i < max_threads; i++) {
		if(next_queues[i]->size() < average_size/2) {
			flag_merge = true;
		}
	}
}

/*** next_queues swap to current_queues and clear ***/
void swap_to_current_queues(int level) {
	#pragma omp parallel for
	for(int i = 0; i < max_threads; i++) {
		Queue_Base *tmp = current_queues[i];
		current_queues[i] = next_queues[i];
		next_queues[i] = tmp;
		next_queues[i]->clear();
	}
}

/*** next_queues merge to merge_queue and clear ***/
void merge_to_merge_queue(int level) {	
	merge_queue.clear();
	merge_queue.reserve(next_queue_size[level]);
	for(int i = 0; i < max_threads; i++) {
		int queue_size = next_queues[i]->size();
		for(int j = 0; j < queue_size; j++) {
			int node = next_queues[i]->top();
			
			merge_queue.push_back(node);
			
			next_queues[i]->pop();
		}
		next_queues[i]->clear();
	}
}

void pbfs(int start_point/*, double* time_group*/) {
	int level = 1;
	/*Windows_Time *timer_1 = new Windows_Time, *timer_2 = new Windows_Time, *timer_3 = new Windows_Time;
	
	time_group[0] = 0.0;
	time_group[1] = 0.0;
	time_group[2] = 0.0;*/
	
	node_parent.clear();
	node_parent.reserve(N);
	
	visited.clear();
	visited.reserve(N);
	
	merge_queue.clear();
	merge_queue.push_back(start_point);
	
	omp_set_num_threads(max_threads);
	#pragma omp parallel
	{
		#pragma omp for nowait
		for(int i = 0; i < N; i++) {
			node_parent[i] = -1;
		}
		#pragma omp for nowait
		for(int i = 0; i < N; i++) {
			visited[i] = false;
		}
		#pragma omp for nowait
		for(int i = 0; i < max_threads; i++) {
			current_queues[i]->clear();
		}
	}
	
	node_parent[start_point] = start_point;
	next_queue_size[0] = 1;
	flag_merge = true;
	
	while(next_queue_size[level-1] > 0) {
		//timer_1->set_start_time();
		if(flag_merge) {
			find_next_queues_by_merge_queue(level);
		}
		else {
			find_next_queues_by_current_queues(level);
		}
		/*timer_1->set_end_time();
		time_group[0] += timer_1->get_time_interval();*/
			
		//timer_2->set_start_time();
		record_next_queues_size(level);
		/*timer_2->set_end_time();
		time_group[1] += timer_2->get_time_interval();*/
			
		//timer_3->set_start_time();
		if(flag_merge) {
			merge_to_merge_queue(level);
		}
		else {
			swap_to_current_queues(level);
		}
		/*timer_3->set_end_time();
		time_group[2] += timer_3->get_time_interval();*/
			
		level++;
	}
	
	/*printf("BFS level %d\n",level);
	for(int i = 0; i < N; i++) {
		printf("node: %d, parent: %d\n", i, node_parent[i]);
	}*/
	
	/*for(int i = 1; i < level; i++) {
		printf("thread next queue size = %d, merge queue size = %d, dif = %d\n", next_queue_size[i], current_queue_size[i], next_queue_size[i] - current_queue_size[i]);
	}*/
	/*for(int i = 1; i < level; i++) {
		printf("level = %d, queue size = %d\n", i, next_queue_size[i]);
	}*/
}

bool validate(int start_point) {
	if(node_parent[start_point] != start_point) {
		printf("err: -1\n");
		return false;
	}
	
	int current_level = 1;
	node_level = std::vector<int> (N);
	std::vector<int> current_level_list, next_level_list;
	next_level_list.clear();
	next_level_list.reserve(N);
	for(int i = 0; i < N; i++) {
		if(node_parent[i] != -1) {
			node_level[i] = 1;
			if(node_parent[i] != start_point) {
				next_level_list.push_back(node_parent[i]);
			}
		}
		else {
			node_level[i] = -1;
		}
	}
	current_level_list = next_level_list;
	while(!current_level_list.empty()) {
		current_level++;
		if(current_level > N) { // there must be a cycle in the tree
			printf("err: -2\n");
			return false;
		}
		
		next_level_list.clear();
		next_level_list.reserve(N);
		for(int i = 0; i < current_level_list.size(); i++) {
			node_level[current_level_list[i]]++;
			if(node_parent[current_level_list[i]] != start_point) {
				next_level_list.push_back(node_parent[current_level_list[i]]);
			}
		}
		current_level_list = next_level_list;
	}
	
	for(int i = 0; i < edge_list[0].size(); i++) {
		if( ! (	(node_level[edge_list[0][i]] == -1 && node_level[edge_list[1][i]] == -1) // both out of graph
			 || (node_level[edge_list[0][i]] != -1 && node_level[edge_list[1][i]] != -1) // both in graph
			 || (abs(node_level[edge_list[0][i]] - node_level[edge_list[1][i]]) <= 1) // each tree edge connects vertices whose BFS levels differ by exactly one
				)){
			printf("err: -3\n");
			printf("e: %d %d, l: %d %d\n", edge_list[0][i], edge_list[1][i], node_level[edge_list[0][i]], node_level[edge_list[1][i]]);
			return false;
		}
	}
	
	node_level.clear();
	return true;
}

void statistics(double* value, int NBFS, double* result_info) {
	double min, firstquartile, median, thirdquartile, max, mean, stddev;
	
	// using default comparison (operator <)
	std::sort(value, value + NBFS);
	min = value[0];
	if(NBFS%4 == 0) {
		firstquartile = (value[NBFS/4-1]+value[NBFS/4])/2;
		thirdquartile = (value[NBFS/4*3-1]+value[NBFS/4*3])/2;
	}
	else {
		firstquartile = value[NBFS/4];
		thirdquartile = value[NBFS/4*3];
	}
	if(NBFS%2 == 0) {
		median = (value[NBFS/2-1]+value[NBFS/2])/2;
	}
	else {
		median = value[NBFS/2];
	}
	max = value[NBFS-1];
	
	double sum_of_value = 0;
	for(int i = 0; i < NBFS; i++) {
		sum_of_value += value[i];
	}
	mean = sum_of_value/NBFS;
	double sum_of_value_mean_diff_square = 0;
	for(int i = 0; i < NBFS; i++) {
		sum_of_value_mean_diff_square += pow((value[i]-mean), 2);
	}
	stddev = sqrt((sum_of_value_mean_diff_square) / (NBFS-1));
	
	result_info[0] = min;
	result_info[1] = firstquartile;
	result_info[2] = median;
	result_info[3] = thirdquartile;
	result_info[4] = max;
	result_info[5] = mean;
	result_info[6] = stddev;
}

void output(int SCALE, int edgefactor, int NBFS, double kernel_1_time, double* kernel_2_time, double* kernel_2_nedge/*, double* total_time_group*/) {
	char strch_SCALE[10], strch_edgefactor[10];
	sprintf(strch_SCALE, "%d", SCALE);
	sprintf(strch_edgefactor, "%d", edgefactor);
	std::string str_SCALE(strch_SCALE), str_edgefactor(strch_edgefactor);
	std::string file_name = "SCALE=" + str_SCALE + "_edgefactor=" + str_edgefactor;
	freopen(file_name.data(), "w", stdout);
	
	printf("SCALE: %d\n", SCALE);
	printf("edgefactor: %d\n", edgefactor);
	printf("NBFS: %d\n", NBFS);
	printf("construction_time: %20.17e\n", kernel_1_time);
	
	double S[7], kernel_2_mean_time;
	
	statistics(kernel_2_time, NBFS, S);
	printf("min_time: %20.17e\n", S[0]);
	printf("firstquartile_time: %20.17e\n", S[1]);
	printf("median_time: %20.17e\n", S[2]);
	printf("thirdquartile_time: %20.17e\n", S[3]);
	printf("max_time: %20.17e\n", S[4]);
	printf("mean_time: %20.17e\n", S[5]);
	printf("stddev_time: %20.17e\n", S[6]);
	
	kernel_2_mean_time = S[5];
	
	statistics(kernel_2_nedge, NBFS, S);
	printf("min_nedge: %20.17e\n", S[0]);
	printf("firstquartile_nedge: %20.17e\n", S[1]);
	printf("median_nedge: %20.17e\n", S[2]);
	printf("thirdquartile_nedge: %20.17e\n", S[3]);
	printf("max_nedge: %20.17e\n", S[4]);
	printf("mean_nedge: %20.17e\n", S[5]);
	printf("stddev_nedge: %20.17e\n", S[6]);
	
	double TEPS[NBFS], tmp[NBFS];
	for(int i = 0; i < NBFS; i++) {
		TEPS[i] = kernel_2_nedge[i] / kernel_2_time[i];
	}
	statistics(TEPS, NBFS, S);
	double sum_of_TEPS_reciprocal = 0;
	for(int i = 0; i < NBFS; i++) {
		if(TEPS[i] > 0) {
			sum_of_TEPS_reciprocal += 1/TEPS[i];
		}
	}
	S[5] = NBFS/sum_of_TEPS_reciprocal;
	for(int i = 0; i < NBFS; i++) {
		if(TEPS[i] > 0) {
			tmp[i] = 1/TEPS[i];
		}
		else {
			tmp[i] = 0;
		}
	}
	for(int i = 0; i < NBFS; i++) {
		tmp[i] = tmp[i] - 1/S[5];
	}
	double sum_of_tmp_square = 0;
	for(int i = 0; i < NBFS; i++) {
		sum_of_tmp_square += pow(tmp[i], 2);
	}
	S[6] = (sqrt(sum_of_tmp_square) / (NBFS-1)) * pow(S[5], 2);
	
	printf("min_TEPS: %20.17e\n", S[0]);
	printf("firstquartile_TEPS: %20.17e\n", S[1]);
	printf("median_TEPS: %20.17e\n", S[2]);
	printf("thirdquartile_TEPS: %20.17e\n", S[3]);
	printf("max_TEPS: %20.17e\n", S[4]);
	printf("harmonic_mean_TEPS: %20.17e\n", S[5]);
	printf("harmonic_stddev_TEPS: %20.17e\n", S[6]);
	
	/*printf("\naverage_time_group[0]: %f (%f%%)\naverage_time_group[1]: %f (%f%%)\naverage_time_group[2]: %f (%f%%)\n"
		, total_time_group[0]/NBFS, (total_time_group[0]/NBFS)/kernel_2_mean_time*100
		, total_time_group[1]/NBFS, (total_time_group[1]/NBFS)/kernel_2_mean_time*100
		, total_time_group[2]/NBFS, (total_time_group[2]/NBFS)/kernel_2_mean_time*100);*/
}

int main() {
	int SCALE, edgefactor;
	int NBFS = 64;
	Windows_Time *timer_kernel_1 = new Windows_Time, *timer_kernel_2 = new Windows_Time;
	double kernel_1_time, kernel_2_time[NBFS], kernel_2_nedge[NBFS];
	
	printf("Please enter SCALE and edgefactor:\n");
	scanf("%d%d", &SCALE, &edgefactor);
	kronecker_generator(SCALE, edgefactor);
	
	//read_edges();
	//read_lists();
	//read_graph();
	//SCALE = log2(N);
	//edgefactor = M/N;
	bool flag_gen_edge_list = true;
	
	timer_kernel_1->set_start_time();
	construct_graph(flag_gen_edge_list);
	timer_kernel_1->set_end_time();
	kernel_1_time = timer_kernel_1->get_time_interval();
	//printf("kernel_1_time: %f\n", kernel_1_time);
	
	rand_search_key = std::vector<int> (N);
	randperm(N, rand_search_key);
	for(int i = 0; i < N; i++) {
		if(adjacency_list[rand_search_key[i]].size() > 0) {
			search_key.push_back(rand_search_key[i]);
		}
	}
	rand_search_key.clear();
	if(search_key.size() > NBFS) {
		search_key.erase(search_key.end()-(N-NBFS), search_key.end());
	}
	else {
		NBFS = search_key.size();
	}
	
	//double time_group[3], total_time_group[3] = {0.0};
	for(int i = 0; i < NBFS; i++) {
		//printf("search key: %d\n", search_key[i]);		
		timer_kernel_2->set_start_time();
		pbfs(search_key[i]/*, time_group*/);	
		timer_kernel_2->set_end_time();
		kernel_2_time[i] = timer_kernel_2->get_time_interval();
		//printf("kernel_2_time[%d]: %f\n", i, kernel_2_time[i]);
		/*for(int j = 0; j < 3; j++) {
			total_time_group[j] += time_group[j];
		}*/
		
		if(!validate(search_key[i])) {
			printf("ERROR! round: %d, search key: %d\n", i, search_key[i]);
			break;
		}
		
		kernel_2_nedge[i] = 0;
		for(int j = 0; j < N; j++) {
			if(node_parent[j] != -1) {
				kernel_2_nedge[i] += in_deg[j];
			}
		}
		kernel_2_nedge[i] /= 2;
	}
	
	output(SCALE, edgefactor, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge/*, total_time_group*/);
	return 0;
}