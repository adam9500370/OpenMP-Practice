#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include "queue.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "time.h"
#include <omp.h>
#include <math.h>

// Seed with a real random value, if available
std::random_device r;

int N, M;

std::vector<double> ii_arr, jj_arr, norm_arr;
std::vector<char> ii_bit, jj_bit;
std::vector<int> p_v, p_e;

struct Edge {
	int x, y;
};
std::vector<Edge> edge_list, tmp_edge_list; // x for start vertex, y for end vertex
std::vector<int> in_deg;
std::vector<int> rand_search_key, search_key;
std::vector<int> merge_queue;
std::vector<Queue_Base*> next_queues;
std::vector<char> visited;
std::vector<std::vector<int> > adjacency_list;
std::vector<int> node_parent;

std::vector<int> next_queue_size, current_queue_size;

int max_threads/* = omp_get_max_threads()*/;

void initialize() {	
	omp_set_num_threads(max_threads);
	
	in_deg = std::vector<int> (N);
	adjacency_list = std::vector<std::vector<int> > (N);
	
	next_queues = std::vector<Queue_Base*> (max_threads);
	for(int i = 0; i < next_queues.size(); i++) {
		next_queues[i] = new Circular_Queue;
		next_queues[i]->clear();
		next_queues[i]->reserve(N + 1);
	}
	
	next_queue_size = std::vector<int> (N);
}

void randperm(int n, std::vector<int> &perm) {
	int j, tmp;
	for(int i = 0; i < n; i++) {
		perm[i] = i;
	}
	
	std::default_random_engine generator(r());
	for(int i = 0; i < n; i++) {	
		std::uniform_int_distribution<int> distribution(0, n-i-1);
		j = distribution(generator) + i;
		// swap i, j
		tmp = perm[j];
		perm[j] = perm[i];
		perm[i] = tmp;
	}
}

void rand_arr(int len, std::vector<double> &arr) {
	std::default_random_engine generator(r());
	std::uniform_int_distribution<int> distribution(0, 10000);
	#pragma omp parallel for
	for(int i = 0; i < len; i++) {		
		arr[i] = (double)distribution(generator) / 10000;
	}
}

void cmp_gt(int len, const std::vector<double> &arr1, const std::vector<double> &arr2, std::vector<char> &result_arr) { // compare whether arr1 > arr2 or not
	#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		result_arr[i] = (arr1[i] > arr2[i]);
	}
}

void cmp_gt(int len, const std::vector<double> &arr1, double value, std::vector<char> &result_arr) { // compare whether arr1 > value or not
	#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		result_arr[i] = (arr1[i] > value);
	}
}

void calculate_norm_arr(double c_norm, double a_norm) {
	#pragma omp parallel for
	for(int k = 0; k < M; k++) {
		norm_arr[k] = c_norm * ii_bit[k] + a_norm * (!ii_bit[k]);
	}
}

void sum_tmp_edge_list(int ib) {
	#pragma omp parallel for
	for(int k = 0; k < M; k++) {
		tmp_edge_list[k].x += pow(2, (ib-1)) * ii_bit[k];
		tmp_edge_list[k].y += pow(2, (ib-1)) * jj_bit[k];
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
	tmp_edge_list = std::vector<Edge> (M);
	
	// Loop over each order of bit.
	double ab, c_norm, a_norm;
	ab = A + B;
	c_norm = C/(1 - (A + B));
	a_norm = A/(A + B);
	
	//printf("N: %d, M: %d, ab: %f, c_norm: %f, a_norm: %f\n", N, M, ab, c_norm, a_norm);
	
	ii_arr.reserve(M);
	jj_arr.reserve(M);
	norm_arr.reserve(M);
	ii_bit.reserve(M);
	jj_bit.reserve(M);
	for(int ib = 1; ib <= SCALE; ib++) {
		// Compare with probabilities and set bits of indices.
		rand_arr(M, ii_arr);
		cmp_gt(M, ii_arr, ab, ii_bit);
		rand_arr(M, jj_arr);
		calculate_norm_arr(c_norm, a_norm);
		cmp_gt(M, jj_arr, norm_arr, jj_bit);
		sum_tmp_edge_list(ib);
	}
	std::vector<double> ().swap(ii_arr);
	std::vector<double> ().swap(jj_arr);
	std::vector<double> ().swap(norm_arr);
	std::vector<char> ().swap(ii_bit);
	std::vector<char> ().swap(jj_bit);
	
	// Permute vertex labels
	p_v.reserve(N);
	randperm(N, p_v);
	#pragma omp parallel for
	for(int k = 0; k < M; k++) {
		tmp_edge_list[k].x = p_v[tmp_edge_list[k].x];
		tmp_edge_list[k].y = p_v[tmp_edge_list[k].y];
	}
	std::vector<int> ().swap(p_v);
	
	// Permute the edge list
	edge_list = std::vector<Edge> (M);
	p_e.reserve(M);
	randperm(M, p_e);
	#pragma omp parallel for
	for(int k = 0; k < M; k++) {
		edge_list[k].x = tmp_edge_list[p_e[k]].x;
		edge_list[k].y = tmp_edge_list[p_e[k]].y;
	}
	std::vector<int> ().swap(p_e);
	std::vector<Edge> ().swap(tmp_edge_list);
}

bool cmp_edge(const Edge &a, const Edge &b) {
    return (a.x < b.x) || (a.x == b.x && a.y < b.y);
}

bool test_self_edge(int index) {
	return (edge_list[index].x == edge_list[index].y);
}

bool test_duplicate_edge(int index) {
	return ((index > 0) && (edge_list[index].x == edge_list[index-1].x) && (edge_list[index].y == edge_list[index-1].y));
}

/*** transform edge list into adjacency list ***/
void construct_graph(bool flag_bidirection) {
	sort(edge_list.begin(), edge_list.end(), cmp_edge);
	for(int i = 0; i < edge_list.size(); i++) {
		if(!test_self_edge(i)) { // no self edge, no duplicate edge, undirected graph
			in_deg[edge_list[i].y]++;				
			if(flag_bidirection) {
				in_deg[edge_list[i].x]++;
			}
			
			if(!test_duplicate_edge(i)) {
				adjacency_list[edge_list[i].x].push_back(edge_list[i].y);
				if(flag_bidirection) {
					adjacency_list[edge_list[i].y].push_back(edge_list[i].x);
				}
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

/*** record each next queue size ***/
void record_next_queues_size(int level) {
	next_queue_size[level] = 0;
	//printf("level %d next_queues size: ", level);
	for(int i = 0; i < max_threads; i++) {
		next_queue_size[level] += next_queues[i]->size();
		//printf("%d ",next_queues[i]->size());
	}
	//printf(", sum = %d\n", next_queue_size[level]);
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

void pbfs(int start_point/*, double* time_group, int group_num*/) {
	int level = 1;
	/*
	std::vector<Windows_Time*> timer.reserve(group_num);
	
	for(int i = 0; i < group_num; i++) {
		timer[i] = new Windows_Time;
		time_group[i] = 0.0;
	}
	*/
	//timer[0]->set_start_time();
	node_parent.clear();
	node_parent.reserve(N);
	
	visited.clear();
	visited.reserve(N);
	
	merge_queue.clear();
	merge_queue.push_back(start_point);
	
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
	}
	
	node_parent[start_point] = start_point;
	next_queue_size[0] = 1;
	/*timer[0]->set_end_time();
	time_group[0] += timer[0]->get_time_interval();*/
	
	while(next_queue_size[level-1] > 0) {
		//timer[1]->set_start_time();
		find_next_queues_by_merge_queue(level);
		/*timer[1]->set_end_time();
		time_group[1] += timer[1]->get_time_interval();*/
		
		//timer[2]->set_start_time();
		record_next_queues_size(level);
		/*timer[2]->set_end_time();
		time_group[2] += timer[2]->get_time_interval();*/
		
		//timer[3]->set_start_time();
		merge_to_merge_queue(level);
		/*timer[3]->set_end_time();
		time_group[3] += timer[3]->get_time_interval();*/
		
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
	std::vector<int> node_level, current_level_list, next_level_list;
	node_level.reserve(N);
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
	std::vector<int> ().swap(next_level_list);
	
	for(int i = 0; i < edge_list.size(); i++) {
		if( ! (	(node_level[edge_list[i].x] == -1 && node_level[edge_list[i].y] == -1) // both out of graph
			 || (node_level[edge_list[i].x] != -1 && node_level[edge_list[i].y] != -1) // both in graph
			 || (abs(node_level[edge_list[i].x] - node_level[edge_list[i].y]) <= 1) // each tree edge connects vertices whose BFS levels differ by exactly one
				)){
			printf("err: -3\n");
			printf("e: %d %d, l: %d %d\n", edge_list[i].x, edge_list[i].y, node_level[edge_list[i].x], node_level[edge_list[i].y]);
			return false;
		}
	}
	std::vector<int> ().swap(node_level);
	
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

void output(int SCALE, int edgefactor, int NBFS, double kernel_1_time, double* kernel_2_time, double* kernel_2_nedge/*, double* total_time_group, int group_num*/) {
	char strch_max_threads[10], strch_SCALE[10], strch_edgefactor[10];
	sprintf(strch_max_threads, "%d", max_threads);
	sprintf(strch_SCALE, "%d", SCALE);
	sprintf(strch_edgefactor, "%d", edgefactor);
	std::string str_max_threads(strch_max_threads), str_SCALE(strch_SCALE), str_edgefactor(strch_edgefactor);
	std::string file_name = "Test_OpenMP_" + str_max_threads + "threads_Kronecker/" + "SCALE=" + str_SCALE + "_edgefactor=" + str_edgefactor; // store in Test_OpenMP_Xthreads_Kronecker dir
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
	/*
	printf("\n");
	for(int i = 0; i < group_num; i++) {
		printf("average_time_group[%d]: %f (%f%%)\n", i, total_time_group[i]/NBFS, (total_time_group[i]/NBFS)/kernel_2_mean_time*100);
	}
	*/
	/*
	printf("\n");
	for(int i = 0; i < N; i++) {
		if(in_deg[i])
			printf("in_deg[%6d] = %6d\n", i, in_deg[i]);
	}
	*/
	/*
	printf("\n");
	for(int i = 0; i < N; i++) {
		if(adjacency_list[i].size())
			printf("adjacency_list[%6d].size = %6d\n", i, adjacency_list[i].size());
	}
	*/
}

int main(int argc, char *argv[]) {
	int SCALE, edgefactor;
	int NBFS = 64;
	Windows_Time *timer_kernel_1 = new Windows_Time, *timer_kernel_2 = new Windows_Time;
	double kernel_1_time, kernel_2_time[NBFS] = {0.0}, kernel_2_nedge[NBFS] = {0.0};
	
	if(argc == 4) { // argv[0] = pbfs
		max_threads = atoi(argv[1]);
		SCALE = atoi(argv[2]);
		edgefactor = atoi(argv[3]);
	}
	else {
		exit(EXIT_FAILURE);
	}
	
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
	
	rand_search_key.reserve(N);
	randperm(N, rand_search_key);
	for(int i = 0; i < N; i++) {
		if(adjacency_list[rand_search_key[i]].size() > 0) {
			search_key.push_back(rand_search_key[i]);
		}
	}
	std::vector<int> ().swap(rand_search_key);
	if(search_key.size() > NBFS) {
		search_key.erase(search_key.end()-(N-NBFS), search_key.end());
	}
	else {
		NBFS = search_key.size();
	}
	/*
	int group_num = 4;
	double time_group[group_num], total_time_group[group_num] = {0.0};
	*/
	for(int i = 0; i < NBFS; i++) {
		//printf("search key: %d\n", search_key[i]);		
		timer_kernel_2->set_start_time();
		pbfs(search_key[i]/*, time_group, group_num*/);	
		timer_kernel_2->set_end_time();
		kernel_2_time[i] = timer_kernel_2->get_time_interval();
		//printf("kernel_2_time[%d]: %f\n", i, kernel_2_time[i]);
		/*for(int j = 0; j < group_num; j++) {
			total_time_group[j] += time_group[j];
		}*/
		
		if(!validate(search_key[i])) {
			printf("ERROR! round: %d, search key: %d\n", i, search_key[i]);
			break;
		}
		
		kernel_2_nedge[i] = 0.0;
		for(int j = 0; j < N; j++) {
			if(node_parent[j] != -1) {
				kernel_2_nedge[i] += in_deg[j];
			}
		}
		kernel_2_nedge[i] /= 2;
	}
	
	output(SCALE, edgefactor, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge/*, total_time_group, group_num*/);
	return 0;
}