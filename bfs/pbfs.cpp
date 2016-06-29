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

std::vector<int> merge_queue;
std::vector<Queue_Base*> current_queues, next_queues;
std::vector<bool> visited;
std::vector<std::vector<int> > adjacency_list;
std::vector<int> node_level;

std::vector<int> next_queue_size, current_queue_size;

int max_threads = omp_get_max_threads();
bool flag_merge;

void initialize(int N) {	
	visited = std::vector<bool> (N);
	adjacency_list = std::vector<std::vector<int> > (N);
	node_level = std::vector<int> (N, -1);
	
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

void read_edges() {
	freopen("in","r",stdin);
	
	int N, node1, node2;
	scanf("%d", &N);
	initialize(N);	
	while(scanf("%d%d", &node1, &node2) != EOF) {
		adjacency_list[node1].push_back(node2);
		adjacency_list[node2].push_back(node1);
	}
}

void read_lists() {
	freopen("in","r",stdin);
	
	int N, size, node2;
	scanf("%d", &N);
	initialize(N);
	for(int node1 = 0; node1 < N; node1++) {
		scanf("%d", &size);
		for(int i = 0; i < size; i++) {
			scanf("%d", &node2);
			adjacency_list[node1].push_back(node2);
			adjacency_list[node2].push_back(node1);
		}
	}
}

void read_graph() {
	freopen("in.graph","r",stdin);
	
	int N, M, node2;
	std::string line;	
	scanf("%d%d%c", &N, &M);
	initialize(N);
	for(int node1 = 0; node1 < N; node1++) {
		getline(std::cin, line);
		std::istringstream in(line);
		while(in >> node2) {
			adjacency_list[node1].push_back(node2 - 1);
		}
	}
}

void output_datas() {
	for(int i = 0; i < node_level.size(); i++) {
		printf("node: %d, level: %d\n", i, node_level[i]);
	}
}

/*** find next_queues elements, each thread maintain their own next queue, directly divide merge queue ***/
void find_next_queues_by_merge_queue(int level) {
	#pragma omp parallel for
	for(int i = 0; i < merge_queue.size(); i++) {
		if(node_level[merge_queue[i]] < 0) {
			node_level[merge_queue[i]] = level-1;
		}
		for(int j = 0; j < adjacency_list[merge_queue[i]].size(); j++) {
			if(!visited[adjacency_list[merge_queue[i]][j]]) {
				visited[adjacency_list[merge_queue[i]][j]] = true;
				next_queues[ omp_get_thread_num() ]->push(adjacency_list[merge_queue[i]][j]);
			}
		}
	}
}

/*** find next_queues elements, each thread maintain their own current and next queue ***/
void find_next_queues_by_current_queues(int level) {
	#pragma omp parallel for
	for(int i = 0; i < max_threads; i++) {
		int queue_size = current_queues[i]->size();
		for(int j = 0; j < queue_size; j++) {
			int node = current_queues[i]->top();
			
			if(node_level[node] < 0) {
				node_level[node] = level-1;
			}
			for(int j = 0; j < adjacency_list[node].size(); j++) {
				if(!visited[adjacency_list[node][j]]) {
					visited[adjacency_list[node][j]] = true;
					next_queues[i]->push(adjacency_list[node][j]);
				}
			}
			
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

void pbfs(int start_point) {
	int level = 1;
	visited[start_point] = true;
	merge_queue.push_back(start_point);
	next_queue_size[0] = 1;
	flag_merge = true;
	
	while(next_queue_size[level-1] > 0) {
		if(flag_merge) {
			find_next_queues_by_merge_queue(level);
		}
		else {
			find_next_queues_by_current_queues(level);
		}
		
		record_next_queues_size(level);
		
		if(flag_merge) {
			merge_to_merge_queue(level);
		}
		else {
			swap_to_current_queues(level);
		}
		
		level++;
	}
	
	/*for(int i = 1; i < level; i++) {
		printf("thread next queue size = %d, merge queue size = %d, dif = %d\n", next_queue_size[i], current_queue_size[i], next_queue_size[i] - current_queue_size[i]);
	}*/
	/*for(int i = 1; i < level; i++) {
		printf("level = %d, queue size = %d\n", i, next_queue_size[i]);
	}*/
}

int main() {
	freopen("out_p","w",stdout);
	
	//read_edges();
	//read_lists();
	read_graph();
	
	Windows_Time* test_time = new Windows_Time;
	
	test_time->set_start_time();
	pbfs(0);	
	test_time->set_end_time();
	
	printf("time: %f\n", test_time->get_time_interval());
	
	output_datas();
	return 0;
}