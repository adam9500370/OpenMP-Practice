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

Queue_Base* current_queue;
std::vector<Queue_Base*> next_queues;
std::vector<bool> thread_visited, visited;
std::vector<std::vector<int> > adjacency_list;
std::vector<int> node_level;

std::vector<int> next_queue_size, current_queue_size;

int max_threads = omp_get_max_threads();

void initialize(int N) {	
	visited = std::vector<bool> (N);
	thread_visited = std::vector<bool> (N);
	adjacency_list = std::vector<std::vector<int> > (N);
	node_level = std::vector<int> (N, -1);
	
	current_queue = new Circular_Queue;
	current_queue->clear();
	current_queue->reserve(N + 1);
	/*for(int i = 0; i < current_queue.size(); i++) {
		current_queue[i] = new Circular_Queue;
		current_queue[i]->clear();
		current_queue[i]->reserve(N + 1);		
	}*/
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

/*** first current_queues top to array ***/
void create_current_queue_array(std::vector<int> &node_array, int level) {
	node_array.clear();
	node_array.reserve(next_queue_size[level-1]);
	//for(int i = 0; i < max_threads; i++) {
		while(!current_queue->empty()) {
			/* current_queue top & pop */
			node_array.push_back(current_queue->top());
			current_queue->pop();
		}
	//}
}

/* find next_queue elements */
/*** each thread maintain their own queue ***/
void find_next_queue(const std::vector<int> &node_array, int level) {
	#pragma omp parallel for
	for(int i = 0; i < node_array.size(); i++) {
		for(int j = 0; j < adjacency_list[node_array[i]].size(); j++) {
			if(!visited[adjacency_list[node_array[i]][j]]) {
			//if(node_level[adjacency_list[node_array[i]][j]] == -1) {
				visited[adjacency_list[node_array[i]][j]] = true;
				//node_level[adjacency_list[node_array[i]][j]] = level;
				int tid = omp_get_thread_num();
				next_queues[tid]->push(adjacency_list[node_array[i]][j]);
			}
		}
	}
}

/*** next_queues swap to current_queues and clear ***/
/*void swap_to_current_queue(int level) {
	next_queue_size[level] = 0;
	for(int i = 0; i < max_threads; i++) {
		next_queue_size[level] += next_queues[i]->size();
		
		Queue_Base *tmp = current_queue[i];
		current_queue[i] = next_queues[i];
		next_queues[i] = tmp;
		next_queues[i]->clear();
	}
}*/

/*** next_queues merge to current_queue and clear ***/
void merge_to_current_queue(int level) {
	next_queue_size[level] = 0;
	for(int i = 0; i < max_threads; i++) {
		next_queue_size[level] += next_queues[i]->size();
		while(!next_queues[i]->empty()) {
			int node = next_queues[i]->top();
			if(!thread_visited[node]) {
			//if(node_level[node] == -1) {
				thread_visited[node] = true;
				node_level[node] = level;
				current_queue->push(node);
			}
			next_queues[i]->pop();
		}
		next_queues[i]->clear();
	}
	current_queue_size[level] = current_queue->size();
}

/*bool is_all_current_queue_empty() {
	for(int i = 0; i < max_threads; i++) {
		if(!current_queue[i]->empty())
			return false;
	}
	return true;
}*/

void pbfs(int start_point) {
	int level = 0;
	visited[start_point] = true;
	thread_visited[start_point] = true;
	current_queue->push(start_point);
	next_queue_size[0] = 1;
	node_level[start_point] = level++;
	
	std::vector<int> node_array;
	while(!current_queue->empty()) {
		create_current_queue_array(node_array, level);
		find_next_queue(node_array, level);
		merge_to_current_queue(level);
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