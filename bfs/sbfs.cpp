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

Queue_Base* current_queue;
Queue_Base* next_queue;
std::vector<bool> visited;
std::vector<std::vector<int> > adjacency_list;
std::vector<int> node_level;

void initialize(int N) {
	visited = std::vector<bool> (N);
	adjacency_list = std::vector<std::vector<int> > (N);
	node_level = std::vector<int> (N, -1);
	
	current_queue = new Circular_Queue;
	current_queue->clear();
	current_queue->reserve(N + 1);
	next_queue = new Circular_Queue;
	next_queue->clear();
	next_queue->reserve(N + 1);	
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

/* find next_queue elements */
void find_next_queue(int level) {
	while(!current_queue->empty()) {
		/* current_queue top & pop */
		int node = current_queue->top();
		current_queue->pop();
		
		for(int i = 0; i < adjacency_list[node].size(); i++) {
			if(!visited[adjacency_list[node][i]]) {
			//if(node_level[adjacency_list[node][i]] == -1) {
				visited[adjacency_list[node][i]] = true;
				node_level[adjacency_list[node][i]] = level;
				next_queue->push(adjacency_list[node][i]);
			}
		}
	}
}

/* next_queue swap to current_queue and clear */
void swap_to_current_queue() {
	Queue_Base *tmp = current_queue;
	current_queue = next_queue;
	next_queue = tmp;
	next_queue->clear();
}

void sbfs(int start_point) {
	int level = 0;
	visited[start_point] = true;
	current_queue->push(start_point);
	node_level[start_point] = level++;
	
	while(!current_queue->empty()) {
		find_next_queue(level);
		swap_to_current_queue();
		level++;
	}
}

int main() {
	freopen("out_s","w",stdout);
	
	//read_edges();
	//read_lists();
	read_graph();
	
	Windows_Time* test_time = new Windows_Time;
	
	test_time->set_start_time();
	sbfs(0);	
	test_time->set_end_time();
	
	printf("time: %f\n", test_time->get_time_interval());
	
	output_datas();
	return 0;
}