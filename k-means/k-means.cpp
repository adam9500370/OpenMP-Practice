#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dataset.h"
#include <vector>
#include <algorithm>
#include <windows.h>

#define k (30)

std::vector<DataSet> datas;
std::vector<DataSet*> clusters[k];
std::vector<int> centers;
double center_sum_distance[k];
bool flag_finish;

void read_datas() {
	freopen("in","r",stdin);
	
	double value, vector[DIM];
	int count_dim = 0;
	for(int id = 0; scanf("%lf",&value) != EOF;) {
		count_dim++;
		vector[count_dim - 1] = value;
		if(count_dim == DIM) {
			DataSet data;
			for(int i = 0; i < DIM; i++) {
				data.vector[i] = vector[i];
			}
			data.id = id;
			datas.push_back(data);
			count_dim = 0;
			id++;
		}
	}
}

void rand_datas() {
	int N = rand() % 100000;
	//printf("N: %d\n",N);
	for(int i = 0; i < N; i++) {
		DataSet data;
		for(int j = 0; j < DIM; j++) {
			data.vector[j] = rand() % (1<<16);
		}
		data.id = i;
		datas.push_back(data);		
	}
}

void output_datas() {
	std::vector<int> id_arr;
	for(int i = 0; i < k; i++) {
		id_arr.clear();
		id_arr.reserve(clusters[i].size());
		for(int j = 0; j < clusters[i].size(); j++) {
			id_arr.push_back(clusters[i][j]->id);
		}
		sort(id_arr.begin(), id_arr.end());
		for(int j = 0; j < id_arr.size(); j++) {
			printf("dn: %d  c: %d\n",id_arr[j],i);
		}
	}
}

/* 初始隨機找k個中心 */
void find_first_k_centers() {
	bool choose[datas.size()] = {false};
	for(int i = 0; i < k; i++) {
		while(true) {
			int random_number = rand() % datas.size();
			if(!choose[random_number]) {
				centers.push_back(random_number);
				choose[random_number] = true;
				break;
			}
		}
	}
	sort(centers.begin(), centers.end());
}

/* 找中心，取距離最小者 */
int find_center(DataSet &b) {
	int center_number = 0;
	double min_distance = -1;
	for(int i = 0; i < k; i++) {
		double distance = b.distance(datas[centers[i]]);
		if(distance < min_distance || min_distance == -1) {
			min_distance = distance;
			center_number = i;
		}
	}
	return center_number;
}

/* 分類所有資料屬於哪個中心(計算距離) */
void classify_datas() {
	for(int i = 0; i < datas.size(); i++) {
		int center_number = find_center(datas[i]);
		clusters[center_number].push_back(&datas[i]);
	}
}

/* 重新計算叢集中心 */
void calculate_centers() {
	for(int i = 0; i < k; i++) {
		center_sum_distance[i] = 0;
		for(int j = 0; j < clusters[i].size(); j++) {
			center_sum_distance[i] += clusters[i][j]->distance(datas[centers[i]]);
		}
	}
	
	flag_finish = true;
	//printf("center pos: ");
	for(int i = 0; i < k; i++) {
		//printf("%d ",centers[i]);
		for(int j = 0; j < clusters[i].size(); j++) {
			double sum_distance = 0;
			for(int m = 0; m < clusters[i].size(); m++) {
				sum_distance += clusters[i][j]->distance(clusters[i][m]);
			}
			if(sum_distance < center_sum_distance[i]) { /* 確認是否與上次計算的中心位置相同 */
				centers[i] = clusters[i][j]->id;
				center_sum_distance[i] = sum_distance;
				flag_finish = false;
			}
		}		
	}
	//printf("\n");
	if(!flag_finish) {
		sort(centers.begin(), centers.end());
		for(int i = 0; i < k; i++) {
			clusters[i].clear();
		}
	}
}

int main() {
	//srand(time(NULL));
	freopen("out","w",stdout);
	
	read_datas();
	//rand_datas();
	
	LARGE_INTEGER frequency, start, end;
	double interval, interval1 = 0.0, interval2 = 0.0;
		
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	
	find_first_k_centers();
	while(true) {
		LARGE_INTEGER start1, end1, start2, end2;
		
		QueryPerformanceCounter(&start1);
		classify_datas();
		QueryPerformanceCounter(&end1);
		interval1 += (double) (end1.QuadPart - start1.QuadPart) / frequency.QuadPart;
		
		QueryPerformanceCounter(&start2);
		calculate_centers();
		QueryPerformanceCounter(&end2);
		interval2 += (double) (end2.QuadPart - start2.QuadPart) / frequency.QuadPart;
		
		if(flag_finish) {
			break;
		}
	}
	
	QueryPerformanceCounter(&end);
	interval = (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;
	
	printf("time1: %f\ntime2: %f\n",interval1, interval2);
	printf("time: %f\n",interval);
	
	output_datas();
	return 0;
}