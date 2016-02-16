#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dataset.h"
#include <vector>
#include <algorithm>
#include <windows.h>
#include <omp.h>

#define k (30)

std::vector<DataSet> datas, clusters[k];
std::vector<int> centers, data_number[k];
double center_sum_distance[k];
bool flag_finish;

void read_datas() {
	freopen("in","r",stdin);
	
	double value, vector[DIM];
	int count_dim = 0;
	while(scanf("%lf",&value) != EOF) {
		count_dim++;
		vector[count_dim - 1] = value;
		if(count_dim == DIM) {
			DataSet data;
			for(int i = 0; i < DIM; i++) {
				data.vector[i] = vector[i];
			}
			datas.push_back(data);
			count_dim = 0;
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
		datas.push_back(data);
	}
}

void output_datas() {
	for(int i = 0; i < k; i++) {
		sort(data_number[i].begin(), data_number[i].end());
		//printf("\ncenter[%d]:\n",i);
		for(int j = 0; j < data_number[i].size(); j++) {
			//printf("%d: (%lf, %lf)\n",data_number[i][j],clusters[i][j].vector[0],clusters[i][j].vector[1]);
			//printf("%lf %lf %d\n",clusters[i][j].vector[0],clusters[i][j].vector[1],i);
			printf("dn: %d  c: %d\n",data_number[i][j],i);
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
	/*** omp_lock ***/
	omp_lock_t lock;
	omp_init_lock(&lock);
	#pragma omp parallel for
	for(int i = 0; i < datas.size(); i++) {
		int center_number = find_center(datas[i]);
		
		omp_set_lock(&lock);
		clusters[center_number].push_back(datas[i]);
		data_number[center_number].push_back(i);
		omp_unset_lock(&lock);
	}
	omp_destroy_lock(&lock);
}

/* 重新計算叢集中心 */
void calculate_centers() {
	for(int i = 0; i < k; i++) {
		center_sum_distance[i] = 0;
		for(int j = 0; j < clusters[i].size(); j++) {
			center_sum_distance[i] += datas[data_number[i][j]].distance(datas[centers[i]]);
		}
	}
	
	flag_finish = true;
	//printf("center pos: ");
	for(int i = 0; i < k; i++) {
		//printf("%d ",centers[i]);
		for(int j = 0; j < clusters[i].size(); j++) {
			double sum_distance = 0;
			for(int m = 0; m < clusters[i].size(); m++) {
				sum_distance += datas[data_number[i][j]].distance(datas[data_number[i][m]]);
			}
			if(sum_distance < center_sum_distance[i]) { /* 確認是否與上次計算的中心位置相同 */
				centers[i] = data_number[i][j];
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
			data_number[i].clear();
		}
	}
}

int main() {
	//srand(time(NULL));
	freopen("out1","w",stdout);
	
	read_datas();
	//rand_datas();
	
	LARGE_INTEGER frequency, start, end;
	double interval;
		
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	
	find_first_k_centers();
	while(true) {
		classify_datas();
		calculate_centers();
		if(flag_finish) {
			break;
		}
	}
	
	QueryPerformanceCounter(&end);
	interval = (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;
	printf("time: %f\n",interval);
	
	output_datas();
	return 0;
}