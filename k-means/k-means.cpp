#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dataset.h"
#include <vector>
#include <algorithm>
#include <windows.h>
#include <omp.h>

#define k (4)

std::vector<DataSet> datas, clusters[k];
std::vector<int> centers, data_number[k];
double center_sum_distance[k];
bool flag_finish;

void read_datas() {
	double vector0, vector1;
	while(scanf("%lf%lf",&vector0,&vector1) != EOF) {
		DataSet data;
		data.vector[0] = vector0;
		data.vector[1] = vector1;
		datas.push_back(data);
	}
}

void rand_datas() {
	int N = rand() % 10000;
	printf("N: %d\n",N);
	for(int i = 0; i < N; i++) {
		DataSet data;
		data.vector[0] = rand() % (1<<16);
		data.vector[1] = rand() % (1<<16);
		datas.push_back(data);
	}
}

void output_datas() {
	for(int i = 0; i < k; i++) {
		printf("\ncenter[%d]:\n",i);
		for(int j = 0; j < clusters[i].size(); j++) {
			printf("%d: (%lf, %lf)\n",data_number[i][j],clusters[i][j].vector[0],clusters[i][j].vector[1]);
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
	double min_distance = (1<<30);
	for(int i = 0; i < k; i++) {
		double distance = b.distance(datas[centers[i]]);
		if(distance < min_distance) {
			min_distance = distance;
			center_number = i;
		}
	}
	center_sum_distance[center_number] += min_distance;
	return center_number;
}

/* 分類所有資料屬於哪個中心(計算距離) */
void classify_datas() {
	for(int i = 0; i < k; i++) {
		center_sum_distance[i] = 0;
	}
	
	for(int i = 0; i < datas.size(); i++) {
		int center_number = find_center(datas[i]);
		clusters[center_number].push_back(datas[i]);
		data_number[center_number].push_back(i);
	}
}

/* 重新計算叢集中心 */
void calculate_centers() {
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
	srand(time(NULL));
	freopen("in","r",stdin);
	freopen("out","w",stdout);
	
	//read_datas();
	rand_datas();
	find_first_k_centers();
	while(true) {
		classify_datas();
		calculate_centers();
		if(flag_finish) {
			break;
		}
	}
	output_datas();
	return 0;
}

/*
LARGE_INTEGER frequency, start, end;
double interval;
	
QueryPerformanceFrequency(&frequency);
QueryPerformanceCounter(&start);

QueryPerformanceCounter(&end);
interval = (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;
printf("%f\n",interval);
*/