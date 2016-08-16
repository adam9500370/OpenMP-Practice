#include "time.h"

void Windows_Time::set_start_time() {
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
}

void Windows_Time::set_end_time() {
	QueryPerformanceCounter(&end);
}

double Windows_Time::get_time_interval() {
	return (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;
}


void Linux_Time::set_start_time() {
	gettimeofday(&start, NULL);
}

void Linux_Time::set_end_time() {
	gettimeofday(&end, NULL);
}

double Linux_Time::get_time_interval() {
	return (double) ( (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0 );
}