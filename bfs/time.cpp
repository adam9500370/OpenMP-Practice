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