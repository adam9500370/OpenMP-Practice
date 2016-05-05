#ifndef __TIME_H__
#define __TIME_H__

#include "time_base.h"
#include <windows.h>

class Windows_Time : public Time_Base {	
	LARGE_INTEGER frequency, start, end;
	double interval;
	
public:
	void set_start_time();
	void set_end_time();
	double get_time_interval();
};

#endif