#ifndef __TIME_BASE_H__
#define __TIME_BASE_H__

class Time_Base {
public:
	virtual void set_start_time() = 0;
	virtual void set_end_time() = 0;
	virtual double get_time_interval() = 0;
};

#endif