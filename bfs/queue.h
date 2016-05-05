#ifndef __QUEUE_H__
#define __QUEUE_H__

#include <vector>
#include "queue_base.h"

class Circular_Queue : public Queue_Base {	
	std::vector<int> queue;
	int front, back;
	
public:
	void push(int value);
	void pop();
	int top();
	int size();
	bool empty();
	void clear();
	void reserve(int size);
};

#endif