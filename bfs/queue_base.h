#ifndef __QUEUE_BASE_H__
#define __QUEUE_BASE_H__

class Queue_Base {
public:
	virtual void push(int value) = 0;
	virtual void pop() = 0;
	virtual int top() = 0;
	virtual int size() = 0;
	virtual bool empty() = 0;
	virtual void clear() = 0;
	virtual void reserve(int size) = 0;
};

#endif