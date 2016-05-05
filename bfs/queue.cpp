#include "queue.h"

void Circular_Queue::push(int value) {
	queue[back] = value;
	back = (back + 1) % queue.capacity();
}

void Circular_Queue::pop() {
	front = (front + 1) % queue.capacity();
}

int Circular_Queue::top() {
	if(!empty())
		return queue[front % queue.capacity()];
}

int Circular_Queue::size() {
	if(front <= back) {
		return back - front;
	}
	// else
	return queue.capacity() + (back - front);
}

bool Circular_Queue::empty() {
	return (front == back);
}

void Circular_Queue::clear() {
	front = back = 0;	
}

void Circular_Queue::reserve(int size) {
	queue.reserve(size);
}