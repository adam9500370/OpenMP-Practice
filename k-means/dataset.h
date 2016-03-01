#ifndef __DATASET_H__
#define __DATASET_H__
#define DIM (200)
class DataSet {
public:
	double vector[DIM];
	int id;
	double distance(const DataSet &b) const {
		double sum = 0;
		for(int i = 0; i < DIM; i++) {
			double d = vector[i] - b.vector[i];
			sum += d*d; // 2-norm
		}
		return sum;
	}
	double distance(const DataSet *b) const {
		return distance(*b);
	}
};
#endif