#ifndef __DATASET_H__
#define __DATASET_H__
class DataSet {
public:
	double vector[2];
	double distance(const DataSet &b) const {
		double sum = 0;
		for(int i = 0; i < 2; i++) {
			double d = vector[i] - b.vector[i];
			sum += d*d; // 2-norm
		}
		return sum;
	}
};
#endif