#ifndef _HASH_SET_H_
#define _HASH_SET_H_

#include "Utility.h"

class HashSet {
private:
	unsigned int* head_;
	unsigned int* val_;
	unsigned int* next_;
	int cur_, capacity_, buckets_;

public:
	HashSet(int size);
	~HashSet();

	void insert(int v);
	int find(int v);
};
#endif
