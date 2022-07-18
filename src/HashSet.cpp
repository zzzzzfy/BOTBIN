#include "HashSet.h"


HashSet::HashSet(int size) {
	head_ = val_ = next_ = NULL;
	buckets_ = size * HASH_SCALE;
	cur_ = 0;
	capacity_ = size;

	head_ = new unsigned int[buckets_];
	val_ = new unsigned int[size];
	next_ = new unsigned int[size];

	for (unsigned int i = 0; i < buckets_; i++) head_[i] = -1;

}

HashSet::~HashSet() {
	if (head_ != NULL) delete[] head_;
	if (val_ != NULL) delete[] val_;
	if (next_ != NULL) delete[] next_;
}

void HashSet::insert(int v) {
#ifdef _DEBUG_
	if (cur_ == capacity_) {
		printf("HashSet is already full.\n");
		return;
	}
#endif

	int target_bucket = v % buckets_;
	val_[cur_] = v;
	next_[cur_] = head_[target_bucket];
	head_[target_bucket] = cur_++;
}

int HashSet::find(int v) {
	int target_bucket = v % buckets_;
	for (unsigned int i = head_[target_bucket]; i != -1; i = next_[i]) {
		if (val_[i] == v) return 1;
	}
	return 0;
}
