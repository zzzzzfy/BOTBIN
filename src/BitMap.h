#ifndef _BIT_MAP_H_
#define _BIT_MAP_H_

#include "Utility.h"

class BitMap {

	int capacity_;
	char* bm_;
	int unit_bits_;

public:
	BitMap();
	BitMap(int size);
	~BitMap();
	int set(int pos, int val);
	int get(int pos);

};


#endif