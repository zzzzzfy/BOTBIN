#include "BitMap.h"

BitMap::BitMap() {
	capacity_ = 0;
	unit_bits_ = 0;
	bm_ = NULL;
}

BitMap::BitMap(int size) {

	capacity_ = size;
	unit_bits_ = sizeof(char);
	int windows_size = (size - 1) / unit_bits_ + 1;
	bm_ = new char[windows_size];
	memset(bm_, 0, unit_bits_ * windows_size);
}

BitMap::~BitMap() {
	if (bm_ != NULL) delete[] bm_;
}

int BitMap::set(int pos, int val) {
	if (pos >= capacity_) {
		printf("Position overflow\n");
		return -1;
	}
	int addr = pos / unit_bits_;
	int addr_offset = pos % unit_bits_;
	int reverse_move = unit_bits_ - addr_offset - 1;


	unsigned char tmp = 0x1 << reverse_move;

	if (val) bm_[addr] |= tmp;
	else {

		bm_[addr] &= ~tmp;
	}

	return 1;
}

int BitMap::get(int pos) {
	if (pos >= capacity_) {
		printf("Position overflow\n");
		return -1;
	}
	int addr = pos / unit_bits_;
	int addr_offset = pos % unit_bits_;
	int reverse_move = unit_bits_ - addr_offset - 1;

	unsigned char tmp = 0x1 << reverse_move;

	return (bm_[addr] & tmp) >> reverse_move;

}
