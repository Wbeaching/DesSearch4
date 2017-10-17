#include "Types.h"
#include "Tables.h"
#include <stdio.h>

//从input中的第index个bool值开始向后取6个
void pick6(u8* output, u32 input, int index){
	if (index > 26){
		int carry = 32 - index;
		*output = ((input << (index - 26)) ^ (input >> (58 - index))) & 0x3f;
	}
	else{
		*output = (input >> (26 - index)) & 0x3f;
	}
}

//扩展置换E
void Expansion(u8* output, u32 input){
	for (int i = 0; i < 8; i++){
		pick6(output + i, input, (4 * i - 1) % 32);
	}
}

//32比特串转换成大小为32的bool数组
void word2bool(bool* output, u32 input){
	for (int i = 0; i < 32; i++){
		*(output + i) = (input >> (31 - i)) & 0x1;
	}
}

//大小为32的bool数组转换成32比特串
void bool2word(u32* output, bool* input){
	u32 temp = 0;
	for (int i = 0; i <31 ; i++){
		temp += *(input + i);
		temp <<= 1;
	}
	temp += *(input + 31);
	*output = temp;
}

//置换P
void Permutation(u32* output, u32 input){
	bool op[32], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 32; i++){
		op[i] = ip[PTable[i]];
	}
	bool2word(output, op);
}

unsigned long long Expansion1(unsigned int x){
	unsigned long long y, z = x, t = 0xfc0;
	y = z >> 31 | z << 1 & 0x3e;
	for (int i = 0, j = 3; i<36; i = i + 6, j = j + 2){
		y |= z << j&t << i;
	}
	y |= z << 15 & 0x7c0000000000 | z << 47 & 0x800000000000;
	return y;
}

unsigned int Permutation1(unsigned int x){
	unsigned int y = 0;
	unsigned int t = 0x80000000;
	for (int i = 0; i<32; i++, t = t >> 1){
		if (PTable[i] - i>0) y |= x << PTable[i] - i&t;
		else y |= x >> i - PTable[i] & t;
	}
	return y;
}
