#include "Types.h"
#include "Tables.h"
#include <stdio.h>

//��input�еĵ�index��boolֵ��ʼ���ȡ6��
void pick6(u8* output, u32 input, int index){
	if (index > 26){
		int carry = 32 - index;
		*output = ((input << (index - 26)) ^ (input >> (58 - index))) & 0x3f;
	}
	else{
		*output = (input >> (26 - index)) & 0x3f;
	}
}

//��չ�û�E
void Expansion(u8* output, u32 input){
	for (int i = 0; i < 8; i++){
		pick6(output + i, input, (4 * i - 1) % 32);
	}
}

//32���ش�ת���ɴ�СΪ32��bool����
void word2bool(bool* output, u32 input){
	for (int i = 0; i < 32; i++){
		*(output + i) = (input >> (31 - i)) & 0x1;
	}
}

//��СΪ32��bool����ת����32���ش�
void bool2word(u32* output, bool* input){
	u32 temp = 0;
	for (int i = 0; i <31 ; i++){
		temp += *(input + i);
		temp <<= 1;
	}
	temp += *(input + 31);
	*output = temp;
}

/*
//64���ش�ת���ɴ�СΪ48��bool����
void word2bool48(bool* output, u64 input){
	for (int i = 0; i < 48; i++){
		*(output + i) = (input >> (47 - i)) & 0x1;
	}
}

//��СΪ48��bool����ת����64���ش�
void bool2word48(u64* output, bool* input){
	u64 temp = 0;
	for (int i = 0; i <47 ; i++){
		temp += *(input + i);
		temp <<= 1;
	}
	temp += *(input + 47);
	*output = temp;
}

//�û�E
void Expansion(u64* output, u32 input){
	bool op[48], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 48; i++){
		op[i] = ip[ETable[i]];
	}
	bool2word48(output, op);
}*/


//�û�P
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

void Substitution(u8* output,u8 input,int index){
	u8 x=input,x1,x2;
	x1=(x&0x1)|((x>>4)&0x2);
	x2=(x>>1)&0xf;
	*output=S[index][x1][x2];
}
