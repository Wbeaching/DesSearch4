#include "Types.h"
#include "Tables.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*-----------------���������----------------------------*/
//����32���������
u32 genRandom32(){
	u32 x;
	srand((unsigned long)time(NULL));
	for(int i=0;i<3;i++){
		x=x<<15|rand();
	}
	return x;
}
//����64���������
u64 genRandom64(){
	u64 x;
	srand((unsigned long)time(NULL));
	for(int i=0;i<5;i++){
		x=x<<15|rand();
	}
	return x;
}

/*-----------------����ת������----------------------------*/
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

//8��6��������ת����64���ش�
void SboxInput2word(u64* output, u8* input){
	u64 temp=0;
	for(int i =0;i<7;i++){
		temp+=*(input+i);
		temp<<=6;
	}
	temp+=*(input+7);
	*output=temp;
}

//8��4��������ת����32���ش�
void SboxOutput2word(u32* output, u8* input){
	u32 temp=0;
	for(int i=0;i<7;i++){
		temp+=*(input+i);
		temp<<=4;
	}
	temp+=*(input+7);
	*output=temp;
}

//64���ش�ת����8��6��������
void word642SboxInput(u8* output, u64 input){
	for(int i=0;i<8;i++){
		*(output+i)=(input>>(42-6*i))&0x3f;
	}
}

//32���ش�ת����8��4��������
void word322SboxOutput(u8* output, u32 input){
	for(int i=0;i<8;i++){
		*(output+i)=(input>>(28-4*i))&0x3f;
	}
}


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
/*
//�û�E
void Expansion(u64* output, u32 input){
	bool op[48], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 48; i++){
		op[i] = ip[ETable[i]];
	}
	bool2word48(output, op);
}*/

/*-----------------��չ�û�----------------------------*/
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
		pick6(output + i, input, ((4 * i - 1) % 32+32)%32);
	}
}

void ExpansionUsingTable(u64* output, u32 input){
	bool op[48], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 48; i++){
		op[i] = ip[ETable[i]];
	}
	bool2word48(output, op);
}

void ExpansionConvUsingShift(u32* output, u64 input){
	u64 x=input;
	u64 y=0;
	y=(x<<31)&0x8000000;
	u64 mask=0x78000000;
	for(int i=15;i>2;i-=2){
		y^=((x>>i)&mask);
		mask>>=4;
	}
	y^=((x>>1)&0x00000007);
	*output=y;
}

//��չ�û�E����֮һ
void ExpansionConv1(u32* output, u64 input){
	bool op[32], ip[48];
	word2bool48(ip, input);
	for (int i = 0; i < 32; i++){
		op[i] = ip[unETable1[i]];
	}
	bool2word(output, op);
}

//��չ�û�E����֮��
void ExpansionConv2(u32* output, u64 input){
	bool op[32], ip[48];
	word2bool48(ip, input);
	for (int i = 0; i < 32; i++){
		op[i] = ip[unETable2[i]];
	}
	bool2word(output, op);
}

bool ExpansionConvExist(u64 input){
	u32 y,z;
	ExpansionConv1(&y,input);
	ExpansionConv2(&z,input);
	return y==z;
}


/*-----------------�û�P----------------------------*/
//�û�P
void Permutation(u32* output, u32 input){
	bool op[32], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 32; i++){
		op[i] = ip[PTable[i]];
	}
	bool2word(output, op);
}

void PermutationConv(u32* output, u32 input){
	bool op[32], ip[32];
	word2bool(ip, input);
	for (int i = 0; i < 32; i++){
		op[i] = ip[unPTable[i]];
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

void SubstitutionDESL(u8* output,u8 input){
	u8 x=input,x1,x2;
	x1=(x&0x1)|((x>>4)&0x2);
	x2=(x>>1)&0xf;
	*output=S_DESL[x1][x2];
}

void print32(const u32 x){
	printf("%x\n",x);
}
void print64(const u64 x){
	printf("%llx\n",x);
}
void print8t8(const u8* y){
	for(int i=0;i<8;i++){
		printf("%x\t",y[i]);
	}
	printf("\n");
}
void fprint8t8(const u8* y,FILE* stream){
	for(int i=0;i<8;i++){
		fprintf(stream,"%x\t",y[i]);
	}
	fprintf(stream,"\n");
}
void fprintnum(const u8* y,FILE* stream){
	int count=0;
	for(int i=0;i<8;i++){
		if(y[i]!=0) count++;
	}
	fprintf(stream,"#active Sboxes:%d\t",count);
}

void fprintTab(int j,FILE* stream){
	for(int i=0;i<j;i++){
		fprintf(stream,"\t");
	}
}