#include "Types.h"
#include "DesFunc.h"
#include "DiffDistribution.h"
#include <math.h>
#include <stdio.h>

unsigned int DDT_int[8][64][16]={0};
double DDT[8][64][16]={0};
u8 DDT_SearchInOrderX[8][9][512]={0};
u8 DDT_SearchInOrderY[8][9][512]={0};
int DDT_SearchInOrderLength[8][9]={0};

u8 DDT_SearchInOrderXWithPrefix[4][8][9][512]={0};
u8 DDT_SearchInOrderYWithPrefix[4][8][9][512]={0};
int DDT_SearchInOrderWithPrefixLength[4][8][9]={0};

u8 DDT_SearchInOrderXWithBifix[4][4][8][9][512]={0};
u8 DDT_SearchInOrderYWithBifix[4][4][8][9][512]={0};
int DDT_SearchInOrderWithBifixLength[4][4][8][9]={0};

double DDT_MaxOutput[8][64];
u8 DDT_MaxOutput_Index[8][64];

u8 DDT_SearchInOrderWithFixedX[8][9][64][16]={0};
int DDT_SearchInOrderWithFixedXLength[8][9][64]={0};

void GenDiffDistributionTable(){
	double frequency;
	for(int Si=0;Si<8;Si++){
		for(u8 Input1=0;Input1<64;Input1++){
			for(u8 Input2=0;Input2<64;Input2++){
				u8 Output1,Output2;
				Substitution(&Output1,Input1,Si);
				Substitution(&Output2,Input2,Si);
				DDT_int[Si][Input1^Input2][Output1^Output2]++;
			}
		}

		for(int i=0;i<64;i++){
			for(int j=0;j<16;j++){
				frequency=DDT_int[Si][i][j];
				if(frequency!=0){
					DDT[Si][i][j]=log(frequency)/log(2.0)-6.0;
				}
			}
		}
	}
}

void GenSearchInOrderWithFixedX(){
	int frequency,index;
	for(int Si=0;Si<8;Si++){
		for(u8 x=0;x<64;x++){
			for(u8 y=0;y<16;y++){
				frequency=DDT_int[Si][x][y]/2;
				if(frequency!=0&&frequency!=32){
					index=DDT_SearchInOrderWithFixedXLength[Si][frequency][x];
					DDT_SearchInOrderWithFixedX[Si][frequency][x][index]=y;
					DDT_SearchInOrderWithFixedXLength[Si][frequency][x]++;
				}
			}
		}
	}
}

void GenSearchInOrder(){
	int frequency,index,prefix,suffix;
	for(int Si=0;Si<8;Si++){
		for(u8 x=0;x<64;x++){
			prefix=x>>4;
			suffix=x&0x3;
			for(u8 y=0;y<16;y++){
				frequency=DDT_int[Si][x][y]/2;
				if(frequency!=0&&frequency!=32){
					index=DDT_SearchInOrderLength[Si][frequency];
					DDT_SearchInOrderX[Si][frequency][index]=x;
					DDT_SearchInOrderY[Si][frequency][index]=y;
					DDT_SearchInOrderLength[Si][frequency]++;
					
					index=DDT_SearchInOrderWithPrefixLength[prefix][Si][frequency];
					DDT_SearchInOrderXWithPrefix[prefix][Si][frequency][index]=x;
					DDT_SearchInOrderYWithPrefix[prefix][Si][frequency][index]=y;
					DDT_SearchInOrderWithPrefixLength[prefix][Si][frequency]++;

					index=DDT_SearchInOrderWithBifixLength[prefix][suffix][Si][frequency];
					DDT_SearchInOrderXWithBifix[prefix][suffix][Si][frequency][index]=x;
					DDT_SearchInOrderYWithBifix[prefix][suffix][Si][frequency][index]=y;
					DDT_SearchInOrderWithBifixLength[prefix][suffix][Si][frequency]++;
				}
			}
		}
	}
}

void GenDiffDistributionTableMax(){
	int frequency;
	u8 index;
	for(int Si=0;Si<8;Si++){
		DDT_MaxOutput[Si][0]=0;
		DDT_MaxOutput_Index[Si][0]=0;
		for(u8 x=1;x<64;x++){
			frequency=0;
			for(u8 y=0;y<16;y++){
				if(DDT_int[Si][x][y]>frequency){
					frequency=DDT_int[Si][x][y];
					index=y;
				}
			}
			DDT_MaxOutput[Si][x]=log((double)frequency)/log(2.0)-6.0;
			DDT_MaxOutput_Index[Si][x]=index;
		}
	}
}

void printMaxOutput(){
	printf("===================\nMaxOutput:\n");
	for(int i=0;i<8;i++){
		printf("Si:%d\n",i);
		for(int j=0;j<64;j++){
			printf("%f %x\t",DDT_MaxOutput[i][j],DDT_MaxOutput_Index[i][j]);
		}
		printf("\n");
	}
	printf("===================\n");
}

void printDDT(int Si){
	printf("===================\nDDT:\n");
	for(int i=0;i<64;i++){
		for(int j=0;j<16;j++){
			printf("%f ",DDT[Si][i][j]);
		}
		printf("\n");
	}
	printf("===================\n");
}


void print(int SboxIndex,u8 inputMask){
	for(int freq=8;freq>0;freq--){
		printf("freq:%d\n",freq);
		for(int index=0;index<DDT_SearchInOrderWithFixedXLength[SboxIndex][freq][inputMask];index++){
			printf("%x\t",DDT_SearchInOrderWithFixedX[SboxIndex][freq][inputMask][index]);
		}
		printf("\n");
	}
	printf("\n");

	for(int i=0;i<16;i++){
		printf("%d\t",DDT_int[SboxIndex][inputMask][i]);
	}
	printf("\n");
	printf("===================\n");
}