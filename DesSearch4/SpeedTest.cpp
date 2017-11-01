#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include "Search.h"
#define TEST_NUM 0x0
void print32(const u32 x){
	printf("%x\n",x);
}
void print8t8(const u8* y){
	for(int i=0;i<8;i++){
		printf("%x",y[i]);
	}
	printf("\n");
}

int main(){
	clock_t start,end;
	u32 x=0,z=0;
	u8 y[8]={0};
	srand((unsigned long)time(NULL));
	for(int i=0;i<3;i++){
		x=x<<15|rand();
	}
	//printf("%x\n",x);
	start = clock();
	for(int i=0;i<TEST_NUM;i++){
		//Expansion(y,x);
		//print8t8(y);
		Permutation(&z,x);
		//print32(z);
	}
	end = clock();
	printf("原始time=%f\n",(double)(end-start)/CLK_TCK);
	
	start = clock();
	GenETableLookUp();
	GenPTableLookUp();
	end = clock();
	printf("查表生成time=%f\n",(double)(end-start)/CLK_TCK);

	start = clock();
	GenDiffDistributionTable();
	GenSearchInOrder();
	end = clock();
	printf("差分分布表生成time=%f\n",(double)(end-start)/CLK_TCK);

	start = clock();
	for(int i=0;i<TEST_NUM;i++){
		//ExpansionTL(y,x);
		//print8t8(y);
		PermutationTL(&z,x);
		//print32(z);
	}
	end = clock();
	printf("查表time=%f\n",(double)(end-start)/CLK_TCK);
	

	//print(0);
	Round_2_(1);
	system("pause");
	return 0;
}