#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include "Search.h"
#define TEST_NUM 0x0


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
	GenDiffDistributionTableMax();
	GenSearchInOrder();
	GenSearchInOrderWithFixedX();
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
	
	//printMaxOutput();
	//printDDT(0);
	//print(5,0x8);
	
	/*u32 x1;
	PermutationTL(&x1,0x60000000);
	printf("%x",x1);
	u8 x2[8];
	Expansion(x2,x1);
	print8t8(x2);*/

	Round_1();
	system("pause");
	return 0;
}