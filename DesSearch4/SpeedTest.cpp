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
	
	start = clock();
	GenETableLookUp();
	GenPTableLookUp();
	GenEConvTableLookUp();
	GenSearchTables();
	end = clock();
	printf("查表生成time=%f\n",(double)(end-start)/CLK_TCK);

	start = clock();
	GenDiffDistributionTable();
	GenDDT_int2DDT();
	GenDiffDistributionTableMax();
	GenSearchInOrder();
	GenSearchInOrderWithFixedX();
	end = clock();
	printf("差分分布表生成time=%f\n",(double)(end-start)/CLK_TCK);

	/*start = clock();
	for(int i=0;i<TEST_NUM;i++){
		//ExpansionTL(y,x);
		//print8t8(y);
		PermutationTL(&z,x);
		//print32(z);
	}
	end = clock();
	printf("查表测试time=%f\n",(double)(end-start)/CLK_TCK);*/
	
	double bound;
	for(rounds=16;rounds<17;rounds++){
		//printf("%d",rounds);
		bound=TestB[rounds]-3.0;
		for(int i=0;i<1;i++){
			B_n_bar=bound;
			trailCount=0;
			characterPr=0;
			start = clock();
			Round_1();
			end = clock();
			printf("%d轮搜索，概率下界为%f，结果为%f，搜索时间=%f\n",rounds,bound,B_n_bar,(double)(end-start)/CLK_TCK);
			printf("迹条数为%d,总概率为%f\n",trailCount,characterPr);
			//printf("&%d&%d&%f\\\\\n",(int)(-bound),trailCount,(double)(end-start)/CLK_TCK);
			bound-=1.0;
		}
	}
	
	/*dx[1][1]=0x0;
	dx[1][2]=0x0;
	dx[1][3]=0x0;
	dx[1][4]=0x0;
	dx[1][5]=0x0;
	dx[1][6]=0x0;
	dx[1][7]=0x0;
	dx[1][8]=0x0;
	for(rounds=16;rounds<17;rounds++){
		bound=TestB[rounds]-2.0;
		for(int i=0;i<1;i++){
			B_n_bar=bound;
			trailCount=0;
			characterPr=0;
			start = clock();
			Round_1_Fix();
			end = clock();
			printf("%d轮搜索，概率下界为%f，结果为%f，搜索时间=%f\n",rounds,bound,B_n_bar,(double)(end-start)/CLK_TCK);
			printf("迹条数为%d,总概率为%f\n",trailCount,characterPr);
			bound-=1.0;
		}
	}*/
	
	system("pause");
	return 0;
}