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
	printf("�������time=%f\n",(double)(end-start)/CLK_TCK);

	start = clock();
	GenDiffDistributionTable();
	GenDDT_int2DDT();
	GenDiffDistributionTableMax();
	GenSearchInOrder();
	GenSearchInOrderWithFixedX();
	end = clock();
	printf("��ֲַ�������time=%f\n",(double)(end-start)/CLK_TCK);

	/*start = clock();
	for(int i=0;i<TEST_NUM;i++){
		//ExpansionTL(y,x);
		//print8t8(y);
		PermutationTL(&z,x);
		//print32(z);
	}
	end = clock();
	printf("������time=%f\n",(double)(end-start)/CLK_TCK);*/
	
	double bound;
	for(rounds=5;rounds<6;rounds++){
		//printf("%d",rounds);
		bound=TestB[rounds];
		for(int i=0;i<4;i++){
			B_n_bar=bound;
			trailCount=0;
			characterPr=0;
			start = clock();
			Round_1();
			end = clock();
			printf("%d�������������½�Ϊ%f�����Ϊ%f������ʱ��=%f\n",rounds,bound,B_n_bar,(double)(end-start)/CLK_TCK);
			printf("������Ϊ%d,�ܸ���Ϊ%f\n",trailCount,characterPr);
			//printf("&%d&%d&%f\\\\\n",(int)(-bound),trailCount,(double)(end-start)/CLK_TCK);
			bound-=1.0;
		}
	}
	

	/*DPL=0x001a0000;
	DPR=0x00000000;
	DCL=0x001a0401;
	DCR=0x00000040;
	for(rounds=5;rounds<6;rounds++){
		bound=TestB[rounds]-4.0;
		for(int i=0;i<1;i++){
			B_n_bar=bound;
			trailCount=0;
			characterPr=0;
			start = clock();
			Round_1_Fix();
			end = clock();
			printf("%d�������������½�Ϊ%f�����Ϊ%f������ʱ��=%f\n",rounds,bound,B_n_bar,(double)(end-start)/CLK_TCK);
			printf("������Ϊ%d,�ܸ���Ϊ%f\n",trailCount,characterPr);
			bound-=1.0;
		}
	}*/
	
	system("pause");
	return 0;
}