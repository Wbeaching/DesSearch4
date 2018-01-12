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
	//printf("原始time=%f\n",(double)(end-start)/CLK_TCK);
	
	start = clock();
	GenETableLookUp();
	GenPTableLookUp();
	GenEConvTableLookUp();
	GenSearchTables();
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

	/*u8 dx1[8]={0, 0, 0, 0, 0, 0, 0, 4 },dy2[8]={1, 0, 0, 0, 0, 0, 0, 0},dx3[8];

	u64 x_i_2;
	u32 x_i_2_EConv,y_i_1,x_i,y_i_1_P;
	SboxInput2word(&x_i_2, dx1);
	ExpansionConv1(&x_i_2_EConv,x_i_2);
	SboxOutput2word(&y_i_1, dy2);
	PermutationTL(&y_i_1_P,y_i_1);
	x_i=x_i_2_EConv^y_i_1_P;
	Expansion(dx3,x_i);
	print8t8(dx1);
	print8t8(dy2);
	print8t8(dx3);
	printf("%f\n",DDT[7][0xa][0x2]);
	printf("%f\n",DDT[0][0x28][0x0]);

	*/

	/*u32 x1=0x60000000,y1=0x00808200,dy1;
	u8 dx1[8];
	Expansion(dx1,x1);
	PermutationConv(&dy1,y1);
	print8t8(dx1);
	print32(dy1);
	printf("%f\n",DDT[0][0xc][0xe]);*/
	
	/*u64 x_i_2;
	u8 dx[9]={0};
	u32 x_i_2_EConv,y_i_1=0x30,x_i,y_i_1_P;
	SboxInput2word(&x_i_2, dx+1);
	ExpansionConv1(&x_i_2_EConv,x_i_2);
	
	PermutationTL(&y_i_1_P,y_i_1);
	x_i=x_i_2_EConv^y_i_1_P;
	Expansion(dx+1,x_i);
	print32(x_i);
	print8t8(dx+1);*/
	
	/*u32 dy;
	u32 dx1=genRandom32();
	u64 dx;
	ExpansionUsingTable(&dx,dx1);
	ExpansionConv1(&dy,dx);
	print32(dy);
	ExpansionConvTL(&dy,dx);
	print32(dy);*/
	//printMaxOutput();
	
	/*u32 x1=0x19600000,y1=0x0,dy;
	u8 dx1[8],dy1[8];
	Expansion(dx1,x1);
	print8t8(dx1);
	PermutationConv(&dy,y1);
	word322SboxOutput(dy1,dy);
	print8t8(dy1);
	for(int i=0;i<8;i++){
		printf("%f\t",DDT[i][dx1[i]][dy1[i]]);
	}*/
	double bound;
	rounds=5;
	bound=-13.0;
	for(int i=0;i<10;i++){
		B_n_bar=bound;
		start = clock();
		Round_1();
		end = clock();
		printf("%d轮搜索，概率下界为%f，搜索时间=%f\n",rounds,bound,(double)(end-start)/CLK_TCK);
		bound-=1.0;
	}
	system("pause");
	return 0;
}