#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#define N 8

double B[N+1];
double p[N+1]={0};
double B_n_bar;
/*
void Round_1(){
	u8 temp[8];
	for(u32 x=0;x<=0xffffffff;x++){
		ExpansionTL(temp,x);
		p[1]=0;
		for(int Si=0;Si<8;Si++){
			p[1]+=DDT_MaxOutput[Si][temp[Si]];
		}
		if((p[1]+B[N-1])>=B_n_bar){
			Round_2();
		}
	}
}

void Round_2(){
	for(u32 */