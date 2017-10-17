#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#ifdef RAND_MAX
#undef RAND_MAX
#define RAND_MAX 0xffffffff
#endif
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
	u32 x=0;
	u8 y[8]={0};
	srand((unsigned long)time(NULL));
	for(int i=0;i<3;i++){
		x=x<<15|rand();
		//printf("%x\n",x);
	}
	start = clock();
	for(int i=0;i<0xffffff;i++){
		Expansion(y,x);
		//print8t8(y);
		//Permutation(&x,x);
		//print32(x);
	}
	end = clock();
	printf("原始time=%f\n",(double)(end-start)/CLK_TCK);
	
	start = clock();
	GenETableLookUp();
	end = clock();
	printf("生成time=%f\n",(double)(end-start)/CLK_TCK);

	start = clock();
	for(int i=0;i<0xffffff;i++){
		ExpansionTL(y,x);
		//print8t8(y);
		//Permutation(&x,x);
		//print32(x);
	}
	end = clock();
	printf("查表time=%f\n",(double)(end-start)/CLK_TCK);
	return 0;
}