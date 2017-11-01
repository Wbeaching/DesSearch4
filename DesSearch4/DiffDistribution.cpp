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

void GenSearchInOrder(){
	int frequency,index;
	for(int Si=0;Si<8;Si++){
		for(u8 x=0;x<64;x++){
			for(u8 y=0;y<16;y++){
				frequency=DDT_int[Si][x][y]/2;
				if(frequency!=0&&frequency!=32){
					index=DDT_SearchInOrderLength[Si][frequency];
					DDT_SearchInOrderX[Si][frequency][index]=x;
					DDT_SearchInOrderY[Si][frequency][index]=y;
					DDT_SearchInOrderLength[Si][frequency]++;
				}
			}
		}
	}
}

void print(int k){
	printf("===================\n");
	for(int i=8;i>0;i--){
		for(int j=0;j<DDT_SearchInOrderLength[k][i];j++){
			printf("%x %x: %f\t",DDT_SearchInOrderX[k][i][j],DDT_SearchInOrderY[k][i][j],DDT[k][DDT_SearchInOrderX[k][i][j]][DDT_SearchInOrderY[k][i][j]]);
		}
		printf("\n");
	}
	printf("===================\n");
}