#pragma once
#include "Types.h"
#include "DesFunc.h"
#include "DiffDistribution.h"
#include <math.h>
#include <stdio.h>

unsigned int DDT_int[8][64][16]={0};
double DDT[8][64][16]={0};
double DDT_MaxOutput[8][64]={0};
unsigned int DDT_MaxOutput_Index[8][64]={0};
u8 DDT_SearchInOrderX[8][9][512]={0};
u8 DDT_SearchInOrderY[8][9][512]={0};
int DDT_SearchInOrderLength[8][9]={0};

void GenDiffDistributionTable(){
	double frequency;
	for(int Si=0;Si<8;Si++){
		//printf("S=%d\n",Si);
		for(u8 Input1=0;Input1<64;Input1++){
			for(u8 Input2=0;Input2<64;Input2++){
				u8 Output1,Output2;
				Substitution(&Output1,Input1,Si);
				Substitution(&Output2,Input2,Si);
				DDT_int[Si][Input1^Input2][Output1^Output2]++;
			}
		}
		/*for(int i=0;i<64;i++){
			printf("deltx=%d ",i);
			for(int j=0;j<16;j++){
				printf("%d ",DDT_int[Si][i][j]);
			}
			printf("\n");
		}*/

		for(int i=0;i<64;i++){
			for(int j=0;j<16;j++){
				frequency=DDT_int[Si][i][j];
				if(frequency==0){
					DDT[Si][i][j]=0;
				}else{
					DDT[Si][i][j]=log(frequency)/log(2.0);
				}
			}
		}

	}
}

void GenSearchInOrder(){
	int count,index;
	for(int Si=0;Si<8;Si++){
		printf("S=%d\n",Si);
		for(u8 x=1;x<64;x++){
			for(u8 y=1;y<16;y++){
				count=DDT_int[Si][x][y]/2;
				index=DDT_SearchInOrderLength[Si][count];
				DDT_SearchInOrderX[Si][count][index]=x;
				DDT_SearchInOrderY[Si][count][index]=y;
				DDT_SearchInOrderLength[Si][count]++;
			}
		}
		/*for(int count=8;count>0;count--){
			printf("count=%d  ",count);
			for(int index=0;index<50;index++){
				printf("%d ",DDT_SearchInOrderX[Si][count][index]);
			}
			printf("\n");
		}*/
	}
	

}

void GenDiffDistributionTableMax(){
	double temp=0;
	u8 temp_y;
	for(int Si=0;Si<8;Si++){
		for(int x=0;x<64;x++){
			for(int y=0;y<16;y++){
				if(DDT[Si][x][y]>temp){
					temp=DDT[Si][x][y];
					temp_y=y;
				}
			}
			DDT_MaxOutput[Si][x]=temp;
			DDT_MaxOutput_Index[Si][x]=temp_y;
		}
	}
}