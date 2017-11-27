#include "Search.h"
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#define N 2
#define WB -6.1

double B[N+1]={0,-2.0,0};

double B_n_bar=-6.1;

int a[9]={0};
double p[N+1]={0};
double p_[N+1][9];

u8 dx[N+1][9]={0};
u8 dy[N+1][9]={0};

FILE* stream;	

void Round_3(){
	fprintf(stream,"dx1:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dx[1][i]);
	}
	fprintf(stream,"\ndy1:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dy[1][i]);
	}
	fprintf(stream,"\ndx2:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dx[2][i]);
	}
	fprintf(stream,"\ndy2:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dy[2][i]);
	}
	/*fprintf(stream,"\na_i:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%d ",a[i]);
	}*/
	fprintf(stream,"\np1:%f\tp2:%f\n==============\n",p[1],p[2]);
}

void ResetCharacter(int k,int l){
	for(int i=k+1;i<=8;i++){
		if(i!=l){
			dx[2][i]=0;
			dy[2][i]=0;
		}
	}
}

void AddWeight(int j){
	p[2]=0;
	for(int k=1;k<=j;k++){
		p[2]+=p_[2][k];
	}
}

void Round_2_(int j){
	for(a[j]=a[j-1]+1;a[j]<=8;a[j]++){
		
		//ÍË³öÌõ¼þ
		if(a[j]==8){
			if(j!=1 && (dx[2][a[j-1]]&0x3)!=0){
				if(a[j-1]!=7){
					break;
				}
			}
			dx[2][8]=0;
			ResetCharacter(a[j-1],8);
			if( 0==(dx[2][7]&0x3) && 0==(dx[2][1]&0x30) ){
				dy[2][8]=0;
				p_[2][j]=0;
				AddWeight(j);
				if(p[2]+p[1]>WB){
					Round_3();
				}
			}
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[j]-1][frequency];index++){
					dx[2][8]=DDT_SearchInOrderX[7][frequency][index];
					ResetCharacter(a[j-1],8);
					if( (dx[2][8]&0x30)==((dx[2][7]&0x3)<<4) && (dx[2][8]&0x3)==((dx[2][1]&0x30)>>4) ){
						dy[2][8]=DDT_SearchInOrderY[7][frequency][index];
						p_[2][j]=DDT[7][dx[2][8]][dy[2][8]];
						AddWeight(j);
						if(p[2]+p[1]>WB){
							Round_3();
						}
					}
				}
			}
		}else if(a[j]==1){
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[j]-1][frequency];index++){
					dx[2][1]=DDT_SearchInOrderX[0][frequency][index];
					ResetCharacter(0,1);
					dy[2][1]=DDT_SearchInOrderY[0][frequency][index];
					p_[2][j]=DDT[0][dx[2][1]][dy[2][1]];
					AddWeight(j);
					if(p[2]>WB){
						Round_2_(j+1);
					}
				}
			}
		}else{
			if(j!=1 && (dx[2][a[j-1]]&0x3)!=0){
				if(a[j]!=(a[j-1]+1)){
					break;
				}
			}
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[j]-1][frequency];index++){
					dx[2][a[j]]=DDT_SearchInOrderX[a[j]-1][frequency][index];
					ResetCharacter(a[j-1],a[j]);
					if( (dx[2][a[j]]&0x30) == ((dx[2][a[j]-1]&0x3)<<4) ){
						dy[2][a[j]]=DDT_SearchInOrderY[a[j]-1][frequency][index];
						p_[2][j]=DDT[a[j]-1][dx[2][a[j]]][dy[2][a[j]]];
						AddWeight(j);
						if(p[2]>WB){
							Round_2_(j+1);
						}
					}
				}
			}
		}
	}
}

void Round_2(){
	Round_2_(1);	
}


void Round_1(){
	stream = fopen( "fprintf.txt", "w" );
	for(u32 x=1;x<0xffffffff;x++){
		ExpansionTL(dx[1]+1,x);
		p[1]=0;
		for(int Si=0;Si<8;Si++){
			p[1]+=DDT_MaxOutput[Si][dx[1][Si+1]];
			dy[1][Si+1]=DDT_MaxOutput_Index[Si][dx[1][Si+1]];
		}
		if((p[1]+B[N-1])>=B_n_bar){
			Round_2();
			printf("##:%x\n",x);
		}
	}
	fclose(stream);
}
