#include "Search.h"
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#define N 4

double B[N]={0,0,-2.0,-4.0};//,-9.60768};

static double B_n_bar=-10.0;

int a[N+1][9]={0};
double p[N+1];
double p_[N+1][9];

u8 dx[N+1][9]={0};
u8 dy[N+1][9]={0};

bool flag=0;
FILE* stream;

void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<=8;i++){
		if(i!=l){
			dx[round][i]=0;
			dy[round][i]=0;
		}
	}
}

void AddWeight(int j,int round){
	p[round]=0;
	for(int k=1;k<=j;k++){
		p[round]+=p_[round][k];
	}
}

double sumWeight(int m){
	double temp=0;
	for(int i=1;i<=m;i++){
		temp+=p[i];
	}
	return temp;
}

void printAndSetBound(){
	
	for(int r=1;r<=N;r++){
		fprintf(stream,"dx%d:",r);
		for(int i=1;i<=8;i++){
			fprintf(stream,"%x ",dx[r][i]);
		}
		fprintf(stream,"\ndy%d:",r);
		for(int i=1;i<=8;i++){
			fprintf(stream,"%x ",dy[r][i]);
		}
		fprintf(stream,"\tp%d:%f\n",r,p[r]);
	}
	B_n_bar=sumWeight(N);
	fprintf(stream,"B_n_bar:%f\n==============\n",B_n_bar);
}

void Round_(int i);

void Round__(int i,int j){
	if(dx[i][j]==0){
		dy[i][j]=0;
		p_[i][j]=0;
		if(j==8){
			Round_(i+1);
		}else{
			Round__(i,j+1);
		}
	}else{
		for(int frequency=8;frequency>0;frequency--){
			for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[i][j]];index++){
				dy[i][j]=DDT_SearchInOrderWithFixedX[j-1][frequency][dx[i][j]][index];
				p_[i][j]=DDT[j-1][dx[i][j]][dy[i][j]];
				AddWeight(j,i);

				//fprintf(stream,"The %d round till the %d sbox:%f\n",i,j,sumWeight(i));

				if((sumWeight(i)+B[N-i])>B_n_bar){
					
					if(j==8){
						Round_(i+1);
					}else{
						Round__(i,j+1);
					}
				}
			}
		}
	}
}

void Round_N_(int j){
	if(dx[N][j]==0){
		dy[N][j]=0;
		p_[N][j]=0;
		if(j==8){
			printAndSetBound();
		}else{
			Round_N_(j+1);
		}
	}else{
		p_[N][j]=DDT_MaxOutput[j-1][dx[N][j]];
		AddWeight(j,N);
		dy[N][j]=DDT_MaxOutput_Index[j-1][dx[N][j]];

		//fprintf(stream,"The N Round till the %d sbox:%f\n",j,sumWeight(N));
		//fprintf(stream,"sumWeight(N)>B_n_bar:%d\tB_n_bar:%f\n",sumWeight(N)>B_n_bar,B_n_bar);

		if(sumWeight(N)>B_n_bar){
			
			if(j==8){
				printAndSetBound();
			}else{
				Round_N_(j+1);
			}
		}
	}
}

void Round_(int i){
	u64 x_i_2;
	u32 x_i_2_EConv,y_i_1,x_i,y_i_1_P;
	SboxInput2word(&x_i_2, dx[i-2]+1);
	ExpansionConv1(&x_i_2_EConv,x_i_2);
	SboxOutput2word(&y_i_1, dy[i-1]+1);
	PermutationTL(&y_i_1_P,y_i_1);
	x_i=x_i_2_EConv^y_i_1_P;
	Expansion(dx[i]+1,x_i);

	//fprintf(stream,"Differential of the %d round:\t",i);
	//fprintnum(dx[i]+1,stream);
	//fprint8t8(dx[i]+1,stream);

	if(i==N){
		p[N]=0;
		Round_N_(1);
	}else{
		p[i]=0;
		Round__(i,1);
	}
}

void Round_2_(int j){
	for(a[2][j]=a[2][j-1]+1;a[2][j]<=8;a[2][j]++){
		
		//ÍË³öÌõ¼þ
		if(a[2][j]==8){
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[2][j-1]!=7){
					break;
				}
			}
			if(j!=1||flag==1){
				dx[2][8]=0;
				ResetCharacter(a[2][j-1],8,2);
				if( 0==(dx[2][7]&0x3) && 0==(dx[2][1]&0x30) ){
					dy[2][8]=0;
					p_[2][j]=0;
					AddWeight(j,2);
					if((p[2]+p[1]+B[N-2])>=B_n_bar){
						Round_(3);
					}
				}
			}
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][8]=DDT_SearchInOrderX[7][frequency][index];
					ResetCharacter(a[2][j-1],8,2);
					if( (dx[2][8]&0x30)==((dx[2][7]&0x3)<<4) && (dx[2][8]&0x3)==((dx[2][1]&0x30)>>4) ){
						dy[2][8]=DDT_SearchInOrderY[7][frequency][index];
						p_[2][j]=DDT[7][dx[2][8]][dy[2][8]];
						AddWeight(j,2);
						if((p[2]+p[1]+B[N-2])>=B_n_bar){
							Round_(3);
						}
					}
				}
			}
		}else if(a[2][j]==1){
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][1]=DDT_SearchInOrderX[0][frequency][index];
					ResetCharacter(0,1,2);
					dy[2][1]=DDT_SearchInOrderY[0][frequency][index];
					p_[2][j]=DDT[0][dx[2][1]][dy[2][1]];
					AddWeight(j,2);
					if((p[2]+p[1]+B[N-2])>=B_n_bar){
						Round_2_(j+1);
					}
				}
			}
		}else{
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[j]!=(a[j-1]+1)){
					break;
				}
			}
			for(int frequency=8;frequency>0;frequency--){
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][a[2][j]]=DDT_SearchInOrderX[a[2][j]-1][frequency][index];
					ResetCharacter(a[2][j-1],a[2][j],2);
					if( (dx[2][a[2][j]]&0x30) == ((dx[2][a[2][j]-1]&0x3)<<4) ){
						dy[2][a[2][j]]=DDT_SearchInOrderY[a[2][j]-1][frequency][index];
						p_[2][j]=DDT[a[2][j]-1][dx[2][a[2][j]]][dy[2][a[2][j]]];
						AddWeight(j,2);
						if((p[2]+p[1]+B[N-2])>=B_n_bar){
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
			print8t8(dx[1]+1);
		}
		fflush(stream);
	}
	fclose(stream);
}