#include "Search.h"
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include <stdio.h>
#include <fstream>
#include <iostream>

/*double B[N]={0,
	0,-2.0,-4.0,-9.607684,-13.215366,
	-19.963827,-23.612152,-30.482869,-31.482869,-38.353586,
	-39.353586,-46.224303,-47.224303,-54.095020,-55.095020,
	-61.965737};*/
double B[N]={0,
	0,-2.0+0.000001,-4.0+0.000001,-9.607682,-13.215364,
	-19.963825,-23.612150,-30.482867,-31.482867,-38.353584,
	-39.353584,-46.224301,-47.224301,-54.095018,-55.095018,
	-61.965735};

double TestB[N]={0,
	0,-2.0,-4.0,-9.0,-13.0,
	-19.0,-23.0,-30.0,-31.0,-38.0,
	-39.0,-46.0,-47.0,-54.0,-55.0,
	-61.0};

int rounds;
double B_n_bar;
int num_a[N+1]={0};
int a[N+1][9]={0};
double p[N+1];
double p_[N+1][9]={0};

u8 dx[N+1][9]={0};
u8 dy[N+1][9]={0};

bool activeflag=1;
FILE* stream;

int freq1[8]={0};
int count1[8]={0};
int freqN[8]={0};
/*void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<=8;i++){
		if(i!=l){
			dx[round][i]=0;
			dy[round][i]=0;
		}
	}
}*/

void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<l;i++){
		dx[round][i]=0;
		dy[round][i]=0;
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
	//B_n_bar=sumWeight(rounds);
	fprintf(stream,"dx1:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dx[1][i]);
	}
	fprintf(stream,"\tp1:%f\n",p[1]);
	for(int r=2;r<=rounds;r++){
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
			
	fprintf(stream,"B_n_bar:%f\n==============\n",B_n_bar);
}

void Round_(int i);

void Round__(int i,int j){
	bool breakflag=0;
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
			p_[i][j]=DDT_int2DDT[frequency];
			for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[i][j]];index++){
				dy[i][j]=DDT_SearchInOrderWithFixedX[j-1][frequency][dx[i][j]][index];
				//p_[i][j]=DDT[j-1][dx[i][j]][dy[i][j]];
				AddWeight(j,i);

				//fprintf(stream,"The %d round till the %d sbox:%f\n",i,j,sumWeight(i));

				if((sumWeight(i)+B[rounds-i])>=B_n_bar){
				//if((sumWeight(i)+B[rounds-i])>(B_n_bar-1e-2)){
					if(j==8){
						Round_(i+1);
					}else{
						Round__(i,j+1);
					}
				}
				else{
					breakflag=1;
				}
			}
			if(breakflag==1) break;
		}
	}
}

void Round_N_(int j){
	if(dx[rounds][j]==0){
		dy[rounds][j]=0;
		p_[rounds][j]=0;
		if(j==8){
			printAndSetBound();
		}else{
			Round_N_(j+1);
		}
	}else{
		p_[rounds][j]=DDT_MaxOutput[j-1][dx[rounds][j]];
		AddWeight(j,rounds);
		dy[rounds][j]=DDT_MaxOutput_Index[j-1][dx[rounds][j]];

		//fprintf(stream,"The N Round till the %d sbox:%f\n",j,sumWeight(N));
		//fprintf(stream,"sumWeight(N)>B_n_bar:%d\tB_n_bar:%f\n",sumWeight(N)>B_n_bar,B_n_bar);

		if(sumWeight(rounds)>=B_n_bar){
			
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
	ExpansionConvTL(&x_i_2_EConv,x_i_2);
	SboxOutput2word(&y_i_1, dy[i-1]+1);
	PermutationTL(&y_i_1_P,y_i_1);
	x_i=x_i_2_EConv^y_i_1_P;
	
	Expansion(dx[i]+1,x_i);

	//fprintf(stream,"Differential of the %d round:\t",i);
	//fprintnum(dx[i]+1,stream);
	//fprint8t8(dx[i]+1,stream);

	if(i==rounds){
		p[rounds]=0;
		Round_N_(1);
	}else{
		p[i]=0;
		Round__(i,1);
	}
}

void Round_2_(int j){

	for(a[2][j]=a[2][j-1]+1;a[2][j]<=8;a[2][j]++){
		ResetCharacter(a[2][j-1],a[2][j],2);
		bool breakflag=0;
		//退出条件
		if(a[2][j]==8){
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[2][j-1]!=7){
					break;
				}
			}
			if(j!=1||activeflag==1){
				dx[2][8]=0;
				
				if( 0==(dx[2][7]&0x3) && 0==(dx[2][1]&0x30) ){
					dy[2][8]=0;
					p_[2][j]=0;
					AddWeight(j,2);
					if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
						Round_(3);
					}
				}
			}
			//int prefix=dx[2][7]&0x3,suffix=dx[2][1]>>4;
			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][8]=DDT_SearchInOrderX[7][frequency][index];
					//ResetCharacter(a[2][j-1],8,2);
					if( (dx[2][8]&0x30)==((dx[2][7]&0x3)<<4) && (dx[2][8]&0x3)==((dx[2][1]&0x30)>>4) ){
						dy[2][8]=DDT_SearchInOrderY[7][frequency][index];
				/*for(int index=0;index<DDT_SearchInOrderWithBifixLength[prefix][suffix][7][frequency];index++){
					dx[2][8]=DDT_SearchInOrderXWithBifix[prefix][suffix][7][frequency][index];
					dy[2][8]=DDT_SearchInOrderYWithBifix[prefix][suffix][7][frequency][index];
					ResetCharacter(a[2][j-1],8,2);*/

						//p_[2][j]=DDT[7][dx[2][8]][dy[2][8]];
						AddWeight(j,2);
						if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
							Round_(3);
						}else{
							breakflag=1;
						}
					}
					if(breakflag==1) break;
				}
			}
		}else if(a[2][j]==1){
			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][1]=DDT_SearchInOrderX[0][frequency][index];
					//ResetCharacter(0,1,2);
					dy[2][1]=DDT_SearchInOrderY[0][frequency][index];
					//p_[2][j]=DDT[0][dx[2][1]][dy[2][1]];
					AddWeight(j,2);
					if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
						Round_2_(j+1);
					}else{
						breakflag=1;
					}
				}
				if(breakflag==1) break;
			}
		}else{
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[2][j]!=(a[2][j-1]+1)){
					break;
				}
			}
			//int prefix=dx[2][a[2][j]-1]&0x3;
			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
					dx[2][a[2][j]]=DDT_SearchInOrderX[a[2][j]-1][frequency][index];
					//ResetCharacter(a[2][j-1],a[2][j],2);
					if( (dx[2][a[2][j]]&0x30) == ((dx[2][a[2][j]-1]&0x3)<<4) ){
						dy[2][a[2][j]]=DDT_SearchInOrderY[a[2][j]-1][frequency][index];
				/*for(int index=0;index<DDT_SearchInOrderWithPrefixLength[prefix][a[2][j]-1][frequency];index++){
					dx[2][a[2][j]]=DDT_SearchInOrderXWithPrefix[prefix][a[2][j]-1][frequency][index];
					dy[2][a[2][j]]=DDT_SearchInOrderYWithPrefix[prefix][a[2][j]-1][frequency][index];
					ResetCharacter(a[2][j-1],a[2][j],2);*/

						//p_[2][j]=DDT[a[2][j]-1][dx[2][a[2][j]]][dy[2][a[2][j]]];
						AddWeight(j,2);
						if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
							Round_2_(j+1);
						}else{
							breakflag=1;
						}
					}
				}
				if(breakflag==1) break;
			}
		}
	}
}

void Round_2(){
	Round_2_(1);
}

//------------------------------
//第一轮搜索
//------------------------------
void Round_1_(int j){
	num_a[1]=j;
	for(a[1][j]=a[1][j-1]+1;a[1][j]<=8;a[1][j]++){
		ResetCharacter(a[1][j-1],a[1][j],1);
//******************************
//将第a[1][j-1]个S盒至第a[1][j]个S盒之间的S盒输入输出全部置为0。
//******************************

		if(a[1][j]==8){
			if(j!=1 && (dx[1][a[1][j-1]]&0x3)!=0){
				if(a[1][j-1]!=7){
					break;
				}
			}
//******************************
//a[1][j-1]后两个比特不为0，那么a[1][j]只能取a[1][j-1]+1；但若j=1，不存在这个问题。
//a[1][j]==8时，满足上述条件，a[1][j-1]!=7。
//******************************

			if(j==1){
				activeflag=0;
			}
			dx[1][8]=0;
//******************************
//j==1且a[1][j]==8且dx[1][8]==0时，第一轮差分为0，activeflag置为0，这时第二轮必须活跃
//******************************

			if( 0==(dx[1][7]&0x3) && 0==(dx[1][1]&0x30) ){
				dy[1][8]=0;
				p_[1][j]=0;
				AddWeight(j,1);
				if((p[1]+B[rounds-1])>=B_n_bar){
					Round_2();
				}
			}
//******************************
//a[1][j]==8且dx[1][8]==0时，dx[1][7]的后两个比特和dx[1][1]的前两个比特必须为0。
//这里j不必为1，前面j==1用于判断第一轮差分是否为0。
//累加概率，剪枝，因a[1][j]==8，通过剪枝则进入下一轮。
//******************************
			
			activeflag=1;
			for(u8 dx_1_8=1;dx_1_8<64;dx_1_8++){
				dx[1][8]=dx_1_8;
				if( (dx_1_8&0x30)==((dx[1][7]&0x3)<<4) && (dx_1_8&0x3)==((dx[1][1]&0x30)>>4) ){
					for(int freq=8,count;freq>0;freq--){
						count=DDT_SearchInOrderWithFixedXLength[7][freq][dx_1_8];
						if(count!=0){
							freq1[7]=freq;
							count1[7]=count;
							p_[1][j]=DDT_int2DDT[freq];
							AddWeight(j,1);
							if((p[1]+B[rounds-1])>=B_n_bar){
								Round_2();
							}
						}
						break;
					}
				}
			}
//******************************
//a[1][j]==8且dx[1][8]!=0时，遍历dx[1][8]的可能，其实这里dx[1][8]的自由度只有两个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，因a[1][j]==8，通过剪枝则进入下一轮。
//******************************


		}else if(a[1][j]==1){
			for(u8 dx_1_1=1;dx_1_1<64;dx_1_1++){
				dx[1][1]=dx_1_1;
				for(int freq=8,count;freq>0;freq--){
					count=DDT_SearchInOrderWithFixedXLength[0][freq][dx_1_1];
					if(count!=0){
						freq1[0]=freq;
						count1[0]=count;
						p_[1][j]=DDT_int2DDT[freq];
						AddWeight(j,1);
						if((p[1]+B[rounds-1])>=B_n_bar){
							Round_1_(j+1);
						}
						break;
					}
				}
			}
//******************************
//a[1][j]==1时，此时也有j=1，遍历dx[1][1]的非零可能，这里dx[1][1]的自由度有六个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

		}else{
			if(j!=1 && (dx[1][a[1][j-1]]&0x3)!=0){
				if(a[1][j]!=(a[1][j-1]+1)){
					break;
				}
			}
//******************************
//a[1][j-1]后两个比特不为0，那么a[1][j]只能取a[1][j-1]+1；但若j=1，不存在这个问题。
//******************************

			for(u8 dx_1_a1j=1;dx_1_a1j<64;dx_1_a1j++){
				dx[1][a[1][j]]=dx_1_a1j;
				if( (dx[1][a[1][j]]&0x30) == ((dx[1][a[1][j]-1]&0x3)<<4) ){
					for(int freq=8,count;freq>0;freq--){
						count=DDT_SearchInOrderWithFixedXLength[a[1][j]-1][freq][dx_1_a1j];
						if(count!=0){
							freq1[a[1][j]-1]=freq;
							count1[a[1][j]-1]=count;
							p_[1][j]=DDT_int2DDT[freq];
							AddWeight(j,1);
							if((p[1]+B[rounds-1])>=B_n_bar){
								Round_1_(j+1);
							}
						}
						break;
					}

				}
			}
		}
//******************************
//a[1][j]在2至7之间时，遍历dx[1][1]的非零可能，这里dx[1][a[1][j]]的自由度有四个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

	}
}

void Round_1(){
	stream = fopen( "fprintf.txt", "w" );
	Round_1_(1);
	fclose(stream);
}