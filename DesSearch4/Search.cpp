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
	0,-2.0+1e-10,-4.0+1e-10,-9.607682,-13.215364,
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

inline void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<l;i++){
		dx[round][i]=0;
		dy[round][i]=0;
	}
}

inline void AddWeight(int j,int round){
	p[round]=0;
	for(int k=1;k<=j;k++){
		p[round]+=p_[round][k];
	}
}



inline double sumWeight(int m){
	double temp=0;
	for(int i=1;i<=m;i++){
		temp+=p[i];
	}
	return temp;
}

//------------------------------
//输出新找到的最佳特征，并设置概率下界
//------------------------------
void printAndSetBound(){
	B_n_bar=sumWeight(rounds);
//******************************
//找到新的最佳概率，则将概率下界设为它
//******************************

//------------------------------
//打印第1轮特征
//------------------------------
	fprintf(stream,"dx1:");
	for(int i=1;i<=8;i++){
		fprintf(stream,"%x ",dx[1][i]);
	}
	fprintf(stream,"\tp1:%f\n",p[1]);
//******************************
//
//******************************

//------------------------------
//打印第2~N-1轮特征
//------------------------------
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

//------------------------------
//打印第N轮特征
//------------------------------

//------------------------------
//打印目前的最佳概率
//------------------------------
	fprintf(stream,"B_n_bar:%f\n==============\n",B_n_bar);
}

void Round_(int i);

void Round__(int i,int j,double pr,double pr_round){
	if(dx[i][j]==0){
		dy[i][j]=0;
		p_[i][j]=0;
		if(j==8){
			p[i]=pr_round;
			Round_(i+1);
		}else{
			Round__(i,j+1,pr,pr_round);
		}
	}else{
		double prob;
		for(int frequency=8;frequency>0;frequency--){
			p_[i][j]=DDT_int2DDT[frequency];
			prob=p_[i][j]+pr_round;
			//AddWeight(j,i);
			if((pr+prob+B[rounds-i])>=B_n_bar){
				for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[i][j]];index++){
					dy[i][j]=DDT_SearchInOrderWithFixedX[j-1][frequency][dx[i][j]][index];
					if(j==8){
						p[i]=prob;
						Round_(i+1);
					}else{
						Round__(i,j+1,pr,prob);
					}
				}
			}else{
				break;
			}
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
	if(i==rounds){
		p[rounds]=0;
		Round_N_(1);
	}else{
		Round__(i,1,sumWeight(i-1),0);
	}
}

//------------------------------
//第二轮搜索子函数
//------------------------------
void Round_2_(int j){

	for(a[2][j]=a[2][j-1]+1;a[2][j]<=8;a[2][j]++){
		ResetCharacter(a[2][j-1],a[2][j],2);
//******************************
//将第a[2][j-1]个S盒至第a[2][j]个S盒之间的S盒输入输出全部置为0。
//******************************

		//bool breakflag=0;
//******************************
//每调用一次子函数时，新建breakflag，后面用于跳出双层循环
//******************************

		
//------------------------------
//a[2][j]==8
//------------------------------
		if(a[2][j]==8){
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[2][j-1]!=7){
					break;
				}
			}
//******************************
//a[2][j-1]后两个比特不为0，那么a[2][j]只能取a[2][j-1]+1；但若j=1，不存在这个问题。
//a[2][j]==8时，满足上述条件，a[2][j-1]!=7。
//******************************

			if(j!=1||activeflag==1){
				dx[2][8]=0;
//******************************
//当第一轮活跃时，第二轮才能不活跃。
//j不为1时，说明第二轮必不活跃；activeflag==1时，说明第一轮是活跃的。
//******************************

				if( 0==(dx[2][7]&0x3) && 0==(dx[2][1]&0x30) ){
					dy[2][8]=0;
					p_[2][j]=0;
					AddWeight(j,2);
					if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
						Round_(3);
					}
				}
			}
//******************************
//a[2][j]==8且dx[2][8]==0时，dx[1][7]的后两个比特和dx[1][1]的前两个比特必须为0。
//累加概率，剪枝，因a[2][j]==8，通过剪枝则进入下一轮。
//******************************

			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				AddWeight(j,2);
				if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][8]=DDT_SearchInOrderX[7][frequency][index];
						if( (dx[2][8]&0x30)==((dx[2][7]&0x3)<<4) && (dx[2][8]&0x3)==((dx[2][1]&0x30)>>4) ){
							dy[2][8]=DDT_SearchInOrderY[7][frequency][index];
							Round_(3);
						}
					}
				}else{
					break;
				}
			}
//******************************
//a[2][j]==8且dx[2][8]!=0时，遍历dx[1][8]的可能，其实这里dx[2][8]的自由度只有两个比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，因a[2][j]==8，通过剪枝则进入下一轮。
//******************************

//------------------------------
//a[2][j]==1
//------------------------------
		}else if(a[2][j]==1){
			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				AddWeight(j,2);
				if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][1]=DDT_SearchInOrderX[0][frequency][index];
						dy[2][1]=DDT_SearchInOrderY[0][frequency][index];
						Round_2_(j+1);
					}
				}else{
					break;
				}
			}
//******************************
//a[2][j]==1时，此时也有j=1，遍历dx[2][1]的非零可能，这里dx[2][1]的自由度有六个比特。
//注：第1个S盒自由度6比特，第2~7个S盒自由度4比特，第8个S盒自由度2比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，因a[2][j]!=8，通过剪枝则进入下一个S盒。
//******************************

//------------------------------
//a[2][j]==2~7
//------------------------------
		}else{
			if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
				if(a[2][j]!=(a[2][j-1]+1)){
					break;
				}
			}
//******************************
//a[2][j-1]后两个比特不为0，那么a[2][j]只能取a[2][j-1]+1；但若j=1，不存在这个问题。
//******************************

			for(int frequency=8;frequency>0;frequency--){
				p_[2][j]=DDT_int2DDT[frequency];
				AddWeight(j,2);
				if((p[2]+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][a[2][j]]=DDT_SearchInOrderX[a[2][j]-1][frequency][index];
						if( (dx[2][a[2][j]]&0x30) == ((dx[2][a[2][j]-1]&0x3)<<4) ){
							dy[2][a[2][j]]=DDT_SearchInOrderY[a[2][j]-1][frequency][index];
							Round_2_(j+1);
						}
					}					
				}else{
					break;
				}
			}
		}
//******************************
//a[2][j]在2至7之间时，遍历dx[2][a[2][j]]的非零可能，这里dx[2][a[2][j]]的自由度有四个比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************
	}
}

//------------------------------
//第二轮搜索
//------------------------------
void Round_2(){
	Round_2_(1);
}

//------------------------------
//第一轮搜索子函数
//------------------------------
void Round_1_(int j){
	for(a[1][j]=a[1][j-1]+1;a[1][j]<=8;a[1][j]++){

		ResetCharacter(a[1][j-1],a[1][j],1);
//******************************
//将第a[1][j-1]个S盒至第a[1][j]个S盒之间的S盒输入输出全部置为0。
//******************************

//------------------------------
//a[1][j]==8
//------------------------------
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
			for(u8 x=1;x<64;x++){
 				dx[1][8]=x;
 				if( (dx[1][8]&0x30)==((dx[1][7]&0x3)<<4) && (dx[1][8]&0x3)==((dx[1][1]&0x30)>>4) ){
 					p_[1][j]=DDT_MaxOutput[7][dx[1][8]];
 					AddWeight(j,1);
 					if((p[1]+B[rounds-1])>=B_n_bar){
 						Round_2();
					}
				}
			}
//******************************
//a[1][j]==8且dx[1][8]!=0时，遍历dx[1][8]的可能，其实这里dx[1][8]的自由度只有两个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，因a[1][j]==8，通过剪枝则进入下一轮。
//******************************

//------------------------------
//a[1][j]==1
//------------------------------
		}else if(a[1][j]==1){
			for(u8 x=1;x<64;x++){
 				dx[1][1]=x;
 				p_[1][j]=DDT_MaxOutput[0][dx[1][1]];
 				AddWeight(j,1);
 				if((p[1]+B[rounds-1])>=B_n_bar){
 					Round_1_(j+1);
				}
			}
//******************************
//a[1][j]==1时，此时也有j=1，遍历dx[1][1]的非零可能，这里dx[1][1]的自由度有六个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

//------------------------------
//a[1][j]==2~7
//------------------------------
		}else{
			if(j!=1 && (dx[1][a[1][j-1]]&0x3)!=0){
				if(a[1][j]!=(a[1][j-1]+1)){
					break;
				}
			}
//******************************
//a[1][j-1]后两个比特不为0，那么a[1][j]只能取a[1][j-1]+1；但若j=1，不存在这个问题。
//******************************

			for(u8 x=1;x<64;x++){
 				dx[1][a[1][j]]=x;
				if( (dx[1][a[1][j]]&0x30) == ((dx[1][a[1][j]-1]&0x3)<<4) ){
 					p_[1][j]=DDT_MaxOutput[a[1][j]-1][dx[1][a[1][j]]];
 					AddWeight(j,1);
 					if((p[1]+B[rounds-1])>=B_n_bar){
 						Round_1_(j+1);
					}
				}
			}
		}
//******************************
//a[1][j]在2至7之间时，遍历dx[1][a[1][j]]的非零可能，这里dx[1][a[1][j]]的自由度有四个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

	}
}

//------------------------------
//第一轮搜索
//------------------------------
void Round_1(){
	stream = fopen( "fprintf.txt", "w" );
	Round_1_(1);
	fclose(stream);
}