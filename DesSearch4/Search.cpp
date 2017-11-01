#include "Search.h"
#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"
#include "DiffDistribution.h"
#include <stdio.h>
#define N 3

double B[N+1];
double B_n_bar;

int a[9]={0};
double p[N+1]={0};
double p_[N+1][9];

u8 dx[N+1][9]={0};
u8 dy[N+1][9]={0};

void Round_3(){
	printf("dx:");
	for(int i=1;i<=8;i++){
		printf("%x ",dx[2][i]);
	}
	printf("\ndy:");
	for(int i=1;i<=8;i++){
		printf("%x ",dy[2][i]);
	}
	printf("\np:%f\n==============\n",p[2]);
}

void Round_2_(int j){
	for(a[j]=a[j-1]+1;a[j]<=8;a[j]++){
		if(j==1){
			if(a[j]==1){
				for(int count=8;count>0;count--){
					for(int index=0;index<DDT_SearchInOrderLength[a[j]-1][count];index++){
						dx[2][a[j]]=DDT_SearchInOrderX[a[j]-1][count][index];
						dy[2][a[j]]=DDT_SearchInOrderY[a[j]-1][count][index];
						p_[2][j]=DDT[a[j]-1][dx[2][a[j]]][dy[2][a[j]]];
						p[2]=0;
						for(int k=1;k<=j;k++){
							p[2]+=p_[2][k];
						}
						//if((p[1]+p[2]+B[N-2])>=B_n_bar){
						if(p[2]>-3){
							Round_3();
							//Round_2_(j+1);
						}
					}
				}
			}
			else{
				for(int count=8;count>0;count--){
					for(int index=0;index<DDT_SearchInOrderLength[a[j]-1][count];index++){
						dx[2][a[j]]=DDT_SearchInOrderX[a[j]-1][count][index];
						if((dx[2][a[j]]&0x30)!=0x00){
							continue;
						}
						for(int i=a[j-1]+1;i<a[j];i++){
							dx[2][i]=0;
							dy[2][i]=0;
						}
						dy[2][a[j]]=DDT_SearchInOrderY[a[j]-1][count][index];
						p_[2][j]=DDT[a[j]-1][dx[2][a[j]]][dy[2][a[j]]];
						p[2]=0;
						for(int k=1;k<=j;k++){
							p[2]+=p_[2][k];
						}
						//if((p[1]+p[2]+B[N-2])>=B_n_bar){
						if(p[2]>-3){
							Round_3();
						}
					}
				}
			}
		}
	}
}
							

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
	void Round_2_(1,0,0);
}

/*

flag=2：与前后S盒无关
flag=0：前一个S盒不活跃，prefix=00****
flag=1：前一个S盒活跃，prefix=??****
*/
/*
u8 prefix[N+1];
u8 suffix[N+1];

void Round_2_(int j,u8 prefixThisRound,u8 suffixThisRound){
	double p2[8];
	
	for(a[j]=a[j-1]+1;a[j]<=8;a[j]++){//可以从8开始减
		if(j==1){
			for(u8 x=0x1;x<=0xf;x){
				dx[2][a[j]]=prefixThisRound|x;
				for(dy[2][a[j]]=0x1;dy[2][a[j]]<=0xf;dy[2][a[j]]++){
					p_[2][j-1]=DDT [a[j]-1] [dx[2][a[j]]] [dy[2][a[j]]];
					
					p[2]=0;
					for(int i=0;i<j;i++){
						p[2]+=p_[2][i];
					}//p_2求和

					if(p[1]+p[2]+B[N-2]>=B_n_bar){
						prefix[2]=(dx[2][a[j]]<<4)&0x30;
						if(a[j]==1){
							suffix[2]=(dx[2][a[j]]>>4)&0x3;
						}else{
							suffix[2]=0&0x3;
						}
						Round_2_(2,prefix[2],0);
					}

					Round_(3);
				}
			}
		}
		else if(j!=8){
			for(dx[2][a[j]]=0x4;dx[2][a[j]]<=0x3f;dx[2][a[j]]++){
				for(dy[2][a[j]]=0x1;dy[2][a[j]]<=0xf;dy[2][a[j]]++){
					p_[2][j-1]=DDT [a[j]-1] [dx[2][a[j]]] [dy[2][a[j]]];
					
					p[2]=0;
					for(int i=0;i<j;i++){
						p[2]+=p_[2][i];
					}//p_2求和

					if(p[1]+p[2]+B[N-2]>=B_n_bar){
						prefix[2]=(dx[2][a[j]]<<4)&0x30;
						if(a[j]==1){
							suffix[2]=(dx[2][a[j]]>>4)&0x3;
						}else{
							suffix[2]=0&0x3;
						}
						Round_2_(2,prefix[2],0);
					}

					Round_3();
				}
			}
		}
	}


		for(u8 x=0;x<=0x3f;x++){
			for(u8 y=0;y<=0xf;y++){
				p_[2][j-1]=DDT[a[j]-1][x][y];
				p[2]=0;
				for(int i=0;i<j;i++){
					p[2]+=p_[2][i];
				}
				if(p[1]+p[2]+B[N-2]>=B_n_bar&&j!=8){
					prefix=(x<<4)&0x30;
					if(j==1){
						if(a[j]-a[j-1]==1){
							suffix=(x>>4)&0x3;
						}
						else{
							suffix=0&0x3;
						}
					}

					if(j==8)
*/					
/*
//遍历8个s盒的输入差分，它们之间存在相互制约
double Prob[8]={64,64,64,64,64,64,64,64};    //记录的是频数的以二为底的真数值,初始值设为最大值6
double P2=0;                              //对 round_two_j中用到的数组进行声明
int deltX[8]={0},deltY[8]={0};
void round_three(){
	printf("(%f,%f,%f,%f,%f,%f,%f,%f)\n",Prob[0],Prob[1],Prob[2],Prob[3],Prob[4],Prob[5],Prob[6],Prob[7]);
}
void round_two_j(int a_count){  
	int aj=(a_count)+1;
	//从第aj个s盒搜到第八个s盒
	for(aj;aj<9;aj++){                                                    
		printf("aj=%x\n",aj);                                                                  //限制条件，相邻S盒之间的输入差分存在相互制约
		if(((aj-a_count)>1)&&(((deltX[a_count-1])&0x3)!=0))return;         
		//如果前一轮搜的是S4，下一轮递归搜S6，S7，S8，那么deltX[S3]的最后两位必须为0，否则不能搜。举例说明，如果S2活跃，搜S3活跃的情况失败了，要搜S4活跃的情况，deltX[1]的后两位不为零，那必须返回上一层递归，找一个deltX[1]的后两位等于零的值，所以应该是return
		for(int i=(a_count)+1;i<aj;i++){      //第a_count+1到第aj-1个s盒都不活跃，输入输出差分为零。
			Prob[i-1]=64;
			deltX[i-1]=deltY[i-1]=0;
		}

		int Si=aj-1;
		if((aj==8)&&((deltX[6]&0x3)==0)){                             //首先考虑S8不活跃的情况,S7的最后两位必须为0
			double temp=Prob[7];
			Prob[7]=64;
			deltX[Si]=0;
			deltY[Si]=0;
			round_three();
			Prob[7]=temp;                  //如果round_three搜不到结果，返回，继续搜S8活跃的情况	
		}
		for(int count=8;count>0;count--){                           //count的取值为0-8
			printf("count=%d\n",count);
			//判断概率是否满足剪枝条件
			for(int i=aj;i<9;i++){     
			Prob[aj]=64;
			deltX[aj]=deltY[aj]=0;
		}
			double temp=Prob[Si];
			Prob[Si]=2*count;
			P2=Prob[0]+Prob[1]+Prob[2]+Prob[3]+Prob[4]   +Prob[5]+Prob[6]+Prob[7];
			if(P2<462) {
				Prob[Si]=temp;
				break;                            //跳出这个count循环，进行下一个aj的搜索
			}

			for(int index=0;index<256;index++){
				printf("index=%x,DDT_SearchInOrderX=%x\n",index,DDT_SearchInOrderX[Si][count][index]);
				if(DDT_SearchInOrderX[Si][count][index]==0) break;     //本行搜索完，进入下一个count的搜索
				//限制条件，相邻S盒之间的输入差分存在相互制约
				int deltx=DDT_SearchInOrderX[Si][count][index];
				if((Si==7)&&((((deltx>>4)&0x3)!=(deltX[Si-1]&0x3))||((deltx&0x3)!=((deltX[0]>>4)&0x3)))) continue;    //S8最后两bit等于S1前两bit，S8前2bit=S7后2bit，若不满足，则进入下一个index搜索
				else if(((deltx>>4)&0x3)!=(deltX[(Si+7)%8]&0x3))continue;                                      //后一个S盒前2bit=前一个S盒后2bit，若不满足，则进入下一个index搜索

				deltX[Si]=DDT_SearchInOrderX[Si][count][index];
				deltY[Si]=DDT_SearchInOrderY[Si][count][index];

				//本轮遍历的s盒不是第八个，递归进入下一个 round_two_j
				if(aj<8)round_two_j(aj);
				//本轮遍历的s盒是第八个，且概率乘积足够大，进入第三轮的搜索
				//但是没有包含进第八个S盒不活跃的情况
				else round_three();
			}
		}

	}
}
*/
