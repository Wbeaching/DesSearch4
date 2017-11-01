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

flag=2����ǰ��S���޹�
flag=0��ǰһ��S�в���Ծ��prefix=00****
flag=1��ǰһ��S�л�Ծ��prefix=??****
*/
/*
u8 prefix[N+1];
u8 suffix[N+1];

void Round_2_(int j,u8 prefixThisRound,u8 suffixThisRound){
	double p2[8];
	
	for(a[j]=a[j-1]+1;a[j]<=8;a[j]++){//���Դ�8��ʼ��
		if(j==1){
			for(u8 x=0x1;x<=0xf;x){
				dx[2][a[j]]=prefixThisRound|x;
				for(dy[2][a[j]]=0x1;dy[2][a[j]]<=0xf;dy[2][a[j]]++){
					p_[2][j-1]=DDT [a[j]-1] [dx[2][a[j]]] [dy[2][a[j]]];
					
					p[2]=0;
					for(int i=0;i<j;i++){
						p[2]+=p_[2][i];
					}//p_2���

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
					}//p_2���

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
//����8��s�е������֣�����֮������໥��Լ
double Prob[8]={64,64,64,64,64,64,64,64};    //��¼����Ƶ�����Զ�Ϊ�׵�����ֵ,��ʼֵ��Ϊ���ֵ6
double P2=0;                              //�� round_two_j���õ��������������
int deltX[8]={0},deltY[8]={0};
void round_three(){
	printf("(%f,%f,%f,%f,%f,%f,%f,%f)\n",Prob[0],Prob[1],Prob[2],Prob[3],Prob[4],Prob[5],Prob[6],Prob[7]);
}
void round_two_j(int a_count){  
	int aj=(a_count)+1;
	//�ӵ�aj��s���ѵ��ڰ˸�s��
	for(aj;aj<9;aj++){                                                    
		printf("aj=%x\n",aj);                                                                  //��������������S��֮��������ִ����໥��Լ
		if(((aj-a_count)>1)&&(((deltX[a_count-1])&0x3)!=0))return;         
		//���ǰһ���ѵ���S4����һ�ֵݹ���S6��S7��S8����ôdeltX[S3]�������λ����Ϊ0���������ѡ�����˵�������S2��Ծ����S3��Ծ�����ʧ���ˣ�Ҫ��S4��Ծ�������deltX[1]�ĺ���λ��Ϊ�㣬�Ǳ��뷵����һ��ݹ飬��һ��deltX[1]�ĺ���λ�������ֵ������Ӧ����return
		for(int i=(a_count)+1;i<aj;i++){      //��a_count+1����aj-1��s�ж�����Ծ������������Ϊ�㡣
			Prob[i-1]=64;
			deltX[i-1]=deltY[i-1]=0;
		}

		int Si=aj-1;
		if((aj==8)&&((deltX[6]&0x3)==0)){                             //���ȿ���S8����Ծ�����,S7�������λ����Ϊ0
			double temp=Prob[7];
			Prob[7]=64;
			deltX[Si]=0;
			deltY[Si]=0;
			round_three();
			Prob[7]=temp;                  //���round_three�Ѳ�����������أ�������S8��Ծ�����	
		}
		for(int count=8;count>0;count--){                           //count��ȡֵΪ0-8
			printf("count=%d\n",count);
			//�жϸ����Ƿ������֦����
			for(int i=aj;i<9;i++){     
			Prob[aj]=64;
			deltX[aj]=deltY[aj]=0;
		}
			double temp=Prob[Si];
			Prob[Si]=2*count;
			P2=Prob[0]+Prob[1]+Prob[2]+Prob[3]+Prob[4]   +Prob[5]+Prob[6]+Prob[7];
			if(P2<462) {
				Prob[Si]=temp;
				break;                            //�������countѭ����������һ��aj������
			}

			for(int index=0;index<256;index++){
				printf("index=%x,DDT_SearchInOrderX=%x\n",index,DDT_SearchInOrderX[Si][count][index]);
				if(DDT_SearchInOrderX[Si][count][index]==0) break;     //���������꣬������һ��count������
				//��������������S��֮��������ִ����໥��Լ
				int deltx=DDT_SearchInOrderX[Si][count][index];
				if((Si==7)&&((((deltx>>4)&0x3)!=(deltX[Si-1]&0x3))||((deltx&0x3)!=((deltX[0]>>4)&0x3)))) continue;    //S8�����bit����S1ǰ��bit��S8ǰ2bit=S7��2bit���������㣬�������һ��index����
				else if(((deltx>>4)&0x3)!=(deltX[(Si+7)%8]&0x3))continue;                                      //��һ��S��ǰ2bit=ǰһ��S�к�2bit���������㣬�������һ��index����

				deltX[Si]=DDT_SearchInOrderX[Si][count][index];
				deltY[Si]=DDT_SearchInOrderY[Si][count][index];

				//���ֱ�����s�в��ǵڰ˸����ݹ������һ�� round_two_j
				if(aj<8)round_two_j(aj);
				//���ֱ�����s���ǵڰ˸����Ҹ��ʳ˻��㹻�󣬽�������ֵ�����
				//����û�а������ڰ˸�S�в���Ծ�����
				else round_three();
			}
		}

	}
}
*/
