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

int a[N+1][9]={0};
double p[N+1];

u8 dx[N+1][9]={0};
u8 dy[N+1][9]={0};

bool activeflag=1;
FILE* stream;

double pr_cut[N+1][9]={0};
int freq1[9]={0};
double pr_whole;

/*void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<=8;i++){
		if(i!=l){
			dx[round][i]=0;
			dy[round][i]=0;
		}
	}
}*/

double addPr(double pr1,double pr2){
	double i,j,t;
	if(pr1<=pr2){
		i=pr1;
		j=pr2;
	}else{
		j=pr1;
		i=pr2;
	}
	t=pow(2,i-j);
	t=log(1+t)/log(2.0);
	return t+j;
}

inline void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<l;i++){
		dx[round][i]=0;
		dy[round][i]=0;
	}
}

/*inline void AddWeight(int j,int round){
	p[round]=0;
	for(int k=1;k<=j;k++){
		p[round]+=p_[round][k];
	}
}*/



/*inline double sumWeight(int m){
	double temp=0;
	for(int i=1;i<=m;i++){
		temp+=p[i];
	}
	return temp;
}*/

//------------------------------
//������ҵ�����������������ø����½�
//------------------------------
int trailCount=0;
double characterPr=0;



void printAndSetBound(){
	trailCount++;
	if(trailCount==1){
		characterPr=pr_whole;
	}else{
		characterPr=addPr(pr_whole,characterPr);
	}

	//B_n_bar=pr_whole;//sumWeight(rounds);
//******************************
//�ҵ��µ���Ѹ��ʣ��򽫸����½���Ϊ��
//******************************

	u64 dx1,dx2;
	u32 dpR,dpL,dy1,dy1AfterP,dx2BeforeE;
	SboxInput2word(&dx1, dx[1]+1);
	SboxInput2word(&dx2, dx[2]+1);
	ExpansionConvTL(&dpR,dx1);
	ExpansionConvTL(&dx2BeforeE,dx2);
	SboxOutput2word(&dy1, dy[1]+1);
	PermutationTL(&dy1AfterP,dy1);
	dpL=dy1AfterP^dx2BeforeE;

	u64 dxN,dxN_1;
	u32 dcR,dcL,dyN,dyNAfterP,dxN_1BeforeE;
	SboxInput2word(&dxN, dx[rounds]+1);
	SboxInput2word(&dxN_1, dx[rounds-1]+1);
	ExpansionConvTL(&dcR,dxN);
	ExpansionConvTL(&dxN_1BeforeE,dxN_1);
	SboxOutput2word(&dyN, dy[rounds]+1);
	PermutationTL(&dyNAfterP,dyN);
	dcL=dyNAfterP^dxN_1BeforeE;

	fprintf(stream,"plaintext diff:%x %x\n",dpL,dpR);
	fprintf(stream,"ciphertext diff:%x %x\n",dcL,dcR);
//******************************
//�������Ĳ�������Ĳ��
//******************************

	for(int r=1;r<=rounds;r++){
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
	fprintf(stream,"B_n_bar:%f\n==============\n",pr_whole);
}


void setdy1(int j){
	if(dx[1][j]==0){
		dy[1][j]=0;
		if(j==8){
			printAndSetBound();
		}else{
			setdy1(j+1);
		}
	}else{
		for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][freq1[j]][dx[1][j]];index++){
			dy[1][j]=DDT_SearchInOrderWithFixedX[j-1][freq1[j]][dx[1][j]][index];
			if(j==8){
				printAndSetBound();
			}else{
				setdy1(j+1);
			}
		}
	}
}

void Setdy1(){
	setdy1(1);
}

void Round_(int i,double pr_former);

void Round__(int i,int j,double pr,double pr_round){
	if(dx[i][j]==0){
		dy[i][j]=0;
		if(j==8){
			p[i]=pr_round;
			Round_(i+1,p[i]+pr);
		}else{
			Round__(i,j+1,pr,pr_round);
		}
	}else{
		double prob;
		for(int frequency=8;frequency>0;frequency--){
			prob=DDT_int2DDT[frequency]+pr_round;
			if((pr+prob+pr_cut[i][j]+B[rounds-i])>=B_n_bar){
				for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[i][j]];index++){
					dy[i][j]=DDT_SearchInOrderWithFixedX[j-1][frequency][dx[i][j]][index];
					if(j==8){
						p[i]=prob;
						Round_(i+1,p[i]+pr);
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

void Round_N_(int j,double pr,double pr_round){
	if(dx[rounds][j]==0){
		dy[rounds][j]=0;
		if(j==8){
			p[rounds]=pr_round;
			pr_whole=p[rounds]+pr;
			Setdy1();
			//printAndSetBound(p[rounds]+pr);
		}else{
			Round_N_(j+1,pr,pr_round);
		}
	}else{
		double prob;
		for(int frequency=8;frequency>0;frequency--){
			prob=DDT_int2DDT[frequency]+pr_round;
		//prob=DDT_MaxOutput[j-1][dx[rounds][j]]+pr_round;
			if(prob+pr+pr_cut[rounds][j]>=B_n_bar){
				for(int index=0;index<DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[rounds][j]];index++){
					dy[rounds][j]=DDT_SearchInOrderWithFixedX[j-1][frequency][dx[rounds][j]][index];
					if(j==8){
						p[rounds]=prob;
						pr_whole=p[rounds]+pr;
						Setdy1();
						//printAndSetBound(p[rounds]+pr);
					}else{
						Round_N_(j+1,pr,prob);
					}
				}
			}
		}
	}
}

void Round_(int i,double pr_former){
	u64 x_i_2;
	u32 x_i_2_EConv,y_i_1,x_i,y_i_1_P;
	SboxInput2word(&x_i_2, dx[i-2]+1);
	ExpansionConvTL(&x_i_2_EConv,x_i_2);
	SboxOutput2word(&y_i_1, dy[i-1]+1);
	PermutationTL(&y_i_1_P,y_i_1);
	x_i=x_i_2_EConv^y_i_1_P;
	Expansion(dx[i]+1,x_i);

	for(int k=7;k>=0;k--){
		if(dx[i][k+1]!=0){
			pr_cut[i][k]=pr_cut[i][k+1]+DDT_MaxOutput[k][dx[i][k+1]];
		}
		else{
			pr_cut[i][k]=pr_cut[i][k+1];
		}
	}

	if((pr_former+pr_cut[i][0]+B[rounds-i])<B_n_bar) return;

	if(i==rounds){
		p[rounds]=0;
		Round_N_(1,pr_former,0);
	}else{
		Round__(i,1,pr_former,0);
	}
}

//------------------------------
//�ڶ��������Ӻ���
//------------------------------
void Round_2_(int j,double pr_round){
	double prob;
	for(a[2][j]=a[2][j-1]+1;a[2][j]<=8;a[2][j]++){
		if(j!=1 && (dx[2][a[2][j-1]]&0x3)!=0){
			if(a[2][j]!=(a[2][j-1]+1)){
				break;
			}
		}
		ResetCharacter(a[2][j-1],a[2][j],2);
//******************************
//����a[2][j-1]��S������a[2][j]��S��֮���S���������ȫ����Ϊ0��
//a[2][j-1]���������ز�Ϊ0����ôa[2][j]ֻ��ȡa[2][j-1]+1������j=1��������������⡣
//******************************
		
//------------------------------
//a[2][j]==8
//------------------------------
		if(a[2][j]==8){
			if(j!=1||activeflag==1){
				dx[2][8]=0;
//******************************
//����һ�ֻ�Ծʱ���ڶ��ֲ��ܲ���Ծ��
//j��Ϊ1ʱ��˵���ڶ��ֱز���Ծ��activeflag==1ʱ��˵����һ���ǻ�Ծ�ġ�
//******************************

				if( 0==(dx[2][7]&0x3) && 0==(dx[2][1]&0x30) ){
					dy[2][8]=0;
					if((pr_round+p[1]+B[rounds-2])>=B_n_bar){
						p[2]=pr_round;
						Round_(3,p[1]+p[2]);
					}
				}
			}
//******************************
//a[2][j]==8��dx[2][8]==0ʱ��dx[1][7]�ĺ��������غ�dx[1][1]��ǰ�������ر���Ϊ0��
//�ۼӸ��ʣ���֦����a[2][j]==8��ͨ����֦�������һ�֡�
//******************************

			for(int frequency=8;frequency>0;frequency--){
				prob=pr_round+DDT_int2DDT[frequency];
				if((prob+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][8]=DDT_SearchInOrderX[7][frequency][index];
						if( (dx[2][8]&0x30)==((dx[2][7]&0x3)<<4) && (dx[2][8]&0x3)==((dx[2][1]&0x30)>>4) ){
							dy[2][8]=DDT_SearchInOrderY[7][frequency][index];
							p[2]=prob;
							Round_(3,p[1]+p[2]);
						}
					}
				}else{
					break;
				}
			}
//******************************
//a[2][j]==8��dx[2][8]!=0ʱ������dx[1][8]�Ŀ��ܣ���ʵ����dx[2][8]�����ɶ�ֻ���������ء�
//���ո��ʴӴ�С�����������롢������ֵ�ԡ�һ��������֦���������������Ϊ֮��ĸ���һ����С��
//�ۼӸ��ʣ���֦����a[2][j]==8��ͨ����֦�������һ�֡�
//******************************

//------------------------------
//a[2][j]==1
//------------------------------
		}else if(a[2][j]==1){
			for(int frequency=8;frequency>0;frequency--){
				prob=pr_round+DDT_int2DDT[frequency];
				if((prob+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][1]=DDT_SearchInOrderX[0][frequency][index];
						dy[2][1]=DDT_SearchInOrderY[0][frequency][index];
						Round_2_(j+1,prob);
					}
				}else{
					break;
				}
			}
//******************************
//a[2][j]==1ʱ����ʱҲ��j=1������dx[2][1]�ķ�����ܣ�����dx[2][1]�����ɶ����������ء�
//ע����1��S�����ɶ�6���أ���2~7��S�����ɶ�4���أ���8��S�����ɶ�2���ء�
//���ո��ʴӴ�С�����������롢������ֵ�ԡ�һ��������֦���������������Ϊ֮��ĸ���һ����С��
//�ۼӸ��ʣ���֦����a[2][j]!=8��ͨ����֦�������һ��S�С�
//******************************

//------------------------------
//a[2][j]==2~7
//------------------------------
		}else{
			for(int frequency=8;frequency>0;frequency--){
				prob=pr_round+DDT_int2DDT[frequency];
				if((prob+p[1]+B[rounds-2])>=B_n_bar){
					for(int index=0;index<DDT_SearchInOrderLength[a[2][j]-1][frequency];index++){
						dx[2][a[2][j]]=DDT_SearchInOrderX[a[2][j]-1][frequency][index];
						if( (dx[2][a[2][j]]&0x30) == ((dx[2][a[2][j]-1]&0x3)<<4) ){
							dy[2][a[2][j]]=DDT_SearchInOrderY[a[2][j]-1][frequency][index];
							Round_2_(j+1,prob);
						}
					}					
				}else{
					break;
				}
			}
		}
//******************************
//a[2][j]��2��7֮��ʱ������dx[2][a[2][j]]�ķ�����ܣ�����dx[2][a[2][j]]�����ɶ����ĸ����ء�
//���ո��ʴӴ�С�����������롢������ֵ�ԡ�һ��������֦���������������Ϊ֮��ĸ���һ����С��
//�ۼӸ��ʣ���֦��ͨ����֦�������һ��S�С�
//******************************
	}
}

//------------------------------
//�ڶ�������
//------------------------------
void Round_2(){
	//printf("2\n");
	Round_2_(1,0);
}

//------------------------------
//��һ�������Ӻ���
//------------------------------
void Round_1_(int j,double pr_round){
	double prob;
	for(a[1][j]=a[1][j-1]+1;a[1][j]<=8;a[1][j]++){

		
		if(j!=1 && (dx[1][a[1][j-1]]&0x3)!=0){
			if(a[1][j]!=(a[1][j-1]+1)){
				break;
			}
		}
		ResetCharacter(a[1][j-1],a[1][j],1);
//******************************
//����a[1][j-1]��S������a[1][j]��S��֮���S���������ȫ����Ϊ0��
//a[1][j-1]���������ز�Ϊ0����ôa[1][j]ֻ��ȡa[1][j-1]+1������j==1��������������⡣
//ע��j!=1ʱa[1][j]!=1��
//******************************

//------------------------------
//a[1][j]==8
//------------------------------
		if(a[1][j]==8){
			dx[1][8]=0;
			if(j==1){
				activeflag=0;
				p[1]=pr_round;
				Round_2();
			}else{
				if( 0==(dx[1][7]&0x3) && 0==(dx[1][1]&0x30) ){
					p[1]=pr_round;
					Round_2();
				}
			}
//******************************
//j==1��a[1][j]==8��dx[1][8]==0ʱ����һ�ֲ��Ϊ0��activeflag��Ϊ0����ʱ�ڶ��ֱ����Ծ
//a[1][j]==8��dx[1][8]==0ʱ��dx[1][7]�ĺ��������غ�dx[1][1]��ǰ�������ر���Ϊ0��
//����j����Ϊ1��ǰ��j==1�����жϵ�һ�ֲ���Ƿ�Ϊ0��
//�ۼӸ��ʣ���֦����a[1][j]==8��ͨ����֦�������һ�֡�
//******************************
			
			activeflag=1;
			for(u8 x=1;x<64;x++){
 				dx[1][8]=x;
 				if( (x&0x30)==((dx[1][7]&0x3)<<4) && (x&0x3)==((dx[1][1]&0x30)>>4) ){
 					for(int frequency=DDT_int_MaxOutput[7][x];frequency>0;frequency--){
						if(DDT_SearchInOrderWithFixedXLength[a[1][j]-1][frequency][x]==0) continue;
 						prob=DDT_int2DDT[frequency]+pr_round;
 						if((prob+B[rounds-1])>=B_n_bar){
 							p[1]=prob;
							freq1[8]=frequency;
							Round_2();
						}else break;
					}
				}
			}
//******************************
//a[1][j]==8��dx[1][8]!=0ʱ������dx[1][8]�Ŀ��ܣ���ʵ����dx[1][8]�����ɶ�ֻ���������ء�
//�̶��������֣�������ȡ������������ֵ��
//�ۼӸ��ʣ���֦����a[1][j]==8��ͨ����֦�������һ�֡�
//******************************

//------------------------------
//a[1][j]==1
//------------------------------
		}else if(a[1][j]==1){
			for(u8 x=1;x<64;x++){
 				dx[1][1]=x;
				for(int frequency=DDT_int_MaxOutput[0][x];frequency>0;frequency--){
					if(DDT_SearchInOrderWithFixedXLength[a[1][j]-1][frequency][x]==0) continue;
					prob=DDT_int2DDT[frequency]+pr_round;
					if((prob+B[rounds-1])>=B_n_bar){
						freq1[1]=frequency;
						Round_1_(j+1,prob);
					}else break;
				}
			}
//******************************
//a[1][j]==1ʱ����ʱҲ��j=1������dx[1][1]�ķ�����ܣ�����dx[1][1]�����ɶ����������ء�
//�̶��������֣�������ȡ������������ֵ��
//�ۼӸ��ʣ���֦��ͨ����֦�������һ��S�С�
//******************************

//------------------------------
//a[1][j]==2~7
//------------------------------
		}else{
			for(u8 x=1;x<64;x++){
 				dx[1][a[1][j]]=x;
				if( (x&0x30) == ((dx[1][a[1][j]-1]&0x3)<<4) ){
					for(int frequency=DDT_int_MaxOutput[a[1][j]-1][x];frequency>0;frequency--){
						if(DDT_SearchInOrderWithFixedXLength[a[1][j]-1][frequency][x]==0) continue;
						prob=DDT_int2DDT[frequency]+pr_round;
						if((prob+B[rounds-1])>=B_n_bar){
							freq1[a[1][j]]=frequency;
 							Round_1_(j+1,prob);
						}else break;
					}
				}
			}
		}
//******************************
//a[1][j]��2��7֮��ʱ������dx[1][a[1][j]]�ķ�����ܣ�����dx[1][a[1][j]]�����ɶ����ĸ����ء�
//�̶��������֣�������ȡ������������ֵ��
//�ۼӸ��ʣ���֦��ͨ����֦�������һ��S�С�
//******************************

	}
}

//------------------------------
//��һ������
//------------------------------
void Round_1(){
	errno_t err;
	err = fopen_s(&stream, "fprintf.txt", "w" );
	Round_1_(1,0);
	fclose(stream);
}

void Round_1_Fix_(int j,double pr_round){
	//printf("j:%d\t",j);
	if(dx[1][j]==0){
		if(j==8){
			//printf("0-8\n");
			p[1]=pr_round;
			Round_2();
		}else{
			//printf("0!8\n");
			Round_1_Fix_(j+1,pr_round);
		}
	}else{
		//printf("!0\n");
		double prob;
		for(int frequency=DDT_int_MaxOutput[j-1][dx[1][j]];frequency>0;frequency--){
			if(DDT_SearchInOrderWithFixedXLength[j-1][frequency][dx[1][j]]==0) continue;
			prob=DDT_int2DDT[frequency]+pr_round;
			if((prob+pr_cut[1][j]+B[rounds-1])>=B_n_bar){
				freq1[j]=frequency;
				if(j==8){
					p[1]=prob;
					Round_2();
				}else{
					Round_1_Fix_(j+1,prob);
				}				
			}else break;
		}
	}
}

void Round_1_Fix(){
	errno_t err;
	err = fopen_s(&stream, "fprintf.txt", "w" );
	
	for(int k=7;k>=0;k--){
		if(dx[1][k+1]!=0){
			pr_cut[1][k]=pr_cut[1][k+1]+DDT_MaxOutput[k][dx[1][k+1]];
		}
		else{
			pr_cut[1][k]=pr_cut[1][k+1];
		}
	}

	if((pr_cut[1][0]+B[rounds-1])<B_n_bar) return;
	Round_1_Fix_(1,0);

	fclose(stream);
}