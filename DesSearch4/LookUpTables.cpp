#include "Types.h"
#include "DesFunc.h"
#include "LookUpTables.h"

u8 ETableLookUp[4][256][8] = { 0 };
u32 PTableLookUp[4][256] = { 0 };
u32 EConvTableLookUp[6][256]={0};
u8 SearchTable1[4][16]={0};//前两位，遍历中间两位
u8 SearchTable2[4][4][4]={0};//前两位，后两位，遍历中间两位

void set(u8* y,const u8* x){
	for(int i=0;i<8;i++){
		y[i]=x[i];
	}
}

void GenETableLookUp(){
	u32 x=0;
	u8 y[8]={0};
	for(int i=0;i<4;i++){
		for(u32 j=0;j<256;j++){
			x=j<<(8*i);
			Expansion(y,x);
			set(ETableLookUp[i][j],y);
		}
	}
}

void ExpansionTL(u8* output, u32 input){
	u32 x=input;
	u8 y[4][8]={0};
	for(int i=0;i<4;i++){
		set(y[i],ETableLookUp[i][(x>>(8*i))&0xff]);
	}
	for(int i=0;i<8;i++){
		output[i]=y[0][i]|y[1][i]|y[2][i]|y[3][i];
	}
}

void GenEConvTableLookUp(){
	u64 x=0;
	u32 y=0;
	for(int i=0;i<6;i++){
		for(u32 j=0;j<256;j++){
			x=j<<(8*i);
			ExpansionConv1(&y,x);
			EConvTableLookUp[i][j]=y;
		}
	}
}

void ExpansionConvTL(u32* output, u64 input){
	u64 x=input;
	u32 y[6]={0};
	for(int i=0;i<6;i++){
		y[i]=EConvTableLookUp[i][(x>>(8*i))&0xff];
	}
	*output=y[0]^y[1]^y[2]^y[3]^y[4]^y[5];
}

void GenPTableLookUp(){
	u32 x=0,y=0;
	for(int i=0;i<4;i++){
		for(u32 j=0;j<256;j++){
			x=j<<(8*i);
			Permutation(&y,x);
			PTableLookUp[i][j]=y;
		}
	}
}

void PermutationTL(u32* output, u32 input){
	u32 x=input,y[4]={0};
	for(int i=0;i<4;i++){
		y[i]=PTableLookUp[i][(x>>(8*i))&0xff];
	}
	*output=y[0]^y[1]^y[2]^y[3];
}

void GenSearchTables(){
	u8 temp;
	for(u8 i=0;i<4;i++){
		temp=i<<4;
		for(u8 j=0;j<16;j++){
			SearchTable1[i][j]=j|temp;
		}
		for(u8 j=0;j<4;j++){
			for(u8 k=0;k<4;k++){
				SearchTable2[i][j][k]=temp|(k<<2)|j;
			}
		}
	}
}