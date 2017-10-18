#include "Types.h"
#include "DesFunc.h"
#include "DiffDistribution.h"

unsigned int DDT[8][64][16]={0};

void GenDiffDistributionTable(){
	for(int Si=0;Si<8;Si++){
		for(u8 Input1=0;Input1<64;Input1++){
			for(u8 Input2=0;Input2<64;Input2++){
				u8 Output1,Output2;
				Substitution(&Output1,Input1,Si);
				Substitution(&Output2,Input2,Si);
				DDT[Si][Input1^Input2][Output1^Output2]++;
			}
		}
	}
}