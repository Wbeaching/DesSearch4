#pragma once
void GenDiffDistributionTable();
void GenSearchInOrder();
extern double DDT[8][64][16];
extern double DDT_int2DDT[9];
void GenDDT_int2DDT();
extern u8 DDT_SearchInOrderX[8][9][512];
extern u8 DDT_SearchInOrderY[8][9][512];
extern int DDT_SearchInOrderLength[8][9];

extern u8 DDT_SearchInOrderXWithPrefix[4][8][9][512];
extern u8 DDT_SearchInOrderYWithPrefix[4][8][9][512];
extern int DDT_SearchInOrderWithPrefixLength[4][8][9];
extern u8 DDT_SearchInOrderXWithBifix[4][4][8][9][512];
extern u8 DDT_SearchInOrderYWithBifix[4][4][8][9][512];
extern int DDT_SearchInOrderWithBifixLength[4][4][8][9];

void printMaxOutput();
void printDDT(int Si);

extern int DDT_int_MaxOutput[8][64];
extern double DDT_MaxOutput[8][64];
extern u8 DDT_MaxOutput_Index[8][64];
extern u8 DDT_MaxOutputs[8][64][16];
extern int DDT_MaxOutputsLength[8][64];

void GenDiffDistributionTableMax();

extern u8 DDT_SearchInOrderWithFixedX[8][9][64][16];
extern int DDT_SearchInOrderWithFixedXLength[8][9][64];
void GenSearchInOrderWithFixedX();
void print(int SboxIndex,u8 inputMask);

extern double DDTDESL[64][16];
extern u8 DDTDESL_SearchInOrderX[9][512];
extern u8 DDTDESL_SearchInOrderY[9][512];
extern int DDTDESL_SearchInOrderLength[9];
extern u8 DDTDESL_SearchInOrderWithFixedX[9][64][16];
extern int DDTDESL_SearchInOrderWithFixedXLength[9][64];
extern int DDTDESL_int_MaxOutput[64];
extern double DDTDESL_MaxOutput[64];
extern u8 DDTDESL_MaxOutputs[64][16];
extern int DDTDESL_MaxOutputsLength[64];