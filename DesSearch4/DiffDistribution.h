#pragma once
void GenDiffDistributionTable();
void GenSearchInOrder();
extern double DDT[8][64][16];
extern u8 DDT_SearchInOrderX[8][9][512];
extern u8 DDT_SearchInOrderY[8][9][512];
extern int DDT_SearchInOrderLength[8][9];
void printMaxOutput();
void printDDT(int Si);
extern double DDT_MaxOutput[8][64];
extern u8 DDT_MaxOutput_Index[8][64];
void GenDiffDistributionTableMax();

extern u8 DDT_SearchInOrderWithFixedX[8][9][64][16];
extern int DDT_SearchInOrderWithFixedXLength[8][9][64];
void GenSearchInOrderWithFixedX();
void print(int SboxIndex,u8 inputMask);