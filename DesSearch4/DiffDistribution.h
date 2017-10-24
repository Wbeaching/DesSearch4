#pragma once
void GenDiffDistributionTable();
void GenDiffDistributionTableMax();
void GenSearchInOrder();
extern double DDT[8][64][16];
extern double DDT_MaxOutput[8][64];
extern unsigned int DDT_MaxOutput_Index[8][64];
extern u8 DDT_SearchInOrderX[8][9][512];
extern u8 DDT_SearchInOrderY[8][9][512];
extern int DDT_SearchInOrderLength[8][9];