#pragma once
void GenDiffDistributionTable();
void GenSearchInOrder();
extern double DDT[8][64][16];
extern u8 DDT_SearchInOrderX[8][9][512];
extern u8 DDT_SearchInOrderY[8][9][512];
extern int DDT_SearchInOrderLength[8][9];
void print(int k);