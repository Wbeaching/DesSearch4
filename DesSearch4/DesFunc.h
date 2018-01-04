#include <fstream>

u32 genRandom32();
u64 genRandom64();

void Expansion(u8* output, u32 input);
void Permutation(u32* output, u32 input);
unsigned long long Expansion1(unsigned int x);
unsigned int Permutation1(unsigned int x);
void Substitution(u8* output,u8 input,int index);
void print32(const u32 x);
void print64(const u64 x);
void print8t8(const u8* y);
void pick6(u8* output, u32 input, int index);

void SboxInput2word(u64* output, u8* input);
void SboxOutput2word(u32* output, u8* input);
void word642SboxInput(u8* output, u64 input);
void word322SboxInput(u8* output, u32 input);
void word2bool48(bool* output, u64 input);
void bool2word48(u64* output, bool* input);

void PermutationConv(u32* output, u32 input);
void ExpansionConv1(u32* output, u64 input);
void ExpansionConv2(u32* output, u64 input);
void ExpansionUsingTable(u64* output, u32 input);

void fprint8t8(const u8* y,FILE* stream);
void fprintnum(const u8* y,FILE* stream);