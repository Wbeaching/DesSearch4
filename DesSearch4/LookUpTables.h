void GenETableLookUp();
void ExpansionTL(u8* output, u32 input);
void GenPTableLookUp();
void PermutationTL(u32* output, u32 input);
void GenEConvTableLookUp();
void ExpansionConvTL(u32* output, u64 input);

extern u32 SearchTable1[4][16];//ǰ��λ�������м���λ
extern u32 SearchTable2[4][4][4];//ǰ��λ������λ�������м���λ
void GenSearchTables();