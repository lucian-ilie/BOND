#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdint.h>
#include <omp.h>

int OLIG_LENGTH = 50;
int seqSimilarity = 75;
int consecMatch = 15;
double oligc = 0.000001, saltc = 0.075;
int tminterval = 10;
int tmmin=1000, tmmax=-1;
int GC_CONT_MAX= 70;
int GC_CONT_MIN = 30;
int dimerLen = 15, dimerStg = 86;
int hpStem = 6, hpminLoop = 1, hpmaxLoop = 3;

unsigned int PRIME_MM;
unsigned int PRIME_FH = 3277283;
unsigned int PRIME_IH = 850027;
int POS_SIZE = 8;
int BLOCK = 64;
uint16_t preNOR[65536];
int MAX_NUM_THREAD;
int geneno, olig_total;
int m;
int nbytes, gSize;
int gseqOff;
bool secStr = false;

using namespace std;

struct entry {
	uint64_t key;
	int *pos;
	int size;
};

void print_Time(){
	struct tm *current;
	time_t now;
	
	time(&now);
	current = localtime(&now);
	
	printf("Time is %i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
}

uint64_t get64bit(int rpos, uint8_t *gCoded){
	int lpos = rpos - 7;
	if (lpos < 0 && rpos >= 0) 
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (int i = lpos; i <= rpos; ++i) 
		genomeChunk = (uint64_t)(genomeChunk << 8 | gCoded[i]);
	return genomeChunk;	
}

uint64_t getvarbit(int byteno, long rpos, uint8_t *Coded){
	if (rpos < 0){
		cout << "error!!! in get64bit\n";
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	long lpos = rpos - byteno + 1;
	if (lpos < 0 && rpos >= 0)
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (long i = lpos; i <= rpos; ++i)
		genomeChunk = (uint64_t)(genomeChunk << 8 | Coded[i]);
	return genomeChunk;
}

uint64_t get64seed(const char *seedSeq){
	int seedWeight=0;
	int seedLength = strlen(seedSeq);
	uint64_t seed=0; 
	for(int j = 0; j < seedLength; ++j) {
		if(seedSeq[j] == '1')
		{
			++seedWeight;
			seed = (uint64_t) (seed << 2 | 3);
		}
		else
			seed = (uint64_t) (seed << 2 | 0);
	}
	return seed;
}

int insert(uint64_t k, int locn, entry *T){
	int i=0, j;
	do {
		j = (k + i) % PRIME_IH;
		if (T[j].pos == NULL){
			T[j].key = k;
			T[j].size = POS_SIZE;
			T[j].pos = (int *) malloc(T[j].size * sizeof(int));
			T[j].pos[0] = locn;
			for (int y = 1; y < T[j].size; ++y)
				T[j].pos[y] = -1;
			return j;
		}
		else if (T[j].key == k)
		{
			int lindex = 0;
			while (lindex < T[j].size && T[j].pos[lindex] != -1)
				++lindex;
			if (lindex < T[j].size)
				T[j].pos[lindex] = locn;
			else {
				T[j].size *= 2;
				T[j].pos = (int *) realloc(T[j].pos, T[j].size * sizeof(int));
				T[j].pos[lindex] = locn;
				for (int y = lindex + 1; y < T[j].size; ++y)
					T[j].pos[y] = -1;
			}
			return j;
		}			
		else				
			++i;
	}while (i!=PRIME_IH);
	
	T[0].pos = NULL;
	T[0].key = k;
	T[0].size = POS_SIZE;
	T[0].pos = (int *) malloc(T[0].size * sizeof(int));
	T[0].pos[0] = locn;
	for (int y = 1; y < T[0].size; ++y)
		T[0].pos[y] = -1;
	return 0;
}

void create_Hash(entry *hashTable, uint8_t *gCoded, uint64_t seed){
	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	int location = -1;
	
	for (int i=0; i < PRIME_IH; ++i) {
		hashTable[i].key = 0;
		hashTable[i].pos = NULL;
		hashTable[i].size = 0;
	}
	
	window = get64bit(nbytes-1, gCoded);
	for (int i = nbytes-1; i >= 12; --i) {		
		for (int pass=3; pass >=0 ; --pass) {	
			hashKey = (uint64_t) (window & seed);  
			location = i*4 + pass; 
			insert(hashKey, location, hashTable);
			window >>=  2;			
		}
		//if (i>=8)
		window |= (((uint64_t) gCoded[i-8]) << (BLOCK - 8));
	}
}

int search_Hash(uint64_t k, entry *T){
	int i=0, j;
	do {
		j = (int)((k + i) % PRIME_IH);
		if (T[j].key == k)
			return j;
		++i;
	}while (i!=PRIME_IH && T[j].pos!=NULL);
	return -1;
}

int check_2Bfactor(uint8_t *gCoded, int pose,  int posc){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
	chitr = get64bit(cbyte, gCoded); chitl = get64bit(cbyte - 8, gCoded);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	chitr >>= 2*(3-coff);
	chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	int kane75;
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity.
	
	if (kane75 > seqSimilarity)
		return 1;
	
	return 0;
}

bool second_Check2B(uint64_t seed, int seedLen, entry *hashTable, uint8_t *gCoded, int pose){
	////////////////////////////////////////////////////////////////Computing Hash value
	uint64_t hashKey=0, ehitr=0, ehitl=0;
	int w; int pos;
	int ebyte, eoff;;
	////////////////////////////////////////////////////////////////////Checking Positions
	ebyte = pose/4; eoff = pose%4;
	ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
		uint8_t elbyte = gCoded[ebyte-16];
		ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	}
	////////////////////////////////////////////////////////////////////	
	int offpos = 0;
	for (int pass=1; pass <= OLIG_LENGTH - seedLen + 1; ++pass){
		hashKey = (uint64_t) (ehitr & seed);  
		pos = search_Hash(hashKey, hashTable);											
		w=0;
		while (w < hashTable[pos].size && hashTable[pos].pos[w] !=-1 ){
			if (hashTable[pos].pos[w] > OLIG_LENGTH - 1 && hashTable[pos].pos[w] + offpos != pose){
				if (check_2Bfactor(gCoded, pose, hashTable[pos].pos[w] + offpos) == 0)		
					++w;
				else
					return false;
				
			}
			else
				++w;
		}
		
		ehitr >>= 2;
		ehitr |= (ehitl & (uint64_t) (3)) << (BLOCK - 2);
		ehitl >>= 2;
		
		++offpos;
	}
	return true;
}

int check_3Bfactor(uint8_t *gCoded, int pose,  int posc){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl=0, ehitm=0, ehitr=0, chitl=0, chitm=0, chitr=0, nolxnl=0, nolxnm=0, nolxnr=0;
	
	ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
	chitr = get64bit(cbyte, gCoded); chitm = get64bit(cbyte - 8, gCoded); chitl = get64bit(cbyte - 16, gCoded);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitm >>= 2*(3-eoff);
	ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	chitr >>= 2*(3-coff);
	chitr |= (chitm & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitm >>= 2*(3-coff);
	chitm |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	int kane75;
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnm = ~(ehitm ^ chitm);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnm >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 8; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity.
	
	if (kane75 > seqSimilarity)
		return 1;
	
	return 0;
}

bool second_Check3B(uint64_t seed, int seedLen, entry *hashTable, uint8_t *gCoded, int pose){
	////////////////////////////////////////////////////////////////Computing Hash value
	uint64_t hashKey=0, ehitr=0, ehitm=0, ehitl=0;
	int w; int pos;
	int ebyte, eoff;;
	////////////////////////////////////////////////////////////////////Checking Positions
	ebyte = pose/4; eoff = pose%4;
	ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitm >>= 2*(3-eoff);
	ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	////////////////////////////////////////////////////////////////////	
	int offpos = 0;
	for (int pass=1; pass <= OLIG_LENGTH - seedLen + 1; ++pass){
		hashKey = (uint64_t) (ehitr & seed);  
		pos = search_Hash(hashKey, hashTable);											
		w=0;
		while (w < hashTable[pos].size && hashTable[pos].pos[w] !=-1 ){
			if (hashTable[pos].pos[w] > OLIG_LENGTH - 1 && hashTable[pos].pos[w] + offpos != pose){
				if (check_3Bfactor(gCoded, pose, hashTable[pos].pos[w] + offpos) == 0)		
					++w;
				else
					return false;
				
			}
			else
				++w;
		}
		
		ehitr >>= 2;
		ehitr |= (ehitm & (uint64_t) (3)) << (BLOCK - 2);
		ehitm >>= 2;
		ehitm |= (ehitl & (uint64_t) (3)) << (BLOCK - 2);
		ehitl >>= 2;
		
		++offpos;
	}
	return true;
}

void kane_75l(uint8_t *gCoded, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded); //ehitl = get64bit(ebyte - 16, gCoded);
	chitr = get64bit(cbyte, gCoded); chitl = get64bit(cbyte - 8, gCoded); //chitl = get64bit(cbyte - 16, gCoded);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
		uint8_t elbyte = gCoded[ebyte-16];
		ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	}
	
	chitr >>= 2*(3-coff);
	chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	if (OLIG_LENGTH > 64 - (3-coff) && cbyte-16 >=0){
		uint8_t clbyte = gCoded[cbyte-16];
		chitl |= ((uint64_t)clbyte & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	}
	
	int kane75;
	
	
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity.
	
	int shiftno = 64 - OLIG_LENGTH - (eoff > coff? eoff: coff);
	//kane75 = nolxnor.count();
	//cout << nole << "\n" << nolc << "\n" << nolxnor << "\n" << kane75 << "\n\n";
	
	while(kane75 > seqSimilarity -1 && posc > OLIG_LENGTH - 2){
		simarray[posc] = simarray[pose] = 1;
		--posc; --pose;
		
		if (shiftno > 0)
		{
			ehitr >>= 2;
			ehitr |= (ehitl & (uint64_t) (3)) << (BLOCK - 2);
			ehitl >>= 2;
			
			chitr >>= 2;
			chitr |= (chitl & (uint64_t) (3)) << (BLOCK - 2);
			chitl >>= 2;
			
			--shiftno;
		}
		else {
			////////////////////////////////// BEGIN: Reading next 50-mer
			ebyte = pose/4; eoff = pose%4; cbyte = posc/4; coff = posc%4;
			ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
			chitr = get64bit(cbyte, gCoded); chitl = get64bit(cbyte - 8, gCoded);
			
			ehitr >>= 2*(3-eoff);
			ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			ehitl >>= 2*(3-eoff);
			
			if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
				uint8_t elbyte = gCoded[ebyte-16];
				ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			}
			
			chitr >>= 2*(3-coff);
			chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
			chitl >>= 2*(3-coff);
			
			if (OLIG_LENGTH > 64 - (3-coff) && cbyte-16 >=0){
				uint8_t clbyte = gCoded[cbyte-16];
				chitl |= ((uint64_t)clbyte & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
			}		
			////////////////////////////////// END: Reading next 50-mer		
			shiftno = 64 - OLIG_LENGTH - (eoff > coff? eoff: coff);
		}
		
		
		////////////////////////////////// BEGIN: Checking for 75% similarity
		kane75 = 0;
		nolxnr = ~(ehitr ^ chitr);
		nolxnl = ~(ehitl ^ chitl);
		int i;
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
		for (i=0; i < (OLIG_LENGTH/8) - 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
		if (OLIG_LENGTH % 8 != 0)
			kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
		kane75 *= (100.0/OLIG_LENGTH);
		////////////////////////////////// END: Checking for 75% similarity.		
	}
}

void kane_75r(uint8_t *gCoded, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
	chitr = get64bit(cbyte, gCoded); chitl = get64bit(cbyte - 8, gCoded);
	
	
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
		uint8_t elbyte = gCoded[ebyte-16];
		ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	}
	
	chitr >>= 2*(3-coff);
	chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	if (OLIG_LENGTH > 64 - (3-coff) && cbyte-16 >=0){
		uint8_t clbyte = gCoded[cbyte-16];
		chitl |= ((uint64_t)clbyte & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	}
	
	int kane75;
	
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity.
	
	while(kane75 > seqSimilarity -1 && pose <= gSize - 1 ){
		simarray[posc] = simarray[pose] = 1;
		++posc; ++pose;
		
		////////////////////////////////// BEGIN: Reading next 50-mer
		ebyte = pose/4; eoff = pose%4; cbyte = posc/4; coff = posc%4;
		ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
		chitr = get64bit(cbyte, gCoded); chitl = get64bit(cbyte - 8, gCoded);
		
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);
		
		if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
			uint8_t elbyte = gCoded[ebyte-16];
			ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		}
		
		chitr >>= 2*(3-coff);
		chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
		chitl >>= 2*(3-coff);
		
		if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
			uint8_t clbyte = gCoded[cbyte-16];
			chitl |= ((uint64_t)clbyte & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
		}		
		////////////////////////////////// END: Reading next 50-mer		
		
		
		////////////////////////////////// BEGIN: Checking for 75% similarity
		kane75 = 0;
		nolxnr = ~(ehitr ^ chitr);
		nolxnl = ~(ehitl ^ chitl);
		int i;
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
		for (i=0; i < (OLIG_LENGTH/8) - 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
		if (OLIG_LENGTH % 8 != 0)
			kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
		kane75 *= (100.0/OLIG_LENGTH);
		////////////////////////////////// END: Checking for 75% similarity.		
	}
	
	
}

void kane_75l3(uint8_t *gCoded, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl=0, ehitm=0, ehitr=0, chitl=0, chitm=0, chitr=0, nolxnl=0, nolxnm=0, nolxnr=0;
	
	ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
	chitr = get64bit(cbyte, gCoded); chitm = get64bit(cbyte - 8, gCoded); chitl = get64bit(cbyte - 16, gCoded);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitm >>= 2*(3-eoff);
	ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	chitr >>= 2*(3-coff);
	chitr |= (chitm & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitm >>= 2*(3-coff);
	chitm |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	int kane75;
	
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnm = ~(ehitm ^ chitm);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnm >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 8; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity.
	
	int shiftno = 96 - OLIG_LENGTH - (eoff > coff? eoff: coff);
	
	while(kane75 > seqSimilarity -1 && posc > OLIG_LENGTH - 2){
		simarray[posc] = simarray[pose] = 1;
		--posc; --pose;
		
		if (shiftno > 0)
		{
			ehitr >>= 2;
			ehitr |= (ehitl & (uint64_t) (3)) << (BLOCK - 2);
			ehitl >>= 2;
			
			chitr >>= 2;
			chitr |= (chitl & (uint64_t) (3)) << (BLOCK - 2);
			chitl >>= 2;
			
			--shiftno;
		}
		else {
			////////////////////////////////// BEGIN: Reading next OLIG-LENGTH-mer
			ebyte = pose/4; eoff = pose%4; cbyte = posc/4; coff = posc%4;
			ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
			chitr = get64bit(cbyte, gCoded); chitm = get64bit(cbyte - 8, gCoded); chitl = get64bit(cbyte - 16, gCoded);
			
			ehitr >>= 2*(3-eoff);
			ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			ehitm >>= 2*(3-eoff);
			ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			ehitl >>= 2*(3-eoff);
			
			chitr >>= 2*(3-coff);
			chitr |= (chitm & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
			chitm >>= 2*(3-coff);
			chitm |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
			chitl >>= 2*(3-coff);		
			////////////////////////////////// END: Reading next OLIG-LENGTH-mer		
			shiftno = 96 - OLIG_LENGTH - (eoff > coff? eoff: coff);
		}
		
		////////////////////////////////// BEGIN: Checking for 75% similarity
		kane75 = 0;
		nolxnr = ~(ehitr ^ chitr);
		nolxnm = ~(ehitm ^ chitm);
		nolxnl = ~(ehitl ^ chitl);
		int i;
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnm >> 16*i))];
		for (i=0; i < (OLIG_LENGTH/8) - 8; ++i)
			kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
		if (OLIG_LENGTH % 8 != 0)
			kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
		kane75 *= (100.0/OLIG_LENGTH);		
		////////////////////////////////// END: Checking for 75% similarity.		
	}
	
	
}

void kane_75r3(uint8_t *gCoded, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl=0, ehitm=0, ehitr=0, chitl=0, chitm=0, chitr=0, nolxnl=0, nolxnm=0, nolxnr=0;
	
	ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
	chitr = get64bit(cbyte, gCoded); chitm = get64bit(cbyte - 8, gCoded); chitl = get64bit(cbyte - 16, gCoded);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitm >>= 2*(3-eoff);
	ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	chitr >>= 2*(3-coff);
	chitr |= (chitm & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitm >>= 2*(3-coff);
	chitm |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	int kane75;
	
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	nolxnm = ~(ehitm ^ chitm);
	nolxnl = ~(ehitl ^ chitl);
	int i;
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
	for (i=0; i < 4; ++i)
		kane75 += preNOR[((uint16_t) (nolxnm >> 16*i))];
	for (i=0; i < (OLIG_LENGTH/8) - 8; ++i)
		kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
	if (OLIG_LENGTH % 8 != 0)
		kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
	kane75 *= (100.0/OLIG_LENGTH);
	////////////////////////////////// END: Checking for 75% similarity
	
	while(kane75 > seqSimilarity -1 && pose <= gSize - 1 ){
		simarray[posc] = simarray[pose] = 1;
		++posc; ++pose;
		
		////////////////////////////////// BEGIN: Reading next OLIG-LENGTH-mer
		ebyte = pose/4; eoff = pose%4; cbyte = posc/4; coff = posc%4;
		ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
		chitr = get64bit(cbyte, gCoded); chitm = get64bit(cbyte - 8, gCoded); chitl = get64bit(cbyte - 16, gCoded);
		
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitm >>= 2*(3-eoff);
		ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);
		
		chitr >>= 2*(3-coff);
		chitr |= (chitm & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
		chitm >>= 2*(3-coff);
		chitm |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
		chitl >>= 2*(3-coff);		
		////////////////////////////////// END: Reading next OLIG-LENGTH-mer		
		
		////////////////////////////////// BEGIN: Checking for 75% similarity
		kane75 = 0;
		nolxnr = ~(ehitr ^ chitr);
		nolxnm = ~(ehitm ^ chitm);
		nolxnl = ~(ehitl ^ chitl);
		int i;
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnr >> 16*i))];
		for (i=0; i < 4; ++i)
			kane75 += preNOR[((uint16_t) (nolxnm >> 16*i))];
		for (i=0; i < (OLIG_LENGTH/8) - 8; ++i)
			kane75 += preNOR[((uint16_t) (nolxnl >> 16*i))];
		if (OLIG_LENGTH % 8 != 0)
			kane75 += preNOR[(uint16_t)( ((uint16_t) (nolxnl >> 16*i)) & (uint16_t) ( (1 << 2*(OLIG_LENGTH % 8)) - 1) )];
		kane75 *= (100.0/OLIG_LENGTH);		
		////////////////////////////////// END: Checking for 75% similarity.		
	}
}

void kane_15(uint8_t *gCoded, int pose,  int posc, uint8_t *simarray, uint8_t *elimTable){
	int rangelc = posc, rangerc = posc + OLIG_LENGTH - consecMatch < gSize? posc + OLIG_LENGTH - consecMatch : gSize -1;
	for (int i=rangelc; i <= rangerc ; ++i)
		simarray[i] = 1;
	elimTable[posc] = 1;
	if(!elimTable[pose]){
		int rangele = pose, rangere = pose + OLIG_LENGTH - consecMatch < gSize? pose + OLIG_LENGTH - consecMatch : gSize -1;
		for (int i=rangele; i <= rangere ; ++i)
			simarray[i] = 1;
		elimTable[pose] = 1;
	}
}

void gc_Content(uint8_t *gCoded, uint8_t *simarray){
	int byte = nbytes -1 , off, pos = gSize -1;
	uint64_t halfl, halfr, nolxnr, nolxnl; 
	halfr = get64bit(byte, gCoded); halfl = get64bit(byte - 8, gCoded);
	int gccount;
	
	////////////////////////////////// BEGIN: Checking for GC
	gccount = 0;
	nolxnr = halfr;
	nolxnl = halfl;
	
	int rmn = 0, lmn = 0;
	
	nolxnl = halfl;
	for (int i=1; i<=32; ++i) {
		if ((nolxnr & 3) == 1 || (nolxnr & 3) == 2){
			if (i==1)
				rmn = 2;
			gccount+=2;
		}
		nolxnr >>= 2;
	}
	for (int i=1; i<=OLIG_LENGTH - 32; ++i) {
		if ((nolxnl & 3) == 1 || (nolxnl & 3) == 2)
			gccount+=2;		
		nolxnl >>= 2;
	}
	////////////////////////////////// END: Checking for GC.
	
	int shiftno = 0;
	
	//cout << nole << "\n" << nolc << "\n" << nolxnor << "\n" << kane75 << "\n\n";
	
	while(pos > OLIG_LENGTH - 2){		
		if (gccount <= GC_CONT_MIN || gccount >= GC_CONT_MAX) 
			simarray[pos] = 1;		
		--pos;
		
		if (shiftno <= 10)
		{
			halfr >>= 2;
			halfr |= (halfl & (uint64_t) (3)) << (BLOCK - 2);
			halfl >>= 2;
			++shiftno;
		}
		else {
			////////////////////////////////// BEGIN: Reading next 50-mer
			byte = pos/4; off = pos%4; 
			halfr = get64bit(byte, gCoded); halfl = get64bit(byte - 8, gCoded);
			
			halfr >>= 2*(3-off);
			halfr |= (halfl & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));
			halfl >>= 2*(3-off);
			
			
			////////////////////////////////// END: Reading next 50-mer		
			shiftno = 0;
		}
		
		
		////////////////////////////////// BEGIN: Checking for GC	
		nolxnr = halfr;
		nolxnl = halfl;
		
		nolxnl >>= 34;
		lmn = ((nolxnl & 3) == 1 || (nolxnl & 3) == 2)? 2 : 0;
		
		gccount += lmn - rmn;
		
		rmn = ((nolxnr & 3) == 1 || (nolxnr & 3) == 2)?	2 : 0;
		
		
		////////////////////////////////// END: Checking for GC.		
	}
}

void melting_Temp(uint8_t *gCoded, uint8_t *simarray, int* geneLoc, float *tmarray){
	
	int byte, off;
	uint64_t halfl2=0, halfl=0, halfr=0, nolxnr=0, nolxnl=0, nolxnl2=0; 		
	
	const float enthalpy[] = {-7.9,-8.4,-7.8,-7.2,-8.5,-8.0,-10.6,-7.8,-8.2,-9.8,-8.0,-8.4,-7.2,-8.2,-8.5,-7.9};
	const float  entropy[] = {-22.2,-22.4,-21.0,-20.4,-22.7,-19.9,-27.2,-21.0,-22.2,-24.4,-19.9,-22.4,-21.3,-22.2,-22.7,-22.2};
	
	float enth, entr, Tm;	
	int length = OLIG_LENGTH;
	
	//while (pos > OLIG_LENGTH - 2) {
	for (int qos=0; qos < geneno; ++qos){
		for(int pos = geneLoc[qos] + OLIG_LENGTH - 1; pos < geneLoc[qos+1]; ++pos){
			if (simarray[pos] == 0){
				byte = pos/4; off = pos%4; 
				halfr = get64bit(byte, gCoded); halfl = get64bit(byte - 8, gCoded); 
				if (byte - 16 >= 0)
					halfl2 = get64bit(byte - 16, gCoded);
				
				halfr >>= 2*(3-off);
				halfr |= (halfl & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));			
				halfl >>= 2*(3-off);	
				halfl |= (halfl2 & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));			
				halfl2 >>= 2*(3-off);
				
				////////////////////////////////// BEGIN: Computing Tm		
				enth = 0; entr = 0;	
				nolxnr = halfr;
				nolxnl = halfl;
				nolxnl2 = halfl2;
				
				if ((nolxnr & 3) == 0 || (nolxnr & 3) == 3){enth+=2.3; entr+=4.1;}
				if ((nolxnr & 3) == 1 || (nolxnr & 3) == 2){enth+=0.1; entr+=-2.8;}
				
				for (int i=1; i<=32; ++i) {
					if (i !=32){
						enth+=enthalpy[nolxnr & 15];
						entr+=entropy[nolxnr & 15];
					}
					else{
						enth+=enthalpy[((nolxnl & 3) << 2) | (nolxnr & 3)];
						entr+=entropy[((nolxnl & 3) << 2) | (nolxnr & 3)];
					}					
					nolxnr >>= 2;
				}	
				int up1 = (length > 64? 64:length);
				for (int i=33; i <= up1; ++i) {
					if (i  != up1 ){
						enth+=enthalpy[nolxnl & 15];
						entr+=entropy[nolxnl & 15];							
					}
					else{				
						if (length == up1){
							if ((nolxnl & 3) == 0 || (nolxnl & 3) == 3){enth+=2.3; entr+=4.1;}
							if ((nolxnl & 3) == 1 || (nolxnl & 3) == 2){enth+=0.1; entr+=-2.8;}
						}
						
						if (length > up1){
							enth+=enthalpy[((nolxnl2 & 3) << 2) | (nolxnl & 3)];
							entr+=entropy[((nolxnl2 & 3) << 2) | (nolxnl & 3)];										
						}
					}
					nolxnl >>= 2;
				}	
				for (int i=65; i <= length; ++i) {		
					if (i != length){
						enth+=enthalpy[nolxnl2 & 15];
						entr+=entropy[nolxnl2 & 15];
					}		
					else{
						if ((nolxnl2 & 3) == 0 || (nolxnl2 & 3) == 3){enth+=2.3; entr+=4.1;}
						if ((nolxnl2 & 3) == 1 || (nolxnl2 & 3) == 2){enth+=0.1; entr+=-2.8;}
					}
					
					nolxnl2 >>= 2;
				}
				
				
				
				/*enth*=1000;
				 Tm = enth / (entr + 1.987*log(oligc/4)) + 12.0 * log10(saltc) -273.15;
				 */
				
				enth=enth*1000;
				Tm = enth / (entr + 1.987*log(oligc/4)) + 12.0 * log10(saltc) -273.15;
				
				
				////////////////////////////////// END: Computing Tm
				tmarray[pos] = Tm;
			}
		}
	}
	
	if (tmmax == -1 || tmmin == 1000){
		int count = 0;	double sumTm = 0.0, maxTm=0.0, minTm=100.0;
		//for (int i = OLIG_LENGTH - 1; i < gSize; ++i)
		for (int i=0; i < geneno; ++i)
			for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
				if (simarray[j] == 0){
					++count;
					sumTm += tmarray[j]; 
					if (tmarray[j]>maxTm) maxTm=tmarray[j];
					if (tmarray[j]<minTm) minTm=tmarray[j];
				}		
		float avtm = sumTm / count;
		printf("max_Tm: %.3f min_Tm: %.3f aveg_Tm: %.3f\n",maxTm,minTm,avtm);
		
		int uppTm,lowTm; uppTm=(int) (100*maxTm); lowTm=(int) (100*minTm); 
		//printf("uppTm: %d \t lowTm: %.d \n",uppTm,lowTm);	
		// Computing Upper Bound Considering the Tm Range 
		bool **tmFreq = (bool **) malloc(geneno * sizeof(bool *));
		for (int i=0; i < geneno; ++i)
			tmFreq[i] = (bool *) malloc ((uppTm+1)*sizeof(bool));
		int *sumPm = (int *) malloc (geneno * sizeof(int));
		
		for (int i=0; i < geneno; ++i){
			for (int j=0; j <= uppTm; ++j)
				tmFreq[i][j]=false;
			sumPm[i] = 0;
		}
		
		for (int i=0; i < geneno; ++i)
			for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
				if (simarray[j] == 0)
					tmFreq[i][(int)(100*tmarray[j])]=true;
		
		int tmRange = 100*tminterval, upperBound, maxUbound = 0, maxLindex = -1, maxRindex = -1;
		
		upperBound = 0;
		for (int k=0; k < geneno; ++k){
			for (int j=uppTm; j > uppTm - tmRange; --j)		
				if (tmFreq[k][j]) ++sumPm[k];
			if(sumPm[k] > 0) 
				++upperBound;
		}
		//printf("Upper Bound in the range [%d,%d) is: %d\n", uppTm - tmRange+1, uppTm+1, upperBound);
		if (upperBound > maxUbound) {maxUbound = upperBound; maxLindex = uppTm - tmRange+1; maxRindex = uppTm;}
		
		for (int i=uppTm-1; i >= lowTm + tmRange - 1; --i){
			for (int k=0; k < geneno; ++k){
				int deltaPm = 0;
				if(tmFreq[k][i+1]) --deltaPm;
				if(tmFreq[k][i-tmRange+1]) ++deltaPm;
				if(deltaPm==1){
					if(sumPm[k]==0) ++upperBound;
				}
				else if (deltaPm==-1){
					if(sumPm[k]==1)	--upperBound;
				}
				sumPm[k]+=deltaPm;
			}
			//printf("Upper Bound in the range [%d,%d) is: %d\n", i - tmRange+1, i+1, upperBound);
			if (upperBound > maxUbound) {maxUbound = upperBound; maxLindex = i - tmRange+1; maxRindex = i;}
		}
		printf("Optimal Tm range of length %d is [%.2f, %.2f)\n\n", tminterval, maxLindex/100.0, (maxRindex+1)/100.0);
		
		for (int i=0; i < geneno; ++i)
			for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
				if (simarray[j] == 0) 
					if(tmarray[j] >=  1.0*maxRindex/100.0+0.01  || tmarray[j] < 1.0*maxLindex/100.0 )
						simarray[j] = 1;
		
		for(int i=0 ; i < geneno ; ++i)
			free(tmFreq[i]);//delete tmFreq[i];
		free(tmFreq);//delete [] tmFreq;
		free(sumPm);//delete [] sumPm;
	}
	else{
		for (int i = OLIG_LENGTH - 1; i < gSize; ++i)
			if ((simarray[i] == 0) && (tmarray[i] > tmmax || tmarray[i] < tmmin))
				simarray[i] = 1;	
	}
}

int hash_Search(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_FH;
		if (T[j] == k)
			return j;
		++i;
	}while (i!=PRIME_FH && T[j]!=s+1);
	return -1;
}

int hash_Insert(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_FH;
		if (T[j] == s+1){
			T[j] = k;
			return j;
		}
		else
			++i;
	}while (i!=PRIME_FH);
	return 0;
}			

void similarity_Check(uint8_t *gCoded, uint8_t *simarray, const char *seedSeq){
	uint64_t seed=0; 
	seed = get64seed(seedSeq);	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	uint64_t *hashTable = (uint64_t *) malloc (PRIME_FH*sizeof(uint64_t));
	int *posTable = (int *) malloc (PRIME_FH*sizeof(int));
	for (int i=0; i < PRIME_FH; ++i) {
		hashTable[i] = seed + 1;
		posTable[i] = -1;
	}
	window = get64bit(nbytes-1, gCoded);
	int shiftc=0;
	
	if (OLIG_LENGTH <= 64)
		for (int i = gSize-1; i >= OLIG_LENGTH-1; --i) {				
			hashKey = (uint64_t) (window & seed);  
			int searchRes = hash_Search(hashKey, hashTable, seed);
			if (( searchRes != -1)){
				if (simarray[i] == 0 || simarray[posTable[searchRes]] == 0){
					kane_75l(gCoded, posTable[searchRes],  i, simarray);
					kane_75r(gCoded, posTable[searchRes],  i, simarray);
				}				
				posTable[searchRes] = i;
			}
			else
				posTable[hash_Insert(hashKey, hashTable, seed)] = i;
			window >>=  2; ++shiftc; 
			if (shiftc % 4 == 0)
				window |= (((uint64_t) gCoded[i/4-8]) << (BLOCK - 8));
		}
	else
		for (int i = gSize-1; i >= OLIG_LENGTH-1; --i) {				
			hashKey = (uint64_t) (window & seed);  
			int searchRes = hash_Search(hashKey, hashTable, seed);
			if (( searchRes != -1)){
				if (simarray[i] == 0 || simarray[posTable[searchRes]] == 0){
					kane_75l3(gCoded, posTable[searchRes],  i, simarray);
					kane_75r3(gCoded, posTable[searchRes],  i, simarray);
				}				
				posTable[searchRes] = i;
			}
			else
				posTable[hash_Insert(hashKey, hashTable, seed)] = i;
			window >>=  2; ++shiftc; 
			if (shiftc % 4 == 0)
				window |= (((uint64_t) gCoded[i/4-8]) << (BLOCK - 8));
		}
	
	free(hashTable);//delete [] hashTable;
	free(posTable);//delete [] posTable;
}

int hash_Search_IDEN(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_MM;
		if (T[j] == k)
			return j;
		++i;
	}while (i!=PRIME_MM && T[j]!=s+1);
	return -1;
}

int hash_Insert_IDEN(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_MM;
		if (T[j] == s+1){
			T[j] = k;
			return j;
		}
		else
			++i;
	}while (i!=PRIME_MM);
	return 0;
}

void identity_Check(uint8_t *gCoded, uint8_t *simarray){
	
	//Begin Coding "111111111111111"
	uint64_t seed=0; 
	for(int j = 0; j < consecMatch; ++j) 
		seed = (uint64_t) (seed << 2 | 3);
	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	uint8_t nexwin=0;
	uint64_t *hashTable = (uint64_t *) malloc(PRIME_MM*sizeof(uint64_t));
	int *posTable = (int *) malloc (PRIME_MM*sizeof(int));
	uint8_t *elimTable = (uint8_t *) malloc (PRIME_MM*sizeof(uint8_t));
	for (int i=0; i < PRIME_MM; ++i) {
		hashTable[i] = seed + 1;
		posTable[i] = -1;
		elimTable[i] = 0;
	}
	
	int SEED_BYTE = consecMatch % 4 == 0 ?  consecMatch / 4 : consecMatch / 4  + 1;
	int rpos = nbytes-1;
	window = getvarbit(SEED_BYTE, rpos, gCoded);
	rpos -= SEED_BYTE;
	int nexwinSize =   SEED_BYTE * 4 - consecMatch;
	if (nexwinSize == 0){
		nexwin = gCoded[rpos];
		--rpos;
		nexwinSize = 4;
	}
	else
		nexwin = (window >> 2*consecMatch) & ((1 << 2*nexwinSize)-1); 
	
	hashKey = (uint64_t) (window & seed);
	
	
	for (int location = gSize-1; location >= consecMatch -1; --location) {				
		int searchRes = hash_Search_IDEN(hashKey, hashTable, seed);
		if (( searchRes != -1)) 
			kane_15(gCoded, posTable[searchRes], location, simarray, elimTable);	
		else
			posTable[hash_Insert_IDEN(hashKey, hashTable, seed)] = location;
		
		hashKey = (hashKey >> 2) | ( ((uint64_t)(nexwin & 3)) << (consecMatch*2 - 2) );
		nexwin >>= 2;  --nexwinSize;
		if (nexwinSize == 0){
			nexwin = gCoded[rpos];
			--rpos;
			nexwinSize = 4;
		}
	}
	
	free (hashTable);//delete [] hashTable;
	free(posTable);//delete [] posTable;
	free(elimTable);//delete [] elimTable;
}

bool secondary_Structure(uint8_t *gCoded, int pose){
	int ebyte = pose/4, eoff = pose%4;
	uint64_t ehitl, ehitm, ehitr, dcwindl=0, rcwindl=0, dcwindm=0, rcwindm=0, dcwindr=0, rcwindr=0, tpwind=0;	
	//OLIG_LENGTH = 70;
	if(OLIG_LENGTH > 64){
		ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitm >>= 2*(3-eoff);
		ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);
		
		dcwindr = ehitr;
		tpwind = ~dcwindr;
		for(int i=1; i <= 32; ++i){
			rcwindl = (rcwindl << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindm = ehitm;
		tpwind = ~dcwindm;
		for(int i=33; i <= 64; ++i){
			rcwindm = (rcwindm << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindl = ehitl;
		tpwind = ~dcwindl;
		for(int i=65; i <= OLIG_LENGTH; ++i){
			rcwindr = (rcwindr << 2) | (tpwind & 3);
			tpwind >>= 2;
		}
		
		rcwindr |= ((uint64_t)(rcwindm & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindm >>= 2*(32-(OLIG_LENGTH % 32));
		
		rcwindm |= ((uint64_t)(rcwindl & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindl >>= 2*(32-(OLIG_LENGTH % 32));
		
		dcwindl &= (( (uint64_t)1 << 2*(OLIG_LENGTH % 32))-1);
		/*printf("dcwind: %16llX \t %16llX \t %16llX \n", dcwindl, dcwindm, dcwindr);
		 printf("rcwind: %16llX \t %16llX \t %16llX \n", rcwindl, rcwindm, rcwindr);*/
		
		// start checking for DIMER
		uint64_t dimerWin=0, dimerTmp=0, dimerSeed=0, dimerKey=0;
		dimerWin = dcwindr;
		dimerTmp = dcwindm;
		dimerSeed = ((uint64_t)1 << 2*dimerLen) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *dimerHash = (uint64_t *) malloc((OLIG_LENGTH-dimerLen+1)*sizeof(uint64_t));
		int dtmpCount=0;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			dimerHash[i] = dimerKey;
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) dimerTmp = dcwindl;
		}
		/*for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i)
		 printf("dimerHash[%d]: %16llX \n",i, dimerHash[i]);*/
		dimerWin = rcwindr;
		dimerTmp = rcwindm;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			//printf("dimerKey: %16llX \n", dimerKey);
			// start core
			for (int j=0; j <= OLIG_LENGTH-dimerLen; ++j){
				int strGncy = 0;
				uint64_t xnorWins = ~(dimerKey ^ dimerHash[j]);
				xnorWins &= dimerSeed;
				for (int k=0; k < (dimerLen-1)/8+1; ++k)
					strGncy += preNOR[((uint16_t) (xnorWins >> 16*k))];
				//cout << (float)strGncy / (float)(dimerLen) << "\n";
				if ( (float)strGncy / (float)(dimerLen) > (dimerStg*1.0)/100.0) {
					/*printf("secondary structure found!\n");
					 printf("dcwind: %16llX \t %16llX \t %16llX %16llX \n", dcwindl, dcwindm, dcwindr, dimerHash[j]);
					 printf("rcwind: %16llX \t %16llX \t %16llX %16llX \n", rcwindl, rcwindm, rcwindr, dimerKey);*/
					free(dimerHash);//delete [] dimerHash;
					return false;
				}
			}
			// end core
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) dimerTmp = rcwindl;
		}
		free(dimerHash);//delete [] dimerHash;
		// end checking for DIMER
		
		// start checking for HAIRPIN
		
		uint64_t hpWin=0, hpTmp=0, hpSeed=0, hpKey=0;
		hpWin = dcwindr;
		hpTmp = dcwindm;
		hpSeed = ((uint64_t)1 << 2*hpStem) - 1;
		uint64_t *hpHash = (uint64_t *) malloc ((OLIG_LENGTH-hpStem+1)*sizeof(uint64_t));
		dtmpCount=0;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			hpHash[i] = hpKey;
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) hpTmp = dcwindl;
		}
		//for (int i=0; i <= OLIG_LENGTH-hpStem; ++i)
		//printf("hpHash[%d]: %16llX \n",i, hpHash[i]);
		hpWin = rcwindr;
		hpTmp = rcwindm;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			for (int j=0; j <= OLIG_LENGTH-hpStem; ++j)
				if (hpKey==hpHash[j] && abs(i-j)>=hpminLoop && abs(i-j)<=hpmaxLoop){
					//printf("Hairpin found!\n");
					//printf("dcwind: %16llX \t %16llX \t %16llX %16llX \n", dcwindl, dcwindm, dcwindr, hpHash[j]);
					//printf("rcwind: %16llX \t %16llX \t %16llX %16llX \n", rcwindl, rcwindm, rcwindr, hpKey);
					free(hpHash);//delete [] hpHash;
					return false;
				}
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) hpTmp = rcwindl;
		}
		free(hpHash);//delete [] hpHash;
		// end checking for HAIRPIN
		
	}
	else{
		ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);// Remember ehitl is not a correct 32Mer (left side). it needs eoff char from ehitl2
		if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
			uint8_t elbyte = gCoded[ebyte-16];
			ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		}
		
		dcwindr = ehitr;
		tpwind = ~dcwindr;
		for(int i=1; i <= 32; ++i){
			rcwindl = (rcwindl << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindl = ehitl;
		tpwind = ~dcwindl;
		for(int i=1; i <= 18; ++i){
			rcwindr = (rcwindr << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}		
		rcwindr |= ((uint64_t)(rcwindl & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindl >>= 2*(32-(OLIG_LENGTH % 32));
		
		dcwindl &= (( (uint64_t)1 << 2*(OLIG_LENGTH % 32))-1);	
		//printf("dcwind: %16llX \t %16llX \n", dcwindl, dcwindr);
		//printf("rcwind: %16llX \t %16llX \n", rcwindl, rcwindr);
		
		// start checking for DIMER
		uint64_t dimerWin=0, dimerTmp=0, dimerSeed=0, dimerKey=0;
		dimerWin = dcwindr;
		dimerTmp = dcwindl;
		dimerSeed = ((uint64_t)1 << 2*dimerLen) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *dimerHash = (uint64_t *) malloc((OLIG_LENGTH-dimerLen+1)*sizeof(uint64_t));
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			dimerHash[i] = dimerKey;
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2;
		}
		/*for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i)
		 printf("dimerHash[%d]: %16llX \n",i, dimerHash[i]);*/
		dimerWin = rcwindr;
		dimerTmp = rcwindl;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			//printf("dimerKey: %16llX \n", dimerKey);
			// start core
			for (int j=0; j <= OLIG_LENGTH-dimerLen; ++j){
				int strGncy = 0;
				uint64_t xnorWins = ~(dimerKey ^ dimerHash[j]);
				xnorWins &= dimerSeed;
				for (int k=0; k < (dimerLen-1)/8+1; ++k)
					strGncy += preNOR[((uint16_t) (xnorWins >> 16*k))];
				//cout << (float)strGncy / (float)(dimerLen) << "\n";
				if ( (float)strGncy / (float)(dimerLen) > (dimerStg*1.0)/100.0) {
					//printf("secondary structure found!\n");
					//printf("dcwind: %16llX \t %16llX \t %16llX \n", dcwindl, dcwindr, dimerHash[j]);
					//printf("rcwind: %16llX \t %16llX \t %16llX \n", rcwindl, rcwindr, dimerKey);
					free(dimerHash);//delete [] dimerHash;
					return false;
				}
			}
			// end core
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2;
		}
		free(dimerHash);//delete [] dimerHash;
		// end checking for DIMER
		
		// start checking for HAIRPIN
		uint64_t hpWin=0, hpTmp=0, hpSeed=0, hpKey=0;
		hpWin = dcwindr;
		hpTmp = dcwindl;
		hpSeed = ((uint64_t)1 << 2*hpStem) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *hpHash = (uint64_t *) malloc ((OLIG_LENGTH-hpStem+1)*sizeof(uint64_t));
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			hpHash[i] = hpKey;
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2;
		}
		hpWin = rcwindr;
		hpTmp = rcwindl;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			for (int j=0; j <= OLIG_LENGTH-hpStem; ++j)
				if (hpKey==hpHash[j] && abs(i-j)>=hpminLoop && abs(i-j)<=hpmaxLoop){
					free(hpHash);//delete [] hpHash;
					return false;
				}
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2;
		}
		free(hpHash);//delete [] hpHash;
		// end checking for HAIRPIN		
		
	}
	return true;
}

unsigned long print_Oligo_S8W9(uint8_t *gCoded, uint8_t *simarray, unsigned char *gSequence, int *geneLoc, 
							   FILE *oligos, float *tmarray){
	unsigned long memoryUsed=0;
	
	const char * seed_Phase2 [] = {
		"111010001101011",
		"1100010010000100010111",
		"1011001010001000010011",
		"1010110000110000001011",
		"11010000001000101000111",
		"110001010000000100101011",
		"110100000100100000010000111",
		"1110010000001000010000010011"
	};
	
	entry *hashTable1 = (entry *) malloc(PRIME_IH * sizeof(entry));
   	uint64_t seedSeq1 = get64seed(seed_Phase2[0]);
	uint64_t seedLen1 = strlen(seed_Phase2[0]);
	
    entry *hashTable2 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq2 = get64seed(seed_Phase2[1]);
	uint64_t seedLen2 = strlen(seed_Phase2[1]);
	
	entry *hashTable3 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq3 = get64seed(seed_Phase2[2]);
	uint64_t seedLen3 = strlen(seed_Phase2[2]);
	
	entry *hashTable4 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq4 = get64seed(seed_Phase2[3]);
	uint64_t seedLen4 = strlen(seed_Phase2[3]);
	
	entry *hashTable5 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq5 = get64seed(seed_Phase2[4]);
	uint64_t seedLen5 = strlen(seed_Phase2[4]);
	
	entry *hashTable6 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq6 = get64seed(seed_Phase2[5]);
	uint64_t seedLen6 = strlen(seed_Phase2[5]);
	
	entry *hashTable7 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq7 = get64seed(seed_Phase2[6]);
	uint64_t seedLen7 = strlen(seed_Phase2[6]);
	
	entry *hashTable8 = (entry *) malloc(PRIME_IH * sizeof(entry));
	uint64_t seedSeq8 = get64seed(seed_Phase2[7]);
	uint64_t seedLen8 = strlen(seed_Phase2[7]);
	
	int nthreads, tid;
	if (MAX_NUM_THREAD >= 8)	
		omp_set_num_threads(8);
	else 
		omp_set_num_threads(MAX_NUM_THREAD);
	
	
#pragma omp parallel shared(gCoded, nthreads) private(tid)
	{		
#pragma omp sections
		{
#pragma omp section
			{
				create_Hash(hashTable1, gCoded, seedSeq1);
			}
#pragma omp section
			{
				create_Hash(hashTable2, gCoded, seedSeq2);
			}
#pragma omp section
			{
				create_Hash(hashTable3, gCoded, seedSeq3);
			}
#pragma omp section
			{
				create_Hash(hashTable4, gCoded, seedSeq4);
			}
#pragma omp section
			{
				create_Hash(hashTable5, gCoded, seedSeq5);
			}
#pragma omp section
			{
				create_Hash(hashTable6, gCoded, seedSeq6);
			}
#pragma omp section
			{
				create_Hash(hashTable7, gCoded, seedSeq7);
			}
#pragma omp section
			{
				create_Hash(hashTable8, gCoded, seedSeq8);
			}
		}
	}	
	
	olig_total = 0; 
	int olig_total_p1 = 0;
	for (int i=0; i < geneno; ++i){
		int j = geneLoc[i] + OLIG_LENGTH -1;
		while (j < geneLoc[i+1] && simarray[j] != 0)
			++j;
		if (j < geneLoc[i+1])
			++olig_total_p1;
	}
	
	
	//int olig_total = 0;
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if(simarray[j] == 0){++olig_total;break;}
	
	
	//cout << "Total number of unique oligos are: " << olig_total << " for "  << geneno  << " genes.\n";
	
	
	
	int *effect_eval = (int *) malloc(geneno*sizeof(int));
	memoryUsed += geneno * sizeof (int);
	
	int *Gene = (int *) malloc (geneno*sizeof(int));
	memoryUsed += geneno * sizeof (int);
	
	for (int i=0; i< geneno; ++i){
		effect_eval[i] = 0;
		Gene[i] = -1;
	}
	
	//printf("Max. number of threads = %d\n", MAX_NUM_THREAD);
	int chunk = 50, i;
	omp_set_num_threads(MAX_NUM_THREAD);	
    
	if(OLIG_LENGTH <= 64){
#pragma omp parallel for shared(gCoded,chunk) private(i) schedule(static,chunk)
		for (i=0; i < geneno; ++i){
			int seenZero, bestZero, seenPlace, bestPlace, canPlace, selPlace;
			bool flag_ok = false;
			
			while(!flag_ok){
				seenZero = 0, bestZero = 0, seenPlace = -1, bestPlace = -1, canPlace=-1, selPlace=-1;
				for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j){
					//if(j < gSize){	// to avoid segmentation fault error!!!		
					if (simarray[j] == 0){
						++seenZero;
						seenPlace = j;
						if (seenZero > bestZero) {bestZero = seenZero; bestPlace = seenPlace;}
					}
					else 
						seenZero = 0;
					//}
				}
				if (bestPlace != -1){				
					canPlace = bestPlace - bestZero/2;
					if (secStr)
						flag_ok = secondary_Structure(gCoded, canPlace);
					else 
						flag_ok = true;
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq1, seedLen1, hashTable1, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq2, seedLen2, hashTable2, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq3, seedLen3, hashTable3, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq4, seedLen4, hashTable4, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq5, seedLen5, hashTable5, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq6, seedLen6, hashTable6, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq7, seedLen7, hashTable7, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check2B(seedSeq8, seedLen8, hashTable8, gCoded, canPlace);
					
					if (flag_ok){
						selPlace = canPlace;
						//++olig_total;
						Gene[i] = canPlace+gseqOff;//fprintf(oligos, "Gene[%d]= %lu\n",i,canPlace+gseqOff);
						/*for (int k = canPlace - OLIG_LENGTH +1; k <= canPlace; ++k)
						 fprintf(oligos, "%c", gSequence[k + gseqOff]);
						 fprintf(oligos, "\n");				*/
						//break;
					}
					else{
						simarray[canPlace]=1;
						effect_eval[i] = 1;
					}
					
				}
				else 
					flag_ok = true;
				
			}
		}
	}
	else{
#pragma omp parallel for shared(gCoded,chunk) private(i) schedule(static,chunk)
		for (i=0; i < geneno; ++i){
			int seenZero, bestZero, seenPlace, bestPlace, canPlace, selPlace;
			bool flag_ok = false;
			
			while(!flag_ok){
				seenZero = 0, bestZero = 0, seenPlace = -1, bestPlace = -1, canPlace=-1, selPlace=-1;
				for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j){
					//if(j < gSize){	// to avoid segmentation fault error!!!		
					if (simarray[j] == 0){
						++seenZero;
						seenPlace = j;
						if (seenZero > bestZero) {bestZero = seenZero; bestPlace = seenPlace;}
					}
					else 
						seenZero = 0;
					//}
				}
				if (bestPlace != -1){				
					canPlace = bestPlace - bestZero/2;
					if (secStr)
						flag_ok = secondary_Structure(gCoded, canPlace);
					else 
						flag_ok = true;
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq1, seedLen1, hashTable1, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq2, seedLen2, hashTable2, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq3, seedLen3, hashTable3, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq4, seedLen4, hashTable4, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq5, seedLen5, hashTable5, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq6, seedLen6, hashTable6, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq7, seedLen7, hashTable7, gCoded, canPlace);
					if (flag_ok)
						flag_ok = second_Check3B(seedSeq8, seedLen8, hashTable8, gCoded, canPlace);
					
					if (flag_ok){
						selPlace = canPlace;
						//++olig_total;
						Gene[i] = canPlace+gseqOff;//fprintf(oligos, "Gene[%d]= %lu\n",i,canPlace+gseqOff);
						/*for (int k = canPlace - OLIG_LENGTH +1; k <= canPlace; ++k)
						 fprintf(oligos, "%c", gSequence[k + gseqOff]);
						 fprintf(oligos, "\n");				*/
						//break;
					}
					else{
						simarray[canPlace]=1;
						effect_eval[i] = 1;
					}
					
				}
				else 
					flag_ok = true;
				
			}
		}
	}
	
	//FILE *tmpview;
	//tmpview = fopen("places","w");
	olig_total = 0;
	for (int j=0; j < geneno; ++j){
		//fprintf(tmpview, "%d\t%d\n", j, Gene[j]);
		if(Gene[j] != -1){
			fprintf(oligos, "Target sequence:%d \t", j+1);
			for (int k = Gene[j]- OLIG_LENGTH + 1; k <= Gene[j]; ++k)
				fprintf(oligos, "%c", gSequence[k]);
			fprintf(oligos, "\t Length=%d", OLIG_LENGTH); 
			fprintf(oligos, "\t Tm=%.2f", tmarray[Gene[j]-gseqOff]); 
			fprintf(oligos, "\t Distance(from 5'):%d", Gene[j]-gseqOff-geneLoc[j]);
			fprintf(oligos, "\t Distance(from 3'):%d", geneLoc[j+1]-(Gene[j]-gseqOff));
			fprintf(oligos, "\n");
			++olig_total;
		}
	}	
	//fclose(tmpview);
	int effectUnique = 0;
	for (int i=0; i< geneno; ++i)
		effectUnique += effect_eval[i];
	
	//cout << "Number of oligos after phase 1 are: " << olig_total_p1 << " for "  << geneno  << " genes.\n";
	//cout << "Number of oligos affected by phase 2 are: " << effectUnique << " for "  << geneno  << " genes.\n";
	
	unsigned long sumPos = 0;
	for (int i=0; i < PRIME_IH; ++i) {
		sumPos += hashTable1[i].size + hashTable2[i].size + 
		hashTable3[i].size + hashTable4[i].size +
		hashTable5[i].size + hashTable6[i].size +
		hashTable7[i].size + hashTable8[i].size ;
	}
	//cout << "sumPos: " << sumPos << endl;
	memoryUsed = sumPos * sizeof (int) + 8 * PRIME_IH * sizeof (entry);
	
	free(effect_eval);//delete [] effect_eval;
	free(Gene);//delete [] Gene;
	
	free(hashTable1->pos);//delete [] hashTable1->pos;
	free(hashTable1);//delete [] hashTable1;
	
	free(hashTable2->pos);//delete [] hashTable2->pos;
	free(hashTable2);//delete [] hashTable2;
	
	free(hashTable3->pos);//delete [] hashTable3->pos;
	free(hashTable3);//delete [] hashTable3;
	
	free(hashTable4->pos);//delete [] hashTable4->pos;
	free(hashTable4);//delete [] hashTable4;
	
	free(hashTable5->pos);//delete [] hashTable5->pos;
	free(hashTable5);//delete [] hashTable5;
	
	free(hashTable6->pos);//delete [] hashTable6->pos;
	free(hashTable6);//delete [] hashTable6;
	
	free(hashTable7->pos);//delete [] hashTable7->pos;
	free(hashTable7);//delete [] hashTable7;
	
	free(hashTable8->pos);//delete [] hashTable8->pos;
	free(hashTable8);//delete [] hashTable8;
	
	return memoryUsed;
}

void pre_Process(){
	for (int i=0; i < 65536;++i)
		preNOR[i]=0;
	for (int i=0; i < 65536;++i){
		uint16_t tmp=i;
		for (int pass=0; pass<8;++pass){
			if ((tmp & 3) == 3) 
				++preNOR[i];
			tmp >>= 2;
		}
	}
	
	/*for (int i=65000; i < 65536;++i)	
	 printf("%X\t%d\n",i,preNOR[i]);
	 
	 exit(0);*/
}

unsigned int pickPrime(){
	if(gSize <  30000000)
		return (34040383);
	else if (gSize < 113493747)
		return (113493747);
	const unsigned long primeTable[] = {
		1114523, 1180043, 1245227, 1310759, 1376447, 1442087, 1507379,
		1573667, 1638899, 1704023, 1769627, 1835027, 1900667, 1966127,
		2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767,
		3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143,
		3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239,
		5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639,
		7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767,
		11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543,
		14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227,
		20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047,
		28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983,
		39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067,
		54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567,
		75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439,
		104858387, 109052183, 113246699, 117440699, 121635467, 125829239,
		130023683, 142606379, 150994979, 159383759, 167772239, 176160779,
		184549559, 192938003, 201327359, 209715719, 218104427, 226493747,
		234882239, 243269639, 251659139, 260047367, 285215507, 301989959,
		318767927, 335544323, 352321643, 369100463, 385876703, 402654059,
		419432243, 436208447, 452986103, 469762067, 486539519, 503316623,
		520094747, 570425399, 603979919, 637534763, 671089283, 704643287,
		738198347, 771752363, 805307963, 838861103, 872415239, 905971007,
		939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679,
		1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119,
		1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967,
		2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539,
		2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823,
		3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783,
		5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119,
		6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599,
		8321499203, 9126806147, 9663676523, 10200548819, 10737418883,
		11274289319, 11811160139, 12348031523, 12884902223, 13421772839,
		13958645543, 14495515943, 15032386163, 15569257247, 16106127887,
		16642998803, 18253612127, 19327353083, 20401094843, 21474837719,
		22548578579, 23622320927, 24696062387, 25769803799, 26843546243,
		27917287907, 28991030759, 30064772327, 31138513067, 32212254947,
		33285996803, 36507222923, 38654706323, 40802189423, 42949673423,
		45097157927, 47244640319, 49392124247, 51539607599, 53687092307,
		55834576979, 57982058579, 60129542339, 62277026327, 64424509847,
		66571993199, 73014444299, 77309412407, 81604379243, 85899346727,
		90194314103, 94489281203, 98784255863, 103079215439, 107374183703,
		111669150239, 115964117999, 120259085183, 124554051983, 128849019059,
		133143986399, 146028888179, 154618823603, 163208757527, 171798693719,
		180388628579, 188978561207, 197568495647, 206158430447, 214748365067,
		223338303719, 231928234787, 240518168603, 249108103547, 257698038539,
		266287975727, 292057776239, 309237645803, 326417515547, 343597385507,
		360777253763, 377957124803, 395136991499, 412316861267, 429496730879,
		446676599987, 463856468987, 481036337207, 498216206387, 515396078039,
		532575944723, 584115552323, 618475290887, 652835029643, 687194768879,
		721554506879, 755914244627, 790273985219, 824633721383, 858993459587,
		893353198763, 927712936643, 962072674643, 996432414899, 1030792152539,
		1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587,
		1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839,
		1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323,
		2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239,
		2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983,
		3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483,
		4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979,
		5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003,
		7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659,
		8521215117407, 9345848837267, 9895604651243, 10445360463947,
		10995116279639, 11544872100683, 12094627906847, 12644383722779,
		13194139536659, 13743895350023, 14293651161443, 14843406975659,
		15393162789503, 15942918604343, 16492674420863, 17042430234443,
		18691697672867, 19791209300867, 20890720927823, 21990232555703,
		23089744183799, 24189255814847, 25288767440099, 26388279068903,
		27487790694887, 28587302323787, 29686813951463, 30786325577867,
		31885837205567, 32985348833687, 34084860462083, 37383395344739,
		39582418600883, 41781441856823, 43980465111383, 46179488367203,
		48378511622303, 50577534878987, 52776558134423, 54975581392583,
		57174604644503, 59373627900407, 61572651156383, 63771674412287,
		65970697666967, 68169720924167, 74766790688867, 79164837200927,
		83562883712027, 87960930223163, 92358976733483, 96757023247427,
		101155069756823, 105553116266999, 109951162779203, 114349209290003,
		118747255800179, 123145302311783, 127543348823027, 131941395333479,
		136339441846019, 149533581378263, 158329674402959, 167125767424739,
		175921860444599, 184717953466703, 193514046490343, 202310139514283,
		211106232536699, 219902325558107, 228698418578879, 237494511600287,
		246290604623279, 255086697645023, 263882790666959, 272678883689987,
		299067162755363, 316659348799919, 334251534845303, 351843720890723,
		369435906934019, 387028092977819, 404620279022447, 422212465067447,
		439804651111103, 457396837157483, 474989023199423, 492581209246163,
		510173395291199, 527765581341227, 545357767379483, 598134325510343,
		633318697599023, 668503069688723, 703687441776707, 738871813866287,
		774056185954967, 809240558043419, 844424930134187, 879609302222207,
		914793674313899, 949978046398607, 985162418489267, 1020346790579903,
		1055531162666507, 1090715534754863
	};
	int i=0;
	while (primeTable[i] < 2*gSize) ++i;
	return primeTable[i];
}

int main (int argc, char * const argv[]) {
	cout << "------------------------------BEGIN------------------------------------\n";
	print_Time();

	time_t seconds, secondf, p1s,p1f,p2s,p2f;	
	seconds = time (NULL);
	clock_t startTotal, finishTotal, startP1, finishP1, startP2, finishP2, start_iden, finish_iden;
	
	
		
	unsigned long memoryUsed = 0, peakUsed = 0;
	
	startTotal = clock();
	startP1 = clock();
	
	pre_Process();
	
	FILE *gf, *oligos;
	unsigned char *tempSequence, *gSequence;
	uint8_t *gCoded;
	uint8_t *simarray;
	int tempSize;
	
	const char *gfname = argv[1];//"/Users/mohamah/Desktop/oligo/data_fasta/ecoli.fsa";// 
	/* Open genome file for reading. */
	if((gf = fopen(gfname, "rb")) == NULL) {
		cout << "Error! cannot open file: "<< gfname << endl;
		exit(EXIT_FAILURE);
	}
	
	const char *oligosname = argv[2];//"/Users/mohamah/Desktop/out_test_1";// 
	if (argv[2][0]=='-'){
		cout << "Error! no output file. \n";
		exit(EXIT_FAILURE);
	}
	
	for (int i=3; i < argc; i+=2){		
		string param = argv[i];
		if (param == "-length") OLIG_LENGTH = atoi(argv[i+1]);
		else if (param == "-seqSim") seqSimilarity = atoi(argv[i+1]); 
		else if (param == "-maxMatch") consecMatch = atoi(argv[i+1]); 
		else if (param == "-maxGC") GC_CONT_MAX = atoi(argv[i+1]); 
		else if (param == "-minGC") GC_CONT_MIN = atoi(argv[i+1]); 
		else if (param == "-dimerLen") dimerLen = atoi(argv[i+1]); 
		else if (param == "-dimerSim") dimerStg = atoi(argv[i+1]); 
		else if (param == "-hairpinLen") hpStem = atoi(argv[i+1]); 				
		else if (param == "-minhpLoop") hpminLoop = atoi(argv[i+1]); 
		else if (param == "-maxhpLoop") hpmaxLoop = atoi(argv[i+1]); 				
		else if (param == "-minTm") tmmin = atoi(argv[i+1]); 
		else if (param == "-maxTm") tmmax = atoi(argv[i+1]);
		else if (param == "-rangeTm") tminterval = atoi(argv[i+1]); 
		else if (param == "-oligCon") oligc = (double)(atoi(argv[i+1]))/1000000000.0; 
		else if (param == "-saltCon") saltc = (double)(atoi(argv[i+1]))/1000.0;	
		else if (param == "-secStr") secStr = true;	
		else {
			cout << "Error! in parameters!\n";
			exit(EXIT_FAILURE);
		}
	}
	
	
	oligos = fopen(oligosname, "w");
	
		
	/* Get the genome file size. */
	if(fseek(gf, 0, SEEK_END) == 0) {
		tempSize = ftell(gf);
		rewind(gf);
		if(tempSize < 0) { 
			cout << "Error! Cannot ftell: " << gfname;
			perror(NULL);
			exit(EXIT_FAILURE);
		}
	} else {
		cout << "Error! Cannot fseek: " << gfname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	/* Allocate memory to genome sequence array. */
	tempSequence = (unsigned char *) malloc(tempSize * sizeof(unsigned char));
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	gSequence = (unsigned char *) malloc(tempSize * sizeof(unsigned char));
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	simarray = (uint8_t *) malloc (tempSize * sizeof(uint8_t));
	memoryUsed += (tempSize) * sizeof(uint8_t);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (long i=0; i<tempSize; ++i)
		simarray[i] = 0;
	
	if((tempSequence == NULL) || (gSequence == NULL)) {
		cout << "Error! Cannot allocate memory.";
		exit(EXIT_FAILURE);
	}
	
	/* Read tempSize bytes of data. */
	if(fread(tempSequence, sizeof(unsigned char), tempSize, gf) != tempSize) {
		cout << "Cannot read from: "<< gfname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	fclose(gf);
	
	
	int gene_num = 0;
	int locSize = 8;
	int blanks=0;
	geneno=0;
	int errorno = 0;
	int i = 0;
	int *geneLoc = (int *) malloc(locSize * sizeof(int));
	while (i < tempSize) {
		switch (tempSequence[i]) {
			case '>':
				if (geneno >= locSize){
					locSize *= 2;
					geneLoc = (int *) realloc(geneLoc, locSize * sizeof(int));
				}
				geneLoc[geneno] = gene_num;
				while (tempSequence[i] != '\n' && tempSequence[i] != '\r') ++i;
				++geneno;
				break;
			case 'a':
				gSequence[gene_num] = 'A'; //goutf << 'A';
				++gene_num;
				break;
			case 'A':
				gSequence[gene_num] = 'A'; //goutf << 'A';
				++gene_num;
				break;
			case 'c':
				gSequence[gene_num] = 'C'; //goutf << 'C';
				++gene_num;
				break;
			case 'C':
				gSequence[gene_num] = 'C'; //goutf << 'C';
				++gene_num;
				break;
			case 'g':
				gSequence[gene_num] = 'G'; //goutf << 'G';
				++gene_num;
				break;
			case 'G':
				gSequence[gene_num] = 'G'; //goutf << 'G';
				++gene_num;
				break;
			case 't':
				gSequence[gene_num] = 'T'; //goutf << 'T';
				++gene_num;
				break;
			case 'T':
				gSequence[gene_num] = 'T'; //goutf << 'T';
				++gene_num;
				break;
			case '\n':
				++blanks;
				break;
			case '\r':
				++blanks;
				break;
			default:
				gSequence[gene_num] = 'A'; 
				int ranger = (gene_num + OLIG_LENGTH <= tempSize)? gene_num + OLIG_LENGTH : tempSize;
				int rangel = (gene_num - OLIG_LENGTH + 1 >= 0)? gene_num - OLIG_LENGTH + 1 : 0;
				for (int j=gene_num; j< ranger; ++j)
					simarray[j] = 1;
				for (int j=rangel; j<= gene_num; ++j)
					simarray[j] = 1;
				++errorno;
				++gene_num;
				break;
		}
		++i;
	}
	free(tempSequence);//delete [] tempSequence;
	memoryUsed -= (tempSize) * sizeof(unsigned char);
	
	gSize = gene_num;
	if (gSize < tempSize){
		gSequence = (unsigned char *) realloc(gSequence, gSize * sizeof(unsigned char));
		simarray = (uint8_t *) realloc(simarray, gSize * sizeof(uint8_t));	
	}
	
	memoryUsed -= (2*tempSize) * sizeof(unsigned char);
	memoryUsed += (2*gSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	if (geneno < locSize)
		geneLoc = (int *) realloc(geneLoc, geneno * sizeof(int));
	
	memoryUsed += (geneno) * sizeof(int);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	
	gseqOff = gSize % 4;
	
	if(gseqOff != 0){
		for(int i=1; i < geneno; ++i)
			geneLoc[i] -= gseqOff;
		for(int i=0; i < gSize - gseqOff; ++i)
			simarray[i] = simarray[i+gseqOff];
	}
	
	/*          Start Encoding the genome 
	 *************************************************** */
	nbytes = gSize/4;
	gCoded = (uint8_t *) malloc (nbytes * sizeof(uint8_t));
	memoryUsed += nbytes * sizeof(uint8_t);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	uint8_t temp = 0;	
	int stpoint = gSize % 4;
	for(int i = stpoint; i < gSize; i += 4) {
		for(int j = 0; j < 4; j++) 
			switch (gSequence[j+i]) {
				case 'A':
					temp = (uint8_t) (temp << 2 | 0);
					break;
				case 'C':
					temp = (uint8_t) (temp << 2 | 1);
					break;
				case 'G':
					temp = (uint8_t) (temp << 2 | 2);
					break;
				case 'T':
					temp = (uint8_t) (temp << 2 | 3);
					break;
				default:
					cout << "Error! " << gSequence[j+i] << "is not a valid nucleotide.\n";
					temp = (uint8_t) (temp << 2 | 0);					
			}	
		gCoded[(i-stpoint)/4] = temp;
	}
	/*          END Encoding the genome 
	 *************************************************** */
	gSize -= gseqOff;
	simarray = (uint8_t *) realloc(simarray, gSize * sizeof(uint8_t));
	////////	////////	////////	//////// END create genome	
	
	cout << "\nGenome size is: " << (gSize + gseqOff) << ", and number of target sequences are: "<< geneno <<"\n\n";
	
	PRIME_MM = pickPrime();
	
	for (int i = 0; i < OLIG_LENGTH - 1; ++i) 
		simarray[i] = 1;
	
	// Kane Identity Check 	
	memoryUsed += PRIME_MM * (sizeof(uint64_t) + sizeof(int));
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	//cout << "Mem. used for identity_check: " << ((memoryUsed/1024)/1024) << endl;
	
	p1s = time(NULL);	
	//start_iden = clock();
	identity_Check(gCoded, simarray);	
	//finish_iden = clock();
	//printf("identity_check time: %.1f sec\n\n", (double)((finish_iden - start_iden)) / (double)CLOCKS_PER_SEC);
	
	memoryUsed -= PRIME_MM * (sizeof(uint64_t) + sizeof(int));
	
	
	// Kane Similarity Check
	const char * seed_Phase1 [] = {
		"111010011010111",
		"1100101001000100011011",
		"1111000100010000110101",
		"1100110000101000101101",
		"11100100010010010000111",
		"11010010000000101001000111",
		"1101000100001000001000100111",
		"1101010000010000010000011011"
	};
	
 	MAX_NUM_THREAD = omp_get_max_threads();	
	//printf("Max. number of threads = %d\n", MAX_NUM_THREAD);
	//printf("Phase 1 is running with seed set: S8W10\n");
	
	if (MAX_NUM_THREAD >= 8)
		omp_set_num_threads(8);
	else
		omp_set_num_threads(MAX_NUM_THREAD);
	
	memoryUsed += MAX_NUM_THREAD*PRIME_FH * (sizeof(uint64_t) + sizeof(int));
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	int tid, chunk=1;
	#pragma omp parallel for shared(gCoded, simarray, seed_Phase1, chunk) private(i, tid) schedule(dynamic,chunk)
	for (i = 0; i < 8; ++i){
		tid = omp_get_thread_num();
		//printf("Thread %d doing section%d\n",tid, i);
		similarity_Check(gCoded, simarray, seed_Phase1[i]);// S8W10L18	
		//printf("Thread %d done section%d.\n",tid, i);
	}		
	
	p1f = time(NULL);
	printf ("Fast homology search done! in %ld sec\n\n", (p1f-p1s));
	//cout << "Peak memory used in MB: " << ((peakUsed/1024)/1024) << endl;
	
	/*int count = 0;
	for (int i = 0; i < gSize; ++i) 
		if (simarray[i] ==1)
			count++;*/
	//cout << "Similar positions in Phase 1: " << count << "\t%" << ( 100 * double(count) / gSize) << "\n\n";
	
	
	memoryUsed -= MAX_NUM_THREAD*PRIME_FH * (sizeof(uint64_t) + sizeof(int));
	
	finishP1 = clock();
	
	//cout << "Mem. used before Phase2: " << ((memoryUsed/1024)/1024) << endl;
	
	
	// G-C Contetnt Check
	gc_Content(gCoded, simarray);
	
	// Tm Check
	p1s = time(NULL);
	float *tmarray = (float *) malloc (gSize*sizeof(float));
	memoryUsed += (gSize) * sizeof(float);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (int i = 0; i < gSize; ++i) 
		tmarray[i] = 0;	
	
	melting_Temp(gCoded, simarray, geneLoc, tmarray);
	
	
	p1f = time(NULL);
	//printf ("Time for Tm: %ld sec\n", (p1f-p1s));
	//memoryUsed -= (gSize) * sizeof(float);
	// Tm Check Done!
	
	startP2 = clock();
	p2s = time(NULL);
	unsigned long memoryEval = 0;
	
	memoryEval = print_Oligo_S8W9(gCoded, simarray, gSequence, geneLoc, oligos, tmarray);
	
	p2f = time(NULL);
	printf ("Intensive homology search done! in %ld sec\n\n", (p2f-p2s));
	finishP2 = clock();
	
	if (memoryUsed + memoryEval > peakUsed) peakUsed = memoryUsed + memoryEval;
	
	
	fclose(oligos);//oligos.close();
	
	finishTotal = clock();	
	
	secondf = time (NULL);
	cout << "Total number of unique oligos are: " << olig_total << " for "  << geneno  << " genes.\n";
	cout << "Peak memory used in MB: " << ((peakUsed/1024)/1024) << endl;
	printf ("Running Time: %ld sec\n\n", (secondf-seconds));	
	free(tmarray);//delete [] tmarray;
	free(gCoded);//delete [] gCoded;
	free(gSequence);//delete [] gSequence;
	free(simarray);	
	free(geneLoc);//delete [] geneLoc;
	cout<<endl;
	print_Time();
	cout << "------------------------------END--------------------------------------\n";
    return 0;
}
