
#include "genome.h"
int badsecStr = 0;
float maxTm=0.0, minTm=100.0;
using namespace std;

struct entry {
	uint64_t key;
	int *pos;
	int size;
};

uint64_t get64bit(int rpos, genome *mygenome){
	int lpos = rpos - 7;
	if (lpos < 0 && rpos >= 0) 
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (int i = lpos; i <= rpos; ++i) 
		genomeChunk = (uint64_t)(genomeChunk << 8 | mygenome->gCoded[i]);
	return genomeChunk;	
}

uint32_t getvarbit(int byteno, int rpos, uint8_t *Coded){
	if (rpos < 0){
		cout << "error!!! in get64bit\n";
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	int lpos = rpos - byteno + 1;
	if (lpos < 0 && rpos >= 0)
		lpos = 0;
	uint32_t genomeChunk = 0;
	for (int i = lpos; i <= rpos; ++i)
		genomeChunk = (uint32_t)(genomeChunk << 8 | Coded[i]);
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
		j = (k + i) % FRIME;
		if (T[j].pos == NULL){
			T[j].key = k;
			T[j].size = POS_SIZE;
			T[j].pos = new int [T[j].size];
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
	}while (i!=FRIME);
	
	T[0].pos = NULL;
	T[0].key = k;
	T[0].size = POS_SIZE;
	T[0].pos = new int [T[0].size];
	T[0].pos[0] = locn;
	for (int y = 1; y < T[0].size; ++y)
		T[0].pos[y] = -1;
	return 0;
}

void create_Hash(entry *hashTable, genome *mygenome, uint64_t seed){
	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	int location = -1;
	
	for (int i=0; i < FRIME; ++i) {
		hashTable[i].key = 0;
		hashTable[i].pos = NULL;
		hashTable[i].size = 0;
	}
	
	window = get64bit(mygenome->nbytes-1, mygenome);
	for (int i = mygenome->nbytes-1; i >= 12; --i) {		
		for (int pass=3; pass >=0 ; --pass) {	
			hashKey = (uint64_t) (window & seed);  
			location = i*4 + pass; 
			insert(hashKey, location, hashTable);
			window >>=  2;			
		}
		//if (i>=8)
		window |= (((uint64_t) mygenome->gCoded[i-8]) << (BLOCK - 8));
	}
}

int search_Hash(uint64_t k, entry *T){
	int i=0, j;
	do {
		j = (int)((k + i) % FRIME);
		if (T[j].key == k)
			return j;
		++i;
	}while (i!=FRIME && T[j].pos!=NULL);
	return -1;
}

int check_50factor(genome *mygenome, int pose,  int posc){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
	chitr = get64bit(cbyte, mygenome); chitl = get64bit(cbyte - 8, mygenome);
	
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
	
	if (kane75 > KANE_SIMILARITY)
		return 1;
	
	return 0;
}

bool second_Check(uint64_t seed, int seedLen, entry *hashTable, genome *mygenome, int pose){
	////////////////////////////////////////////////////////////////Computing Hash value
	uint64_t hashKey1=0, ehitr=0, ehitl=0;
	int w; int pos1;
	int ebyte, eoff;;
	////////////////////////////////////////////////////////////////////Checking Positions
	ebyte = pose/4; eoff = pose%4;
	ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	////////////////////////////////////////////////////////////////////	
	int offpos = 0;
	for (int pass=1; pass <= OLIG_LENGTH - seedLen + 1; ++pass){
		hashKey1 = (uint64_t) (ehitr & seed);  
		pos1 = search_Hash(hashKey1, hashTable);											
		w=0;
		while (w < hashTable[pos1].size && hashTable[pos1].pos[w] !=-1 ){
			if (hashTable[pos1].pos[w] > OLIG_LENGTH - 1 && hashTable[pos1].pos[w] + offpos != pose){
				if (check_50factor(mygenome, pose, hashTable[pos1].pos[w] + offpos) == 0)		
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
			
void kane_75l(genome *mygenome, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome); //ehitl = get64bit(ebyte - 16, mygenome);
	chitr = get64bit(cbyte, mygenome); chitl = get64bit(cbyte - 8, mygenome); //chitl = get64bit(cbyte - 16, mygenome);
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
		uint8_t elbyte = mygenome->gCoded[ebyte-16];
		ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	}
	
	chitr >>= 2*(3-coff);
	chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	if (OLIG_LENGTH > 64 - (3-coff) && cbyte-16 >=0){
		uint8_t clbyte = mygenome->gCoded[cbyte-16];
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
	
	while(kane75 > KANE_SIMILARITY -1 && posc > OLIG_LENGTH - 2){
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
			ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
			chitr = get64bit(cbyte, mygenome); chitl = get64bit(cbyte - 8, mygenome);
			
			ehitr >>= 2*(3-eoff);
			ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			ehitl >>= 2*(3-eoff);
			
			if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
				uint8_t elbyte = mygenome->gCoded[ebyte-16];
				ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
			}
			
			chitr >>= 2*(3-coff);
			chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
			chitl >>= 2*(3-coff);
			
			if (OLIG_LENGTH > 64 - (3-coff) && cbyte-16 >=0){
				uint8_t clbyte = mygenome->gCoded[cbyte-16];
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

void kane_75r(genome *mygenome, int pose,  int posc, uint8_t *simarray){
	int ebyte = pose/4, eoff = pose%4, cbyte = posc/4, coff = posc%4;
	uint64_t ehitl, ehitr, chitl, chitr, nolxnr, nolxnl;
	
	ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
	chitr = get64bit(cbyte, mygenome); chitl = get64bit(cbyte - 8, mygenome);
	
	
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	
	chitr >>= 2*(3-coff);
	chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
	chitl >>= 2*(3-coff);
	
	//bitset<64> b1=ehitl,b2=ehitr, d1=chitl, d2=chitr; 
	//cout << b1 << "\t" << b2 << "\n" << d1 << "\t" << d2 << "\n";
	//bitset<100> ehit, chit;
	//ehit = ehitl; ehit <<= BLOCK; ehit |= ehitr;
	//chit = chitl; chit <<= BLOCK; chit |= chitr;
	
	
	//bitset<100> nole,nolc, nolxnor;
	//nole = get100bit(pose, mygenome);
	//nolc = get100bit(posc, mygenome);
	//nolxnor = ~(nole ^ nolc);
	
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
	int shiftno = 0;
	//kane75 = nolxnor.count();
	//cout << nole << "\n" << nolc << "\n" << nolxnor << "\n" << kane75 << "\n\n";
	
	while(kane75 > KANE_SIMILARITY -1 && pose <= mygenome->gSize - 1 ){
		simarray[posc] = simarray[pose] = 1;
		++posc; ++pose;
		
		/*if (shiftno <= 10)
		 {
		 ehitr >>= 2;
		 ehitr |= (ehitl & (uint64_t) (3)) << (BLOCK - 2);
		 ehitl >>= 2;
		 
		 chitr >>= 2;
		 chitr |= (chitl & (uint64_t) (3)) << (BLOCK - 2);
		 chitl >>= 2;
		 
		 ++shiftno;
		 }*/
		//else {
		////////////////////////////////// BEGIN: Reading next 50-mer
		ebyte = pose/4; eoff = pose%4; cbyte = posc/4; coff = posc%4;
		ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
		chitr = get64bit(cbyte, mygenome); chitl = get64bit(cbyte - 8, mygenome);
		
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);
		
		chitr >>= 2*(3-coff);
		chitr |= (chitl & (uint64_t) ( (1 << 2*(3-coff)) - 1)) << (BLOCK - 2*(3-coff));
		chitl >>= 2*(3-coff);		
		////////////////////////////////// END: Reading next 50-mer		
		shiftno = 0;
		//}
		
		
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

void kane_15(genome *mygenome, int pose,  int posc, uint8_t *simarray, uint8_t *elimTable){
	int rangelc = posc, rangerc = posc + OLIG_LENGTH - KANE_IDENTITY < mygenome->gSize? posc + OLIG_LENGTH - KANE_IDENTITY : mygenome->gSize -1;
	for (int i=rangelc; i <= rangerc ; ++i)
		simarray[i] = 1;
	elimTable[posc] = 1;
	if(!elimTable[pose]){
		int rangele = pose, rangere = pose + OLIG_LENGTH - KANE_IDENTITY < mygenome->gSize? pose + OLIG_LENGTH - KANE_IDENTITY : mygenome->gSize -1;
		for (int i=rangele; i <= rangere ; ++i)
			simarray[i] = 1;
		elimTable[pose] = 1;
	}
}

void gc_Content(genome *mygenome, uint8_t *simarray){
	int byte = mygenome->nbytes -1 , off, pos = mygenome->gSize -1;
	uint64_t halfl, halfr, nolxnr, nolxnl; 
	halfr = get64bit(byte, mygenome); halfl = get64bit(byte - 8, mygenome);
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
			halfr = get64bit(byte, mygenome); halfl = get64bit(byte - 8, mygenome);
			
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

void melting_Temp(genome *mygenome, uint8_t *simarray, int* geneLoc, float *tmarray){
	
	int byte, off;
	uint64_t halfl2=0, halfl=0, halfr=0, nolxnr=0, nolxnl=0, nolxnl2=0; 		
	
	const float enthalpy[] = {-7.9,-8.4,-7.8,-7.2,-8.5,-8.0,-10.6,-7.8,-8.2,-9.8,-8.0,-8.4,-7.2,-8.2,-8.5,-7.9};
	const float  entropy[] = {-22.2,-22.4,-21.0,-20.4,-22.7,-19.9,-27.2,-21.0,-22.2,-24.4,-19.9,-22.4,-21.3,-22.2,-22.7,-22.2};
	
	float enth, entr, Tm, oligc = 0.000001, saltc = 0.075;	
	int length = OLIG_LENGTH;
	
	//while (pos > OLIG_LENGTH - 2) {
	for (int qos=0; qos < geneno; ++qos){
		for(int pos = geneLoc[qos] + OLIG_LENGTH - 1; pos < geneLoc[qos+1]; ++pos){
			if (simarray[pos] == 0){
				byte = pos/4; off = pos%4; 
				halfr = get64bit(byte, mygenome); halfl = get64bit(byte - 8, mygenome); 
				if (byte - 16 >= 0)
					halfl2 = get64bit(byte - 16, mygenome);
				
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
				
				
				/*	enth=enth*4.18*1000;
				 entr=entr*4.18;
				 Tm = enth/(entr+0.368*(length-1)*log(saltc)*4.18+1.9865*4.18*log(oligc/4))-273.15;*/
				
				/// Tm checking and elimination.
				
				////////////////////////////////// END: Computing Tm
				tmarray[pos] = Tm;
			}
		}
	}
	
	int count = 0;	double sumTm = 0.0;
	//for (int i = OLIG_LENGTH - 1; i < mygenome->gSize; ++i)
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j] == 0){
				++count;
				sumTm += tmarray[j]; 
				if (tmarray[j]>maxTm) maxTm=tmarray[j];
				if (tmarray[j]<minTm) minTm=tmarray[j];
			}		
	float avtm = sumTm / count;
	//printf("max_Tm: %.3f \t min_Tm: %.3f \t aveg_Tm: %.3f\n",maxTm,minTm,avtm);
	
	int uppTm=100*maxTm, lowTm=100*minTm; 
	//printf("uppTm: %d \t lowTm: %.d \n",uppTm,lowTm);	
	// Computing Upper Bound Considering the Tm Range 
	bool **tmFreq = new bool * [geneno];
	for (int i=0; i < geneno; ++i)
		tmFreq[i] = new bool [uppTm+1];
	int *sumPm = new int [geneno];
	
	for (int i=0; i < geneno; ++i){
		for (int j=0; j <= uppTm; ++j)
			tmFreq[i][j]=false;
		sumPm[i] = 0;
	}
	
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j] == 0)
				tmFreq[i][(int)(100*tmarray[j])]=true;
	
	int tmRange = 1000, upperBound, maxUbound = 0, maxLindex = -1, maxRindex = -1;
	
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
	//printf("Max. Upper Bound in the range [%.2f, %.2f) is: %d\n\n", maxLindex/100.0, (maxRindex+1)/100.0, maxUbound);
	
	for (int i=0; i < geneno; ++i)
	 	for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j] == 0) 
				if(tmarray[j] >= 82.56   || tmarray[j] < 64.36 )
					simarray[j] = 1;
	
	for(int i=0 ; i < geneno ; ++i)
		delete tmFreq[i];
	delete [] tmFreq;
	delete [] sumPm;
}

int hash_Search(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME;
		if (T[j] == k)
			return j;
		++i;
	}while (i!=PRIME && T[j]!=s+1);
	return -1;
}

int hash_Insert(uint64_t k, uint64_t *T, uint64_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME;
		if (T[j] == s+1){
			T[j] = k;
			return j;
		}
		else
			++i;
	}while (i!=PRIME);
	return 0;
}			

void similarity_Check(genome *mygenome, uint8_t *simarray, const char *seedSeq){
	uint64_t seed=0; 
	seed = get64seed(seedSeq);	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	uint64_t *hashTable = new uint64_t[PRIME];
	int *posTable = new int[PRIME];
	for (int i=0; i < PRIME; ++i) {
		hashTable[i] = seed + 1;
		posTable[i] = -1;
	}
	window = get64bit(mygenome->nbytes-1, mygenome);
	for (int i = mygenome->nbytes-1; i >= 8; --i) {				
		for (int pass=3; pass >=0 ; --pass) {	
			hashKey = (uint64_t) (window & seed);  
			int searchRes = hash_Search(hashKey, hashTable, seed);
			if (( searchRes != -1)){// && gene_find(posTable[searchRes], FindGene)!=gene_find(i*4 + pass, FindGene) ){ 
				if (simarray[i*4+pass] == 0 || simarray[posTable[searchRes]] == 0){
					//++hitTotalP1;
					kane_75l(mygenome, posTable[searchRes],  i*4 + pass, simarray);
					kane_75r(mygenome, posTable[searchRes],  i*4 + pass, simarray);
					//if(simarray[i*4+pass] == 1 && simarray[posTable[searchRes]] == 1)
					//	++hitGoodP1;
				}				
				posTable[searchRes] = i*4 + pass;
			}
			else
				posTable[hash_Insert(hashKey, hashTable, seed)] = i*4 + pass;
			
			window >>=  2;			
		}
		window |= (((uint64_t) mygenome->gCoded[i-8]) << (BLOCK - 8));
	}
	//printf("total hit: %ld\t good hit: %ld\n", hitTotalP1, hitGoodP1);
	delete [] hashTable;
	delete [] posTable;
}

int hash_Search_IDEN(uint32_t k, uint32_t *T, uint32_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_IDEN;
		if (T[j] == k)
			return j;
		++i;
	}while (i!=PRIME_IDEN && T[j]!=s+1);
	return -1;
}

int hash_Insert_IDEN(uint32_t k, uint32_t *T, uint32_t s){
	int i=0, j;
	do {
		j = (k + i) % PRIME_IDEN;
		if (T[j] == s+1){
			T[j] = k;
			return j;
		}
		else
			++i;
	}while (i!=PRIME_IDEN);
	return 0;
}

void identity_Check(genome *mygenome, uint8_t *simarray){
	
	//Begin Coding "111111111111111"
	uint32_t seed=0;
	for(int j = 0; j < KANE_IDENTITY; ++j)
		seed = (uint32_t) (seed << 2 | 3);
	
	//Begin Hashing
	uint32_t window = 0, hashKey=0;
	uint8_t nexwin=0;
	uint32_t *hashTable = new uint32_t[PRIME_IDEN];
	int *posTable = new int[PRIME_IDEN];
	for (int i=0; i < PRIME_IDEN; ++i) {
		hashTable[i] = seed + 1;
		posTable[i] = -1;
	}
	uint8_t *elimTable = new uint8_t [mygenome->gSize];
	for (int i=0; i < mygenome->gSize; ++i) 
		elimTable[i] = 0;		

	int SEED_BYTE = KANE_IDENTITY % 4 == 0 ?  KANE_IDENTITY / 4 : KANE_IDENTITY / 4  + 1;
	int rpos = mygenome->nbytes-1;
	window = getvarbit(SEED_BYTE, rpos, mygenome->gCoded);
	rpos -= SEED_BYTE;
	int nexwinSize =   SEED_BYTE * 4 - KANE_IDENTITY;
	if (nexwinSize == 0){
		nexwin = mygenome->gCoded[rpos];
		--rpos;
		nexwinSize = 4;
	}
	else
		nexwin = (window >> 2*KANE_IDENTITY) & ((1 << 2*nexwinSize)-1);
	
	hashKey = (uint32_t) (window & seed);
	
	
	for (int location = mygenome->gSize-1; location >= KANE_IDENTITY -1; --location) {
		int searchRes = hash_Search_IDEN(hashKey, hashTable, seed);
		if (( searchRes != -1))
			kane_15(mygenome, posTable[searchRes], location, simarray, elimTable);
		else
			posTable[hash_Insert_IDEN(hashKey, hashTable, seed)] = location;
		
		hashKey = (hashKey >> 2) | ( ((uint32_t)(nexwin & 3)) << (KANE_IDENTITY*2 - 2) );
		nexwin >>= 2;  --nexwinSize;
		if (nexwinSize == 0){
			nexwin = mygenome->gCoded[rpos];
			--rpos;
			nexwinSize = 4;
		}
	}
	
	delete [] hashTable;
	delete [] posTable;
	delete [] elimTable;
}

bool secondary_Structure(genome *mygenome, int pose, const char *seedSeq){
	uint64_t seed=0;
	seed = get64seed(seedSeq);
	int ebyte = pose/4, eoff = pose%4;
	uint64_t ehitl, ehitr, dcwindl=0, rcwindl=0, dcwindr=0, rcwindr=0, tpwind=0;	
	
	ehitr = get64bit(ebyte, mygenome); ehitl = get64bit(ebyte - 8, mygenome);
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);// Remember ehitl is not a correct 32Mer (left side). it needs eoff char from ehitl2
	
	//printf("KMER: %16llX \t %16llX \n", ehitl, ehitr);
	
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
	
	rcwindr |= (((uint64_t)(rcwindl & ((1 << 2*(32-18))-1))) << 2*18);
	rcwindl >>= 2*(32-18);
	
	dcwindl &= ((((uint64_t)1) << 2*(18))-1);
	
	//printf("dcwind: %16llX \t %16llX \n", dcwindl, dcwindr);
	//printf("rcwind: %16llX \t %16llX \n", rcwindl, rcwindr);
	
	
	/*for(int location = pose; location >= pose - OLIG_LENGTH + 1; --location){
		insert_ss(kmer, ssTable);
		dcwind = (dcwind >> 2) | ( ((uint64_t)(ehitl & 3)) << (BLOCK - 2));
		ehitl >>= 2;  --nexwinSize;
		if (nexwinSize == 0){
			nexwin = gCoded[rpos];
			--rpos;
			nexwinSize = 4;
		}
	}*/
	
	
	
	
	int simStr;
	uint64_t nolxnr, nolxnl;
	////////////////////////////////// BEGIN: Checking for 75% similarity
	simStr = 0;
	nolxnr = ~(dcwindr ^ rcwindr);
	nolxnl = ~(dcwindl ^ rcwindl);
	for (int i=1; i<=32; ++i) {
		if ((nolxnr & 3) == 3) 
			simStr+=2;
		nolxnr >>= 2;
	}
	for (int i=1; i<=OLIG_LENGTH-32; ++i) {
		if ((nolxnl & 3) == 3) 
			simStr+=2;
		nolxnl >>= 2;
	}
	
	if(simStr < 40){
		return true;
	}
	++badsecStr;
	return false;
}

unsigned long print_Oligo_S8W9(genome *mygenome, uint8_t *simarray, unsigned char *gSequence, int *geneLoc, FILE* oligos, int gseqOff, float *tmarray){
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
	
	entry *hashTable1 = new entry [FRIME];
   	uint64_t seedSeq1 = get64seed(seed_Phase2[0]);
	uint64_t seedLen1 = strlen(seed_Phase2[0]);
	
    entry *hashTable2 = new entry [FRIME];
	uint64_t seedSeq2 = get64seed(seed_Phase2[1]);
	uint64_t seedLen2 = strlen(seed_Phase2[1]);

	entry *hashTable3 = new entry [FRIME];
	uint64_t seedSeq3 = get64seed(seed_Phase2[2]);
	uint64_t seedLen3 = strlen(seed_Phase2[2]);

	entry *hashTable4 = new entry [FRIME];
	uint64_t seedSeq4 = get64seed(seed_Phase2[3]);
	uint64_t seedLen4 = strlen(seed_Phase2[3]);

	entry *hashTable5 = new entry [FRIME];
	uint64_t seedSeq5 = get64seed(seed_Phase2[4]);
	uint64_t seedLen5 = strlen(seed_Phase2[4]);

	entry *hashTable6 = new entry [FRIME];
	uint64_t seedSeq6 = get64seed(seed_Phase2[5]);
	uint64_t seedLen6 = strlen(seed_Phase2[5]);

	entry *hashTable7 = new entry [FRIME];
	uint64_t seedSeq7 = get64seed(seed_Phase2[6]);
	uint64_t seedLen7 = strlen(seed_Phase2[6]);

	entry *hashTable8 = new entry [FRIME];
	uint64_t seedSeq8 = get64seed(seed_Phase2[7]);
	uint64_t seedLen8 = strlen(seed_Phase2[7]);

	int nthreads, tid;
	if (MAX_NUM_THREAD >= 8)	
		omp_set_num_threads(8);
	else 
		omp_set_num_threads(MAX_NUM_THREAD);
	
	
	#pragma omp parallel shared(mygenome, nthreads) private(tid)
	{		
		tid = omp_get_thread_num();
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			//printf("Number of threads = %d\n", nthreads);
		}
		//printf("Thread %d starting...\n",tid);
		
		#pragma omp sections
		{
			#pragma omp section
			{
				//printf("Thread %d doing section 1\n",tid);
				create_Hash(hashTable1, mygenome, seedSeq1);
				//printf("Thread %d done section 1.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 2\n",tid);
				create_Hash(hashTable2, mygenome, seedSeq2);
				//printf("Thread %d done section 2.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 3\n",tid);
				create_Hash(hashTable3, mygenome, seedSeq3);
				//printf("Thread %d done section 3.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 4\n",tid);
				create_Hash(hashTable4, mygenome, seedSeq4);
				//printf("Thread %d done section 4.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 5\n",tid);
				create_Hash(hashTable5, mygenome, seedSeq5);
				//printf("Thread %d done section 5.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 6\n",tid);
				create_Hash(hashTable6, mygenome, seedSeq6);
				//printf("Thread %d done section 6.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 7\n",tid);
				create_Hash(hashTable7, mygenome, seedSeq7);
				//printf("Thread %d done section 7.\n",tid); 
			}
			#pragma omp section
			{
				//printf("Thread %d doing section 8\n",tid);
				create_Hash(hashTable8, mygenome, seedSeq8);
				//printf("Thread %d done section 8.\n",tid); 
			}
		}
	}	
	
	int count = 0;
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if(simarray[j] == 0) ++count;
	//printf("Left before extensive check: %d.\n", count);
	
	int olig_total = 0;
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if(simarray[j] == 0){++olig_total;break;}
	//cout << "Total number of unique oligos before extensive check are: " << olig_total << " for "  << geneno  << " genes.\n";

	int *Gene = new int [geneno];
	memoryUsed += geneno * sizeof (int);

	for (int i=0; i< geneno; ++i)
		Gene[i] = -1;
	
	//printf("Max. number of threads = %d\n", MAX_NUM_THREAD);
	int chunk = 10, i;
	omp_set_num_threads(MAX_NUM_THREAD);	
    #pragma omp parallel for shared(mygenome,chunk) private(i) schedule(static,chunk)
	for (i=0; i < geneno; ++i){
		bool flag_ok;
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j]==0){				
						flag_ok = second_Check(seedSeq1, seedLen1, hashTable1, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq2, seedLen2, hashTable2, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq3, seedLen3, hashTable3, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq4, seedLen4, hashTable4, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq5, seedLen5, hashTable5, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq6, seedLen6, hashTable6, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq7, seedLen7, hashTable7, mygenome, j);
					if (flag_ok)
						flag_ok = second_Check(seedSeq8, seedLen8, hashTable8, mygenome, j);
					
					if (flag_ok)
						Gene[i] = j+gseqOff;
					else
						simarray[j]=1;
			}
	}

	count = 0;
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if(simarray[j] == 0) ++count;
	//printf("Left after extensive check: %d.\n", count);
	
	olig_total = 0;
	for (int j=0; j < geneno; ++j)
		if(Gene[j] != -1)
			++olig_total;
	//cout << "Total number of unique oligos after extensive check are: " << olig_total << " for "  << geneno  << " genes.\n";

 

	// BEGIN Computing Upper Bound Considering the Tm Range 
	int uppTm=100*maxTm, lowTm=100*minTm; 
	bool **tmFreq = new bool * [geneno];
	for (int i=0; i < geneno; ++i)
		tmFreq[i] = new bool [uppTm+1];
	int *sumPm = new int [geneno];
	
	for (int i=0; i < geneno; ++i){
		for (int j=0; j <= uppTm; ++j)
			tmFreq[i][j]=false;
		sumPm[i] = 0;
	}
	
	for (int i=0; i < geneno; ++i)
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j] == 0)
				tmFreq[i][(int)(100*tmarray[j])]=true;
	
	int tmRange = 1000, upperBound, maxUbound = 0, maxLindex = -1, maxRindex = -1;
	
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
	fprintf(oligos, "Estimated maximum number of oligos is: %d\n", maxUbound);
	
	for (int i=0; i < geneno; ++i)
	 	for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)
			if (simarray[j] == 0) 
				if(tmarray[j] >=  1.0*maxRindex/100.0+0.01  || tmarray[j] < 1.0*maxLindex/100.0 )
					simarray[j] = 1;
	
	// END Computing Upper Bound Considering the Tm Range 
	
	for (int i=0; i < geneno; ++i){
		int seenZero, bestZero, seenPlace, bestPlace, canPlace, selPlace;
		seenZero = 0, bestZero = 0, seenPlace = -1, bestPlace = -1, canPlace=-1, selPlace=-1;
		for(int j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j){
			if (simarray[j] == 0){
				++seenZero;
				seenPlace = j;
				if (seenZero > bestZero) {bestZero = seenZero; bestPlace = seenPlace;}
			}
			else 
				seenZero = 0;
		}
		if (bestPlace != -1){				
			canPlace = bestPlace - bestZero/2;
			Gene[i] = canPlace+gseqOff;
		}		
		else 
			Gene[i] = -1;
	}

	olig_total = 0;
	/*for (int j=0; j < geneno; ++j)
	if(Gene[j] != -1){
		for (int k = Gene[j]- OLIG_LENGTH + 1; k <= Gene[j]; ++k)
			fprintf(oligos, "%c", gSequence[k]);
		fprintf(oligos, "\n");
		++olig_total;
	}*/
	//cout << "Total number of unique oligos copied into the output file are: " << olig_total << " for "  << geneno  << " genes.\n";

	for(int i=0 ; i < geneno ; ++i)
		delete tmFreq[i];
	delete [] tmFreq;
	delete [] sumPm;

	unsigned long sumPos = 0;
	for (int i=0; i < FRIME; ++i) {
		sumPos += hashTable1[i].size + hashTable2[i].size + 
		hashTable3[i].size + hashTable4[i].size +
		hashTable5[i].size + hashTable6[i].size +
		hashTable7[i].size + hashTable8[i].size ;
	}
	memoryUsed = sumPos * sizeof (int) + 8 * FRIME * sizeof (entry);
	
	delete [] Gene;
	
	delete [] hashTable1->pos;
	delete [] hashTable1;
	
	delete [] hashTable2->pos;
	delete [] hashTable2;
	
	delete [] hashTable3->pos;
	delete [] hashTable3;
	
	delete [] hashTable4->pos;
	delete [] hashTable4;
	
	delete [] hashTable5->pos;
	delete [] hashTable5;
	
	delete [] hashTable6->pos;
	delete [] hashTable6;
	
	delete [] hashTable7->pos;
	delete [] hashTable7;
	
	delete [] hashTable8->pos;
	delete [] hashTable8;
	
	return memoryUsed;
}

void destroy_Genome(genome* mygenome){
	delete [] mygenome->gCoded;
	delete [] mygenome;
}

void print_Time(){
	struct tm *current;
	time_t now;
	
	time(&now);
	current = localtime(&now);
	
	printf("Time is %i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
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
