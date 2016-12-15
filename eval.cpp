#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <bitset>
#include <cmath>
#include <stdint.h>

#define PRIME 200003
#define POS_SIZE 8
#define BLOCK 64
#define MAX_OLIG_LENGTH 70
#define KANE_SIMILARITY 75
#define KANE_IDENTITY 15
#define OLIG_LENGTH 50
int gSize, nbytes, gseqOff;
int geneno, m, oligcur;
int DELTATM = 5;

using namespace std;

struct entry {
	uint64_t key;
	int *pos;
	int size;
};

struct poslist {
	int *pos;
	int *sim;
	int *gen;
	int size;
};

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

uint64_t get64bit(int rpos, uint8_t *Coded){
	if (rpos < 0){
		cout << "error!!! in get64bit\n";
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	int lpos = rpos - 7;
	if (lpos < 0 && rpos >= 0) 
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (int i = lpos; i <= rpos; ++i) 
		genomeChunk = (uint64_t)(genomeChunk << 8 | Coded[i]);
	return genomeChunk;	
}

int gene_find(int p, int *FindGene){
	int gene;
	int i = p >> m;
	int k = FindGene[i];
	if (k<=0)
		gene = - k;
	else{
		int j=i+1;
		int x = FindGene[j];
		while (x > 0){	
			j=j+1;
			x = FindGene[j]; 
		}
		if (p <= k)
			gene = -x-j+i; 
		else
			gene = -x-j+i+1;
	}
	return gene;
}

int insert(uint64_t k, int locn, entry *T){
	int i=0, j;
	do {
		j = (k + i) % PRIME;
		if (T[j].pos == NULL){
			T[j].key = k;
			T[j].size = POS_SIZE;
			T[j].pos = (int *) malloc(T[j].size * sizeof(int));
			//T[j].pos = new int [T[j].size];
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
	}while (i!=PRIME);
	
	T[0].pos = NULL;
	T[0].key = k;
	T[0].size = POS_SIZE;
	T[0].pos = (int *) malloc(T[0].size * sizeof(int));
	//T[0].pos = new int [T[0].size];
	T[0].pos[0] = locn;
	for (int y = 1; y < T[0].size; ++y)
		T[0].pos[y] = -1;
	return 0;
}

void create_Hash(entry *hashTable, uint8_t *gCoded, int nbytes, const char *seedSeq){
	//Begin coding seed
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
	//End coding seed
	
	//Begin Hashing
	uint64_t window = 0, hashKey=0;
	int location = -1;
	
	for (int i=0; i < PRIME; ++i) {
		hashTable[i].key = 0;
		hashTable[i].pos = NULL;
		hashTable[i].size = 0;
	}
	if(nbytes-1 <0) cout<<"create_Hash\n";
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
		j = (int)((k + i) % PRIME);
		if (T[j].key == k)
			return j;
		++i;
	}while (i!=PRIME && T[j].pos!=NULL);
	return -1;
}

void check_50factor(int pose, uint8_t *gCoded, int posc, uint8_t *oCoded, int olig_no, 
					poslist *posRep, int *FindGene){
	
	
	int ebyte = pose/4, eoff = pose%4, cbyte = posc;
	uint64_t ehitl2, ehitl, ehitr, chitl2, chitl, chitr, nolxnr, nolxnl, nolxnl2;
	
	if(ebyte - 8 < 0 || cbyte - 8 < 0)
		cout << "Error in get64bit in check_50factor! pose and posc are: " << pose << "\t" << posc << "\n";
	
	ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded); 
	chitr = get64bit(cbyte, oCoded); chitl = get64bit(cbyte - 8, oCoded); 
	
	if (ebyte <16)
		cout << "here\n";
	ehitl2 = get64bit(ebyte - 16, gCoded);
	
	
	ehitr >>= 2*(3-eoff);
	ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) -1)) << (BLOCK - 2*(3-eoff));
	ehitl >>= 2*(3-eoff);
	ehitl |= (ehitl2 & (uint64_t) ( (1 << 2*(3-eoff)) -1)) << (BLOCK - 2*(3-eoff));
	
	
	
	int kane75;
	
	////////////////////////////////// BEGIN: Checking for 75% similarity
	kane75 = 0;
	nolxnr = ~(ehitr ^ chitr);
	
	int stcCounter = 0, bestStc = 0;
	
	for (int i=1; i<=32; ++i) {
		if ((nolxnr & 3) == 3){
			kane75+=2;
			++ stcCounter;
			if (stcCounter > bestStc) bestStc = stcCounter;
			
		}
		else 
			stcCounter = 0;

		nolxnr >>= 2;
	}
	if (oligcur <= 64){
		nolxnl = ~(ehitl ^ chitl);
		for (int i=1; i<=oligcur-32; ++i) {
			if ((nolxnl & 3) == 3){
				kane75+=2;
				++ stcCounter;
				if (stcCounter > bestStc) bestStc = stcCounter;
			}
			else 
				stcCounter = 0;
			
			nolxnl >>= 2;
		}
	}
	else{
		
		
		if (ebyte <16 || cbyte <16)
			cout << "here in else check_50_factor\n";
		chitl2 = get64bit(cbyte - 16, oCoded);
		
		ehitl2 >>= 2*(3-eoff);
		
		nolxnl = ~(ehitl ^ chitl);
		nolxnl2 = ~(ehitl2 ^ chitl2);
		
		for (int i=1; i<=32; ++i) {
			if ((nolxnl & 3) == 3){ 
				kane75+=2;
				++ stcCounter;
				if (stcCounter > bestStc) bestStc = stcCounter;
			}
			else 
				stcCounter = 0;
			
			nolxnl >>= 2;
		}
		for (int i=1; i<=oligcur-64; ++i) {
			if ((nolxnl2 & 3) == 3){
				kane75+=2;
				++ stcCounter;
				if (stcCounter > bestStc) bestStc = stcCounter;
			}
			else 
				stcCounter = 0;
			
			nolxnl2 >>= 2;
		}
	}
	
	
	
	////////////////////////////////// END: Checking for 75% similarity.
	kane75 *= (50.0 / oligcur); 
	
	if (  (kane75 > KANE_SIMILARITY - 1)  ||  (bestStc >= 15)  ){// && kane75 < 100){
		
		int geneid = gene_find(pose, FindGene);
		
		bool found = false;
		int sindex = 0;
		while (sindex < posRep[olig_no].size && (!found)) {
			if (posRep[olig_no].gen[sindex] == geneid)
				found = true;
			++sindex;
		}
		if (!found){
			/*fout << olig_no << "#:\t" << pose << "\t";
			 for (int i=0; i <= posc; ++i){
			 bs = oCoded[i];
			 fout << bs; 
			 }
			 fout << "\t" << kane75 << "\n";*/	
			int lastpos = 0;
			
			while (posRep[olig_no].pos[lastpos] !=-1 && lastpos < posRep[olig_no].size)
				++lastpos;
			
			if (lastpos < posRep[olig_no].size){
				posRep[olig_no].pos[lastpos] = pose;
				posRep[olig_no].sim[lastpos] = kane75; // new 
				posRep[olig_no].gen[lastpos] = geneid; // new 
			}
			else {
				posRep[olig_no].size *= 2;
				posRep[olig_no].pos = (int *) realloc(posRep[olig_no].pos, posRep[olig_no].size * sizeof(int));
				posRep[olig_no].pos[lastpos] = pose;
				for (int y = lastpos + 1; y < posRep[olig_no].size; ++y)
					posRep[olig_no].pos[y] = -1;
				
				posRep[olig_no].sim = (int *) realloc(posRep[olig_no].sim, posRep[olig_no].size * sizeof(int));	//new
				posRep[olig_no].sim[lastpos] = kane75; // new
				
				posRep[olig_no].gen = (int *) realloc(posRep[olig_no].gen, posRep[olig_no].size * sizeof(int));	//new
				posRep[olig_no].gen[lastpos] = geneid; // new
			}
		}
	}
	//else if (kane75 >= 50 && bestStc >= 15)
		
}

void create_Report(int oSize, unsigned char *oSequence, const char *seedSeq1, entry *hashTable1, 
				   uint8_t *gCoded, poslist *posRep, int * FindGene){
	int oligcul = 0, stbyte, olig_no = 0;
	uint8_t *oCoded;
	while (oligcul < oSize){
		////////////////////////////////////////////////////////////////Computing coded oligo
		////////////////////////////////
		++olig_no;
		oligcur = 0; // length of oligo		
		while (oSequence[oligcul+oligcur] != '\n')
			++oligcur;
		stbyte = oligcur % 4 == 0? 0 : 1;
		oCoded = (uint8_t *) malloc((oligcur/4 + stbyte) * sizeof(uint8_t));
		//oCoded = new uint8_t [oligcur/4 + stbyte];
		uint8_t temp = 0;
		if (oligcur % 4 != 0) { // leftmost byte of coded olig
			for(int j = oligcul; j < oligcul + oligcur % 4; ++j) {
				switch (oSequence[j]) {
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
						//cout << "Error! " << oSequence[j] << "is not a valid nucleotide.\n";
						temp = (uint8_t) (temp << 2 | 0);					
				}	
			}
			oCoded[0] = temp;
		}
		for(int i = oligcul + oligcur % 4; i < oligcul + oligcur; i += 4) {
			for(int j = 0; j < 4; j++) 
				switch (oSequence[j+i]) {
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
						//cout << "Error! " << oSequence[j+i] << "is not a valid nucleotide.\n";
						temp = (uint8_t) (temp << 2 | 0);					
				}	
			oCoded[(i-oligcul)/4+ stbyte] = temp;
		}
		oligcul += oligcur +1; 
		////////////////////////////////
		////////////////////////////////////////////////////////////////END Computing coded oligo
		
		
		////////////////////////////////////////////////////////////////Computing Hash value
		////////////////////////////////
		uint64_t window = 0, hashKey1;
		uint64_t seed1=get64seed(seedSeq1); 
		int w; int pos1;
		////////////////////////////////////////////////////////////////////Checking Positions
		if(oligcur/4 + stbyte-1 <0) cout << "In create_Report, oligocur: " << oligcur << "\t oligo_no: " << olig_no << endl;
		int rpos=oligcur/4 + stbyte-1;
		window = get64bit(rpos, oCoded);
		rpos-=8;
		
		int offpos = 0;
		
		uint8_t nexwin=0;
		if (rpos >= 0) nexwin=oCoded[rpos];
		--rpos;
		int nexwinSize = 4;
		
		for (int pass = 1; pass <= oligcur - strlen(seedSeq1) +1 ; ++pass) {		
			hashKey1 = (uint64_t) (window & seed1);  
			pos1 = search_Hash(hashKey1, hashTable1);											
			w=0;
			while (w < hashTable1[pos1].size && hashTable1[pos1].pos[w] !=-1 ){
				if (hashTable1[pos1].pos[w] > 63)
					check_50factor(hashTable1[pos1].pos[w] + offpos , gCoded, 
								   oligcur/4 + stbyte -1 , oCoded, olig_no-1, posRep, FindGene); // new idea offpos
				++w;
			}
			
			window = (window >> 2) | ( ((uint64_t)(nexwin & 3)) << (BLOCK - 2));
			nexwin >>= 2;  --nexwinSize;
			if (nexwinSize == 0){
				if(rpos>=0){
					nexwin = oCoded[rpos];
					--rpos;
				}
				nexwinSize = 4;
			}	

			++offpos;
		}
		
		/*for (int i = oligcur/4 + stbyte -1; i >= 7; --i) {		
			for (int pass=3; pass >=0 ; --pass) {	
				
				hashKey1 = (uint64_t) (window & seed1);  
				pos1 = search_Hash(hashKey1, hashTable1);											
				w=0;
				while (w < hashTable1[pos1].size && hashTable1[pos1].pos[w] !=-1 ){
					if (hashTable1[pos1].pos[w] > 63)
						check_50factor(hashTable1[pos1].pos[w] + offpos , gCoded, 
									   oligcur/4 + stbyte -1 , oCoded, olig_no-1, posRep, FindGene); // new idea offpos
					++w;
				}
				window >>=  2;
				++offpos;
			}
			if (i>=8)
				window |= (((uint64_t) oCoded[i-8]) << (BLOCK - 8));
		}*/
		
		
		////////////////////////////////////////////////////////////////////Checking Positions
		
		////////////////////////////////////////////////////////////////Computing Hash value
		////////////////////////////////
	}
	free(oCoded);
	
}

int main (int argc, char * const argv[]) {
	clock_t start, finish;
	start = clock();
	//////////////////////////////////////////////////creating genome
	FILE *gf, * genef;
	unsigned char *tempSequence, *gSequence;
	int tempSize;
	
	const char *gfname = argv[1];//
	const char *ofname = argv[2];//
	const char *foutfn = argv[3];//
	
	
	
	/* Open genome file for reading. */
	if((gf = fopen(gfname, "rb")) == NULL) {
		cout << "Error! cannot open file: "<< gfname << endl;
		exit(EXIT_FAILURE);
	}
	
	
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
	gSequence = (unsigned char *) malloc(tempSize * sizeof(unsigned char));
	
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
	//int *geneLoc = new int [locSize];
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
				++errorno;
				++gene_num;
				break;
		}
		++i;
	}
	free(tempSequence);//free(tempSequence;
	
	gSize = gene_num;
	if (gSize < tempSize)
		gSequence = (unsigned char *) realloc(gSequence, gSize * sizeof(unsigned char));
	
	if (geneno < locSize)
		geneLoc = (int *) realloc(geneLoc, geneno * sizeof(int));
	
	gseqOff = gSize % 4;
	
	if(gseqOff != 0)
		for(int i=1; i < geneno; ++i)
			geneLoc[i] -= gseqOff;

	
	
	/*          Start Encoding the genome 
	 *************************************************** */
	nbytes = gSize/4;
	uint8_t *gCoded = (uint8_t *) malloc(nbytes * sizeof(uint8_t));
	//uint8_t *gCoded = new uint8_t[nbytes];
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
	gSize = gSize - gSize % 4;
	
	//////////////////////////////////////////////////creating genome
		
	//////////////////////////////////////////////////BEGIN Read Oligo Set
	FILE *oligf;
	unsigned char *oSequence;
	int oSize;
	if((oligf = fopen(ofname, "rb")) == NULL) {
		cout << "Error! cannot open file: " << ofname << endl;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	/* Get the genome file size. */
	if(fseek(oligf, 0, SEEK_END) == 0) {
		oSize = ftell(oligf);
		rewind(oligf);
		if(oSize < 0) {
			cout << "Error! Cannot ftell: " << ofname;
			perror(NULL);
			exit(EXIT_FAILURE);
		}
	} else {
		cout << "Error! Cannot fseek: " << ofname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	/* Allocate memory to genome sequence array. */
	oSequence = (unsigned char *) malloc(oSize * sizeof(unsigned char));
	//oSequence = new unsigned char[oSize];
	if((oSequence == NULL)) {
		cout << "Error! Cannot allocate memory.";
		exit(EXIT_FAILURE);
	}
	
	/* Read gSize bytes of data. */
	if(fread(oSequence, sizeof(unsigned char), oSize, oligf) != oSize) {
		cout << "Cannot read from: "<< ofname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	//////////////////////////////////////////////////END Read Oligo Set
	
	int numolig = 0;
	for (int i=0; i< oSize; ++i) 
		if (oSequence[i] == '\n')
			++numolig;
	int *oligLength = (int *) malloc(numolig * sizeof(int));
	//int *oligLength = new int [numolig];
	numolig = 0;
	int counter = 0;	
	for (int i=0; i< oSize; ++i){		
		++counter;
		if (oSequence[i] == '\n'){
			oligLength[numolig] = counter-1;
			++numolig;
			counter = 0;
		}			
	}
	
	//////////////// BEGIN Eval
	poslist *posRep = (poslist *) malloc(numolig * sizeof(poslist));
	//poslist *posRep = new poslist [numolig];
	for (int i=0; i < numolig; ++i){
		posRep[i].pos = (int *) malloc(POS_SIZE * sizeof(int));
		//posRep[i].pos = new int [POS_SIZE];
		for (int j=0; j < POS_SIZE; ++j)
			posRep[i].pos[j] = -1;
		posRep[i].sim = (int *) malloc(POS_SIZE * sizeof(int));
		//posRep[i].sim = new int [POS_SIZE];
		for (int j=0; j < POS_SIZE; ++j)
			posRep[i].sim[j] = -1;
		posRep[i].gen = (int *) malloc(POS_SIZE * sizeof(int));
		//posRep[i].gen = new int [POS_SIZE];
		for (int j=0; j < POS_SIZE; ++j)
			posRep[i].gen[j] = -1;
		posRep[i].size = POS_SIZE;
	}
	FILE *fout;
	fout = fopen(foutfn, "w");
	
	///
	
	int minLen = gSize - geneLoc[geneno-1] + 1; int minLenInd = geneno-1;
	for (int i=0; i< geneno -1; ++i)
		if (minLen > geneLoc[i+1]-geneLoc[i]){
			minLen = geneLoc[i+1]-geneLoc[i];
			minLenInd = i;
		}
	m = log2(minLen);
	int two2m = 1; 
	two2m <<= m;
	int sizeB2m = gSize >> m;
	int *FindGene = (int *) malloc((sizeB2m + 2) * sizeof(int));
	//int *FindGene = new int [sizeB2m + 2];
	int glind=1;
	for (int i=0; i <= sizeB2m; ++i)
		if (i*two2m+two2m-1 < geneLoc[glind])
			FindGene[i] = -glind;
		else 
			if (glind <= geneno-1){
				FindGene[i]=geneLoc[glind]-1;
				++glind;
			}
			else 
				FindGene[i] = -glind;
	FindGene[sizeB2m+1] = (gSize % two2m == 1) ? -geneno-1: -geneno;	
	///
		
	
	entry *hashTable1 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];//
	const char *seedSeq1 = "11101011";
	create_Hash(hashTable1, gCoded, nbytes, seedSeq1);
	create_Report(oSize, oSequence, seedSeq1, hashTable1, gCoded, posRep, FindGene);
	free(hashTable1->pos);//free(hashTable1->pos);
	free(hashTable1);//free(hashTable1);
	cout << "Hash 1 is done!\n";
	
	entry *hashTable2 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq2 = "1100100111";
	create_Hash(hashTable2, gCoded, nbytes, seedSeq2);
	create_Report(oSize, oSequence, seedSeq2, hashTable2, gCoded, posRep, FindGene);
	free(hashTable2->pos);
	free(hashTable2);
	cout << "Hash 2 is done!\n";
	
	entry *hashTable3 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq3 = "1101000001011";
	create_Hash(hashTable3, gCoded, nbytes, seedSeq3);
	create_Report(oSize, oSequence, seedSeq3, hashTable3, gCoded, posRep, FindGene);
	free(hashTable3->pos);
	free(hashTable3);
	cout << "Hash 3 is done!\n";
	
	entry *hashTable4 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq4 = "10100010000111";
	create_Hash(hashTable4, gCoded, nbytes, seedSeq4);
	create_Report(oSize, oSequence, seedSeq4, hashTable4, gCoded, posRep, FindGene);
	free(hashTable4->pos);
	free(hashTable4);
	cout << "Hash 4 is done!\n";
	
	entry *hashTable5 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq5 = "101001001000101";
	create_Hash(hashTable5, gCoded, nbytes, seedSeq5);
	create_Report(oSize, oSequence, seedSeq5, hashTable5, gCoded, posRep, FindGene);
	free(hashTable5->pos);
	free(hashTable5);
	cout << "Hash 5 is done!\n";
	
	entry *hashTable6 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq6 = "1100010000100101";
	create_Hash(hashTable6, gCoded, nbytes, seedSeq6);
	create_Report(oSize, oSequence, seedSeq6, hashTable6, gCoded, posRep, FindGene);
	free(hashTable6->pos);
	free(hashTable6);
	cout << "Hash 6 is done!\n";
	
	entry *hashTable7 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq7 = "11010000001000101";
	create_Hash(hashTable7, gCoded, nbytes, seedSeq7);
	create_Report(oSize, oSequence, seedSeq7, hashTable7, gCoded, posRep, FindGene);
	free(hashTable7->pos);
	free(hashTable7);
	cout << "Hash 7 is done!\n";
	
	entry *hashTable8 = (entry *) malloc(PRIME * sizeof(entry));//new entry [PRIME];
	const char *seedSeq8 = "110000001000010011";
	create_Hash(hashTable8, gCoded, nbytes, seedSeq8);
	create_Report(oSize, oSequence, seedSeq8, hashTable8, gCoded, posRep, FindGene);
	free(hashTable8->pos);
	free(hashTable8);
	cout << "Hash 8 is done!\n";
	
	
	fprintf(fout, "Evaluation report for: %s \n\n", ofname); //fout << "Evaluation report for: " << ofname << "\n\n";
	fprintf(fout, "Positions: Similarity: Gene# \n\n"); //fout << "Positions: Similarity\n\n";
	int repolig = 0;
	
	int *geneTar = (int *) malloc(geneno * sizeof(int));
	//int *geneTar = new int [geneno];
	for (int i=0; i < geneno; ++i)
		geneTar[i] = 0;
	
	for (int i=0; i < numolig; ++i){
		//if (posRep[i].pos[1] != -1){ 
		
		int freq[101], freqSum=0;
		for (int fc=0; fc <101; ++fc)
			freq[fc] = 0;
		
		int j=0;
		while (j < posRep[i].size && posRep[i].pos[j] !=-1){
			++freq[posRep[i].sim[j]];
			++j;
		}
		
		for (int fc=0; fc <101; ++fc)
			freqSum += freq[fc];
		
		if( freqSum > 1){
				fprintf(fout, "Bad Oligo \n");
				++repolig;
				fprintf(fout, "Oligo#%d :\n", i+1);//fout << "Oligo #" << i << ":\n";
				j=0;
				while (j < posRep[i].size && posRep[i].pos[j] !=-1){
					fprintf(fout, "%d : %d : %d : \t\t", posRep[i].pos[j] + gseqOff, posRep[i].sim[j], posRep[i].gen[j]);
					for (int k = posRep[i].pos[j]-(oligLength[i]-1); k <= posRep[i].pos[j]; ++k)							
						fprintf(fout, "%c", gSequence[k + gseqOff]);//fout << gSequence[k + gseqOff];
					fprintf(fout, "\n");//fout << "\n";
					++j; 							
				}
				
				for (int fc=100; fc >= 0; --fc)
					if (freq[fc] > 0)
						fprintf(fout, "Sim%d: %d  ", fc,freq[fc]);				
				fprintf(fout, "\n");
				fprintf(fout, "--------------------------------------------------------------------------------------------------------");
				fprintf(fout, "\n");
				
		}
		//}
	}
	
	fprintf(fout, "Non-Unique Oligos / Total Oligos: %d / %d \n", repolig, numolig);
	fprintf(fout, "Unique Oligos / Total Oligos: %d / %d \n", numolig-repolig, numolig);
	fprintf(fout, "Percentage of Non-Unique Oligos: %f \n", (repolig * 100.0 / numolig));
	
	fclose(fout);
	fclose(oligf);
	//////////////// END Eval
	
	finish = clock();
	
	cout << "\nGenome size is: " << gSize << "\n";
	printf("Running time: %.4f sec\n\n", (double)((finish - start)) / (double)CLOCKS_PER_SEC);
	

	free(gSequence);
	free(gCoded);
	free(oSequence);
	free(posRep);
	free(geneLoc);
	free(FindGene);
	free(oligLength);
	free(geneTar);
	//End Hashing
	return 0;
}
