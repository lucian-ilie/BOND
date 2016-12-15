#include "genome.h"

unsigned int PRIME_IDEN = 113493747;//34040383;//520094747;//226493747;//113493747;//226493747;//436208447;//
unsigned int PRIME = 3277283;//8126747;//22020227;//75012279;//34040383;//
unsigned int FRIME = 850027;//3277283;//750083;//650227;//1376447;//750077;//200003;//
int TEST_PARAM = 89;
int POS_SIZE = 8;
int BLOCK = 64;
int KANE_SIMILARITY = 75;
int KANE_IDENTITY = 15;
int OLIG_LENGTH = 50;
int DELTATM = 5;
int GC_CONT_MIN = 30;
int GC_CONT_MAX= 70;
int geneno;
int m;
uint16_t preNOR[65536];
int MAX_NUM_THREAD;

using namespace std;

int main (int argc, char * const argv[]) {
	cout << "------------------------------BEGIN------------------------------------\n";
	print_Time();

	time_t seconds, secondf, p1s,p1f,p2s,p2f;	
	seconds = time (NULL);
	clock_t startTotal, finishTotal, startP1, finishP1, startP2, finishP2, start_iden, finish_iden;
	
	
	pre_Process();
		
	unsigned long memoryUsed = 0, peakUsed = 0;
	
	startTotal = clock();
	startP1 = clock();
	
	FILE *gf, *oligos;
	unsigned char *tempSequence, *gSequence;
	genome *mygenome;
	mygenome = new genome;
	uint8_t *simarray;
	int tempSize;
	
	const char *gfname = argv[1];//
	const char *oligosname = argv[2];//

	oligos = fopen(oligosname, "w");
	
	/* Open genome file for reading. */
	if((gf = fopen(gfname, "rb")) == NULL) {
		cout << "Error! cannot open file: "<< gfname;
		perror(NULL);
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
	tempSequence = new unsigned char[tempSize];
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	gSequence = new unsigned char[tempSize];
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	simarray = new uint8_t [tempSize];
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
	int *geneLoc = new int [locSize];
	while (i < tempSize) {
		switch (tempSequence[i]) {
			case '>':
				if (geneno >= locSize){
					locSize *= 2;
					geneLoc = (int *) realloc(geneLoc, locSize * sizeof(int));
				}
				geneLoc[geneno] = gene_num;
				while (tempSequence[i] != '\n') ++i;
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
	delete [] tempSequence;
	memoryUsed -= (tempSize) * sizeof(unsigned char);
	
	mygenome->gSize = gene_num;
	if (mygenome->gSize < tempSize){
		gSequence = (unsigned char *) realloc(gSequence, mygenome->gSize * sizeof(unsigned char));
		simarray = (uint8_t *) realloc(simarray, mygenome->gSize * sizeof(uint8_t));	
	}
	
	memoryUsed -= (2*tempSize) * sizeof(unsigned char);
	memoryUsed += (2*mygenome->gSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	if (geneno < locSize)
		geneLoc = (int *) realloc(geneLoc, geneno * sizeof(int));
	
	memoryUsed += (geneno) * sizeof(int);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	
	int gseqOff = mygenome->gSize % 4;
	
	if(gseqOff != 0){
		for(int i=1; i < geneno; ++i)
			geneLoc[i] -= gseqOff;
		for(int i=0; i < mygenome->gSize - gseqOff; ++i)
			simarray[i] = simarray[i+gseqOff];
	}
	
	/*          Start Encoding the genome 
	 *************************************************** */
	mygenome->nbytes = mygenome->gSize/4;
	mygenome->gCoded = new uint8_t[mygenome->nbytes];
	memoryUsed += mygenome->nbytes * sizeof(uint8_t);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	uint8_t temp = 0;	
	int stpoint = mygenome->gSize % 4;
	for(int i = stpoint; i < mygenome->gSize; i += 4) {
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
		mygenome->gCoded[(i-stpoint)/4] = temp;
	}
	/*          END Encoding the genome 
	 *************************************************** */
	mygenome->gSize -= gseqOff;
	simarray = (uint8_t *) realloc(simarray, mygenome->gSize * sizeof(uint8_t));
	////////	////////	////////	//////// END create genome	
	
	cout << "\nGenome size is: " << (mygenome->gSize + gseqOff) << ", and number of target sequences are: "<< geneno <<"\n\n";
	
	//uint8_t *simarray = new uint8_t [mygenome->gSize];
	//memoryUsed += (mygenome->gSize) * sizeof(uint8_t);
	//if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (int i = 0; i < OLIG_LENGTH - 1; ++i) 
		simarray[i] = 1;
	
	// Kane Identity Check 	
	memoryUsed += PRIME_IDEN * (sizeof(uint32_t) + sizeof(int));
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	//cout << "Mem. used for identity_check: " << ((memoryUsed/1024)/1024) << endl;
	
	start_iden = clock();
	identity_Check(mygenome, simarray);	
	finish_iden = clock();
	//printf("identity_check time: %.1f sec\n\n", (double)((finish_iden - start_iden)) / (double)CLOCKS_PER_SEC);
	
	memoryUsed -= PRIME_IDEN * (sizeof(uint32_t) + sizeof(int));
	
	
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
	
	p1s = time(NULL);	
 	MAX_NUM_THREAD = omp_get_max_threads();	
	//printf("Max. number of threads = %d\n", MAX_NUM_THREAD);
	//printf("Phase 1 is running with seed set: S8W10\n");
	
	if (MAX_NUM_THREAD >= 8)
		omp_set_num_threads(8);
	else
		omp_set_num_threads(MAX_NUM_THREAD);
	
	memoryUsed += MAX_NUM_THREAD*PRIME * (sizeof(uint64_t) + sizeof(int));
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	int tid, chunk=1;
	#pragma omp parallel for shared(mygenome, simarray, seed_Phase1, chunk) private(i, tid) schedule(dynamic,chunk)
	for (i = 0; i < 8; ++i){
		tid = omp_get_thread_num();
		//printf("Thread %d doing section%d\n",tid, i);
		similarity_Check(mygenome, simarray, seed_Phase1[i]);// S8W10L18	
		//printf("Thread %d done section%d.\n",tid, i);
	}		
	
	p1f = time(NULL);
	//printf ("Time for phase 1: %ld sec\n", (p1f-p1s));
	//cout << "Peak memory used in MB: " << ((peakUsed/1024)/1024) << endl;
	
	int count = 0;
	for (int i = 0; i < mygenome->gSize; ++i) 
		if (simarray[i] ==1)
			count++;
	//cout << "Similar positions in Phase 1: " << count << "\t%" << ( 100 * double(count) / mygenome->gSize) << "\n\n";
	
	
	memoryUsed -= MAX_NUM_THREAD*PRIME * (sizeof(uint64_t) + sizeof(int));
	
	finishP1 = clock();
	
	//cout << "Mem. used before Phase2: " << ((memoryUsed/1024)/1024) << endl;
	
	
	// G-C Contetnt Check
	gc_Content(mygenome, simarray);
	
	// Tm Check
	float *tmarray = new float [mygenome->gSize];
	memoryUsed += (mygenome->gSize) * sizeof(float);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (int i = 0; i < mygenome->gSize; ++i) 
		tmarray[i] = 0;	
	
	melting_Temp(mygenome, simarray, geneLoc, tmarray);
	
	// Tm Check Done!
	
	startP2 = clock();
	p2s = time(NULL);
	unsigned long memoryEval = 0;
	
	if (TEST_PARAM == 89){
		//printf("Phase 2 is running with seed set: S8W9\n");
		memoryEval = print_Oligo_S8W9(mygenome, simarray, gSequence, geneLoc, oligos, gseqOff, tmarray);
	}
		
	p2f = time(NULL);
	//printf ("Time for phase 2: %ld sec\n\n", (p2f-p2s));
	finishP2 = clock();
	
	delete [] tmarray;
	memoryUsed -= (mygenome->gSize) * sizeof(float);
	
	if (memoryUsed + memoryEval > peakUsed) peakUsed = memoryUsed + memoryEval;
	
	//cout << "Peak memory used in MB: " << ((peakUsed/1024)/1024) << endl;
	
	fclose(oligos);//oligos.close();
	
	finishTotal = clock();	
	
	secondf = time (NULL);
	printf ("Running Time: %ld sec\n\n", (secondf-seconds));	
	destroy_Genome(mygenome);
	delete [] gSequence;
	delete [] simarray;	
	delete [] geneLoc;
	cout<<endl;
	print_Time();
	cout << "------------------------------END--------------------------------------\n";
    return 0;
}
