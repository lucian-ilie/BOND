#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <stdint.h>
using namespace std;

int main (int argc, char * const argv[]) {	
	FILE *gf;
	ofstream goutf;
	const char *gfname = argv[1];//"../../PICKY/Result/ecoli.picky";//
	  goutf.open(argv[2]);
	unsigned char *gSequence;
	uint64_t gSize;
	
	/* Open genome file for reading. */
	if((gf = fopen(gfname, "rb")) == NULL) {
		cout << "Error! cannot open file: "<< gfname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	/* Get the genome file size. */
	if(fseek(gf, 0, SEEK_END) == 0) {
		gSize = ftell(gf);
		rewind(gf);
		if(gSize < 0) {
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
	gSequence = new unsigned char[gSize];
	if((gSequence == NULL)) {
		cout << "Error! Cannot allocate memory.";
		exit(EXIT_FAILURE);
	}
	
	/* Read gSize bytes of data. */
	if(fread(gSequence, sizeof(unsigned char), gSize, gf) != gSize) {
		cout << "Cannot read from: "<< gfname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	
	
	
	long incur = 0;
	/*while(gSequence[incur] != '\n')
		++incur;*/
	//++incur;
	while (incur < gSize) {
		while(gSequence[incur] != '\t' && incur < gSize)
			++incur;
		++incur;
		while(gSequence[incur] != '\t' && incur < gSize){			
			goutf << gSequence[incur];
			++incur;
		}
		while(gSequence[incur] != '\n') ++incur;
		++incur;
		goutf << "\n";
	}
	//goutf << "\n";
	fclose(gf);
	delete [] gSequence;
	goutf.close();
	cout<<"done!\n";
    return 0;

}


