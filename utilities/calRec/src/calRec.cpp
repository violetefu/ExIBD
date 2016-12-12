/***
 Converts physical position (GRCh37, hg19) to genetic position
 
 07/28/2016 By Wenqing Fu
 ***/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

typedef unsigned long int ulnt;

int main(int argc, char *argv[]){
	ulnt Lnum, Rnum, *pos, *posi, prer;
	int c, chr;
	double *recrate, *recmap, reco;
	char str[256];
	FILE *inp, *out;
	
	static const struct option long_opt[] = {
		{"help", no_argument, 0, 'h'},
		{"input", required_argument, 0, 'i'},
		{"numMarker", required_argument, 0, 'l'},
		{"chr", required_argument, 0, 'c'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	
	while ((c = getopt_long(argc, argv, "hi:l:c:o:", long_opt, NULL)) != -1) {
		switch (c) {
			case 'h':
				printf("Usage: %s [OPTIONS]\n",argv[0]);
				printf("-h, --help:         print this help and exit\n");
				printf("-i, --input:        name for input file\n");
				printf("-l, --numMarker:    number of markers\n");
				printf("-c, --chr:          No. chromosome\n");
				printf("-o, --output:       name for output file\n");
				return 0;
			case 'i': inp = fopen(optarg, "rt"); break;
			case 'l': Lnum = atol(optarg); break;
			case 'c': chr = atol(optarg); break;
			case 'o': out = fopen(optarg, "wt"); break;
		}
	}
	
	pos = new ulnt[Lnum];
	for(ulnt l=0;l<Lnum;l++) fscanf(inp,"%lu",&pos[l]);
	fclose(inp);
	
	sprintf(str,"genetic_map_GRCh37_chr%d.txt",chr);
	inp=fopen(str,"rt");
	fscanf(inp,"%lu",&Rnum);
	posi=new ulnt[Rnum];recrate=new double[Rnum];recmap=new double[Rnum];
	for(ulnt r=0;r<Rnum;r++) fscanf(inp,"%*s%lu%lf%lf",&posi[r],&recrate[r],&recmap[r]);
	fclose(inp);
	
	prer = 0;
	for(ulnt l = 0; l < Lnum; l++){
		if(pos[l] < posi[0]) reco=0;
		else{
			ulnt r;
			for(r = prer; r < (Rnum-1); r++){
				if(pos[l] >= posi[r] && pos[l] < posi[r+1]){
					reco=recmap[r]+(pos[l]-posi[r])*recrate[r]/1000000;
					prer=r; break;
				}
			}
			if(r == Rnum-1) reco=recmap[r-1]+(pos[l]-posi[r-1])*recrate[r-1]/1000000;
		}
		fprintf(out,"%lu\t%g\n",pos[l], reco);
	}
	fclose(out);
	delete []pos;delete []posi;delete []recrate;delete []recmap;
	return 1;
}

