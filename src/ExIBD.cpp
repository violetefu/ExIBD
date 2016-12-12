/***
 ExIBD: to detect IBD segments in exome sequencing data through a pipeline
 
 07/25/2016 By Wenqing Fu
***/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <omp.h>

typedef unsigned long int ulnt;

int main(int argc, char *argv[]){
	ulnt Lnum, Cnum;
	int c, threads=1, FastIBDnum=10, IBDnum=5, Enum=0, chrom;
	double fastibdthreshold=1e-10, ibd2nonibd=0.01, nonibd2ibd=0.0001, ibderror=0.005, ibdscale=2.0, cutoff=0.5, fdr=0.1;
	char *beagle=0, *file=0, *refname=strdup("./reference/hg19_FDRref.txt"), code[1024], ch;
	
	static const struct option long_opt[] = {
		{"help", no_argument, 0, 'h'},
		{"file", required_argument, 0, 'f'},
		{"numMarker", required_argument, 0, 'l'},
		{"chr", required_argument, 0, 'c'},
		{"roundFastIBD", required_argument, 0, 'r'},
		{"roundIBD", required_argument, 0, 'i'},
		{"extend", required_argument, 0, 'e'},
		{"beagle", required_argument, 0, 'b'},
		{"fastibdthreshold", required_argument, 0, 'd'},
		{"ibd2nonibd", required_argument, 0, 'u'},
		{"nonibd2ibd", required_argument, 0, 'v'},
		{"ibderror", required_argument, 0, 'w'},
		{"ibdscale", required_argument, 0, 's'},
		{"cutoff", required_argument, 0, 'p'},
		{"reference", required_argument, 0, 'a'},
		{"fdr", required_argument, 0, 'g'},
		{"threads", required_argument, 0, 't'},
		{0, 0, 0, 0}
	};
	
	while ((c = getopt_long(argc, argv, "hf:l:c:r:i:e:b:d:u:v:w:s:p:a:g:t:", long_opt, NULL)) != -1) {
		switch (c) {
			case 'h':
				printf("Usage: %s [OPTIONS]\n",argv[0]);
				printf("-h, --help:             print this help and exit\n");
				printf("-f, --file:             name for input files\n");
				printf("-l, --numMarker:        number of markers\n");
				printf("-c, --chr:              No. chromosome\n");
				printf("-r, --roundFastIBD:     number of Beagle fastIBD rounds (default: 10)\n");
				printf("-i, --roundIBD:         number of Beagle IBD rounds (default: 5)\n");
				printf("-e, --extend:           number of extended variants from candidate IBD breakpoints (default: 0)\n");
				printf("-b, --beagle:           basic commands for BEAGLE IBD\n");
				printf("-d, --fastibdthreshold: fastIBD score threshold (default: 1e-10)\n");
				printf("-u, --ibd2nonibd:       IBD to non-IBD transition rate (default: 0.01)\n");
				printf("-v, --nonibd2ibd:       non-IBD to IBD transition rate (default: 0.0001)\n");
				printf("-w, --ibderror:         genotype error rate (default: 0.005)\n");
				printf("-s, --ibdscale:         IBD and fastIBD tuning parameter (default: 2.0)\n");
				printf("-p, --cutoff:           cutoff of IBD probability for each position (default: 0.5)\n");
				printf("-a, --reference:        file for pre-assigned FDRs (default: ./reference/hg19_FDRref.txt)\n");
				printf("-g, --fdr:              cutoff of FDR (default: 0.1)\n");
				printf("-t, --threads:          number of threads (default: 1)\n");
				return(0);
			case 'f': file = strdup(optarg); break;
			case 'l': Lnum = atol(optarg); break;
			case 'c': chrom = atoi(optarg); break;
			case 'r': FastIBDnum = atoi(optarg); break;
			case 'i': IBDnum = atoi(optarg); break;
			case 'e': Enum = atoi(optarg); break;
			case 'b': beagle = strdup(optarg); break;
			case 'd': fastibdthreshold = atof(optarg); break;
			case 'u': ibd2nonibd = atof(optarg); break;
			case 'v': nonibd2ibd = atof(optarg); break;
			case 'w': ibderror = atof(optarg); break;
			case 's': ibdscale = atof(optarg); break;
			case 'p': cutoff = atof(optarg); break;
			case 'a': refname = strdup(optarg); break;
			case 'g': fdr = atof(optarg); break;
			case 't': threads = atoi(optarg); break;
		}
	}
	
	srand(time(NULL));
	ulnt unique = rand();
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(int re = 0; re < FastIBDnum; re++){
		char str[1024];
		sprintf(str,"%s out=%d fastibd=true fastibdthreshold=%g ibdscale=%g seed=%d", beagle, re, fastibdthreshold, ibdscale, rand());
		system(str);
	}
	
	sprintf(code, "./ExIBD_Candidate -f %s -r %d", file, FastIBDnum);
	system(code);
	
	Cnum = 0;
	sprintf(code, "%s.candidate", file);
	FILE *inp = fopen(code,"rt");
	while(!feof(inp)){
		ch = fgetc(inp);
		if(ch == '\n') Cnum++;
	}
	fclose(inp);
	
	sprintf(code, "./ExIBD_Refined -f %s -l %lu -n %lu -r %d -e %d -b \"%s\" -u %g -v %g -w %g -s %g -p %g -t %d", file, Lnum, Cnum, IBDnum, Enum, beagle, ibd2nonibd, nonibd2ibd, ibderror, ibdscale, cutoff, threads);
	system(code);
	
	sprintf(code, "./ExIBD_Filtered -f %s -c %d -r %s -t %g", file, chrom, refname, fdr);
	system(code);
	
	return 1;
}
