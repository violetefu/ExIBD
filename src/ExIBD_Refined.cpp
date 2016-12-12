/***
 ExIBD_Refined: the 2nd step for ExIBD;
 to refine the breakpoints of candidate IBD segments by Beagle IBD
 
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
	ulnt Lnum, Cnum, *start, *end;
	int c, threads=1, Rnum=5, Enum=0, *pairnum;
	double ibd2nonibd=0.01, nonibd2ibd=0.0001, ibderror=0.005, ibdscale=2.0, cutoff=0.5;
	char *beagle=0, *file=0, filename[1024], **marker, ***pairs;
	FILE *inp, **out;
	
	static const struct option long_opt[] = {
		{"help", no_argument, 0, 'h'},
		{"file", required_argument, 0, 'f'},
		{"numMarker", required_argument, 0, 'l'},
		{"numFIBD", required_argument, 0, 'n'},
		{"round", required_argument, 0, 'r'},
		{"extend", required_argument, 0, 'e'},
		{"beagle", required_argument, 0, 'b'},
		{"ibd2nonibd", required_argument, 0, 'u'},
		{"nonibd2ibd", required_argument, 0, 'v'},
		{"ibderror", required_argument, 0, 'w'},
		{"ibdscale", required_argument, 0, 's'},
		{"cutoff", required_argument, 0, 'p'},
		{"threads", required_argument, 0, 't'},
		{0, 0, 0, 0}
	};
	
	while ((c = getopt_long(argc, argv, "hf:l:n:r:e:b:u:v:w:s:p:t:", long_opt, NULL)) != -1) {
		switch (c) {
			case 'h':
				printf("Usage: %s [OPTIONS]\n",argv[0]);
				printf("-h, --help:         print this help and exit\n");
				printf("-f, --file:         name for input files\n");
				printf("-l, --numMarker:    number of markers\n");
				printf("-n, --numFIBD:      row number of the candidate IBD file\n");
				printf("-r, --round:        number of Beagle IBD rounds (default: 5)\n");
				printf("-e, --extend:       number of extended variants from candidate IBD breakpoints (default: 0)\n");
				printf("-b, --beagle:       basic commands for BEAGLE IBD\n");
				printf("-u, --ibd2nonibd:   IBD to non-IBD transition rate (default: 0.01)\n");
				printf("-v, --nonibd2ibd:   non-IBD to IBD transition rate (default: 0.0001)\n");
				printf("-w, --ibderror:     genotype error rate (default: 0.005)\n");
				printf("-s, --ibdscale:     IBD and fastIBD tuning parameter (default: 2.0)\n");
				printf("-p, --cutoff:       cutoff of IBD probability for each position (default: 0.5)\n");
				printf("-t, --threads:      number of threads (default: 1)\n");
				return(0);
			case 'f': file = strdup(optarg); break;
			case 'l': Lnum=atol(optarg); break;
			case 'n': Cnum = atol(optarg); break;
			case 'r': Rnum = atoi(optarg); break;
			case 'e': Enum = atoi(optarg); break;
			case 'b': beagle = strdup(optarg); break;
			case 'u': ibd2nonibd = atof(optarg); break;
			case 'v': nonibd2ibd = atof(optarg); break;
			case 'w': ibderror = atof(optarg); break;
			case 's': ibdscale = atof(optarg); break;
			case 'p': cutoff = atof(optarg); break;
			case 't': threads = atoi(optarg); break;
		}
	}
	
	marker = new char *[Lnum];
	for(ulnt l = 0; l < Lnum; l++) marker[l] = new char[20];
	start = new ulnt[Cnum]; end = new ulnt[Cnum]; pairnum = new int[Cnum]; pairs = new char **[Cnum];
	out = new FILE *[threads];
	
	// Read Marker Identifiers
	sprintf(filename,"%s.info", file);
	inp = fopen(filename, "rt");
	for(ulnt l = 0; l < Lnum; l++) fscanf(inp,"%s%*s%*s%*s", marker[l]);
	fclose(inp);
	
	// Read candidate IBD segments
	sprintf(filename,"%s.candidate", file);
	inp = fopen(filename, "rt");
	for(ulnt i = 0; i < Cnum; i++){
		fscanf(inp,"%lu%lu%d",&start[i],&end[i],&pairnum[i]);
		start[i] = (start[i] < Enum) ? 0 : (start[i] - Enum);
		end[i] = ((end[i]+Enum) >= (Lnum-1)) ? (Lnum-1) : (end[i]+Enum);
		pairs[i] = new char *[2*pairnum[i]];
		for(int p = 0; p < 2*pairnum[i]; p++) pairs[i][p] = new char[256];
		for(int p = 0; p < pairnum[i]; p++) fscanf(inp,"%s%s", pairs[i][2*p], pairs[i][2*p+1]);
	}
	fclose(inp);
	
	for(int t = 0; t < threads; t++){
		sprintf(filename, "%s.%d.ribd",file, t);
		out[t] = fopen(filename,"wt");
	}
	
	srand(time(NULL));
	ulnt unique = rand();
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(ulnt i = 0; i < Cnum; i++){
		ulnt sLnum = end[i]-start[i]+1, refinel[2];
		int core = omp_get_thread_num();
		double indpro, **pro;
		char str[1024];
		FILE *tp;
		
		// Create file for IBD pairs
		sprintf(str, "%lu_%lu.pair", unique, i);
		tp = fopen(str, "wt");
		for(int p = 0; p < pairnum[i]; p++) fprintf(tp, "%s %s\n",pairs[i][2*p], pairs[i][2*p+1]);
		fclose(tp);
		
		// Create file for excluded markers
		sprintf(str, "%lu_%lu.exl", unique, i);
		tp = fopen(str, "wt");
		for(ulnt l = 0; l < Lnum; l++) if(l < start[i] || l > end[i]) fprintf(tp, "%s\n", marker[l]);
		fclose(tp);
		
		pro = new double *[pairnum[i]];
		for(int p = 0; p < pairnum[i]; p++) pro[p] = new double[sLnum];
		for(int re = 0; re < Rnum; re++){
			// Run Beagle IBD
			sprintf(str,"%s markers=%s.info ibdpairs=%lu_%lu.pair excludemarkers=%lu_%lu.exl out=%lu_%lu ibd2nonibd=%g nonibd2ibd=%g ibderror=%g ibdscale=%g seed=%d", beagle, file, unique, i, unique, i, unique, i, ibd2nonibd, nonibd2ibd, ibderror, ibdscale, rand());
			system(str);
			
			// Summary of Beagle IBD
			sprintf(str,"%lu_%lu.%s.ibd", unique, i, file);
			tp = fopen(str, "rt");
			fscanf(tp,"%*s");
			for(int p = 0; p < pairnum[i]; p++) fscanf(tp,"%*s%*s%*s%*s");
			fscanf(tp,"%*s");
			for(int p = 0; p < pairnum[i]; p++) fscanf(tp,"%*s%*s%*s%*s");
			for(ulnt l = 0; l < sLnum; l++){
				fscanf(tp,"%*s");
				for(int p = 0; p < pairnum[i]; p++){
					fscanf(tp,"%lf%*s%*s%*s",&indpro);
					if(re==0 || pro[p][l]<indpro) pro[p][l]=indpro;
				}
			}
			fclose(tp);
		}
		
		// refined IBD segments calling
		for(int p = 0; p < pairnum[i]; p++){
			refinel[0] = refinel[1] = 0;
			for(ulnt l = 0; l < sLnum; l++) if(pro[p][l] >= cutoff){
				refinel[0] = refinel[1] = l;
				for(ulnt ll = l+1; ll < sLnum; ll++){
					if(pro[p][ll] >= cutoff) refinel[1] = ll;
					else{
						fprintf(out[core],"%s\t%s\t%lu\t%lu\n",pairs[i][2*p],pairs[i][2*p+1],refinel[0]+start[i],refinel[1]+start[i]);
						l = ll; break;
					}
					if(ll == (sLnum-1) && pro[p][ll] >= cutoff){
						fprintf(out[core],"%s\t%s\t%lu\t%lu\n",pairs[i][2*p],pairs[i][2*p+1],refinel[0]+start[i],refinel[1]+start[i]);
						l = sLnum;				
					}
				}
			}
		}
		for(int p = 0; p < pairnum[i]; p++) delete pro[p];
		delete []pro;
		sprintf(str,"rm ./%lu_%lu.*", unique, i);
		system(str);
	}
	for(int t = 0; t < threads; t++) fclose(out[t]);
	
	sprintf(filename,"cat %s.*.ribd > %s.ribd", file, file);
	system(filename);
	sprintf(filename,"rm %s.*.ribd", file);
	system(filename);
	
	for(ulnt l = 0; l < Lnum; l++) delete marker[l];
	for(ulnt i = 0; i < Cnum; i++){
		for(int p = 0; p < 2*pairnum[i]; p++) delete pairs[i][p];
		delete pairs[i];
	}
	delete []marker; delete []start; delete []end; delete []pairnum; delete []pairs; delete []out;
	return 1;
}
