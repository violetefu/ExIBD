'''
 ExIBD_Filtered: the 3rd step for ExIBD;
 to filter out IBD segment if it spans a genomic region where the corresponding 
 locus specific false discovery rate (FDR) exceeds a desired cutoff.
 
 07/25/2016 By Wenqing Fu
'''
import sys
import getopt
import pandas as pd

def IBDType(x):
    num = 10
    if x < 0.5:
        return 0
    for k in range(1, (num+1), 1):
        if x >= (k-0.5) and x < (k+0.5):
            return k
    return 'L'

def getFDR(center, refpop1, refpop2, IBDsize, reference):
    try:
        val1 = reference[(center >= reference.Gstart) & (center < reference.Gend)][refpop1+'_FDR'+str(IBDsize)].values[0]
        if refpop1 != refpop2:
            val2 = reference[(center >= reference.Gstart) & (center < reference.Gend)][refpop2+'_FDR'+str(IBDsize)].values[0]
            val1 = val1 if val1 > val2 else val2
        return val1
    except KeyError:
        return 1
    
def main(argv):
    refile = './FDRref.txt'
    fdr = 0.1
    
    try:
        opts, args = getopt.getopt(argv,'hf:c:r:t:',['help','file=','chr=','reference=','fdr='])
    except getopt.GetoptError:
        print 'python ./ExIBD_Filtered.py -f <filename> -c <chr> -r <reference file> -t <fdr cutoff>'
        print '-h, --help:       print this help and exit'
        print '-f, --file:       name for input files'
        print '-c, --chr:        No. chromosome'
        print '-r, --reference:  file for pre-assigned FDRs (default: ./FDRref.txt)'
        print '-t, --fdr:        cutoff of FDR (default: 0.1)'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print 'python ./ExIBD_Filtered.py -f <filename> -c <chr> -r <reference file> -t <fdr cutoff>'
            print '-h, --help:       print this help and exit'
            print '-f, --file:       name for input files'
            print '-c, --chr:        No. chromosome'
            print '-r, --reference:  file for pre-assigned FDRs (default: ./FDRref.txt)'
            print '-t, --fdr:        cutoff of FDR (default: 0.1)'
            sys.exit(1)
        elif opt in ('-f', '--file'):
            filename = arg
        elif opt in ('-c', '--chr'):
            chrom = int(arg)
        elif opt in ('-r', '--reference'):
            refile = arg
        elif opt in ('-t', '--fdr'):
            fdr = float(arg)
    
    if refile == './FDRref.txt':
        print 'Note that FDR filtering is based on the default setting of ExIBD (i.e.,'
        print '1) ten independent runs of Beagle fastIBD with the fastIBD score of 1e-10;'
        print '2) five independent runs of Beagle IBD with ibd2nonibd=0.01 and nonibd2ibd=0.0001).'
    
    # Annotate marker information and sample information
    ribd = pd.read_table(filename+'.ribd', header = None, names = ['ID1', 'ID2', 'start', 'end'])
    marker = pd.read_table(filename+'.info', header = None, names = ['pos', 'gpos', 'ref', 'alt'], usecols = ['pos', 'gpos'])
    sample = pd.read_table(filename+'.rpop', header = None, names = ['ID', 'refPop'])
    
    ribd['chr'] = chrom
    ribd['start_pos'] = ribd.start.map(lambda x: marker.loc[x, 'pos'])
    ribd['end_pos'] = ribd.end.map(lambda x: marker.loc[x, 'pos'])
    ribd['start_gpos'] = ribd.start.map(lambda x: marker.loc[x, 'gpos'])
    ribd['end_gpos'] = ribd.end.map(lambda x: marker.loc[x, 'gpos'])
    ribd.start_pos = ribd.start_pos.astype(int)
    ribd.end_pos = ribd.end_pos.astype(int)
    
    ribd = ribd.merge(sample, how='left', left_on='ID1', right_on='ID', left_index=True).rename(columns = {'refPop': 'refPop1'}).drop(['ID'], axis = 1)
    ribd = ribd.merge(sample, how='left', left_on='ID2', right_on='ID', left_index=True).rename(columns = {'refPop': 'refPop2'}).drop(['ID'], axis = 1)
    
    # FDR filtering
    reference = pd.read_table(refile)
    reference = reference[reference.chr == chrom]
    
    ribd['center'] = 0.5 * (ribd.start_gpos + ribd.end_gpos)
    ribd['IBDsize'] = (ribd.end_gpos - ribd.start_gpos).apply(IBDType)
    ribd['fdr'] = ribd.apply(lambda x: getFDR(x.center, x.refPop1, x.refPop2, x.IBDsize, reference), axis = 1)
    ribd = ribd.loc[ribd.fdr < fdr, ['ID1', 'ID2', 'chr', 'start_pos', 'end_pos', 'start_gpos', 'end_gpos']]
    ribd.to_csv(filename+'.exIBD', sep = '\t', header = False, index = False)
      
if __name__ == '__main__':
    main(sys.argv[1:])
         
