'''
 ExIBD_Candidate: part of the 1st step for ExIBD;
 to combine results from multiple runs of Beagle fastIBD
 
 07/25/2016 By Wenqing Fu
'''
import sys
import getopt
import pandas as pd

def main():
    numround = 10
    fibdscore = 1e-10
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hf:r:t:',['help', 'file=','round=','fastibdthreshold='])
    except getopt.GetoptError:   		
        print 'python ./ExIBD_Candidate.py -f <filename> -r <#round> -t <fibdthreshold>'
        print '-h, --help:           print this help and exit'
        print '-f, --file:           name for input files'
        print '-r, --round:          number of Beagle fastIBD rounds (default: 10)'
        print '-t, --fastibdthreshold:  cutoff of fastIBD score (default: 1e-10)'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print 'python ./ExIBD_Candidate.py -f <filename> -r <#round> -t <fibdthreshold>'
            print '-h, --help:           print this help and exit'
            print '-f, --file:           name for input files'
            print '-r, --round:          number of Beagle fastIBD rounds (default: 10)'
            print '-t, --fastibdthreshold:  cutoff of fastIBD score (default: 1e-10)'
            sys.exit(1)
        elif opt in ('-f', '--file'):
            filename = arg
        elif opt in ('-r', '--round'):
            numround = int(arg)
        elif opt in ('-t', '--fastibdthreshold'):
            fibdscore = float(arg)
    
    for k in range(numround):
        if k == 0:
            dat = pd.read_table(str(k)+'.'+filename+'.fibd.gz', compression='gzip', header=None, names = ['ID1', 'ID2', 'start', 'end', 'score'])
        else:
            subdat = pd.read_table(str(k)+'.'+filename+'.fibd.gz', compression='gzip', header=None, names = ['ID1', 'ID2', 'start', 'end', 'score'])
            dat = dat.append(subdat, ignore_index = True)
    del subdat
    
    dat.end = dat.end-1
    dat = dat[(dat.score <= fibdscore)]
    dat.sort_values(by = ['ID1', 'ID2', 'start', 'end'], inplace = True)
    dat.reset_index(drop = True, inplace = True)
    dat['sel'] = 1 
    
    Rnum = dat.shape[0]
    for i in range(Rnum):
        if dat.sel[i] == 1:
            for j in range(i+1, Rnum, 1):
                if dat.ID1[j] == dat.ID1[i] and dat.ID2[j] == dat.ID2[i]:
                    if dat.start[j] <= dat.end[i]:
                        if dat.score[j] < dat.score[i]:
                            dat.loc[i,'sel'] = 0
                            break
                        else:
                            dat.loc[j,'sel'] = 0
                else:
                    break
    
    dat = dat[dat.sel == 1].sort_values(by = ['start', 'end'])
    dat.to_csv(filename+'.fibd', sep = '\t', columns=['ID1', 'ID2', 'start', 'end', 'score'], header = False, index = False)
    dat['pair'] = dat.ID1.map(str) + ' ' + dat.ID2.map(str)
    
    datgr = dat.groupby(['start','end'])
    
    out = open(filename+'.candidate', 'w')
    for key, group in datgr:
        start, end = key
        out.write(str(start) + '\t' + str(end) + '\t' + str(len(group)) + '\t')
        for index, row in group.iterrows():
            out.write(row.pair + '\t')
        out.write('\n')
    out.close()
    return len(datgr)
    
if __name__ == '__main__':
    main()
         
