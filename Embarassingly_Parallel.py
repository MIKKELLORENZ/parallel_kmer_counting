#####################################################
######## WELCOME TO THE KMER FINDER PROGRAM #########
#####################################################

##################### IMPORTS #######################
import os, sys, time
from collections import Counter
from joblib import Parallel, delayed

##################### FUNCTIONS #####################
#### FASTA INDEXING FUNCTION ###
def indexfasta(filename):
    
    # Open file
    infile = open(filename, 'rb')
    
    # Set chunksize
    chunksize = 1024*1024
    filepos = 0
    headstart = list()
    headend = list()
    
    # Read chunck size
    while True:
        content = infile.read(chunksize)
        if len(content) == 0:
            break
        
        # Find Header Start
        chunkpos = 0
        while chunkpos != -1:
            chunkpos = content.find(b'>', chunkpos)
            if chunkpos != -1:
                headstart.append(chunkpos + filepos)
                chunkpos += 1
                
        # Find Header End
        for i in range(len(headend), len(headstart)):
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()
    
    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart)-1, 0, -1):
        if headend[i] == headend[i-1]:
            del headstart[i]
            del headend[i]
    headstart.append(filepos)
    fastaindex = list()
    for i in range(len(headend)):
        seq_start = headend[i]+1
        seq_end = headstart[i+1] - 1
        fastaindex.append((seq_start, seq_end, seq_end-seq_start))
    
    # Load Balancing        
    #fastaindex = sorted(fastaindex, key=lambda x: x[2], reverse=False)
    
    return fastaindex

#### ENTRY INDEXING FUNCTION ####
def indexsequence(seq):
    pointer = 0
    seqindex = list()
    
    while len(seq) > pointer:
        
        # Find start of seq    
        potenstart = [ seq.find(b'a', pointer), seq.find(b't', pointer), seq.find(b'c', pointer), seq.find(b'g', pointer)]
            
        realstart = min(potenstart)
        if realstart == -1:
            # happens rarely, so slow code is ok
            potenstart = [ i for i in potenstart if i > -1 ]
            if len(potenstart) == 0:
                break
            realstart = min(potenstart)
        realend = seq.find(b'N', realstart)
        if realend == -1:
            realend = len(seq)
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex


#### KMER FINDING FUNCITON ####
def find_kmers(fasta,idx):
    
    # Read sequence
    infile = open(fasta, 'rb')
    infile.seek(idx[0])     
    seq = infile.read(idx[1]-idx[0]+1).translate(transtable, b'\r\n\t ')
    infile.close()
    subdict = dict()

    
    # Index sequence
    seqindex = indexsequence(seq)
    
    # Slide through sequences and add to dictionary
    for start,stop in seqindex:
        for i in range(start, stop-kmer_len+1):
            kmer = seq[i:i+kmer_len]
            if kmer not in subdict:
                subdict[kmer] = 1
            else:
                subdict[kmer] += 1 

    return subdict

if __name__ == '__main__':
    
    # SET VARIABLES   
    file = "humantest.fsa"
    kmer_len = 5
    transtable = bytes.maketrans(b'ATCGMRYKVHDBWmrykvhdbxnsw', b'atcgNNNNNNNNNNNNNNNNNNNNN')
    n_worker = os.cpu_count()
    
    # INDEXING
    start = time.time()
    my_indexes = indexfasta(file)
    index_time = time.time()-start
     
    # FIND KMER
    start = time.time()
    results = Parallel(n_jobs=1)(delayed(find_kmers)(file,z) for z in my_indexes)
    search_time = time.time()-start
    
    # MERGE DICTS
    final_dict = dict()
    start = time.time() 
    for r in results:
        for merkey in r.keys():
            if merkey in final_dict:
                final_dict[merkey] += r[merkey]
            else:
                final_dict[merkey] = r[merkey]
    dict_sum_time = time.time()-start
    del results
    del r

    print("")
    print("RESULTS FOR KMER FINDER:")
    print("Indexing time:",round(index_time,3))
    print("Finding kmers time:",round(search_time,3))
    print("Merging dicts time:",round(dict_sum_time,3))
    print("TOTAL:",round(index_time + search_time + dict_sum_time,3))
   
