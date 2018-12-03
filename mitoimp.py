#! /usr/bin/env python
# -*- coding:utf-8 -*-

import time
import cPickle as pickle
import bz2
from collections import Counter
import multiprocessing as mp
import argparse
import numpy
import csv
import os.path
import multiprocessing
from subprocess import check_call
import distutils.spawn
import sys

usage = 'mitoimp.py -i input.fasta <optional>\nVersion: 1.0.0\n'
parser = argparse.ArgumentParser(usage=usage,description='MitoIMP:1.0.0')
parser.add_argument('-i',action='store',dest='fasta',help='Query sequence(.fasta):Query will be automatically aligned with MAFFT.')
parser.add_argument('-p',action='store',dest='panel',help='Panel sequences(.fasta)',default='ALL_panel')
parser.add_argument('-w',action='store',dest='window',help='window-size(default:16569)',default=16569)
parser.add_argument('-k',action='store',dest='k_num',help='K-number(default:5)',default=5)
parser.add_argument('-f',action='store',dest='freq',help='The threshold frequency to determine a genotype. (default:0.7)',default=0.7)
parser.add_argument('-no_aln', action='store_true',dest='non_align',help='Set a switch to non-alignment mode (default:Disable)',default=False)
parser.add_argument('-t',action='store',dest='threads',help='Multiprocessing numbers (default:The max number of CPU-threads available on your system)',default=-1)
parser.add_argument('-v',action='version',version='Version : 1.0.0')

args = parser.parse_args()

print "\n"
print "============================================================================"
print "                                                                            "
print "                             ◇ mitoimp ◇                                    "
print "                                                                            "
print "  　　  An imputation tool for low-coverage human mitochondrial genome       "
print "                   based on shared allele pairwise distances.               "
print "                                                                            "
print "                          > Version : 1.0.0 (beta)                          "
print "                          > Date : 30, Dec., 2018                           "
print "                          > Written by K.I.                                 "
print "                                                                            "
print "============================================================================"

#Allele-Sharing Distance (Malaspinas et al., 2014)
def ASD(s1, s2):
    k = len([1 for i,j in enumerate(s1) if j != 'N' and j != 'n' and s2[i] != 'N' and s2[i] != 'n'])
    if k != 0:
        d = len([1 for i,j in enumerate(s1) if j != s2[i] and j != 'N' and j != 'n' and s2[i] != 'N' and s2[i] != 'n'])
        return float(d)/k
    elif k == 0:
        return 0

#FASTA parser
class FastaParser(object):
    def __init__(self, header, sequence):
        self.header    = header
        self.sequence = sequence

def read_fasta(input_fasta):
    line = input_fasta.readline()
    if not line.startswith(">"):
        raise TypeError("Error:Input file is not a valid (single) FASTA: {0}".format(line))
    header = line[1:].rstrip()

    sequence_lines = []
    while 1:
        line = input_fasta.readline().rstrip()
        if line == "":
            break
        sequence_lines.append(line)
    sequence = "".join(sequence_lines)

    return FastaParser(header, sequence)

#Alignment func(requirement:MAFFT)
def align(seq,mafft_path,proc_num):
    load_refseq = open("./refseq")
    load_query = open(seq)
    #Load a reference sequence
    ref=read_fasta(load_refseq)
    #Load a query sequence
    query=read_fasta(load_query)
    #Output a multiple fasta
    tmp = open("{0}".format(os.path.splitext(os.path.abspath(seq))[0]+'_tmp'),"w")
    print >> tmp , '>{0}\n'.format(ref.header)+ref.sequence+'\n'+'>{0}\n'.format(query.header)+query.sequence
    tmp.close()
    load_refseq.close()
    load_query.close()

    #Perform an alignment
    print "Aligning an input file with MAFFT ...."
    cmd_mafft="{0} --localpair --maxiterate 1000 --thread {1} {2} > {3} 2> .log".format(mafft_path,proc_num,os.path.splitext(os.path.abspath(seq))[0]+'_tmp',os.path.splitext(os.path.abspath(seq))[0]+'_tmp_aligned')
    check_call(cmd_mafft,shell=True)
    print "Alignment complete !"

    #Load an aligned fasta
    aligned_list = open(os.path.splitext(os.path.abspath(seq))[0]+'_tmp_aligned').readlines()

    #Retrive aligned sequences
    for k,v in enumerate(aligned_list):
        if v[0] == '>' and v[0] != '>refseq\n':
            start_index = k
    ref_aligned = ''.join([i.strip().upper() for i in aligned_list[1:start_index]])
    query_aligned=''.join([i.strip().upper() for i in aligned_list[start_index+1:]])

    #Relocate the sequence position of rCRS (not including indels)
    ref_gaps=[k for k,v in enumerate(ref_aligned) if v == '-']
    if len(ref_gaps) != 0:
        query_aligned_list=list(query_aligned)
        for i in ref_gaps:
            query_aligned_list[i] = ''
            new_query=''.join(query_aligned_list)
    else :
        new_query = query_aligned

    query_gaps=[k for k,v in enumerate(new_query) if v == '-']

    if len(query_gaps) != 0:
        query_list2=list(new_query)
        for i in query_gaps:
            query_list2[i]='N'
        new_query=''.join(query_list2)

    os.remove("{0}".format(os.path.splitext(os.path.abspath(seq))[0]+"_tmp"))
    os.remove("{0}".format(os.path.splitext(os.path.abspath(seq))[0]+"_tmp_aligned"))

    return new_query

def imp_run(cov_fasta, proc_num,window_size,k_num,freq_threshold,non_align):

    #Check the enviroment
    if distutils.spawn.find_executable('mafft') == None:
        print "MAFFT could not find on your PATH environment. Please install MAFFT before running mitoimp."
        sys.exit()
    else:
        mafft_path = distutils.spawn.find_executable('mafft')

    sttime = time.time()

    with open('./panel.bz2','rb') as panel:
        print "ALL Haplogroup Panel (ver:0.1)"
        print "Loading Panel Data ...."
        panel_dict=pickle.loads(bz2.decompress(panel.read()))

    output_csv=os.path.splitext(os.path.abspath(cov_fasta))[0]+'_stats.csv'

    f_pre=open(output_csv,"w")
    header_list=['Input','Process','Window_size','K_number','Freq_threshold','Missing_sites_before_imputation','Missing_sites_after_imputation','Coverage_before','Coverage_after','Run_time']
    writer = csv.writer(f_pre, lineterminator='\n')
    writer.writerows([header_list])
    f_pre.close()

    output_fasta=os.path.splitext(os.path.abspath(cov_fasta))[0]+'_imputed.fasta'
    load_fasta = open(cov_fasta)
    fasta_record=read_fasta(load_fasta)
    con_seq_id=fasta_record.header
    load_fasta.close()

    if non_align == False:
        con_seq=align(seq=os.path.abspath(cov_fasta),mafft_path=mafft_path,proc_num=proc_num)
    else:
        con_seq=fasta_record.sequence
        if len(con_seq) != 16569:
            print "Error: The query sequence is not a valid length. Please switch off non-alignment mode."
            sys.exit()

    missing_sites=[]
    non_missing_sites=[]

    for k,v in enumerate(con_seq):
        if v == 'N':
            missing_sites.append(k)
        else:
            non_missing_sites.append(k)

    mitogenome_site=len(con_seq)

    print "Detected missing sites : {0} -> {1}% / mitogenome".format(len(missing_sites),round((float(len(missing_sites))/mitogenome_site)*100,2))

    # These sites are not used for the calculation.
    non_informative_sites = panel_dict['NonInfo']
    deletion_sites=list(set(missing_sites+non_informative_sites))
    deletion_sites.sort()

    #k-mer splitting
    def split_str(s, n):
        length = len(s)
        return [s[i:i+n] for i in range(0, length, n)]

    #Deletion sites are assign to N
    def trans_delete(s):
        s_list=list(s)
        for i in deletion_sites:
            s_list[i]='N'
        new_s="".join(s_list)
        return new_s

    #Divide each sequence of group panel into k-mer array based on the window size.
    k_mer_db2={}

    for k,v in panel_dict['Panel'].items():
        split_list=[]
        split_list=split_str(trans_delete(v['seq']),window_size)
        k_mer_db2[k]=split_list

    #Remove deletion sites from an input sequence.
    consensus_window2=split_str(trans_delete(con_seq),window_size)

    #Calculate the distance between k-mer and an input sequence for each window.
    distance_dict2={}
    missing_sites_window={}
    non_missing_sites_window={}
    for i in xrange(len(consensus_window2)):
        distance_dict2[i]={}
        missing_sites_window[i]=[]
        non_missing_sites_window[i]=[]

    proc=proc_num
    panel_count = len(panel_dict['Panel'])
    L = panel_count
    print "Calculating distances (processes:{0}) ...".format(proc)

    def mp_func(queue,p):
        # Divide ranges
        ini = L * p / proc
        fin = L * (p+1) / proc
        sub_dict={k : {} for k in k_mer_db2.keys()[ini:fin]}
        for k,v in k_mer_db2.items()[ini:fin]:
            for g,h in enumerate(v):
                sub_dict[k]={'ratio':1.0 - ASD(consensus_window2[g],k_mer_db2[k][g])}
        queue.put(sub_dict)

    #Make queues
    queue = mp.Queue()

    #Provide the specified number of processes
    ps=[mp.Process(target=mp_func, args=(queue, i)) for i in xrange(proc)]

    for p in ps:
        p.start()
    for i in range(proc):
        distance_dict2[0].update(queue.get())

    print "Complete imputation !"

    #K-nearest neighbor algorithm

    index_list=[i for i in range(len(trans_delete(con_seq)))]
    index_list_window = split_str(index_list,window_size)

    for i in missing_sites:
        for k,v in enumerate(index_list_window):
            if i in index_list_window[k] : missing_sites_window[k].append(i)

    for i in non_missing_sites:
        for k,v in enumerate(index_list_window):
            if i in index_list_window[k]:non_missing_sites_window[k].append(i)

    impute_list=list(con_seq)

    for k,v in missing_sites_window.items():
        seq_id=[]
        seq_id=[k2 for k2,v2 in sorted(distance_dict2[k].items(),key=lambda x:x[1],reverse=True)]
        seq_id=seq_id[0:k_num]
        for v3 in v:
            if float(Counter([panel_dict['Panel'][i]['seq'][v3] for i in seq_id]).most_common(1)[0][1])/len([panel_dict['Panel'][i]['seq'][v3] for i in seq_id]) > freq_threshold:
                impute_list[v3]=Counter([panel_dict['Panel'][i]['seq'][v3] for i in seq_id]).most_common(1)[0][0]
            else:
                impute_list[v3]='N'

    imputation_seq="".join(impute_list)
    imputation_seq_output=imputation_seq

    obj = open("{0}".format(output_fasta),"w")
    print >> obj , '>{0}_imputed\n'.format(con_seq_id)+imputation_seq_output
    obj.close()

    missing_sites_after=0
    for k,v in enumerate(imputation_seq_output):
        if v == 'N' or v == 'n' : missing_sites_after += 1

    #mitogenome coverage
    def coverage_cal(s):
        s1=s
        mis=Counter(s1)['N']+Counter(s1)['n']
        return int(round((float(len(s1)-mis)/mitogenome_site)*100,1))

    run_time=round((time.time() - sttime),3)

    stat_list=[os.path.basename(cov_fasta),proc_num,window_size,k_num,freq_threshold,len(missing_sites),missing_sites_after,coverage_cal(con_seq),coverage_cal(imputation_seq),run_time]
    csv_list=[stat_list]

    f=open(output_csv,"a")
    writer = csv.writer(f, lineterminator='\n')
    writer.writerows(csv_list)
    f.close()

if __name__ == "__main__":

    #Set a parameter for multiprocessing
    if args.threads != -1:
        t_num = args.threads
    else:
        t_num=multiprocessing.cpu_count()

    imp_run(cov_fasta=args.fasta,proc_num=t_num,window_size=args.window,k_num=args.k_num,freq_threshold=args.freq,non_align=args.non_align)
