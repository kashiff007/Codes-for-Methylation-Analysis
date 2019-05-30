#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio.Alphabet import generic_dna
import sys
import re

fOut = open("C_content_output.txt",'w')
# read names and postions from bed file
positions = defaultdict(list)
with open('output_sample_max.txt') as f:
    for line in f:
        line = line.split("\t")
        num, cov, CG, CHG, CHH = line[3], line[4], line[5], line[6], line[7]
        name,start,stop = line[2],line[0],line[1]
        positions[name].append((int(start), int(stop)))

# parse faste file and turn into dictionary
records = SeqIO.to_dict(SeqIO.parse(open('../../Genome_with_repeat_35791/Seq_with_repeat_35791.fasta'), 'fasta'))


short_seq_records = []
for name in positions:
    for (start, stop) in positions[name]:
        long_seq_record = records[name]
#        print (long_seq_record)
        long_seq = long_seq_record.seq
#        print (long_seq)
        alphabet = long_seq.alphabet
#        print (alphabet)
        short_seq = str(long_seq)[start-1:stop]
        short_seq = short_seq.lower()
        cg = (str(short_seq.count("cg")))
        chg = (int(short_seq.count("cag"))+int(short_seq.count("ctg"))+int(short_seq.count("ccg")))
        chh = (int(short_seq.count("caa"))+int(short_seq.count("cat"))+int(short_seq.count("cac"))+int(short_seq.count("cta"))+int(short_seq.count("ctt"))+int(short_seq.count("ctc"))+int(short_seq.count("cca"))+int(short_seq.count("cct"))+int(short_seq.count("ccc")))
        cg_r = (str(short_seq.count("cg")))
        chg_r = (int(short_seq.count("ctg"))+int(short_seq.count("cag"))+int(short_seq.count("cgg")))
        chh_r = (int(short_seq.count("ttg"))+int(short_seq.count("atg"))+int(short_seq.count("gtg"))+int(short_seq.count("tag"))+int(short_seq.count("aag"))+int(short_seq.count("gag"))+int(short_seq.count("tgg"))+int(short_seq.count("agg"))+int(short_seq.count("ggg")))
 
        fOut.write(str(start)+"_"+str(stop)+"_"+ str(name)+"\t"+str(cg) +"\t"+str(chg) +"\t"+str(chh)+"\t"+str(cg_r)+"\t"+str(chg_r) +"\t"+str(chh_r)+"\n")        
        short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name, description='')
#        print (short_seq_record)
        short_seq_records.append(short_seq_record)
#        print (short_seq_record)

    
fOut.close()

