#!/usr/bin/python

import csv
import operator
import sys
import numpy as np

# dictionary mapping row[0] to a list of (int(row[1]), row) values
report_map = {}
report_map_CG_Minus = {}

with open("../../BS_Forward.fastq_bismark_pe.CX_report.txt", 'r', newline='') as report:
    reader = csv.reader(report, delimiter='\t' )
    for row in reader:
        if row[5] != "CG":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map.setdefault(key, []).append((value, line))
        
with open("../../../Reapeat_analysis/FIMO/fimo_out/fimo.gff", 'r', newline='') as fimo, \
     open("x_CG.txt", 'w', newline='') as fout:
    reader = csv.reader(fimo, delimiter='\t')
    writer = csv.writer(fout, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])
