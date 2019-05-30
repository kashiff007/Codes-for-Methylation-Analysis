#!/usr/bin/python


import sys
import re

f1 = open("../../Reapeat_analysis/FIMO/fimo_out/fimo.gff",'r')
f2 = open("3407_A_run459_ATGAGCAT_S9_L007_R1_001.fastq_bismark_pe.sam",'r')
f3 = open("../bigwig/1_CG_meth.bedgraph",'r')
f4 = open("../bigwig/1_CHG_meth.bedgraph",'r')
f5 = open("../bigwig/1_CHH_meth.bedgraph",'r')

f1.seek(0)
f2.seek(0)
f3.seek(0)
f4.seek(0)
f5.seek(0)

fOut = open("output_sample_x1.txt",'w')
for line1 in f1:
    
    if not line1.startswith("#"):
        line1 = line1.split('\t')
        s1 = int(line1[3])
        s2 = int(line1[4])
        s0 = str(line1[0])
        count = 0
        percent = 0
        lt = []
        for line2 in f2:
            line2 = line2.split("\t")
            t1 = int(line2[3])
            t2 = int(line2[3]) + 150
            
            if (s0 == str(line2[2])):
                if (t2 >= s1 >= t1):
                    percent = percent + ((t2-s1)/(s2-s1))*100
                    count = count + 1
                    if count >= 4:
                        try:        
                            maxi = int(percent)/int(count)
                        except ZeroDivisionError:
                            maxi = 0
                        CG = []    
                        for l3 in f3:
                            l3 = l3.split("\t")
                            if (s0 == str(l3[0])) and (s1 <= int(l3[1]) <= s2):
                                CG.append(float(l3[3]))
                        f3.seek(0)
                        try:        
                            CG_M = sum(CG)/len(CG)
                        except ZeroDivisionError:
                            CG_M = 0    
                        CHG = []
                        for l4 in f4:
                            l4 = l4.split("\t")
                            if (s0 == str(l4[0])) and (s1 <= int(l4[1]) <= s2):
                                CHG.append(float(l4[3]))
                        f4.seek(0)
                        try:        
                            CHG_M = sum(CHG)/len(CHG)
                        except ZeroDivisionError:
                            CHG_M = 0
                        CHH = []
                        for l5 in f5:
                            l5 = l5.split("\t")
                            if (s0 == str(l5[0])) and (s1 <= int(l5[1]) <= s2):
                                CHH.append(float(l5[3]))
                            
                        f5.seek(0)
                        try:        
                            CHH_M = sum(CHH)/len(CHH)
                        except ZeroDivisionError:
                            CHH_M = 0						
                        fOut.write(str(s1) +"\t"+ str(s2)+"\t"+ str(s0)+"\t"+ str(t1)+"\t"+ str(t2)+"\t"+ str(line2[2])+"\t"+str(count)+"\t"+str(maxi)+"\t"+str(CG_M)+"\t"+str(CHG_M)+"\t"+str(CHH_M)+"\n")
                elif (t2 >= s2 >= t1):
                    percent = percent + ((s2-t1)/(s2-s1))*100
                    count = count + 1
                    if count >= 4:
                        try:        
                            maxi = int(percent)/int(count)
                        except ZeroDivisionError:
                            maxi = 0
                        CG = []    
                        for l3 in f3:
                            l3 = l3.split("\t")
                            if (s0 == str(l3[0])) and (s1 <= int(l3[1]) <= s2):
                                CG.append(float(l3[3]))
                        f3.seek(0)
                        try:        
                            CG_M = sum(CG)/len(CG)
                        except ZeroDivisionError:
                            CG_M = 0
                        CHG = []
                        for l4 in f4:
                            l4 = l4.split("\t")
                            if (s0 == str(l4[0])) and (s1 <= int(l4[1]) <= s2):
                                CHG.append(float(l4[3]))
                        f4.seek(0)
                        try:        
                            CHG_M = sum(CHG)/len(CHG)
                        except ZeroDivisionError:
                            CHG_M = 0
                        CHH = []
                        for l5 in f5:
                            l5 = l5.split("\t")
                            if (s0 == str(l5[0])) and (s1 <= int(l5[1]) <= s2):
                                CHH.append(float(l5[3]))
                        f5.seek(0)
                        try:        
                            CHH_M = sum(CHH)/len(CHH)
                        except ZeroDivisionError:
                            CHH_M = 0						
                        fOut.write(str(s1) +"\t"+ str(s2)+"\t"+ str(s0)+"\t"+ str(t1)+"\t"+ str(t2)+"\t"+ str(line2[2])+"\t"+str(count)+"\t"+str(maxi)+"\t"+str(CG_M)+"\t"+str(CHG_M)+"\t"+str(CHH_M)+"\n")
                elif (t1 >= s1 and t2 <= s2):
                    percent = percent + ((t2-t1)/(s2-s1))*100
                    count = count + 1
                    if count  >= 4:
                        try:
                            maxi = int(percent)/int(count)
                        except ZeroDivisionError:
                            maxi = 0
                        CG = []
                        for l3 in f3:
                            l3 = l3.split("\t")
                            if (s0 == str(l3[0])) and (s1 <= int(l3[1]) <= s2):
                                CG.append(float(l3[3]))
                        f3.seek(0)
                        try:        
                            CG_M = sum(CG)/len(CG)
                        except ZeroDivisionError:
                            CG_M = 0
                        CHG = []
                        for l4 in f4:
                            l4 = l4.split("\t")
                            if (s0 == str(l4[0])) and (s1 <= int(l4[1]) <= s2):
                                CHG.append(float(l4[3]))
                        f4.seek(0)
                        try:        
                            CHG_M = sum(CHG)/len(CHG)
                        except ZeroDivisionError:
                            CHG_M = 0
                        CHH = []
                        for l5 in f5:
                            l5 = l5.split("\t")
                            if (s0 == str(l5[0])) and (s1 <= int(l5[1]) <= s2):
                                CHH.append(float(l5[3]))
                        f5.seek(0)
                        try:        
                            CHH_M = sum(CHH)/len(CHH)
                        except ZeroDivisionError:
                            CHH_M = 0						
                        fOut.write(str(s1) +"\t"+ str(s2)+"\t"+ str(s0)+"\t"+ str(t1)+"\t"+ str(t2)+"\t"+ str(line2[2])+"\t"+str(count)+"\t"+str(maxi)+"\t"+str(CG_M)+"\t"+str(CHG_M)+"\t"+str(CHH_M)+"\n")
        f2.seek(0)
fOut.close()          

                                CHG.append(float(l4[3]))
                        f4.seek(0)
                        try:        
                            CHG_M = sum(CHG)/len(CHG)
                        except ZeroDivisionError:
                            CHG_M = 0
                        CHH = []
                        for l5 in f5:
                            l5 = l5.split("\t")
                            if (s0 == str(l5[0])) and (s1 <= int(l5[1]) <= s2):
                                CHH.append(float(l5[3]))
                        f5.seek(0)
                        try:        
                            CHH_M = sum(CHH)/len(CHH)
                        except ZeroDivisionError:
                            CHH_M = 0                       
                        fOut.write(str(s1) +"\t"+ str(s2)+"\t"+ str(s0)+"\t"+ str(t1)+"\t"+ str(t2)+"\t"+ str(line2[2])+"\t"+str(count)+"\t"+str(maxi)+"\t"+str(CG_M)+"\t"+str(CHG_M)+"\t"+str(CHH_M)+"\n")
                elif (t1 >= s1 and t2 <= s2):
                    percent = percent + ((t2-t1)/(s2-s1))*100
                    count = count + 1
                    if count  >= 4:
                        try:
                            maxi = int(percent)/int(count)
                        except ZeroDivisionError:
                            maxi = 0
                        CG = []
                        for l3 in f3:
                            l3 = l3.split("\t")
                            if (s0 == str(l3[0])) and (s1 <= int(l3[1]) <= s2):
                                CG.append(float(l3[3]))
                        f3.seek(0)
                        try:        
                            CG_M = sum(CG)/len(CG)
                        except ZeroDivisionError:
                            CG_M = 0
                        CHG = []
                        for l4 in f4:
                            l4 = l4.split("\t")
                            if (s0 == str(l4[0])) and (s1 <= int(l4[1]) <= s2):
                                CHG.append(float(l4[3]))
                        f4.seek(0)
                        try:        
                            CHG_M = sum(CHG)/len(CHG)
                        except ZeroDivisionError:
                            CHG_M = 0
                        CHH = []
                        for l5 in f5:
                            l5 = l5.split("\t")
                            if (s0 == str(l5[0])) and (s1 <= int(l5[1]) <= s2):
                                CHH.append(float(l5[3]))
                        f5.seek(0)
                        try:        
                            CHH_M = sum(CHH)/len(CHH)
                        except ZeroDivisionError:
                            CHH_M = 0                       
                        fOut.write(str(s1) +"\t"+ str(s2)+"\t"+ str(s0)+"\t"+ str(t1)+"\t"+ str(t2)+"\t"+ str(line2[2])+"\t"+str(count)+"\t"+str(maxi)+"\t"+str(CG_M)+"\t"+str(CHG_M)+"\t"+str(CHH_M)+"\n")
        f2.seek(0)
fOut.close()          
