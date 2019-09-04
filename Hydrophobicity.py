#Modules Needed
import Bio
import numpy as np 
import matplotlib.pyplot as plt
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import matplotlib.pyplot as plt 
import os
import argparse

#The hydropathy scale used
kyle_doolittle = {'A' : 1.8, 'C' : 2.5, 'D' : -3.5, 'E' : -3.5, 'F' : 2.8, 'G' : -0.4, 
                  'H' : -3.2, 'I' : 4.5, 'K' : -3.90, 'L' : 3.8, 'M' : 1.9, 'N' : -3.5,
                  'P' : -1.6, 'Q' : -3.5, 'R' : -4.5, 'S' : -0.8, 'T': -0.7, 'V' : 4.2, 
                  'W' : -0.9, 'Y' : -1.3}

#Various Methods needed
def filter_profiles(profile): #Filters the length of the area of interest
    new_profile = []
    for x in profile:
        if len(x) > args.size:
            new_profile.append(x)
    return new_profile

def generate_region_average(profile): #Generates an average based on the area of interest
    avg_profile = []
    for x in profile:
        yvalues = []
        for y in x[1:]:
            yvalues.append(y[1])
        avg_profile.append([x[0],(sum(yvalues))/(len(yvalues)),[x[1][0],x[(len(x) - 1)][0]]])
    return avg_profile

##([minx1,maxx1],[minx1,maxx2])
def majority_overlap(xrange1,xrange2): #Catches any potential overlaps to keep stuff from repeating
    xray1 = range(xrange1[0],xrange1[1] + 1)
    xray2 = range(xrange2[0],xrange2[1] + 1)
    if xray1 in xray2:
        return True
    elif xray2 in xray1:
        return True
    else:
        lst3 = [value for value in xray1 if value in xray2] 
        if (len(lst3) > args.overlap):
            return True
        else:
            return False

#Sorting Everything
def mysort(col1,col2):
    col1_c = col1
    col2_c = col2
    ncol = []
    while len(ncol) != (len(col1) + len(col2)):
        if len(col1_c) == 0:
            for x in col2_c:
                ncol.append(x)
        elif len(col2_c) == 0:
            for x in col1_c:
                ncol.append(x)
        else:
            if col1_c[0][2][0] < col2_c[0][2][0]:
                ncol.append(col1_c[0])
                col1_c = col1_c[1:]
            else:
                ncol.append(col2_c[0])
                col2_c = col2_c[1:]
    return ncol

def compilation(all_trends): #Analyzes the collection of regions and their averages to produce a min and max for each region
    finished_analysis = []
    region_analyzed = []
    one_trend = []
    count = 0
    while len(all_trends) != count:
        if len(one_trend) == 0:
            region_analyzed = all_trends[0][2]
            one_trend.append(all_trends[count][1])
            count = count + 1
        #Makes sure that a region only contains averages that are either all positive or negative
        elif majority_overlap(all_trends[count][2],region_analyzed) and (all_trends[count][1] > 0) and (one_trend[0] > 0):
            one_trend.append(all_trends[count][1])
            count = count + 1
        elif majority_overlap(all_trends[count][2],region_analyzed) and (all_trends[count][1] < 0) and (one_trend[0] < 0):
            one_trend.append(all_trends[count][1])
            count = count + 1
        else:
            all_trends = all_trends[count:]
            count = 0
            finished_analysis.append([region_analyzed,min(one_trend),max(one_trend)])
            one_trend = []
            region_analyzed = []
    return finished_analysis


def create_cors(h_plot,window_mid): #Sets coordinates in a format to work with
    xcor = window_mid
    cors = []
    for y in h_plot:
        cors.append([xcor,y])
        xcor = xcor + 1
    return cors

def analyze_profile(cors,profile=[]): 
    #Creates a "profile" for each sequence 
    #consisting of the range of positions, values, and hydropathy
    if len(cors) == 0:
        return profile
    a_profile = [cors[0]]
    sign = a_profile[0][1] > 0
    a_profile.insert(0,sign)
    for x in cors[1:]:
        if (sign and x[1] >= 0) or (not sign and x[1] <= 0):
            a_profile.append(x)
        else:
            profile.append(a_profile)
            break
    if len(cors) == (len(a_profile) - 1):
        profile.append(a_profile)
        return profile
    return analyze_profile(cors[(len(a_profile) - 1):],profile)

def Nmaxelements(list1, N): 
    final_list = [] 
  
    for i in range(0, N):  
        max1 = [0,0]
          
        for j in range(len(list1)):      
            if list1[j][1] > max1[1]: 
                max1 = list1[j]; 
                  
        list1.remove(max1); 
        final_list.append(max1) 
          
    return final_list


seq_record = []
seq_names = []
seq_analysis = []

parser = argparse.ArgumentParser(description='HydroPathy Plots')
parser.add_argument('input_reference',type=str,help='Input reference FASTA file')
parser.add_argument('input_test',type=str,help='Input test FASTA file')
parser.add_argument('output',type=str,help='Output FASTA file')
parser.add_argument('-w','--size',type=int,default=7,help='Filtering out the regions lengths to be analyzed')
parser.add_argument('-e','--error',type=float,default=.1,help='The error range')
parser.add_argument('-n','--output_size',type=int,default=5,help='The number of sequences desired')
parser.add_argument('-v','--view_more',type=bool,default=False,help='Toggles additional information about matched areas')
parser.add_argument('-o','--overlap',type=int,default=5,help='The amount of positions defined as overlapping')

args = parser.parse_args()

#Parsing the preestablished sequences
for record in SeqIO.parse(args.input_reference,'fasta'):
    seq_record.append(str(record.seq))
    seq_names.append(str(record.id))

    
for index in range(len(seq_names)):
    seq = ProteinAnalysis(seq_record[index])
    test = seq.protein_scale(kyle_doolittle,7)
    seq_analysis.append(test)


all_trends = []
s_regions = []

#Generates Average for preestablished data
for index in range(len(seq_analysis)):
    a = generate_region_average(filter_profiles(analyze_profile(create_cors(seq_analysis[index],4),[])))
    s_regions.append(a)

for index in range(len(s_regions)):
    all_trends = mysort(s_regions[index],all_trends)

#Lists used to record the data   
new_seq_record = []
new_seq_names = []
new_seq_analysis = []
new_seq_descr = []
new_s_regions = []
all_correlation = []

#parsing for the test sequences
for record in SeqIO.parse(args.input_test,'fasta'):
    new_seq_record.append(str(record.seq))
    new_seq_names.append(str(record.id))
    new_seq_descr.append(str(record.description))
    
for index in range(len(new_seq_names)):
    seq = ProteinAnalysis(new_seq_record[index])
    test = seq.protein_scale(kyle_doolittle,7)
    new_seq_analysis.append(test)
    
for index in range(len(new_seq_analysis)):
    a = generate_region_average(filter_profiles(analyze_profile(create_cors(new_seq_analysis[index],4),[])))
    new_s_regions.append(a)

compared = compilation(all_trends)
for sequence in new_s_regions:
    #Gets each indiviual sequence test data
    seq_correlation = []
    for input_region in sequence:
        #Tries to match it with a region in the preestablished data
        for profile_region in compared:
            if profile_region[0][0] > input_region[2][1]:
                seq_correlation.append('Not Matched')
                break
            if majority_overlap(input_region[2],profile_region[0]) and (input_region[0] == (profile_region[1] > 0)) and (input_region[1] <= (profile_region[2] + args.error)) and (input_region[1] >= (profile_region[1] - args.error )):
                seq_correlation.append('Match')
                break
    all_correlation.append(seq_correlation)
            
seq_correlations = []

#Counts the correlations
for index in range(len(new_seq_names)):
    print(new_seq_names[index])
    seq_correlation = 0
    for index2 in range(len(all_correlation[index])):
        if (all_correlation[index][index2] == 'Match'):
            seq_correlation = seq_correlation + 1
            #Prints out specific region correlations
            if args.view_more:
                print('Match at windows ' + str(new_s_regions[index][index2][2]))
    print('Overall correlation = ' + str(seq_correlation/len(new_s_regions[index])))
    seq_correlations.append(([index,seq_correlation/len(all_correlation[index])]))

#Filters out the best correlated sequences
seq_correlations = Nmaxelements(seq_correlations,args.output_size)
final_output = []
for index in range(len(seq_correlations)):
    #Writes them onto a FASTA file
   final_output.append(SeqRecord(Seq(new_seq_record[seq_correlations[index][0]], generic_protein), 
                                     id = new_seq_names[seq_correlations[index][0]], 
                                     description = new_seq_descr[seq_correlations[index][0]]))
        
SeqIO.write(final_output,args.output,'fasta')




