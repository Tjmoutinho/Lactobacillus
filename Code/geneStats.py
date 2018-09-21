#!/usr/bin/python2.7
'''USAGE: geneStats.py gene_fasta
This script calculates various statistics about the provided fasta file.
'''
import sys
import math

# This function reads in fasta file, appends the length of each sequence to a list and sorts them
def read_lengths(fasta):

        len_lst = []
        seq = ''
        firstLine = 'placeholder'
        while firstLine[0] != '>':
                firstLine = fasta.readline()

        for line in fasta:
                if line[0] == '>': 
                        len_lst.append(len(seq))
                        seq = ''
                else:
                        seq += line.strip()

        len_lst.append(len(seq))
        len_lst.sort()

        return(len_lst)

# Function t0 calculate standard deviation for a list of numbers
def standDev(values):
        x_mean = sum(values) / len(values)
        sd_list = []
        for x in values:
                y = (x - x_mean) ** 2
                sd_list.append(y)
        y_mean = sum(sd_list) / len(sd_list)   
        sd = math.sqrt(y_mean)
        return(sd)


# This function calculates and returns all the printed statistics.
def calc_stats(lengths):

        shortest = lengths[0]
        longest = lengths[-1]
        total_seqs = len(lengths) # Total number of sequences
        len_sum = sum(lengths) # Total number of residues
        mid_pos = int(round(total_seqs/2))

        median_len = lengths[mid_pos] # Median sequence length
        q1 = lengths[0:mid_pos][int(len(lengths[0:mid_pos])/2)]
        q3 = lengths[mid_pos:-1][int(len(lengths[mid_pos:-1])/2)]
        iqr = q3 - q1

        mean_len = round((sum(lengths) / len(lengths)), 2)
        sd = round(standDev(lengths), 2)

        return(total_seqs, len_sum, median_len, q1, q3, iqr, mean_len, sd, shortest, longest)

#----------------------------------------------------------------------------------------#

with open(sys.argv[1], 'r') as contigs:
        gene_lengths = read_lengths(open(sys.argv[1], 'r'))

if len(gene_lengths) < 5:
        print('\nToo few contigs in ' + str(sys.argv[1]) + ' to calculate useful statistics. Exiting.\n')
        sys.exit()

stat_list = list(calc_stats(gene_lengths))

output_string = """
Fasta:          {filename}
Genes:          {genes}
Residues:       {residues}
Median:         {medianLen}
Q1:             {q1Len}
Q3:             {q3Len}
IQR:            {iqrLen}
Mean:           {meanLen}
Std:            {sdLen}
Shortest:       {shortLen}
Longest:        {longLen}
""".format(filename = str(sys.argv[1]).split('/')[-1],  
        genes = stat_list[0], 
        residues = stat_list[1], 
        medianLen = stat_list[2], 
        q1Len = stat_list[3], 
        q3Len = stat_list[4],
        iqrLen = stat_list[5],
        meanLen = stat_list[6],
        sdLen = stat_list[7],
        shortLen = stat_list[8], 
        longLen = stat_list[9])

print(output_string)
