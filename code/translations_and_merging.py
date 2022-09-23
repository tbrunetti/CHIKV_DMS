import Bio
from Bio import SeqIO
import sys
import os
sys.path.append(os.path.join(os.getcwd(), 'modules'))
from NucleicAcids import *
from Codon import *
from funcsForRefs import *
from typing import *
import pandas

codon_library = generate_codon_reference()

refSeq = ''

with open("pCHIKV_AF15561.fasta", 'r') as seq:
    for line in seq:
        if line[0] != ">":
            refSeq = line

# starts at base number 9820 but it is 1-indexed so change to 0-indexed
# total length of seq should be 1219 base pairs
# codon actually starts at 9821, 1-indexed (GAC), the amino acid sequence starting as DVESN

refSeq = refSeq[9819:11038] # truncate to region of interest
refInfo = {}
codonPosStart = 2 # 0-indexed, so 2 here means last position in codon
for pos, nt in enumerate(refSeq, start = 9820):
    refInfo[str(pos)] = { 'base' : nt, 'codonPos' : codonPosStart }   
    if codonPosStart == 2:
        codonPosStart = 0
    else:
        codonPosStart += 1



    
# required since sinple produces files with unequal columns per line
variantInfo = []
maxPad = 0
with open("mutDNA_S15_trimmed_sorted_9820_11038.variants", 'r') as unpadded:
    for line in unpadded:
        if len(line.strip().split('\t')) > maxPad:
            maxPad = len(line.strip().split('\t'))
        variantInfo.append(line.strip().split('\t'))

# pad list to match match columns
for var in variantInfo:
    if len(var) < maxPad:
        reps = maxPad - len(var)
        var.extend([''] * reps)

# pad and add unique names to header
consistentHeader = ['chrom', 'pos']
reps = int((maxPad - 2)/4)

for i in range(0, reps):
    consistentHeader.extend(['variant_{}'.format(i), 'total_reads_{}'.format(i), 'avg_read_qual_{}'.format(i), 'posterior_probability_{}'.format(i)])


# convert to df
variantInfoDf = pandas.DataFrame(variantInfo, columns = consistentHeader)


# merge reference info into df
tmp = pandas.DataFrame.transpose(pandas.DataFrame.from_dict(refInfo))
tmp['pos'] = tmp.index # since index is the rownames which are the base pair positions
results = pandas.DataFrame.merge(variantInfoDf, tmp, how = 'left', on = 'pos')


# rearrange pandas dataframe
results = results[['chrom', 'pos', 'base', 'codonPos'] + list(results.columns.values[2:(len(results.columns.values)-3)])]
results.insert(4, 'total_reads', 0)

# WARNING! Truncating first row since that should not be considered for analysis
results=results.loc[1:,:]

# insert column as specific location df.insert(loc, column, value)
totalReadsCols = [i for i in list(results.columns.values) if i.startswith('total_reads')]
for idx,row in results.iterrows():
    results.loc[idx, 'total_reads'] = sum(row[totalReadsCols].loc[row[totalReadsCols] != ''].astype(int)) # get sums of all reads for that variant

for idx,row in results.iterrows():
    varNum = 0
    while (varNum < reps) and row['variant_{}'.format(varNum)] != '':
        if 'af_{}'.format(varNum) in list(results.columns):
            results.loc[idx, 'af_{}'.format(varNum)] = int(row['total_reads_{}'.format(varNum)])/int(row['total_reads'])
            if row['variant_{}'.format(varNum)] in ['A', 'T', 'G', 'C']:    
                if refInfo[row['pos']]['codonPos'] == 0:
                    pos0 = int(row['pos'])
                    pos1 = int(row['pos']) + 1
                    pos2 = int(row['pos']) + 2
                    codonCheck = row['variant_{}'.format(varNum)] + refInfo[str(pos1)]['base'] + refInfo[str(pos2)]['base']
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 

                elif refInfo[row['pos']]['codonPos'] == 1:
                    pos0 = int(row['pos']) - 1
                    pos1 = int(row['pos'])
                    pos2 = int(row['pos']) + 1
                    codonCheck =  refInfo[str(pos0)]['base'] + row['variant_{}'.format(varNum)] + refInfo[str(pos2)]['base']
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 
                    
                elif refInfo[row['pos']]['codonPos'] == 2:
                    pos0 = int(row['pos']) - 2
                    pos1 = int(row['pos']) - 1
                    pos2 = int(row['pos'])
                    codonCheck =  refInfo[str(pos0)]['base'] + refInfo[str(pos1)]['base'] + row['variant_{}'.format(varNum)]
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 
                    
                else: 
                    print("codon position not found")
            else: 
                print('{} is not a DNA nucleotide, therefore skipping...'.format(row['variant_{}'.format(varNum)]))
        else:
            # determine where to insert new columns
            insertLoc = list(results.columns).index('variant_{}'.format(varNum)) + 1
            # insert new columns
            results.insert(insertLoc, 'af_{}'.format(varNum), 0.0)
            results.insert(insertLoc + 1, 'aa_change_{}'.format(varNum), '')
            results.insert(insertLoc + 2, 'polarity_change_{}'.format(varNum), '')
            results.insert(insertLoc + 3, 'charge_change_{}'.format(varNum), '')
            results.insert(insertLoc + 4, 'hydropathy_change_{}'.format(varNum), '')
            results.insert(insertLoc + 5, 'chemical_change_{}'.format(varNum), '')
            results.insert(insertLoc + 6, 'hydrogen_donor_change_{}'.format(varNum), '')

            results.loc[idx, 'af_{}'.format(varNum)] = int(row['total_reads_{}'.format(varNum)])/int(row['total_reads'])
            
            if row['variant_{}'.format(varNum)] in ['A', 'T', 'G', 'C']:    
                if refInfo[row['pos']]['codonPos'] == 0:
                    pos0 = int(row['pos'])
                    pos1 = int(row['pos']) + 1
                    pos2 = int(row['pos']) + 2
                    codonCheck = row['variant_{}'.format(varNum)] + refInfo[str(pos1)]['base'] + refInfo[str(pos2)]['base']
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 
                    
                elif refInfo[row['pos']]['codonPos'] == 1:
                    pos0 = int(row['pos']) - 1
                    pos1 = int(row['pos'])
                    pos2 = int(row['pos']) + 1
                    codonCheck =  refInfo[str(pos0)]['base'] + row['variant_{}'.format(varNum)] + refInfo[str(pos2)]['base']
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 
                    
                elif refInfo[row['pos']]['codonPos'] == 2:
                    pos0 = int(row['pos']) - 2
                    pos1 = int(row['pos']) - 1
                    pos2 = int(row['pos'])
                    codonCheck =  refInfo[str(pos0)]['base'] + refInfo[str(pos1)]['base'] + row['variant_{}'.format(varNum)]
                    newCodon = Dna(codonCheck).transcribe().sequence
                    results.loc[idx, 'aa_change_{}'.format(varNum)] = '{} ({})'.format(codon_library[newCodon].translate_shortname(), codon_library[newCodon].translate_symbol()) 
                    results.loc[idx, 'polarity_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_polarity()) 
                    results.loc[idx, 'charge_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_charge()) 
                    results.loc[idx, 'hydropathy_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_hydropathy()) 
                    results.loc[idx, 'chemical_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_chemical_class()) 
                    results.loc[idx, 'hydrogen_donor_change_{}'.format(varNum)] = '{}'.format(codon_library[newCodon].get_donor_status()) 
                    
                else: 
                    print("codon position not found")
            else: 
                print('{} is not a DNA nucleotide, therefore skipping...'.format(row['variant_{}'.format(varNum)]))
        varNum += 1


results.to_csv('test.txt', sep = '\t', index = False)


            
            
        

    # calc allele freq
    # get amino acid change for codon


