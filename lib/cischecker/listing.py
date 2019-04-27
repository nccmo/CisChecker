#! /usr/bin/env python2

import pandas as pd
import sys
import itertools

def find_multiple_muts(maf_file, output_file, gene, sample, vaf):

    hout = open(output_file, 'w')

    # skip comment line(s)
    with open(maf_file) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    df = pd.read_table(maf_file, sep='\t', skiprows=skip_rows)

    # select gene
    if gene != "":
        gene_list = gene.split(':')
        df = df[(df['Hugo_Symbol'].isin(gene_list))]    

    # select sample
    if sample != "":
        sample_list = sample.split(':')
        df = df[(df['Tumor_Sample_Barcode'].isin(sample_list))]

    #rename
    df = df.rename(columns={'Hugo Symbol': 'Hugo_Symbol'})
    df = df.rename(columns={'Chr': 'Chromosome'})
    df = df.rename(columns={'HGVSp Short': 'HGVSp_Short'})
    df = df.rename(columns={'amino_acid_change': 'HGVSp_Short'})
    df = df.rename(columns={'aa_change': 'HGVSp_Short'})
    df = df.rename(columns={'aachange': 'HGVSp_Short'})
    df = df.rename(columns={'Protein_Change': 'HGVSp_Short'})
    df = df.rename(columns={'Protein Change': 'HGVSp_Short'})
    df = df.rename(columns={'Start_Position': 'Start_position'})
    df = df.rename(columns={'Start position': 'Start_position'})
    df = df.rename(columns={'Start Position': 'Start_position'})
    df = df.rename(columns={'End_Position': 'End_position'})
    df = df.rename(columns={'End position': 'End_position'})
    df = df.rename(columns={'End Position': 'End_position'})
    df = df.rename(columns={'Tumor Sample Barcode': 'Tumor_Sample_Barcode'})
    df = df.rename(columns={'Reference Allele': 'Reference_Allele'})
    df = df.rename(columns={'Tumor Seq Allele': 'Tumor_Seq_Allele'})
    df = df.rename(columns={'Tumor Seq Allele2': 'Tumor_Seq_Allele2'})

    if 'Tumor_Seq_Allele2' in df.columns.values:
        df = df.rename(columns={'Tumor_Seq_Allele2': 'Tumor_Seq_Allele'})
    elif 'Tumor_Seq_Allele1' in df.columns.values:
        df = df.rename(columns={'Tumor_Seq_Allele1': 'Tumor_Seq_Allele'})

    if vaf:
        df.loc[:,'VAF'] = df['t_alt_count'] / (df['t_ref_count'] + df['t_alt_count'])

    #select samples with multiple mutations in a same gene
    df2 = df[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome']]
    df2.loc[:,'n'] = 1
    df3 = df2.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome']).sum().reset_index()
    multiple_samples = df3[(df3['n'] >= 2)]

    if vaf:
      print >> hout, '\t'.join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Pos1", "Ref1", "Alt1", "Protein_Change1", "Pos2", "Ref2", "Alt2", "Protein_Change2", "VAF1", "VAF2"])
    else:
      print >> hout, '\t'.join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Pos1", "Ref1", "Alt1", "Protein_Change1", "Pos2", "Ref2", "Alt2", "Protein_Change2"])

    # get all two combinations
    for i in range(len(multiple_samples)):
        multiple_mutation_file = df[(df['Hugo_Symbol'] == multiple_samples.iloc[i,0]) & (df['Tumor_Sample_Barcode'] == multiple_samples.iloc[i,1]) & (df['Chromosome'] == multiple_samples.iloc[i,2])]

        # get all combinations
        tmp_list = list(itertools.combinations(range(len(multiple_mutation_file)),2))

        for j in range(len(tmp_list)):

            llist = [str(multiple_samples.iloc[i,0]), str(multiple_samples.iloc[i,1]), str(multiple_samples.iloc[i,2])]
            pos1  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'Start_position']
            ref1  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'Reference_Allele']
            alt1  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'Tumor_Seq_Allele']
            prot1 = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'HGVSp_Short']
            pos2  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'Start_position']
            ref2  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'Reference_Allele']
            alt2  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'Tumor_Seq_Allele']
            prot2 = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'HGVSp_Short']

            llist.append(str(pos1))
            llist.append(str(ref1))
            llist.append(str(alt1))
            llist.append(str(prot1))
            llist.append(str(pos2))
            llist.append(str(ref2))
            llist.append(str(alt2))
            llist.append(str(prot2))

            if vaf:
              vaf1  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'VAF']
              llist.append(str(vaf1))
              vaf2  = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'VAF']
              llist.append(str(vaf2))

            print >> hout, '\t'.join(llist)