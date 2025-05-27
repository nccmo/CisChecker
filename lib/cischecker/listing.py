#! /usr/bin/env python3

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

    # rename columns
    df = df.rename(columns={
        'Hugo Symbol': 'Hugo_Symbol',
        'Chr': 'Chromosome',
        'HGVSp Short': 'HGVSp_Short',
        'amino_acid_change': 'HGVSp_Short',
        'aa_change': 'HGVSp_Short',
        'aachange': 'HGVSp_Short',
        'Protein_Change': 'HGVSp_Short',
        'Protein Change': 'HGVSp_Short',
        'Start_Position': 'Start_position',
        'Start position': 'Start_position',
        'Start Position': 'Start_position',
        'End_Position': 'End_position',
        'End position': 'End_position',
        'End Position': 'End_position',
        'Tumor Sample Barcode': 'Tumor_Sample_Barcode',
        'Reference Allele': 'Reference_Allele',
        'Tumor Seq Allele': 'Tumor_Seq_Allele',
        'Tumor Seq Allele2': 'Tumor_Seq_Allele2'
    })

    if 'Tumor_Seq_Allele2' in df.columns:
        df = df.rename(columns={'Tumor_Seq_Allele2': 'Tumor_Seq_Allele'})
    elif 'Tumor_Seq_Allele1' in df.columns:
        df = df.rename(columns={'Tumor_Seq_Allele1': 'Tumor_Seq_Allele'})

    if vaf:
        df['VAF'] = df['t_alt_count'] / (df['t_ref_count'] + df['t_alt_count'])

    # select samples with multiple mutations in the same gene
    df2 = df[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome']]
    df2['n'] = 1
    df3 = df2.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome']).sum().reset_index()
    multiple_samples = df3[df3['n'] >= 2]

    if vaf:
        print('\t'.join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Pos1", "Ref1", "Alt1", "Protein_Change1", "Pos2", "Ref2", "Alt2", "Protein_Change2", "VAF1", "VAF2"]), file=hout)
    else:
        print('\t'.join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Pos1", "Ref1", "Alt1", "Protein_Change1", "Pos2", "Ref2", "Alt2", "Protein_Change2"]), file=hout)

    # get all two combinations
    for i in range(len(multiple_samples)):
        multiple_mutation_file = df[(df['Hugo_Symbol'] == multiple_samples.iloc[i,0]) & (df['Tumor_Sample_Barcode'] == multiple_samples.iloc[i,1]) & (df['Chromosome'] == multiple_samples.iloc[i,2])]

        # get all combinations
        tmp_list = list(itertools.combinations(range(len(multiple_mutation_file)), 2))

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

            llist.extend([str(pos1), str(ref1), str(alt1), str(prot1), str(pos2), str(ref2), str(alt2), str(prot2)])

            if vaf:
                vaf1 = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][0]], 'VAF']
                vaf2 = multiple_mutation_file.at[multiple_mutation_file.index[tmp_list[j][1]], 'VAF']
                llist.extend([str(vaf1), str(vaf2)])

            print('\t'.join(llist), file=hout)

    hout.close()
