#!/usr/bin/env python3

import pandas as pd
import subprocess, sys
import argparse
from . import listing
from . import cis_check_indel
import numexpr
import re

def is_valid_nuc(nuc):
    return bool(re.match('^[ATGC]+$', nuc))

def listing_main(args):
    listing.find_multiple_muts(args.MAF_file, args.output_file, args.gene, args.sample, args.vaf)

def cis_check_main(args):

    bam_dict = {}
    with open(args.sample_list, 'r') as hin:
        first_line = next(hin)
        skip_rows = 1 if first_line.startswith('#') else 0

    with open(args.sample_list, 'r') as hin:
        if skip_rows == 1:
            header = hin.readline()
        for line in hin:
            F = line.rstrip('\n').split('\t')
            bam_dict[F[0]] = F[1]

    with open(args.output_file, 'w') as hOUT:

        with open(args.multiple_mutation_list, 'r') as hIN:
            header = hIN.readline().rstrip('\n').split('\t')
            print('\t'.join(header) + '\t' + "Ref_Ref" + '\t' + "Alt_Alt" + '\t' + "Ref_Alt" + '\t' + "Alt_Ref" + '\t' + "Type1" + '\t' + "Type2", file=hOUT)

            for line in hIN:
                F = line.rstrip('\n').split('\t')

                # check if bam file exists
                if F[1] not in bam_dict:
                    print(f"bam file not found {F[1]}", file=sys.stderr)
                    print('\t'.join(F) + '\t' + "NO BAM" + '\t' + "NO BAM" + '\t' + "NO BAM" + '\t' + "NO BAM" + '\t' + "NO BAM" + '\t' + "NO BAM", file=hOUT)
                    continue

                first_snv_flag = (len(F[4]) == len(F[5])) and is_valid_nuc(F[4]) and is_valid_nuc(F[5])
                first_del_flag = (len(F[4]) > len(F[5])) or F[5] == "-"
                first_ins_flag = (len(F[4]) < len(F[5])) or F[4] == "-"

                second_snv_flag = (len(F[8]) == len(F[9])) and is_valid_nuc(F[8]) and is_valid_nuc(F[9])
                second_del_flag = (len(F[8]) > len(F[9])) or F[9] == "-"
                second_ins_flag = (len(F[8]) < len(F[9])) or F[8] == "-"

                first = second = ""

                if first_del_flag:
                    if args.mapping_difference:
                        cis_check_indel.del_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", F[2], F[3], F[4], F[5], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    else:
                        cis_check_indel.del_reads_accurate(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", F[2], F[3], F[4], F[5], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    first = "DEL"
                if first_ins_flag:
                    if args.mapping_difference:
                        cis_check_indel.ins_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", F[2], F[3], F[4], F[5], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    else:
                        cis_check_indel.ins_reads_accurate(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", F[2], F[3], F[4], F[5], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    first = "INS"
                if first_snv_flag:
                    cis_check_indel.snv_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", F[2], F[3], F[4], F[5], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    first = "SNV"

                if second_del_flag:
                    if args.mapping_difference:
                        cis_check_indel.del_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", F[2], F[7], F[8], F[9], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    else:
                        cis_check_indel.del_reads_accurate(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", F[2], F[7], F[8], F[9], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    second = "DEL"
                if second_ins_flag:
                    if args.mapping_difference:
                        cis_check_indel.ins_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", F[2], F[7], F[8], F[9], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    else:
                        cis_check_indel.ins_reads_accurate(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", F[2], F[7], F[8], F[9], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    second = "INS"
                if second_snv_flag:
                    cis_check_indel.snv_reads(bam_dict[F[1]], f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", F[2], F[7], F[8], F[9], args.baseq, args.mapq, args.overlaps, args.stepper, args.orphans)
                    second = "SNV"

                cis_check_indel.merge_results(f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp1.txt", 
                                              f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis_tmp2.txt", 
                                              f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis.txt")

                df = pd.read_table(f"{F[0]}_{F[1]}_{F[6]}_{F[10]}_combined_cis.txt", sep='\t')
                ref_num = len(df[df['output'] == "ref:ref"])
                cis_num = len(df[df['output'] == "alt:alt"])
                ref_alt_num = len(df[df['output'] == "ref:alt"])
                alt_ref_num = len(df[df['output'] == "alt:ref"])
                print('\t'.join(F) + '\t' + str(ref_num) + '\t' + str(cis_num) + '\t' + str(ref_alt_num) + '\t' + str(alt_ref_num) + '\t' + first + '\t' + second, file=hOUT)

    if not args.debug:
        subprocess.call("rm -f *combined_cis_tmp1.txt", shell=True)
        subprocess.call("rm -f *combined_cis_tmp2.txt", shell=True)
        subprocess.call("rm -f *combined_cis.txt", shell=True)
