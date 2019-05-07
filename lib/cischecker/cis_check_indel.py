#!/usr/bin/env python2

import pysam
import pandas as pd
import sys
import re


def cis_check_snv_snv(bam_file, output_file, chromosome, pos1, pos2, ref1, alt1, ref2, alt2, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    print >> hout, '\t'.join(["read_ID", "pos1", "pos2", "output1", "output2"])

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    # check
    if int(pos1) == int(pos2):
        print >> sys.stderr, "positions are same.", bam_file, chromosome, pos1, pos2
    
    left_read_dict = {}
    tabixErrorFlag1 = 0
    try:
        records_one = bamfile.pileup(str(chromosome), int(pos1)-1, int(pos1)+1, min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag1 = 1
    if tabixErrorFlag1 == 0:
        for pileupcolumn in records_one:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and (pileupcolumn.pos) == int(pos1)-1:
                    left_read_dict[pileupread.alignment.query_name] = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(ref1)]
    right_read_dict = {}
    tabixErrorFlag2 = 0
    try:
        records_two = bamfile.pileup(str(chromosome), int(pos2)-1, int(pos2)+1, min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag2 = 1
    if tabixErrorFlag2 == 0:
        for pileupcolumn in records_two:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and (int(pileupcolumn.pos) == int(pos2)-1) and (str(pileupread.alignment.query_name) in list(left_read_dict.keys())):
                    right_read_dict[pileupread.alignment.query_name] = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(ref2)]
    
    for key in (set(left_read_dict.keys()) & set(right_read_dict.keys())):
        print >> hout, key + '\t' + pos1 + '\t' + pos2 + '\t' + left_read_dict[key] + '\t' + right_read_dict[key]
    
    bamfile.close()
    hout.close()


def del_reads(bam_file, output_file, chromosome, pos, ref, alt, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    if alt != "-":
        del_len = len(ref) - len(alt)
    else:
        del_len = len(ref)

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    records_one = bamfile.pileup(str(chromosome), int(pos)-1, int(pos), min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)

    for pileupcolumn in records_one:
        if int(pileupcolumn.pos) >= int(pos)-2-3 and int(pileupcolumn.pos) <= int(pos)-2+3:
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel < 0 and pileupread.indel >= -del_len-1 and pileupread.indel <= -del_len+1:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "alt"
                else:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "ref"
    hout.close()


def del_reads_accurate(bam_file, output_file, chromosome, pos, ref, alt, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    if alt != "-":
        del_len = len(ref) - len(alt)
    else:
        del_len = len(ref)

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    records_one = bamfile.pileup(str(chromosome), int(pos)-1, int(pos), min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)

    for pileupcolumn in records_one:
        if int(pileupcolumn.pos) == int(pos)-2:
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel == -del_len:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "alt"
                else:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "ref"
    hout.close()


def ins_reads(bam_file, output_file, chromosome, pos, ref, alt, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    if ref != "-":
        ins_len = len(alt) - len(ref)
    else:
        ins_len = len(alt)

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    records_one = bamfile.pileup(str(chromosome), int(pos)-1, int(pos), min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)
    for pileupcolumn in records_one:
        if int(pileupcolumn.pos) >= int(pos)-1-3 and int(pileupcolumn.pos) <= int(pos)-1+3:
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel > 0 and pileupread.indel >= ins_len-1 and pileupread.indel <= ins_len+1:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "alt"
                else:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "ref"
    hout.close()


def ins_reads_accurate(bam_file, output_file, chromosome, pos, ref, alt, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    if ref != "-":
        ins_len = len(alt) - len(ref)
    else:
        ins_len = len(alt)

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    records_one = bamfile.pileup(str(chromosome), int(pos)-1, int(pos), min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)
    for pileupcolumn in records_one:
        if int(pileupcolumn.pos) == int(pos)-1:
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel == ins_len:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "alt"
                else:
                    print >> hout, str(pileupread.alignment.query_name) + '\t' + "ref"
    hout.close()


def snv_reads(bam_file, output_file, chromosome, pos, ref, alt, baseq, mapq, overlaps, stepper, orphans):

    hout = open(output_file, 'w')
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    records_one = bamfile.pileup(str(chromosome), int(pos)-1, int(pos), min_base_quality=baseq, min_mapping_quality=mapq, ignore_overlaps=overlaps, stepper=stepper, ignore_orphans=orphans)
    for pileupcolumn in records_one:
        for pileupread in pileupcolumn.pileups:
            if int(pileupcolumn.pos) == int(pos)-1:
                if pileupread.query_position != None and pileupread.alignment.query_sequence != None: 
                    tmp = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(ref)]
                    if str(tmp) == str(ref):
                        print >> hout, str(pileupread.alignment.query_name) + '\t' + "ref"
                    elif str(tmp) == str(alt):
                        print >> hout, str(pileupread.alignment.query_name) + '\t' + "alt"
    hout.close()


def merge_results(result1, result2, merged_result):

    hout = open(merged_result, "w")
    print >> hout, '\t'.join(["read_ID", "output"])

    dict1 = {}
    with open(result1, "r") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            dict1[F[0]] = F[1]

    dict2 = {}
    with open(result2, "r") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            dict2[F[0]] = F[1]

    for key in (set(dict1.keys()) & set(dict2.keys())):
        print >> hout, key + '\t' + dict1[key] + ':' + dict2[key]
    hout.close()


def combine_two_results(input_file1, input_file2, output_file):

    hout = open(output_file, 'w')
    print >> hout, '\t'.join(["read_ID", "pos1", "pos2", "merge"])

    dict2 = {}
    with open(input_file2, 'r') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            dict2[F[0]] = F[1]
    hIN.close()

    with open(input_file1, 'r') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            print >> hout, '\t'.join(F) + '\t' + dict2[F[0]] + '\t' + ":".join([F[1], dict2[F[0]]])
    hIN.close()
