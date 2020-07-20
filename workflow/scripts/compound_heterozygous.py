from pysam import VariantFile
from collections import defaultdict
import numpy as np
import argparse
import tempfile, os

#from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("cand_mother_field")
parser.add_argument("cand_father_field")
args = parser.parse_args()

transcripts_mother = defaultdict(list)
transcripts_father = defaultdict(list)

p_transcript_compound = defaultdict(lambda: 0)

def Q(P):
    return -10 * np.log10(P)

def P(Q):
    return np.power(10,-Q/10)

def at_least_one(ps):
    return 1 - np.prod(1 - P(ps))

records = []

with VariantFile(args.vcf) as f:
    for record in f:
        records.append(record)
        p_mother = record.info[args.cand_mother_field]
        p_father = record.info[args.cand_father_field]
        transcripts = set(ann.split("|")[6] for ann in record.info["ANN"])
        for transcript in transcripts:
            transcripts_mother[transcript].extend(p_mother)
            transcripts_father[transcript].extend(p_father)

    for transcript in set(transcripts_mother) & set(transcripts_father):
        p_mother = np.array(transcripts_mother[transcript], dtype=np.float128)
        p_father = np.array(transcripts_father[transcript], dtype=np.float128)
        p = at_least_one(p_mother) * at_least_one(p_father)
        if p > 0:
            p_transcript_compound[transcript] = p
        else:
            p_transcript_compound[transcript] = np.nan

    f.header.info.add("PROB_COMPOUND_HETEROZYGOUS","A","Float","Posterior probability for compound heterozygous variant (PHRED)")
    print(f.header, end="")

    for record in records:
        probs = []
        new_ann = []
        for ann in record.info["ANN"]:
            transcript = ann.split("|")[6]
            p = p_transcript_compound[transcript]
            if p > 0:
                new_ann.append(f"{ann}|{p}")
                probs.append(p)
        record.info["ANN"] = new_ann
        if not new_ann:
            continue

        record.info["PROB_COMPOUND_HETEROZYGOUS"] = float(Q(max(probs)))
        print(record, end="")