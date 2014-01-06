#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue
echo -e "\tBAMRead.get_indels/cigar-varieties.bam:"
PYTHONPATH=$dirname/../lib $dirname/unit.test.py BAMRead.get_indels $dirname/cigar-varieties.bam | diff -s - $dirname/cigar-varieties.get_indels.out
echo -e "\tBAMRead.indel_at/cigar-varieties.bam:"
PYTHONPATH=$dirname/../lib $dirname/unit.test.py BAMRead.indel_at $dirname/cigar-varieties.bam -i $dirname/cigar-varieties.indel_at.tsv | diff -s - $dirname/cigar-varieties.indel_at.out
