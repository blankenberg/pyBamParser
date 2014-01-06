#!/usr/bin/env python
import os
import sys
# hack to prefer local version over installed egg
if os.environ.get('PYTHONPATH'):
  newpath = [os.path.abspath(os.environ.get('PYTHONPATH'))]
  newpath.extend(sys.path)
  sys.path = newpath
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
from optparse import OptionParser

def main():

  FUNCTIONS = {
    'BAMRead.get_indels':BAMRead_get_indels,
    'BAMRead.indel_at':BAMRead_indel_at,
  }

  OPT_DEFAULTS = {'indels_file':'', 'int':0, 'bool':False}
  USAGE = "USAGE: %prog [options] function.to.test reads.bam"
  DESCRIPTION = """Run test on a given function and input BAM and print results.
  Give one of the following function names: """+', '.join(FUNCTIONS)
  EPILOG = """ """

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-i', '--indels-file', dest='indels_file',
    default=OPT_DEFAULTS.get('indels_file'),
    help="""A file containing indels to test for. Required for BAMRead.indel_at.
Format: One indel per line, 3 tab-separated columns: 1. chrom, 2. coordinate
(1-based), "I" or "D" or "ID" for insertion, deletion, or both.""")

  (options, arguments) = parser.parse_args()

  if len(arguments) == 2:
    (function, bamfilename) = arguments
  else:
    parser.print_help()
    fail('Error: Provide a function name and a BAM file.')

  if function not in FUNCTIONS:
    fail('Error: function "'+function+'" not supported. Please pick one from '
      +'the list: '+', '.join(FUNCTIONS))
  if not os.path.exists(bamfilename):
    fail('Error: cannot find BAM file "'+bamfilename+'"')

  bam_reader = Reader(bamfilename)

  FUNCTIONS[function](bam_reader, options)


def BAMRead_get_indels(bam_reader, options):
  for read in bam_reader:
    print "\t".join([read.get_read_name(),str(read.get_position()),read.get_sam_cigar()])
    (insertions, deletions) = read.get_indels()
    if insertions:
      print "\t"+str(insertions)
    if deletions:
      print "\t"+str(deletions)


def BAMRead_indel_at(bam_reader, options):
  if not options.indels_file:
    fail('Error: BAMRead.indel_at requires you to specify an indels file.')

  indels = []
  with open(options.indels_file) as indels_file:
    for line in indels_file:
      line = line.strip()
      if not line:
        continue
      fields = line.split()
      try:          #  coordinate       I or D
        indels.append((int(fields[1]), fields[2]))
      except IndexError, ValueError:
        fail('Error: Bad format in indels file "'+options.indels_file+'"')

  for read in bam_reader:
    read_qname = read.get_read_name()
    read_pos   = read.get_position()
    read_cigar = read.get_sam_cigar()
    read_end   = read.get_end_position()
    print "\t".join([read_qname, str(read_pos), read_cigar])
    for (indel_pos, indel_type) in indels:
      if read_pos <= indel_pos <= read_end:
        has_indel = read.indel_at(
          indel_pos,
          check_insertions=('I' in indel_type),
          check_deletions=('D' in indel_type)
        )
        print "\t".join([str(indel_pos), indel_type, str(has_indel)])


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
