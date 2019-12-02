#!/usr/bin/env python3

import pysam
import os
import sys

genome = sys.argv[1]

fh    = pysam.FastxFile(genome)
f_out = os.path.basename(os.path.splitext(genome)[0]) + "_repeats.bed"

f_o = open(f_out, 'w')

for entry in fh:
    seq            = entry.sequence
    nuc_pos = 0
    sm_start = 0
    sm_end = 0
    for nuc in seq:
        nuc_pos += 1
    ## see if we should begin a softmask region
        if nuc.islower() and sm_start == 0:
            sm_start  = nuc_pos
            ## reset the end position until we hit regular sequence
        elif nuc.islower() and sm_start != 0:
            sm_end = nuc_pos
        else:
            if sm_start == 0:
                continue
            else:
                out_line = "{}\t{}\t{}\n".format(entry.name, sm_start, sm_end)
                f_o.write(out_line)
                sm_start = 0
                sm_end = 0
    ## if we end in a soft masked region, still print
    if sm_start != 0:
        out_line = "{}\t{}\t{}\n".format(entry.name, sm_start, sm_end)
        f_o.write(out_line)


