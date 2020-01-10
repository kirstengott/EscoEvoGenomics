#!/usr/bin/env python3

import re
import os
import sys


def main(argv):
    if len(argv) != 4:
        sys.stderr.write("Usage: %s <id_map> <nwk> <id1> <idx2>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
        
    key_f = sys.argv[1]
    nwk = sys.argv[2]

    nwk_f = open(nwk, 'r')

    id_conv = {}


    key_file = open(key_f, 'r') ## file holding the id conversion


    for line in key_file:
        line = line.rstrip().split()
        if line[0] == '':
            title_line = 1
            continue
        else:
            key = line[int(sys.argv[3])]
            value = line[int(sys.argv[4])]
            id_conv[key] = value


    line = nwk_f.readlines()

    if len(line) == 1:
        line = line[0]
        new_line = ''
        for key in id_conv:
            re_key = re.compile(key + ':')
            replace = id_conv[key] + ":"
            if new_line == '':
                new_line = re.sub(pattern = re_key, repl = replace, string = line)
            else:
                new_line = re.sub(pattern = re_key, repl = replace, string = new_line)

        print(new_line)

    else:
        message('Expected a one line tree file, got many, exit and fix code')
        sys.exit()

        
if __name__ == '__main__':
    main(sys.argv[1:])
