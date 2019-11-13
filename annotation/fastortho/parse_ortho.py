#!/usr/bin/env python3

import sys

ortho = open(sys.argv[1], 'r')

for line in ortho:
    line = line.strip().split(':')
    print(line[0])
