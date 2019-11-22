#!/usr/bin/env python3

import json
import sys
import re

file_in = sys.argv[1]
file_sub = re.sub('_report.json', '', file_in)

with open(file_in) as json_file:
     data = json.load(json_file)
    # print(sys.argv[1], data["summary"]["after_filtering"]["read1_mean_length"])
     string = "~/Analysis/genomescope/genomescope.R {}.trimmed.histo 31 {} genomescope/{}".format(file_sub,
                                                                                           data["summary"]["after_filtering"]["read1_mean_length"],
                                                                                           file_sub)
     print(string)
