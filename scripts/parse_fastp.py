#!/usr/bin/env python3

import json
import sys
import re

file_in = sys.argv[1]
file_sub = re.sub('_report.json', '', file_in).split("_")[0]

with open(file_in) as json_file:
     data = json.load(json_file)
     data_sub = data["summary"]["after_filtering"]
     for i in data_sub:
          print(file_in, i, data_sub[i])
