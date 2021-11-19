#! /usr/bin/env python

import sys
from subprocess import call

in_fn=sys.argv[1]

ascp_cmd = [
        "/beevol/home/riemondy/.aspera/connect/bin/ascp",
        "-QT",
        "-l",
        "100m",
        "-P",
        "33001",
        "-i",
        "/beevol/home/riemondy/.aspera/connect/etc/asperaweb_id_dsa.openssh"]

with open(in_fn) as f:
    for line in f:
        if line.split()[0].startswith("study_accession"):
            continue
        fields = line.split("\t")
        
        time_point = float(fields[10].split("_")[2])

        if time_point < 9.5:
            fqs = fields[9]
            fq_1 = fqs.split(";")[0]
            fq_2 = fqs.split(";")[1]
            
            fq1_cmd = ascp_cmd + ["era-fasp@" + fq_1, "."]
            fq2_cmd = ascp_cmd + ["era-fasp@" + fq_2, "."]
            
            call(fq1_cmd)
            call(fq2_cmd)

