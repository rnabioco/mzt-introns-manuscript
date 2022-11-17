#! /usr/bin/env python3

import os
import sys
import argparse
import urllib.request
import urllib.error
from subprocess import call

""" Utility to download fastq and metadata files from ENA database using
just the project identifier (i.e SRX.... or PRJN....)

Uses either wget (slow) or aspera (fast)
"""

# mdata values to capture
mdata_vals = {"study_accession",
  "secondary_study_accession",
  "sample_accession",
  "experiment_accession",
  "run_accession",
  "tax_id",
  "scientific_name",
  "experiment_title",
  "study_title",
  "fastq_ftp",
  "fastq_aspera",
  "submitted_ftp",
  "submitted_format",
  "sra_ftp",
  "sample_title"}

def format_metadata(mdata):
  "format mdata file into dictionary, run_accession will be key"
  
  d = {} 
  hdr = mdata[0].rstrip().split("\t")
  
  required_cols = ["sample_title", "fastq_ftp", "fastq_aspera"]
  if not all([x in hdr for x in required_cols]):
      out_str = ",".join(required_cols)
      sys.exit(f"unable to find required metadata columns {out_str}")
  
  for line in mdata[1:]:
      vals = line.rstrip().split("\t")
      sample_dict = {}
      for idx,val in enumerate(vals):
          val_type = hdr[idx]
          if val_type == "run_accession":
              run_key = val
          else:
              sample_dict[val_type] = val
      d[run_key] = sample_dict        
  
  return d

    
def get_study_metadata(study_id, logfile, mdata_vals):
  
  baseurl = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
  fields_url = "&result=read_run&fields="
  fields = ",".join(list(mdata_vals))
  full_url = baseurl + study_id + fields_url + fields

  # get study info
  try: 
    fp = urllib.request.urlopen(full_url)
  except:
    pass
  
  try:
    new_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
    full_url = new_url + study_id + fields_url + fields
    fp = urllib.request.urlopen(full_url)
  except urllib.error.HTTPError as e:
      sys.exit("unable to retrieve study from ENA: error {}".format(e.code)) 
  except urllib.error.URLError as e:
      sys.exit("unable to retrieve study from ENA: error {}".format(e.reason)) 
    
  mdata = []
  for line in fp:
    line = line.decode("utf8")
    logfile.write(line)
    mdata.append(line)
  
  fp.close()
  
  if len(mdata) < 2:
      # ENA will return only headers with some input strings
      sys.exit("unable to retrieve study from ENA")

  mdata = format_metadata(mdata)

  return mdata

def download_files(metadata, fq_ids, dl_prog, dl_prog_path, ssh_key,
        output_dir):

  ascp_cmd = [
        dl_prog_path,
        "-QT",  
        "-l",
        "200m",
        "-P",
        "33001",
        "-i",
        ssh_key]

  wget_cmd = [
    dl_prog_path
  ]

  if fq_ids:
    filter_fqs = True
  else:
    filter_fqs = False

  for k,v in metadata.items():

     if filter_fqs:
       if k not in fq_ids:
         continue

     ftp_urls = v["fastq_ftp"].split(";")
     for i,furl in enumerate(ftp_urls):
         if furl.startswith('ftp://'):
             pass
         else:
             ftp_urls[i] =  "ftp://" + furl 
     
     ascp_urls = v["fastq_aspera"].split(";")
     
     if len(ftp_urls) == 2 and len(ascp_urls) == 2:
         libtype = "paired_end"
     elif len(ftp_urls) == 1 and len(ascp_urls) == 1:
         libtype = "single_end"
     else:
         sys.exit("unknown ftp or ascp urls in mdata")
     
     dl_cmds = []
     if dl_prog == 'ascp':
       fqs = [os.path.join(output_dir, os.path.basename(x)) for x in ascp_urls]
          
       for idx,fq in enumerate(fqs):
         dl_cmd = ascp_cmd + ["era-fasp@" + ascp_urls[idx], fq]
         dl_cmds.append(dl_cmd)

     elif dl_prog == 'wget':
       fqs = [os.path.join(output_dir, os.path.basename(x)) for x in ftp_urls]

       for idx,fq in enumerate(fqs):
         dl_cmd = wget_cmd + ["-O", fq, ftp_urls[idx]]
         dl_cmds.append(dl_cmd)

     else:
       sys.exit("unknown download program {}".format(dl_prog))

     
     for cmd in dl_cmds:
       print("downloading sample {} from accession {}".format(k, v["sample_title"]))
       print(" ".join(cmd))
       call(cmd)



def main():
  parser = argparse.ArgumentParser(description = """
    Utility to download fastqs from European Nucleotide Archive
    """)

  parser.add_argument('-s',
                      '--study',
                      help = """
                      Study accessions numbers(ERP, SRP, DRP, PRJ prefixes)
                      e.g. SRX017289 or PRJNA342320
                      """,
                      required = True)

  parser.add_argument('-i',
                      '--ids',
                      help = """
                      Specific fastqs to download, defaults to all fastqs
                      in study
                      """,
                      required = False, nargs = "+")
  
  parser.add_argument('-d',
                      '--downloader',
                      help ="""either 'ascp' or 'wget' (default) """,
                   required = False,
                   default = "wget")
  
  parser.add_argument('-p',
                      '--prog_path',
                      help ="""path to downloader binary, defaults to name of program,
                      i.e. 'ascp', 'wget' '""",
                   required = False)
  
  parser.add_argument('-k',
                      '--sshkey',
                      help ="""path to aspera openssh key file (i.e something like
                      ".aspera/connect/etc/asperaweb_id_dsa.openssh") """,
                   required = False)
  
  parser.add_argument('-o',
                      '--outputdir',
                      help ="""output directory to place fastqs, defaults
                      to '.' """,
                      default = ".",
                   required = False)
  parser.add_argument('-M',
                      '--metadata-vals',
                      help ="""
                      additional attributes to add to metadata table 
                      supply as space deliminated values""",
                      nargs = "+",
                   required = False)
  
  parser.add_argument('-m',
                      '--metadatafile',
                      help ="""
                      preparsed metadata file for downloading
                      """,
                   required = False)

  args=parser.parse_args()
  
  study_id = args.study
  fq_ids = args.ids
  dl_prog = args.downloader
  dl_prog_path = args.prog_path
  ssh_key = args.sshkey
  output_dir = args.outputdir 
   
  if dl_prog == "ascp" and ssh_key is None:
    sys.exit("-k required for using aspera")
  if dl_prog_path is None:
    dl_prog_path = dl_prog
  
  if output_dir:
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
  else:
    output_dir = "."
  
  logfile = os.path.join(output_dir, study_id + "_download_log.txt")

  log_fp = open(logfile, 'w')
  
  if args.metadata_vals:
      for i in args.metadata_vals:
          mdata_vals.add(i)

  if args.metadatafile:
      mdata = []
      with open(args.metadatafile) as f:
          for line in f:
              mdata.append(line)
  else:
      mdata = get_study_metadata(study_id, log_fp, mdata_vals)
  
  download_files(mdata,
          fq_ids,
          dl_prog, 
          dl_prog_path, 
          ssh_key,
          output_dir)

  log_fp.close()

if __name__ == '__main__': main()

