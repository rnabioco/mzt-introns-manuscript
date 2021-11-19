import sys
import os

preamble =  """
track RisslandData
superTrack on show
shortLabel RisslandData
longLabel Rissland Lab Data 

    track rissland_bigwigs
    parent RisslandData
    compositeTrack on show
    shortLabel RNA-seq 
    longLabel Ribo-zero RNA-seq
    type bigWig
    autoScale on
    windowingFunction mean
    visibility full
    maxHeightPixels 100:30:8
    subGroup1 tpt Timepoint {timepoints}
    subGroup2 gtype Genotype {gtypes}
    subGroup3 strand Strand pos=Positive neg=Negative
    dimensions dimX=tpt dimY=gtype dimA=strand
    sortOrder strand=- tpt=+ gtype=-
"""


track =  """
        track {tid}
        parent rissland_bigwigs on
        shortLabel {shortlabel} 
        bigDataUrl {url_str} 
        longLabel  {longlabel}
        subGroups gtype={gtype} tpt={timepoint} strand={strand}
        maxHeightPixels 30:30:10
        color {col}
        type bigWig
"""

Spectral_palette = [
"213,62,79",
"252,141,89",
"254,224,139",
"255,255,191",
"230,245,152",
"153,213,148",
"50,136,189"
]

colors = [Spectral_palette[6], Spectral_palette[0]]

id_col = 5
tpt_col = 8
sample_col = 11


id_dict = {}
tpts = set()
samples = set()
ids = []
with open("mdata.txt") as f:
    for line in f:
       fields = line.rstrip().split("\t")
       id = fields[id_col]
       tpt = fields[tpt_col]
       tpt = tpt[:6]
       tpt = tpt.replace(" ", "_")
       sample = fields[sample_col]
       id_dict[id] = [tpt, sample]


for val in id_dict.values():
  tpts.add(val[0])
  samples.add(val[1])

tpts = [x + "=" + x for x in tpts]
samples = [x + "=" + x for x in samples]

print(preamble.format(timepoints = " ".join(tpts), 
                      gtypes = " ".join(samples)))

url = "http://amc-sandbox.ucdenver.edu/User33/rissland/hub/dm6/rissland/{fn}"

suffixes = ["_fwd.bw", "_rev_neg.bw"]
strands = ["pos", "neg"]

for key,vals in id_dict.items():  
  for idx,bw in enumerate(suffixes):
    tid=key + bw
    shortlabel=vals[0] + ":" + vals[1] + ":" + strands[idx]
    longlabel = shortlabel
    url_str=url.format(fn=tid)
    timepoint = vals[0]
    gtype = vals[1]
    strand = strands[idx]
    
    track_str = track.format(tid =tid,
      shortlabel=shortlabel,
      longlabel=longlabel,
      url_str=url_str,
      timepoint=timepoint,
      gtype=gtype,
      strand=strand,
      col=colors[idx])
    print(track_str)


