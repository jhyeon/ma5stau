import os, sys

#This code is to read SAF file, 'sum of weights'.

ver = sys.argv[1]

base_path = '/home/minerva1993/phenostudy/madanalysis_workshop2020/madanalysis5/tools/PAD/Output/SAF/list_MASS_CHIRAL_taueff/ATLAS_SUSY_2018_04_VER/Cutflows/SR.saf'

files = []
for sr in ['SRlow', 'SRhigh']:
  for mass in ['120', '280']:
    for chiral in ['LH', 'RH']:
      path = base_path.replace('SR',sr).replace('MASS',mass).replace('CHIRAL', chiral).replace('VER', ver)
      files.append(path)

for path in files:
  print path
  with open(path) as f:
    lines = f.readlines()

    for line in lines:
      if 'sum of weights' in line and '^2' not in line:
        pos = line.find('e+')
        evt = float(line[:pos])
        exp = int(line[pos+2:pos+4])
        nevt = evt * 10**exp
        print evt* 10**exp
