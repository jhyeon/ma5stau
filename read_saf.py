import os, sys

#This code is to read SAF file, 'sum of weights'.

ver = sys.argv[1]

base_path = '/home/minerva1993/phenostudy/madanalysis_workshop2020/madanalysis5/tools/PAD/Output/SAF/list_MASS_CHIRAL_taueff/ATLAS_SUSY_2018_04_VER/Cutflows/SR.saf'

files = []
for mass in ['100','120_20','120_30','120_50','140','140_30','140_50','160','180_60','180_80','200_50','200_70','200_80','200_100','220_80','220_100','250_70','250_100','250_120','250_150','300_70','300_100','300_120','300_140','320','320_50','340_100','340_120','360_30','360_50','360_60','360_80','360_100','360_110','380_30','380_50','380_80']:
#for sr in ['SRlow', 'SRhigh']:
  #for mass in ['120', '280']:
  for sr in ['SRlow', 'SRhigh']:
    #for chiral in ['LH', 'RH']:
    for chiral in ['LH']:
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
