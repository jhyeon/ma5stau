import os, sys

#This code is to read SAF file, 'sum of weights'.

ver = sys.argv[1]

base_path = '/home/minerva1993/phenostudy/madanalysis_workshop2020/madanalysis5/tools/PAD/Output/SAF/list_MASS_CHIRAL_taueff/ATLAS_SUSY_2018_04_VER/Histograms/histos.saf'

files = []
for mass in ['120', '280']:
  for chiral in ['LH', 'RH']:
    path = base_path.replace('MASS',mass).replace('CHIRAL', chiral).replace('VER', ver)
    files.append(path)

for path in files:
  print path
  with open(path) as f:
    lines = f.readlines()

    data = []
    print_name = False
    store_data = False

    for line in lines:
      if print_name:
        print_name = False
        print line

      if '<Description>' in line:
        print_name = True

      if store_data:
        if '</Data>' in line:
          store_data = False
          print data
          data = []
        else:
          pos = line.find('e')
          evt = float(line[:pos])
          exp = int(line[pos+1:pos+4])
          nevt = round(evt * 10**exp, 2)
          if 'bin 1' in line:
            data[0] = data[0] + nevt #add to underflow
          elif 'overflow' in line:
            data[len(data)-1] = data[len(data)-1] + nevt #add to overflow
          else: data.append(nevt)

      if '<Data>' in line:
        store_data = True
