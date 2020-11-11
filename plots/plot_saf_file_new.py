#!/usr/bin/env python3

import sys
import os
from operator import add
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

outnames = ['m120_SRlow.pdf', 'm120_SRhigh.pdf', 'm280_SRlow.pdf', 'm280_SRhigh.pdf']

for outfile in outnames:

  pp = PdfPages(outfile)

  nbins = 5
  xmin = 70.
  xmax = 120.0

  plot_title = ""
  ma5_data = []
  comp_data = []
  comp_err = []
  #https://www.hepdata.net/record/ins1765529
  if "m120" in outfile:
      #Weights from yield excel file
      LH_weight = 0.019587
      RH_weight = 0.0210939
      if "SRlow" in outfile:
          plot_title = r"SR-lowMass for $m_{\tilde{\tau}} = 120$ GeV"
          comp_data = [4.73, 2.62, 1.88, 0.44, 0.114]
          comp_err = [0.50, 0.32, 0.27, 0.11, 0.072]
          ma5_data_LH = [175.56, 135.73, 74.24, 33.46, 16.25]
          ma5_data_RH = [109.29, 91.13, 51.3, 28.36, 12.739999999999998]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
      elif "SRhigh" in outfile:
          plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 120$ GeV"
          comp_data = [6.36, 0.82, 0, 0, 0]
          comp_err = [0.59, 0.23, 0, 0, 0]
          ma5_data_LH = [215.11, 44.61, 3.53, 0.44, 0.0]
          ma5_data_RH = [136.05, 33.57, 2.21, 0.44, 0.88]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
  elif "m280" in outfile:
      LH_weight = 0.00129629
      RH_weight = 0.0012234
      if "SRlow" in outfile:
          plot_title = r"SR-lowMass for $m_{\tilde{\tau}} = 280$ GeV"
          comp_data = [0.82, 0.94, 0.446, 0.66, 3.22]
          comp_err = [0.32, 0.24, 0.08, 0.10, 0.67]
          ma5_data_LH = [610.47, 651.9, 661.14, 596.46, 1566.65]
          ma5_data_RH = [114.07, 94.31, 48.43, 26.13, 12.75]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
      elif "SRhigh" in outfile:
          plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 280$ GeV"
          comp_data = [2.23, 2.49, 5.32, 2.47, 1.83]
          comp_err = [0.36, 0.35, 0.80, 0.44, 0.44]
          ma5_data_LH = [2057.93, 2424.11, 2797.8, 2123.31, 1478.85]
          ma5_data_RH = [127.65, 31.8, 0.44, 2.65, 0.44]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]

  data = list( map(add, ma5_data_LH, ma5_data_RH) )

  steps = (xmax - xmin)/nbins
  x = [xmin + steps/2.0 + steps*i for i in range(nbins)]
  thebins = [xmin + steps*i for i in range(nbins+1)]

  print "Plotting ",plot_title
  #Normalization
  norm_data = [ i/sum(data) for i in data ]
  norm_comp_data = [ i/sum(comp_data) for i in comp_data ]
  norm_comp_err = [ i/sum(comp_data) for i in comp_err ]
  print(data)
  print(norm_data)

  print(comp_data)
  plt.figure()
  plt.rc("axes",labelsize=14)
  plt.rc("legend",fontsize=12)
  plt.hist(x, bins=thebins, weights=norm_data, label='MA5', histtype='step', fill=False)
  plt.errorbar(x, norm_comp_data, yerr=norm_comp_err, label='ATLAS', fmt='.k')
  plt.xlabel(r"$m_{T2}$")
  plt.ylabel('Normalized Events')
  plt.legend()
  plt.title(plot_title)
  pp.savefig()
  pp.close()
  print "wrote results to", outfile
