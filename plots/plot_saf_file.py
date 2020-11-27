import sys
import os
from operator import add
from matplotlib import pyplot as plt

outnames = ['m120_SRlow.pdf', 'm120_SRhigh.pdf', 'm280_SRlow.pdf', 'm280_SRhigh.pdf']

for outfile in outnames:

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
          ma5_data_LH = [174.95, 135.26, 73.98, 33.34, 16.19]
          ma5_data_RH = [108.91, 90.81, 51.12, 28.26, 12.71]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
      elif "SRhigh" in outfile:
          plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 120$ GeV"
          comp_data = [6.36, 0.82, 0, 0, 0]
          comp_err = [0.59, 0.23, 0, 0, 0]
          ma5_data_LH = [201.02, 41.69, 3.3, 0.41, 0.0]
          ma5_data_RH = [127.13, 31.37, 2.06, 0.41, 0.83]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
  elif "m280" in outfile:
      LH_weight = 0.00129629
      RH_weight = 0.0012234
      if "SRlow" in outfile:
          plot_title = r"SR-lowMass for $m_{\tilde{\tau}} = 280$ GeV"
          comp_data = [0.82, 0.94, 0.446, 0.66, 3.22]
          comp_err = [0.32, 0.24, 0.08, 0.10, 0.67]
          ma5_data_LH = [608.37, 649.65, 658.85, 594.4, 1561.24]
          ma5_data_RH = [113.67, 93.99, 48.26, 26.04, 12.7]
          ma5_data_LH = [x * LH_weight for x in ma5_data_LH]
          ma5_data_RH = [x * RH_weight for x in ma5_data_RH]
      elif "SRhigh" in outfile:
          plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 280$ GeV"
          comp_data = [2.23, 2.49, 5.32, 2.47, 1.83]
          comp_err = [0.36, 0.35, 0.80, 0.44, 0.44]
          ma5_data_LH = [1923.12, 2265.32, 2614.52, 1984.21, 1381.97]
          ma5_data_RH = [119.29, 29.72, 0.41, 2.48, 0.41]
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
  #print(data)
  #print(norm_data)
  #print(comp_data)
  plt.rc("axes", labelsize=14)
  plt.rc("legend", fontsize=14)
  fig, ax = plt.subplots()
  plt.hist(x, bins=thebins, weights=norm_data, label='MA5', histtype='step', fill=False)
  plt.errorbar(x, norm_comp_data, yerr=norm_comp_err, label='ATLAS', fmt='.k', markersize=10)
  plt.xlabel(r"$m_{T2}$ [GeV]")
  plt.ylabel('Normalized Events')
  plt.title(plot_title, fontsize=20)
  for tick in ax.xaxis.get_major_ticks():
      tick.label.set_fontsize(12)
  for tick in ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(12)
  plt.legend()
  plt.savefig(outfile, bbox_inches='tight')
  print "wrote results to", outfile
