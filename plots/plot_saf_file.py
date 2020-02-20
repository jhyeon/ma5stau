#!/usr/bin/env python3

import sys
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

def extract_runname(textfile):
    """
    Given filepath extract runname.
    """
    runname = textfile
    if 'SMS' in runname:
        rn = 'SMS' + runname.split('SMS')[1]
        rn = rn[:rn.index('/')]
        return rn
    if 'Output/' in textfile:
        runname = textfile[textfile.index('Output') + 7:]
    if 'Histograms/' in textfile: runname = runname.replace('Histograms/', '')
    runname  = runname.replace('/', '_')
    if '.saf' in runname: runname = runname.replace('.saf', '')
    return runname


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if 'histos.saf' in sys.argv[1]:
            textfile = sys.argv[1]
        else:
            textfile = os.path.join(sys.argv[1], "Histograms/histos.saf")
    else:
        print("Usage: python plot_saf_file.py /path/to/ma5/output/runname/analysisname/Histograms/histos.saf [outfile]")
        print("Or: python plot_saf_file.py /path/to/ma5/output/runname/analysisname [outfile]")
        textfile = '/home/jory/bin/madanalysis5/ra2b/Output/SMS-GlGl_onejet_T1bbbb_1500_100/cms_sus_16_033_0/Histograms/histos.saf'
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        runname = extract_runname(textfile)
        outfile = runname + '.pdf'
    pp = PdfPages(outfile)
    fopen = open(textfile, 'r')
    line  = fopen.readline()
    while not '<SAFfooter>' in line:
        line = fopen.readline()

        ## Go to next histogram
        while line != '' and not '<Description>' in line:
            line = fopen.readline()

        ## Now read the title
        if line == '': break
        line = fopen.readline()
        title = line.replace("\n","")
        title = title.replace('"','')

        # Now go on and read binning
        line = fopen.readline()
        line = fopen.readline()
        nbins = int(line.split()[0])
        xmin = float(line.split()[1])
        xmax = float(line.split()[2])
        while not line == '' and not 'nevents' in line:
            line = fopen.readline()
        nentries = line.split()[0]
        while not line == '' and not '<Data>' in line:
            line = fopen.readline()
        line = fopen.readline()

        data = []
        while not line == '' and not '</Data>' in line:
            if not 'flow' in line:
                data.append(float(line.split()[0]))
            elif 'overflow' in line:
                data[nbins-1] += float(line.split()[0])
            line = fopen.readline()
        
        plot_title = ""
        comp_data = []
        comp_err = []
        if "m120" in outfile:
            if "SRlow" in title:
                plot_title = r"SR-lowMass for $m_{\tilde{\tau}} = 120$ GeV"
                comp_data = [4.73, 2.62, 1.88, 0.44, 0.114]
                comp_err = [0.50, 0.32, 0.27, 0.11, 0.072]
            elif "SRhigh" in title:
                plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 120$ GeV"
                comp_data = [6.36, 0.82, 0, 0, 0]
                comp_err = [0.59, 0.23, 0, 0, 0]
        elif "m280" in outfile:
            if "SRlow" in title:
                plot_title = r"SR-lowMass for $m_{\tilde{\tau}} = 280$ GeV"
                comp_data = [0.82, 0.94, 0.446, 0.66, 3.22]
                comp_err = [0.32, 0.24, 0.08, 0.10, 0.67]
            elif "SRhigh" in title:
                plot_title = r"SR-highMass for $m_{\tilde{\tau}} = 280$ GeV"
                comp_data = [2.23, 2.49, 5.32, 2.47, 1.83]
                comp_err = [0.36, 0.35, 0.80, 0.44, 0.44]
        
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
    fopen.close()
    pp.close()
    print "wrote results to", outfile
