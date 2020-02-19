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
        title = line

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
        steps = (xmax - xmin)/nbins
        x = [xmin + steps/2.0 + steps*i for i in range(nbins)]
        thebins = [xmin + steps*i for i in range(nbins+1)]
        plt.figure()
        plt.hist(x, bins=thebins, weights=data, label=title + '\nEntries ' + nentries, histtype='step', fill=False)
        plt.xlabel(title)
        plt.ylabel('number of events')
        plt.legend()
        plt.title(title)
        pp.savefig()
    fopen.close()
    pp.close()
    print "wrote results to", outfile
