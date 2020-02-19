# ma5stau
Stau production analysis at 2nd MadAnalysis 5 workshop on LHC recasting @ Korea

## Analysis being worked on

- [ATLAS-SUSY-2018-04](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-04/)

## Packages in use

- [MA5_v1.8beta.tgz](https://launchpad.net/madanalysis5)

## github How-to
```{.Bash}
git init
git remote add origin git@github.com:jhyeon/ma5stau.git
mv madanalysis5/tools/PAD/Build/SampleAnalyzer/User/Analyzer/analysisList.h ~/
git pull origin master
```

## Plots
```{.Bash}
python plot_saf_file.py /path_to_Output_SAF_Histograms/histos.saf
```
