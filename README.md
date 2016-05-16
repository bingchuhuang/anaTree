# pico anaTree
main script and macro for submitting jobs starting from MuDst: script/AuAu200.xml and macros/makeAnaTree.C
Some options: 
prodMod: 0 for mb, 1 for ht, 2 for mtd (it's different in picoDstMaker, there is no "2" definedor different definition in picoDstMaker). It requires a run list as input. For Au+Au, a recenter file is required to correct event plane.
To save the recenter correction histograms: use macros/makeRecenter.C
All setters can be found in StPicoAnaTreeMaker.h 

submit jobs: script/submitAll.sh

hadd histograms and trees for each run: script/haddHist.sh and script/haddTree.sh
