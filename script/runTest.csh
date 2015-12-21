#!/bin/csh

rm log_test
root4star -b -l << EOF >& log_test
.O2
.x makeAnaTree.C(15094070,"root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_mid_2014/ReversedFullField/P15ic/2014/094/15094070/st_physics_15094070_raw_0000007.MuDst.root",false,0)
.q
EOF
#.x makeAnaTree.C
#.x makeAnaTree.C(15080058,"root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_2014/ReversedFullField/P15ic/2014/080/15080058/st_physics_15080058_raw_1000002.MuDst.root",false,0)
