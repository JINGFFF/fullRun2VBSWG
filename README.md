cmsrel CMSSW_9_4_9_cand2

cd CMSSW_9_4_9_cand2/src

cmsenv

git cms-init

git cms-merge-topic cms-egamma:EgammaID_949

git cms-merge-topic cms-egamma:EgammaPostRecoTools_940

cd $CMSSW_BASE/src

git cms-merge-topic cms-met:METFixEE2017_949_v2

git cms-addpkg RecoMET/METFilters

git clone -b 17 https://github.com/JINGFFF/fullRun2VBSWG.git VAJets
