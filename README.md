cmsrel CMSSW_10_2_18

cd CMSSW_10_2_18/src

cmsenv

git cms-init

git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier

git-cms-addpkg RecoEgamma/PhotonIdentification

git cms-addpkg RecoMET/METFilters

git clone -b 17 https://github.com/JINGFFF/fullRun2VBSWG.git
