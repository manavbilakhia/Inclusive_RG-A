#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <bitset>
#include "clas12reader.h"
#include "jsonFileMerger.h"
#include "clas12databases.h"
#include "QADB.h"

// this file create a root file with ttree with all the information for QA DB from hipo files
const int nToProcess = -1;
void processHipo(TString inputFile);
int PCALFidXY(float x, float y, int cutLevel);


void ana12GeVShortFCQA(){
  TString inputFile;
  int isHipo = -1;
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".dat") || opt.Contains(".txt"))){
      inputFile=opt(5,opt.Sizeof());
      isHipo = 1;
    }
    else if(opt.Contains(".root")){
      inputFile=opt(5,opt.Sizeof());
      isHipo = 0;
    } 
  }
  cout << "isHipo: "<< isHipo << endl;
  if(isHipo < 0)  {
    std::cout << " *** please provide a root or text input file name..." << std::endl;
    exit(0);
  }
  processHipo(inputFile);
}

bool wCut(double e_px, double e_py, double e_pz, double e_E, double Ebeam){
  double lowCut = 0.7;
  double highCut = 2.7;
    
    
  //double lowCut = 0.7;
  //double highCut = 1.4;
    
  TLorentzVector ele(e_px,e_py,e_pz,e_E);
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	TLorentzVector fCM = fGamma + target;
  double W = fCM.M();
	return (W > lowCut) && (W < highCut);
}


bool thetaCut(double e_px, double e_py, double e_pz, double e_E, double Ebeam){

  double lowCut = 10.;
  double highCut = 27.;
  float toRD = 57.2958;
    
  TLorentzVector ele(e_px,e_py,e_pz,e_E);
  double theta = ele.Theta()*toRD;
    
	return (theta > lowCut) && (theta < highCut);
}


bool q2Cut(double e_px, double e_py, double e_pz, double e_E, double Ebeam){
  //return true;
    
  double lowCut = 2.55;
  double highCut = 10.4;
    
  TLorentzVector ele(e_px,e_py,e_pz,e_E);
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	const double Q2 =  -fGamma.M2();
	return (Q2 > lowCut) && (Q2 < highCut);
}

struct calParam{
  int sec;
  int layer;
  int strip;
};

void processHipo(TString inputFile){
    
  double Ebeam = 10.6041;
    
  bool isData = 0;
  bool isEmptyTarget = 0;
  bool isKinemAndFidCut = 0;


    
  vector<int> sectorE;
  vector<int> ele_ndf;
  vector<double> ele_tr_chi2;
  vector<double> p4_ele_px;
  vector<double> p4_ele_py;
  vector<double> p4_ele_pz;
  vector<double> p4_ele_vx;
  vector<double> p4_ele_vy;
  vector<double> p4_ele_vz;
  vector<double> p4_ele_E;
  vector<double> ele_chi2;

  vector<double> htccNPE;
  vector<double> htccX;
  vector<double> htccY;

  vector<double> pcalHX;
  vector<double> pcalHY;
  vector<double> pcalHZ;
    
  vector<double> pcalLu;
  vector<double> pcalLv;
  vector<double> pcalLw;

  vector<double> ecinHX;
  vector<double> ecinHY;
  vector<double> ecinHZ;

  vector<double> ecoutHX;
  vector<double> ecoutHY;
  vector<double> ecoutHZ;
 // FTOF p1b, ECin, and ECou. 

  vector<double> ftofHX;
  vector<double> ftofHY;

  vector<double> dcXR1;
  vector<double> dcYR1;
  vector<double> dcZR1;

  vector<double> dcXR2;
  vector<double> dcYR2;
  vector<double> dcZR2;

  vector<double> dcXR3;
  vector<double> dcYR3;
  vector<double> dcZR3;


  vector<double> pcalE;
  vector<double> ecinE;
  vector<double> ecoutE;

  vector<float>  vMC_px;
  vector<float>  vMC_py;
  vector<float>  vMC_pz;
    
 // Data only scin sec 6 eff fix:
  vector<int> sectorSci;
  vector<int> layerSci;
  vector<int> compSci;
    

  TFile outFile(Form("../../data/outH2R_test/%s_QADBtest_first_electron.root",inputFile.Data()), "recreate");
  TTree out_tree("out_tree","out_tree");
    
  //electrons                                                                                                                                                                                                 
  out_tree.Branch("sectorE", &sectorE);
  out_tree.Branch("ele_ndf", &ele_ndf);
  out_tree.Branch("ele_tr_chi2", &ele_tr_chi2);
    
  out_tree.Branch("p4_ele_px", &p4_ele_px);
  out_tree.Branch("p4_ele_py", &p4_ele_py);
  out_tree.Branch("p4_ele_pz", &p4_ele_pz);
  out_tree.Branch("p4_ele_E", &p4_ele_E);


  out_tree.Branch("p4_ele_vx", &p4_ele_vx);
  out_tree.Branch("p4_ele_vy", &p4_ele_vy);
  out_tree.Branch("p4_ele_vz", &p4_ele_vz);

  out_tree.Branch("ele_chi2", &ele_chi2);

  out_tree.Branch("htccNPE", &htccNPE);
  out_tree.Branch("htccX", &htccX);
  out_tree.Branch("htccY", &htccY);

  out_tree.Branch("pcalHX", &pcalHX);
  out_tree.Branch("pcalHY", &pcalHY);
  out_tree.Branch("pcalHZ", &pcalHZ);
    
  out_tree.Branch("pcalLu", &pcalLu);
  out_tree.Branch("pcalLv", &pcalLv);
  out_tree.Branch("pcalLw", &pcalLw);    

  out_tree.Branch("ecinHX", &ecinHX);
  out_tree.Branch("ecinHY", &ecinHY);
  out_tree.Branch("ecinHZ", &ecinHZ);

  out_tree.Branch("ecoutHX", &ecoutHX);
  out_tree.Branch("ecoutHY", &ecoutHY);
  out_tree.Branch("ecoutHZ", &ecoutHZ);
    
    // FTOF p1b
    
  out_tree.Branch("ftofHX", &ftofHX);
  out_tree.Branch("ftofHY", &ftofHY);


  out_tree.Branch("pcalE", &pcalE);
  out_tree.Branch("ecinE", &ecinE);
  out_tree.Branch("ecoutE", &ecoutE);

  out_tree.Branch("dcXR1", &dcXR1);
  out_tree.Branch("dcYR1", &dcYR1);
  out_tree.Branch("dcZR1", &dcZR1);
    
  out_tree.Branch("dcXR2", &dcXR2);
  out_tree.Branch("dcYR2", &dcYR2);
  out_tree.Branch("dcZR2", &dcZR2);
    
  out_tree.Branch("dcXR3", &dcXR3);
  out_tree.Branch("dcYR3", &dcYR3);
  out_tree.Branch("dcZR3", &dcZR3);

  //simulation                                                                                                                                                                                                
  out_tree.Branch("gen_px", &vMC_px);
  out_tree.Branch("gen_py", &vMC_py);
  out_tree.Branch("gen_pz", &vMC_pz);
       
// Data only scin sec 6 eff fix:
  out_tree.Branch("sectorSci", &sectorSci);
  out_tree.Branch("layerSci", &layerSci);
  out_tree.Branch("compSci", &compSci);

  HipoChain chain;
  TString nextFile;
  ifstream chainIn(inputFile);
  while (chainIn >> nextFile){
    chain.Add(nextFile);
  }
  int counter = 0;
  long long counter_generated = 0;
  auto config_c12=chain.GetC12Reader(); //in case you want to configure, use this

  config_c12->applyQA("pass1");//GETPASSSTRINGHERE="latest", "pass1, "pass2",...
  config_c12->db()->qadb_addQARequirement("MarginalOutlier");
  config_c12->db()->qadb_addQARequirement("TotalOutlier");
  config_c12->db()->qadb_addQARequirement("TerminalOutlier");
  config_c12->db()->qadb_addQARequirement("MarginalOutlier");
  config_c12->db()->qadb_addQARequirement("SectorLoss");
  config_c12->db()->qadb_addQARequirement("LowLiveTime");

  
  /*
   
  if(config_c12->qadb()!=nullptr && isEmptyTarget != 1){
    config_c12->db()->qadb_requireGolden(true);
    config_c12->applyQA();
  }
  */

  if (isData) config_c12->addAtLeastPid(11,1);
    
    
  gBenchmark->Start("db");
  auto& c12=chain.C12ref();
  while(chain.Next()){
    counter++;
    if (counter%1000000 == 0) cout << "processed "<< counter/1000000 << "M events" << endl;
    if (counter == nToProcess) break;

    
    p4_ele_px.clear();
    p4_ele_py.clear();
    p4_ele_pz.clear();
    
    ele_ndf.clear();
    ele_tr_chi2.clear();
    p4_ele_vx.clear();
    p4_ele_vy.clear();
    p4_ele_vz.clear();
    p4_ele_E.clear();
    ele_chi2.clear();
    sectorE.clear();
    htccNPE.clear();
    htccX.clear();
    htccY.clear();
    pcalHX.clear();
    pcalHY.clear();
    pcalHZ.clear();
    pcalLu.clear();
    pcalLv.clear();
    pcalLw.clear();
    
    ecinHX.clear();
    ecinHY.clear();
    ecinHZ.clear();
    ecoutHX.clear();
    ecoutHY.clear();
    ecoutHZ.clear();
    
    ftofHX.clear();
    ftofHY.clear();
    pcalE.clear();
    ecinE.clear();
    ecoutE.clear();
  
    dcXR1.clear();
    dcYR1.clear();
    dcZR1.clear();
    
    dcXR2.clear();
    dcYR2.clear();
    dcZR2.clear();
    
    dcXR3.clear();
    dcYR3.clear();
    dcZR3.clear();
    vMC_px.clear();
    vMC_py.clear();
    vMC_pz.clear();
    
    // Data only scin sec 6 eff fix:
    
    sectorSci.clear();
    layerSci.clear();
    compSci.clear();
    auto eventbank = c12->event();
    auto runconfigbank = c12->runconfig();
    if (isData){
        int nEvent = runconfigbank->getEvent();
        int runNumber = runconfigbank->getRun();
    }
    // auto particles = c12->getDetParticles();
    auto electron=c12->getByID(11);
    
    if (electron.size()==0) continue; // discard events with no electrons
    for (int u = 0; u < electron.size(); u++){
        TVector3 rotVector(electron[u]->cal(PCAL)->getHx(), electron[u]->cal(PCAL)->getHy(), electron[u]->cal(PCAL)->getHz());
        rotVector.RotateZ(-60*(electron[u]->getSector()-1)/57.2958);
        rotVector.RotateY(-25/57.2958);
        int isFid = PCALFidXY(rotVector.X(), rotVector.Y(), 0);
      
      
 	      double e_px = electron[u]->par()->getPx();
 	      double e_py = electron[u]->par()->getPy();
      	double e_pz = electron[u]->par()->getPz();
        double e_E = sqrt(electron[u]->par()->getPx()*electron[u]->par()->getPx() + electron[u]->par()->getPy()*electron[u]->par()->getPy() + electron[u]->par()->getPz()*electron[u]->par()->getPz());
        if (electron[u]->getRegion() == FD
          //&& electron[u]->par()->getVz() > -10 && electron[u]->par()->getVz() < 5 && isFid == 1){
        
          //&& electron[u]->par()->getVz() > -10 && electron[u]->par()->getVz() < 5    
          &&(wCut(e_px, e_py, e_pz, e_E, Ebeam) && q2Cut(e_px, e_py, e_pz, e_E, Ebeam)  )){
      
        
          //cout<<electron[u]->cal(ECOUT)->getHx()<< " -1- " <<electron[u]->par()->getChi2Pid()<<endl;
          //cout<<electron[u]->sci(FTOF1B)->getHX()<< " -2- " <<electron[u]->sci(FTOF1B)->getComponent()<<endl;
	        htccNPE.push_back(electron[u]->che(HTCC)->getNphe());
	        htccX.push_back(electron[u]->traj(HTCC,1)->getX());
	        htccY.push_back(electron[u]->traj(HTCC,1)->getY());
        
	        dcXR1.push_back(electron[u]->traj(DC,DC1)->getX());
	        dcYR1.push_back(electron[u]->traj(DC,DC1)->getY()); 
	        dcZR1.push_back(electron[u]->traj(DC,DC1)->getZ()); 
        
	        dcXR2.push_back(electron[u]->traj(DC,DC3)->getX()); 
	        dcYR2.push_back(electron[u]->traj(DC,DC3)->getY()); 
	        dcZR2.push_back(electron[u]->traj(DC,DC3)->getZ()); 
        
	        dcXR3.push_back(electron[u]->traj(DC,DC6)->getX()); 
	        dcYR3.push_back(electron[u]->traj(DC,DC6)->getY()); 
	        dcZR3.push_back(electron[u]->traj(DC,DC6)->getZ()); 
        
	        sectorE.push_back(electron[u]->getSector());
          ele_ndf.push_back(electron[u]->trk(DC)->getNDF());
          ele_tr_chi2.push_back(electron[u]->trk(DC)->getChi2());
        
	        p4_ele_vx.push_back(electron[u]->par()->getVx());
	        p4_ele_vy.push_back(electron[u]->par()->getVy());
	        p4_ele_vz.push_back(electron[u]->par()->getVz());
	        p4_ele_px.push_back(e_px);
	        p4_ele_py.push_back(e_py);
	        p4_ele_pz.push_back(e_pz);
	        p4_ele_E.push_back(e_E);
        
        
	        pcalHX.push_back(electron[u]->cal(PCAL)->getHx());
	        pcalHY.push_back(electron[u]->cal(PCAL)->getHy());
	        pcalHZ.push_back(electron[u]->cal(PCAL)->getHz());
	        pcalLu.push_back(electron[u]->cal(PCAL)->getLu());
	        pcalLv.push_back(electron[u]->cal(PCAL)->getLv());
	        pcalLw.push_back(electron[u]->cal(PCAL)->getLw());  
	        ecinHX.push_back(electron[u]->cal(ECIN)->getHx());
	        ecinHY.push_back(electron[u]->cal(ECIN)->getHy());
	        ecinHZ.push_back(electron[u]->cal(ECIN)->getHz());
	        ecoutHX.push_back(electron[u]->cal(ECOUT)->getHx());
	        ecoutHY.push_back(electron[u]->cal(ECOUT)->getHy());
	        ecoutHZ.push_back(electron[u]->cal(ECOUT)->getHz());
	        pcalE.push_back(electron[u]->cal(PCAL)->getEnergy());
	        ecinE.push_back(electron[u]->cal(ECIN)->getEnergy());
	        ecoutE.push_back(electron[u]->cal(ECOUT)->getEnergy());
	        ele_chi2.push_back(electron[u]->par()->getChi2Pid());
        
          ftofHX.push_back(electron[u]->sci(FTOF1B)->getHX());
          ftofHY.push_back(electron[u]->sci(FTOF1B)->getHY());
        
          sectorSci.push_back(electron[u]->sci(FTOF1B)->getSector());
          layerSci.push_back(electron[u]->sci(FTOF1B)->getLayer());
          compSci.push_back(electron[u]->sci(FTOF1B)->getComponent());
        }  
      }
    auto mcpbank=c12->mcparts();
    for(Int_t i=0;i<mcpbank->getRows();i++){
      mcpbank->setEntry(i);
      if ( mcpbank->getPid() == 11){
    	vMC_px.push_back(mcpbank->getPx());
    	vMC_py.push_back(mcpbank->getPy());
    	vMC_pz.push_back(mcpbank->getPz());
        cout<<mcpbank->getPx()<<endl;
      }
    }
    counter_generated++;
    out_tree.Fill();
  }
  
  cout <<"Accumulated charge post QA 2: " << chain.TotalBeamCharge()<<" nC"<<endl;
  cout << "Accumulated charge good 2: " << chain.db()->qa()->getAccCharge() << " nC" << endl;

    
  outFile.cd();
  outFile.Write();
  outFile.Close();
}
int PCALFidXY(float x, float y, int cutLevel){
  int fidCut = 0;
  float cutLimit = 0;
  if (cutLevel == 0) cutLimit = 260;
  if (cutLevel == 1) cutLimit = 250;
  if (y > -0.5*(x + cutLimit) && y < +0.5*(x + cutLimit)) fidCut = 1;
  return fidCut;
}