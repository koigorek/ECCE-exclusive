//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in excl_ana.h.
//
// excl_ana(const std::string &name = "excl_ana")
// everything is keyed to excl_ana, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// excl_ana::~excl_ana()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int excl_ana::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int excl_ana::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int excl_ana::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int excl_ana::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int excl_ana::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int excl_ana::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int excl_ana::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void excl_ana::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "excl_ana.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <g4eval/CaloEvalStack.h>
//#include <calobase/RawCluster.h>
#include <g4eval/CaloRawClusterEval.h>
//#include <calobase/RawClusterContainer.h>
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4eval/SvtxEvalStack.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>


/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>



// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

using namespace std;


//____________________________________________________________________________..
//excl_ana::excl_ana(const std::string &name):
// SubsysReco(name)
//{
//  std::cout << "excl_ana::excl_ana(const std::string &name) Calling ctor" << std::endl;
//}


excl_ana::excl_ana(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "Diff_Tagg_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

_caloevalstackFEMC=nullptr;
_caloevalstackEEMC=nullptr;
_caloevalstackCEMC=nullptr;
_svtxEvalStack = nullptr;
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

  initializeVariables();
  initializeTrees();


}




//____________________________________________________________________________..
excl_ana::~excl_ana()
{
  delete tree;
if(_caloevalstackCEMC) delete _caloevalstackCEMC;
if(_caloevalstackFEMC) delete _caloevalstackFEMC;
if(_caloevalstackEEMC) delete _caloevalstackEEMC;
if(_svtxEvalStack) delete _svtxEvalStack;
  gsl_rng_free(m_RandomGenerator);

  std::cout << "excl_ana::~excl_ana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int excl_ana::Init(PHCompositeNode *topNode)
{

  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TNtuple("hitntup", "G4Hits", "x0:y0:z0:x1:y1:z1:edep");

  std::cout << "excl_ana::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;


  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -50, 50, 200, -50, 50);

  h2_ZDC_XY_double = new TH2F("ZDC_XY_double", "ZDC XY Double gamma", 200, -50, 50, 200, -50, 50);

  h1_E_dep = new TH1F("E_dep", "E Dependence", 120, 0.0, 60.0);

  h1_E_dep_smeared = new TH1F("E_dep_smeared", "E Dependence Smeared", 120, 0.0, 60.0);



  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::InitRun(PHCompositeNode *topNode)
{
  std::cout << "excl_ana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::process_event(PHCompositeNode *topNode)
{
//  std::cout << "excl_ana::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  nHits=0;
  hitsEEMC=0;
  hitsFEMC=0;
  hitsCEMC=0;
    RP1 = 0;
    RP2 =0;
    RPhits=0;
    B0hits=0;
  ntr=0;
  cout<<" event = "<<event_itt<<endl;
  ZDC_hit = 0;

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;




//=========================


PHG4TruthInfoContainer* m_TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
PHG4Particle* primary = m_TruthInfoContainer->GetParticle(3);

//PHG4TruthInfoContainer::ConstRange particles = m_TruthInfoContainer->GetPrimaryParticleRange();
Epx = primary->get_px();
Epy = primary->get_py();
Epz = primary->get_pz();

primary = m_TruthInfoContainer->GetParticle(1);
Ppx = primary->get_px();
Ppy = primary->get_py();
Ppz = primary->get_pz();

primary = m_TruthInfoContainer->GetParticle(2);
Gpx = primary->get_px();
Gpy = primary->get_py();
Gpz = primary->get_pz();

//cout<<Epx<<"\t"<<Ppx<<"\t"<<Gpx<<endl;

//========================

if(!_caloevalstackFEMC){
 _caloevalstackFEMC = new CaloEvalStack(topNode, "FEMC");
     _caloevalstackFEMC->set_strict(true);
}
else{
  _caloevalstackFEMC->next_event(topNode);
}

if(!_caloevalstackEEMC){
 _caloevalstackEEMC = new CaloEvalStack(topNode, "EEMC");
     _caloevalstackEEMC->set_strict(true);
}
else{
  _caloevalstackEEMC->next_event(topNode);
}


if(!_caloevalstackCEMC){
 _caloevalstackCEMC = new CaloEvalStack(topNode, "CEMC");
     _caloevalstackCEMC->set_strict(true);
}
else{
 _caloevalstackCEMC->next_event(topNode);
}

  std::string caloName;
  caloName="CLUSTER_FEMC";
  process_ClusterCalo(topNode,caloName);
  caloName="CLUSTER_EEMC";
  process_ClusterCalo(topNode,caloName);
  caloName="CLUSTER_CEMC";
  process_ClusterCalo(topNode,caloName);
    // process_g4hits(topNode);


  process_tracks(topNode);

  process_RomanPots(topNode,1);
  process_B0(topNode);
 // process_RomanPots(topNode,2);
//cout<<"B0 hts  = "<<B0hits<<endl;



  tree->Fill();


cout<<" Event filled, next "<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "excl_ana::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::EndRun(const int runnumber)
{
  std::cout << "excl_ana::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::End(PHCompositeNode *topNode)
{
  std::cout << "excl_ana::End(PHCompositeNode *topNode) This is the End..." << std::endl;

if(_caloevalstackCEMC) delete _caloevalstackCEMC;
if(_caloevalstackFEMC) delete _caloevalstackFEMC;
if(_caloevalstackEEMC) delete _caloevalstackEEMC;
//  h2_ZDC_XY->Write();
//  h2_ZDC_XY_double->Write();
//  
//  h1_E_dep->Write();
//  h1_E_dep_smeared->Write();

  
  outfile->cd();
  tree->Write();
  g4hitntuple->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::Reset(PHCompositeNode *topNode)
{
 std::cout << "excl_ana::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void excl_ana::Print(const std::string &what) const
{
  std::cout << "excl_ana::Print(const std::string &what) const Printing info for " << what << std::endl;
}


//***************************************************
//


int excl_ana::process_g4hits(PHCompositeNode* topNode)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
//  nodename << "G4HIT_" << detector;
//  nodename << "G4HIT_" << "ZDC";
  nodename << "G4HIT_" << "ZDC";
//  nodename << "G4HIT_" << "EEMC";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

//  cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " << endl;


  float smeared_E;


  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

    {
//	HIT_IN_ZDC=true;
	ZDC_hit++;
    }

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {


//	ZDC_hit++;

//      cout << "AAA" << endl;
      // the pointer to the G4Hit is hit_iter->second
      g4hitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_x(1),
                        hit_iter->second->get_y(1),
                        hit_iter->second->get_z(1),
                        hit_iter->second->get_edep());



//      cout << hit_iter->second->get_x(0)-90 << "   " << hit_iter->second->get_y(0) << endl;


      h2_ZDC_XY->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0)); 
//
      smeared_E = EMCAL_Smear(hit_iter->second->get_edep());
//
      if (ZDC_hit == 2 ) {

//      cout << hit_iter->second->get_x(0)-90 << "   " << hit_iter->second->get_y(0) << endl;

        h2_ZDC_XY_double->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0)); 
//      h1_E_dep->Fill(hit_iter->second->get_edep()); 
//
        h1_E_dep->Fill(hit_iter->second->get_edep()); 
        h1_E_dep_smeared->Fill(smeared_E); 
//
      }
//
//
//
////	hit_iter->get_avg_t();

    }
  }

//  cout << "BB" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//***************************************************
// Geting hits from Calo FEMC
int excl_ana::process_ClusterCalo(PHCompositeNode* topNode, string caloName)
{

//PHG4TruthInfoContainer* m_TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
//PHG4Particle* primaryPart = m_TruthInfoContainer->GetParticle(2);


	if(caloName == "CLUSTER_FEMC"){
		CaloRawClusterEval* clusterevalFEMC = _caloevalstackFEMC->get_rawcluster_eval();
		RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
		if(clusters){
			RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
			for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
				if(!cluster_iter->second) continue;
				float clX = cluster_iter->second->get_x();
				float clY = cluster_iter->second->get_y();
				float clZ = cluster_iter->second->get_z();
				float clEn = cluster_iter->second->get_energy();

				hitX[nHits] = clX;
				hitY[nHits] = clY;
				hitZ[nHits] = clZ;
				hitE[nHits] = clEn;
				caloInd[nHits] = 1;
				hitsNtowers[nHits] = cluster_iter->second->getNTowers();
				if(clusterevalFEMC){
					PHG4Particle* primary = clusterevalFEMC->max_truth_primary_particle_by_energy(cluster_iter->second);
					if(primary){
						hitPid[nHits] = primary->get_pid();
					}
				}
				nHits++;
				hitsFEMC++;
			} // for loop
		} // clusters

	} // FEMC

	if(caloName == "CLUSTER_EEMC"){
		CaloRawClusterEval* clusterevalEEMC = _caloevalstackEEMC->get_rawcluster_eval();
		RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
		if(clusters){
			RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
			for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
				if(!cluster_iter->second) continue;				
				float clX = cluster_iter->second->get_x();
				float clY = cluster_iter->second->get_y();
				float clZ = cluster_iter->second->get_z();
				float clEn = cluster_iter->second->get_energy();

				hitX[nHits] = clX;
				hitY[nHits] = clY;
				hitZ[nHits] = clZ;
				hitE[nHits] = clEn;
				caloInd[nHits] = 2;
				hitsNtowers[nHits] = cluster_iter->second->getNTowers();
				
                if(clusterevalEEMC){
					PHG4Particle* primary = clusterevalEEMC->max_truth_primary_particle_by_energy(cluster_iter->second);
					if(primary){	
						hitPid[nHits] = primary->get_pid();			
					}
				}
				nHits++;
				hitsEEMC++;
			} // for loop
		} // clusters

	} // EEMC


	
	if(caloName == "CLUSTER_CEMC"){
		CaloRawClusterEval* clusterevalCEMC = _caloevalstackCEMC->get_rawcluster_eval();
		RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, caloName.c_str());
		if(clusters){
			RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
			for(RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++){
				if(!(cluster_iter->second)) continue;
				float clX = cluster_iter->second->get_x();
				float clY = cluster_iter->second->get_y();
				float clZ = cluster_iter->second->get_z();
				float clEn = cluster_iter->second->get_energy();
				hitX[nHits] = clX;
				hitY[nHits] = clY;
				hitZ[nHits] = clZ;
				hitE[nHits] = clEn;
				caloInd[nHits] = 3;
				hitsNtowers[nHits] = cluster_iter->second->getNTowers();
				
				if(clusterevalCEMC){
					  PHG4Particle* primary = clusterevalCEMC->max_truth_primary_particle_by_energy(cluster_iter->second);

//					int val;

					if(primary)
					{
						hitPid[nHits] = primary->get_pid();
					}
				}
				nHits++;
				hitsCEMC++;
			} // for loop
			
		} // clusters

	} // EEMC

return Fun4AllReturnCodes::EVENT_OK;
}


//***************************************************
// Getting the RomanPots hits

int excl_ana::process_RomanPots(PHCompositeNode* topNode, int rPot)
{


  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
//  nodename << "G4HIT_" << detector;
//  nodename << "G4HIT_" << "ZDC";
//if(rPot==1) {  nodename << "G4HIT_" << "RomanPots_0";}
//else 
//nodename << "G4HIT_" << "rpTruth";
  nodename << "G4HIT_" << "rpTruth";

//cout<<" nodemae = "<<nodename.str()<<endl;
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();

      int counter = 0;
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
/*cout<<" Roman Hit "<<endl;
cout<<hit_iter->second->get_z(0)<<"\t"<<hit_iter->second->get_z(1)<<"\t"<<hit_iter->second->get_z(2)<<endl;
cout<<hit_iter->second->get_x(0)<<"\t"<<hit_iter->second->get_x(1)<<"\t"<<hit_iter->second->get_x(2)<<endl;
cout<<hit_iter->second->get_y(0)<<"\t"<<hit_iter->second->get_y(1)<<"\t"<<hit_iter->second->get_y(2)<<endl;
*/
        RPx[RPhits] = hit_iter->second->get_x(0);
        RPy[RPhits] = hit_iter->second->get_y(0);
        RPz[RPhits] = hit_iter->second->get_z(0);
	if(TMath::Abs(RPz[RPhits] - 2600.0)<50) rPot = 1;
	if(TMath::Abs(RPz[RPhits] - 2800.0)<50) rPot = 2;
        RPind[RPhits ] = rPot;
        RPhits++;
        counter++;

      }
RP1 = counter;
RP2 = counter;

    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//*****************************************************

float excl_ana::EMCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}


//*****************************************************

float excl_ana::HCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.50*.50/E + 0.1*0.1);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float excl_ana::PbWO4_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

//*****************************************************

float excl_ana::Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}
//***********************************************************
void excl_ana::initializeTrees()
{
  tree = new TTree("T", "A tree for DVCS");
  tree->Branch("Epx", &Epx, "Epx/F");
  tree->Branch("Epy", &Epy, "Epy/F");
  tree->Branch("Epz", &Epz, "Epz/F");

  tree->Branch("Ppx", &Ppx, "Ppx/F");
  tree->Branch("Ppy", &Ppy, "Ppy/F");
  tree->Branch("Ppz", &Ppz, "Ppz/F");

  tree->Branch("Gpx", &Gpx, "Gpx/F");
  tree->Branch("Gpy", &Gpy, "Gpy/F");
  tree->Branch("Gpz", &Gpz, "Gpz/F");

  tree->Branch("nHits",&nHits,"nHits/I");
  tree->Branch("caloInd",&caloInd,"caloInd[nHits]/I");
  tree->Branch("hitX",&hitX,"hitX[nHits]/F");
  tree->Branch("hitY",&hitY,"hitY[nHits]/F");
  tree->Branch("hitZ",&hitZ,"hitZ[nHits]/F");
  tree->Branch("hitE",&hitE,"hitE[nHits]/F");
  tree->Branch("hitPid",&hitPid,"hitPid[nHits]/I");
  tree->Branch("hitsNtowers",&hitsNtowers,"hitsNtowers[nHits]/I");

    tree->Branch("RP1",&RP1,"RP1/I");
    tree->Branch("RP2",&RP2,"RP2/I");
    tree->Branch("RPhits",&RPhits,"RPhits/I");
    tree->Branch("RPx",&RPx,"RPx[RPhits]/F");
    tree->Branch("RPy",&RPy,"RPy[RPhits]/F");
    tree->Branch("RPz",&RPz,"RPz[RPhits]/F");
    tree->Branch("RPind",&RPind,"RPind[RPhits]/I");

    tree->Branch("B0hits",&B0hits,"B0hits/I");
    tree->Branch("B0x",&B0x,"B0Px[B0hits]/F");
    tree->Branch("B0y",&B0y,"B0y[B0hits]/F");
    tree->Branch("B0z",&B0z,"B0z[B0hits]/F");
    tree->Branch("B0ind",&B0ind,"B0ind[B0hits]/I");

    
tree->Branch("ntr",&ntr,"ntr/I");
tree->Branch("tr_px",&tr_px,"tr_px[ntr]/F");
tree->Branch("tr_py",&tr_py,"tr_py[ntr]/F");
tree->Branch("tr_pz",&tr_pz,"tr_pz[ntr]/F");
tree->Branch("tr_p",&tr_p,"tr_p[ntr]/F");

tree->Branch("tr_x",&tr_x,"tr_x[ntr]/F");
tree->Branch("tr_y",&tr_y,"tr_y[ntr]/F");
tree->Branch("tr_z",&tr_z,"tr_z[ntr]/F");

tree->Branch("tr_phi",&tr_phi,"tr_phi[ntr]/F");
tree->Branch("tr_eta",&tr_eta,"tr_eta[ntr]/F");
tree->Branch("tr_Pid",&tr_Pid,"tr_Pid[ntr]/I");
tree->Branch("charge",&charge,"charge[ntr]/F");


tree->Branch("hitsEEMC",&hitsEEMC,"hitsEEMC/I");
tree->Branch("hitsFEMC",&hitsFEMC,"hitsFEMC/I");
tree->Branch("hitsCEMC",&hitsCEMC,"hitsCEMC/I");

}
//**************************************88
void excl_ana::initializeVariables()
{
	Epx= -1000;
	Epy= -1000;
	Epz= -1000;

	Ppx= -1000;
	Ppy= -1000;
	Ppz= -1000;
	Gpx= -1000;
	Gpy= -1000;
	Gpz= -1000;

	nHits=0;
	ntr=0;
    RPhits=0;
    B0hits=0;
    RP1=0;
    RP2=0;

	hitsEEMC=0;
	hitsFEMC=0;
	hitsCEMC=0;

	for(int i=0;i<10000;i++){
		hitsNtowers[i]=0;
		caloInd[i]=-1;
		hitX[i]=-1000;
		hitY[i]=-1000;
		hitZ[i]=-1000;
		hitE[i]=-1000;
		hitPid[i]=0;

		tr_px[i]=-1000;
		tr_py[i]=-1000;
		tr_pz[i]=-1000;
		tr_p[i]=-1000;
		tr_phi[i]=-1000;
		tr_eta[i]=-1000;
		charge[i]=-1000;
		tr_x[i]=-1000;
		tr_y[i]=-1000;
		tr_z[i]=-1000;
		tr_Pid[i]=0;
        
        RPx[i]=-1000;
        RPy[i]=-1000;
        RPz[i]=-1000;
        RPind[i]=-1;

        B0x[i]=-1000;
        B0y[i]=-1000;
        B0z[i]=-1000;
        B0ind[i]=-1;

	}
}

//***********************************
int excl_ana::process_tracks(PHCompositeNode *topNode)
{
	SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");

	if (!trackmap)
	    return Fun4AllReturnCodes::EVENT_OK;

 	/// EvalStack for truth track matching
  	if(!_svtxEvalStack)
    	{
      		_svtxEvalStack = new SvtxEvalStack(topNode);
      		_svtxEvalStack->set_verbosity(Verbosity());
    	}
    else{
        _svtxEvalStack->next_event(topNode);
    }

	SvtxTrackEval *trackeval = _svtxEvalStack->get_track_eval();
//	PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

	for (SvtxTrackMap::Iter iter = trackmap->begin();iter != trackmap->end();++iter)
	  {
	    SvtxTrack *track = iter->second;

	    /// Get the reconstructed track info
	    tr_px[ntr] = track->get_px();
	    tr_py[ntr] = track->get_py();
	    tr_pz[ntr] = track->get_pz();
	    tr_p[ntr] = TMath::Sqrt(tr_px[ntr] * tr_px[ntr] + tr_py[ntr] * tr_py[ntr] + tr_pz[ntr] * tr_pz[ntr]);

	    tr_phi[ntr] = track->get_phi();
	    tr_eta[ntr] = track->get_eta();

	    charge[ntr] = track->get_charge();
	    tr_x[ntr] = track->get_x();
	    tr_y[ntr] = track->get_y();
	    tr_z[ntr] = track->get_z();

          /// Get truth track info that matches this reconstructed track
          PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);
          if(truthtrack) tr_Pid[ntr] = truthtrack->get_pid();
          
	    ntr++;
	}




	return Fun4AllReturnCodes::EVENT_OK;
}


//***************************************************
// Getting the B0 hits

int excl_ana::process_B0(PHCompositeNode* topNode)
{
  ostringstream nodename;


//  cout << "Entering Romanpot?" << endl;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << "b0Truth";

 // cout << "Detector: " << nodename.str().c_str() <<" hits = "<<B0hits<< endl;

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  Int_t layer=-1;
  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {

//cout<<" z - "<<hit_iter->second->get_z(0)<<endl;

	B0x[B0hits] = (Float_t)hit_iter->second->get_x(0);
        B0y[B0hits] = (Float_t)hit_iter->second->get_y(0);
        B0z[B0hits] = (Float_t)hit_iter->second->get_z(0);

//cout<<" B0z = "<<B0z[B0hits]<<endl;
	if(TMath::Abs(B0z[B0hits] - 541.0)<5) {layer = 1;}
	if(TMath::Abs(B0z[B0hits] - 565.0)<5) {layer = 2;}
	if(TMath::Abs(B0z[B0hits] - 589.0)<5) {layer = 3;}
	if(TMath::Abs(B0z[B0hits] - 613.0)<5) {layer = 4;}

        B0ind[B0hits ] = layer;
        B0hits++;

      }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


