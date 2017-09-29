// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the lambo class.
//

#include "lambo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace std;

namespace {

    inline int isFromW(tcPPtr p, int ie)  // DEBUG version of isFromW
    {
        if (p->parents().size()!=1 )	cout << "ieve = " << ie << " pid = " << p->id() << " PDGname = " << p->PDGName() << " and p->parents().size() = " <<  p->parents().size() << endl;
	while (p->parents()[0] && p->parents().size() == 1) 
	{
	    cout << "ieve = " << ie << " pid = " << p->id() << " PDGname = " << p->PDGName() << " and p->parents().size() = " <<  p->parents().size() << endl;
	    if(p->id() != 21) {p = p->parents()[0];} else {return 2;} // 2 = propagete info the parton id was not checked
	    if(abs(p->id()) == 24) {cout << "mother found!" << endl; return 1;}
	}
	return 0; 
    } 

  inline int isFromW(tcPPtr p)  
  {
	//while (p->parents()[0] && p->parents().size() == 1) 
	if(p->parents().size() != 1) cout << "parents size is not 1 " << endl;
	while (p->parents()[0] && p->parents().size() == 1) 
	{
	    if(p->id() != 21 ) {p = p->parents()[0];} else {return 2;} // 2 = propagete info the parton id was not checked
	    if (abs(p->id()) == 24) return 1;
	}
    return 0; 
  } 


}

lambo::lambo() {}
lambo::~lambo() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif


void lambo::analyze(tEventPtr event, long ieve, int loop, int state) {
    AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).

  // --- initialize tree variables
    loop_       = 0;
    eveW_       = 0;
    mW_         = 0;
    flag_       = 0;
    nIndex_     = nIndex;
    for(int ii = 0; ii <nIndex; ++ii) totCh_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) chMult_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) nuMult_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) px_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) py_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) pz_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) E_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) cE_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) nE_[ii] = 0;
    for(int ii = 0; ii <nIndex; ++ii) invM_[ii] = 0;


    loop_       = loop;
    eveW_       = event->weight();
    mW_         = mW_;

    // --- find W parton info
    StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
    StepVector::const_iterator send=event->primaryCollision()->steps().end();
    for(;sit!=send;++sit) 
    {
	while(mW_ == 0) // once you find W, stop looping to gain execution speed 
	{
	    ParticleSet part=(**sit).all();
	    ParticleSet::const_iterator iter=part.begin();
	    ParticleSet::const_iterator end =part.end();
	    for( ;iter!=end;++iter) 
	    {
		if ((**iter).id()==ParticleID::Wplus) 
		{ 
	//		std::cout <<  "rW+ found with children = " << (**iter).children().size() << " and m = " << (**iter).mass()/GeV << std::endl;
		    mW_ = (**iter).mass()/GeV;
		    break;
		}
	    }
	}
    }
  
    // --- loop over the stable particles

    set<tcPPtr> particles;
    event->selectFinalState(inserter(particles));

    float tot_px = 0;
    float tot_py = 0;
    float tot_pz = 0;
    float tot_E  = 0;
    float tot_ch = 0; 

    int provanceFlag(0);

    for(set<tcPPtr>::const_iterator it = particles.begin(); it != particles.end(); ++it) 
    {
               
	float px_tmp = float(((**it).momentum()).x()/GeV);
	float py_tmp = float((**it).momentum ().y()/GeV);
       	float pz_tmp = float((**it).momentum ().z()/GeV);
	float E_tmp  = float((**it).momentum ().t()/GeV);
	float ch_tmp = (**it).data().charge()/eplus;

        LorentzVector<double> WP4(px_tmp, py_tmp, pz_tmp, E_tmp);
        float eta_tmp = WP4.eta();
        float pt_tmp  = sqrt(px_tmp*px_tmp + py_tmp*py_tmp);

	totCh_[all] += ch_tmp ;
	if (ch_tmp !=0) {chMult_[all] += eveW_;} else {nuMult_[all] += eveW_;}

	px_[all] += px_tmp;
	py_[all] += py_tmp;
	pz_[all] += pz_tmp;
	E_[all]  += E_tmp;
	if (ch_tmp !=0) {cE_[all] += E_tmp;} else {nE_[all] += E_tmp;}

	provanceFlag = isFromW(*it);
        if( provanceFlag == 1) isFromW(*it, ieve);
	if( provanceFlag == 1) 
	{

	    totCh_[fromW] += ch_tmp ;
	    if (ch_tmp !=0) {chMult_[fromW] += eveW_;} else {nuMult_[fromW] += eveW_;}

	    px_[fromW] += px_tmp;
	    py_[fromW] += py_tmp;
	    pz_[fromW] += pz_tmp;
	    E_[fromW]  += E_tmp;
	    if (ch_tmp !=0) {cE_[fromW] += E_tmp;} else {nE_[fromW] += E_tmp;}
            flag_ = 1;
            
            pt_fromW->Fill(pt_tmp, eveW_);
            eta_fromW->Fill(fabs(eta_tmp), eveW_);
	}
	if( provanceFlag == 2 ) flag_ = 2;

        if( pt_tmp > 1 && fabs(eta_tmp) < 2.5 && ch_tmp !=0)
        {
	    totCh_[det] += ch_tmp ;
	    if (ch_tmp !=0) {chMult_[det] += eveW_;} else {nuMult_[det] += eveW_;}

	    px_[det] += px_tmp;
	    py_[det] += py_tmp;
	    pz_[det] += pz_tmp;
	    E_[det]  += E_tmp;
	    if (ch_tmp !=0) {cE_[det] += E_tmp;} else {nE_[det] += E_tmp;}
        }
    }
 
    if( provanceFlag == 0 || provanceFlag == 2) cout << "didn't find mother W for ieve = " << ieve << " provanceFlag = " << provanceFlag << endl;

    invM_[all] = sqrt(E_[all]*E_[all] - px_[all]*px_[all] - py_[all]*py_[all] - pz_[all]*pz_[all]);
    invM_[fromW] = sqrt(E_[fromW]*E_[fromW] - px_[fromW]*px_[fromW] - py_[fromW]*py_[fromW] - pz_[fromW]*pz_[fromW]);
    invM_[det] = sqrt(E_[det]*E_[det] - px_[det]*px_[det] - py_[det]*py_[det] - pz_[det]*pz_[det]);

//  LorentzVector<double> WP4(tot_px, tot_py, tot_pz, tot_E); 
//  float mass2 = tot_E*tot_E - tot_px*tot_px - tot_py*tot_py - tot_pz*tot_pz;
//  std::cout << "mass2 = " << mass2 << " WP4.m2() = " << WP4.m2()  << std::endl;

//  if(fabs(mW_ - WP4.m()) >0.1){
//  std::cout << "ieve = " << ieve << " tot_px = " << tot_px << " tot_py = " << tot_py << " tot_pz = " << tot_pz << " tot_E = " << tot_E << " tot_ch = ";
//  cout << tot_ch << " P4 stable mass = "<< WP4.m(); 
//  cout << "  mW_ = " << mW_<< std::endl;
//  }

    // --- fill the tree
    events_->Fill(); 
}

LorentzRotation lambo::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void lambo::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void lambo::analyze(tPPtr, double weight) {}

void lambo::dofinish() {
  AnalysisHandler::dofinish();
  // *** ATTENTION *** Normalize and post-process histograms here.
  events_ ->GetCurrentFile();
  events_ ->Write();
  pt_fromW->Write();
  eta_fromW->Write();
  myfile_ ->Close();
}

void lambo::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
  std::cout <<"creating a TTree to be saved in output.root" << std::endl; 
  myfile_ = new TFile("output.root","RECREATE");
  events_ = new TTree("events","events");
  events_ ->Branch("loop", &loop_, "loop/I");  loop_ = 0;
  events_ ->Branch("eveW", &eveW_, "eveW/F");  eveW_ = 0;
  events_ ->Branch("mW", &mW_, "mW/F");  mW_ = 0;
  events_ ->Branch("flag", &flag_, "flag/I");  flag_ = 0;
  events_ ->Branch("nIndex", &nIndex_, "nIndex/I");  nIndex_ = nIndex;
  events_ ->Branch("totCh", totCh_, "totCh[nIndex]/I");  for(int ii = 0; ii <nIndex; ++ii) totCh_[ii] = 0;
  events_ ->Branch("chMult", chMult_, "chMult[nIndex]/I");  for(int ii = 0; ii <nIndex; ++ii) chMult_[ii] = 0;
  events_ ->Branch("nuMult", nuMult_, "nuMult[nIndex]/I");  for(int ii = 0; ii <nIndex; ++ii) nuMult_[ii] = 0;
  events_ ->Branch("px", px_, "px[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) px_[ii] = 0;
  events_ ->Branch("py", py_, "py[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) py_[ii] = 0;
  events_ ->Branch("pz", pz_, "pz[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) pz_[ii] = 0;
  events_ ->Branch("E", E_, "E[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) E_[ii] = 0;
  events_ ->Branch("cE", cE_, "cE[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) cE_[ii] = 0;
  events_ ->Branch("nE", nE_, "nE[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) nE_[ii] = 0;
  events_ ->Branch("invM", invM_, "invM[nIndex]/F");  for(int ii = 0; ii <nIndex; ++ii) invM_[ii] = 0;

 
  pt_fromW = new TH1F("pt_fromW","pt_fromW",100,0,100);
  eta_fromW = new TH1F("eta_fromW","eta_fromW",100,0, 2.4);

}


IBPtr lambo::clone() const {
  return new_ptr(*this);
}

IBPtr lambo::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void lambo::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void lambo::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<lambo,AnalysisHandler>
  describeThePEGlambo("ThePEG::lambo", "lambo.so");

void lambo::Init() {

  static ClassDocumentation<lambo> documentation
    ("There is no documentation for the lambo class");

}

