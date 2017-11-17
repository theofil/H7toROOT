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

lambo::lambo() {}
lambo::~lambo() {}


#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void lambo::findDecayProducts(tcPPtr myPart ,std::vector<tcPPtr> &products)
{
    int familySize = myPart->children().size(); 
    for (int ii = 0 ; ii< familySize ; ++ii) // loop over all daughters 
    {
	tcPPtr daughter = myPart->children()[ii];
	if(daughter->children().size() == 0) // if final state particle 
	{
	    products.push_back(daughter);
            //cout << daughter->PDGName() << endl;
	} 
	else
	{
	    findDecayProducts(daughter, products);
	}
     } 
}


void lambo::analyze(tEventPtr event, long ieve, int loop, int state) {
    AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).

  // --- fill histogram with event statistics
  hstats_->Fill(1, event->weight());

  // --- initialize tree variables
    loop_       = 0;
    eveW_       = 0;
    mW_         = 0;
    nTracks_    = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) pt_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) eta_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) phi_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) dr_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) dp_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) dz_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) ID_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) MID_[ii] = 0;
    for(int ii = 0; ii <nTracksMax; ++ii) GMID_[ii] = 0;
    invM_       = 0;
    etaW_       = 0;
    ptW_        = 0;
    nCh_        = 0;
    nNu_        = 0;


    loop_       = loop;
    eveW_       = event->weight();

    // --- find W parton info
    std::vector<tcPPtr> WbosonProducts;
    std::vector<tcPPtr> ZbosonProducts;

    StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
    StepVector::const_iterator send=event->primaryCollision()->steps().end(); 

    bool isLEP(false);
    for(;sit!=send;++sit) 
    {
        if(mW_>0) break;  // no need to search further
        if(isLEP) break;  

      	ParticleSet part=(**sit).all();
	ParticleSet::const_iterator iter=part.begin();
	ParticleSet::const_iterator end =part.end();
	for( ;iter!=end;++iter) 
	{

            if((**iter).PDGName() == "e-")isLEP=true;
            if((**iter).PDGName() == "e+")isLEP=true;
	    if(isLEP) break;

	    if ((**iter).id()==ParticleID::Wplus) 
	    { 
		mW_ = (**iter).mass()/GeV;
		findDecayProducts(*iter, WbosonProducts);
		break;
    	    }
	}
    }
   

    // --- loop over the stable particles

    set<tcPPtr> particles;
    event->selectFinalState(inserter(particles));

    TLorentzVector recoW(0,0,0,0);
    TLorentzVector recoZ(0,0,0,0);

    for(set<tcPPtr>::const_iterator it = particles.begin(); it != particles.end(); ++it) 
    {
        // --- get P4 of the particle
	float pxTmp = float(((**it).momentum()).x()/GeV);
	float pyTmp = float((**it).momentum ().y()/GeV);
	float pzTmp = float((**it).momentum ().z()/GeV);
	float energyTmp = float((**it).momentum ().t()/GeV);
	float chTmp = (**it).data().charge()/eplus;

        const Lorentz5Momentum particleP4 = (**it).momentum ();
        float ptTmp  = sqrt(pxTmp*pxTmp + pyTmp*pyTmp);
        float etaTmp = particleP4.eta();
        float phiTmp = particleP4.phi();

        TLorentzVector tmpV(pxTmp, pyTmp, pzTmp, energyTmp);

        // --- do something within acceptance
        bool isInAcceptanceAndCharged = (ptTmp > 1 && fabs(etaTmp) < 2.5 && chTmp !=0) ? 1:0;
        bool isFromW(false);

        if (std::find(WbosonProducts.begin(), WbosonProducts.end(), *it) != WbosonProducts.end()) isFromW = true;

        if(isFromW) 
        {
          recoW +=  TLorentzVector(pxTmp, pyTmp, pzTmp, energyTmp);
          if(chTmp != 0) nCh_++;
          if(chTmp == 0) nNu_++;
        }

        if(isLEP) 
        {
          recoZ +=  TLorentzVector(pxTmp, pyTmp, pzTmp, energyTmp);
          if(chTmp != 0) nCh_++;
          if(chTmp == 0) nNu_++;
        }

        if(isInAcceptanceAndCharged || isFromW || isLEP)
        {
            if(nTracks_> nTracksMax) break;
            LorentzPoint vtxLab = (*it)->labVertex();
            float dx = vtxLab.x()/micrometer ;
            float dy = vtxLab.y()/micrometer ;
            float dz = vtxLab.z()/micrometer ;
            float dr = sqrt(dx*dx + dy*dy);
            float dp = vtxLab.phi();
            int MID  = 0;
            int GMID = 0; 

            if((*it)->parents().size() > 0)
            {
		Particle *mother = (*it)->parents()[0];  
		MID = mother->id();

                if(MID!=0 && mother->parents().size() > 0) 
                {
                   Particle *gmother = mother->parents()[0];
                   GMID = gmother->id();
                }
            }
             
            if(dr < 10) dr = 0; // don't save floating point for anything that small        

	    pt_[nTracks_]       = ptTmp;
	    eta_[nTracks_]      = etaTmp;
	    phi_[nTracks_]      = phiTmp;
	    dr_[nTracks_]       = dr;
	    dp_[nTracks_]       = dp;
	    dz_[nTracks_]       = dz;
            ID_[nTracks_]       = (*it)->id(); 
            MID_[nTracks_]      = MID; 
            GMID_[nTracks_]     = GMID; 
            charge_[nTracks_]   = chTmp; 
            isFromW_[nTracks_]  = isFromW; 
            nTracks_++;
        }
    }
 
    float invMass  = 0;
    if(isLEP)  invMass = recoZ.M();
    if(!isLEP) invMass = recoW.M();
    invM_ = invMass ; // it's either W+ or Z0 events, by default recoW.M() and recoZ.M() are initialized 0 so the one is greater than 0
    etaW_ = !isLEP ?  recoW.Rapidity(): 0;
    ptW_  = !isLEP ?  recoW.Pt(): 0;

    // --- fill the tree
    if(nCh_<=21)
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
  hstats_ ->Write();
  myfile_ ->Close();
}

void lambo::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
  hstats_ = new TH1I("hstats","hstats",1,0,2);
  std::cout <<"creating a TTree to be saved in output.root" << std::endl; 
  myfile_ = new TFile("output.root","RECREATE");
  events_ = new TTree("events","events");
  events_ ->Branch("loop",        &loop_,    "loop/I");  
  events_ ->Branch("eveW",        &eveW_,    "eveW/F");  
  events_ ->Branch("mW",          &mW_,      "mW/F");  
  events_ ->Branch("nTracks",     &nTracks_, "nTracks/I");  
  events_ ->Branch("pt",          pt_,       "pt[nTracks]/F");
  events_ ->Branch("eta",         eta_,      "eta[nTracks]/F");
  events_ ->Branch("phi",         phi_,      "phi[nTracks]/F");
  events_ ->Branch("dr",          dr_,       "dr[nTracks]/F");
  events_ ->Branch("dp",          dp_,       "dp[nTracks]/F");
  events_ ->Branch("dz",          dz_,       "dz[nTracks]/F");
  events_ ->Branch("ID",          ID_,       "ID[nTracks]/I");  
  events_ ->Branch("MID",         MID_,      "MID[nTracks]/I");  
  events_ ->Branch("GMID",        GMID_,     "GMID[nTracks]/I");  
  events_ ->Branch("charge",      charge_,   "charge[nTracks]/I");  
  events_ ->Branch("isFromW",     isFromW_,  "isFromW[nTracks]/O");  
  events_ ->Branch("invM",        &invM_,    "invM/F");  
  events_ ->Branch("ptW",         &ptW_,     "ptW/F");  
  events_ ->Branch("etaW",        &etaW_,    "etaW/F");  
  events_ ->Branch("nCh",         &nCh_,     "nCh/I");  
  events_ ->Branch("nNu",         &nNu_,     "nNu/I");  

 

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

