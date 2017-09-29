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


void lambo::analyze(tEventPtr event, long ieve, int loop, int state) {
    AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).

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


    loop_       = loop;
    eveW_       = event->weight();

    // --- find W parton info
    Particle *Wboson;
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
                    Wboson = *iter;
		    break;
		}
	    }
	}
    }
   

    // --- loop over the stable particles

    set<tcPPtr> particles;
    event->selectFinalState(inserter(particles));

    for(set<tcPPtr>::const_iterator it = particles.begin(); it != particles.end(); ++it) 
    {
               
        // --- get P4 of the particle
	float pxTmp = float(((**it).momentum()).x()/GeV);
	float pyTmp = float((**it).momentum ().y()/GeV);
       	float pzTmp = float((**it).momentum ().z()/GeV);
       	float energyTmp = float((**it).momentum ().t()/GeV);
	float chTmp = (**it).data().charge()/eplus;

        LorentzVector<double> WP4(pxTmp, pyTmp, pzTmp, energyTmp);
        float ptTmp  = sqrt(pxTmp*pxTmp + pyTmp*pyTmp);
        float etaTmp = WP4.eta();
        float phiTmp = WP4.phi();


        // --- do something

        if( ptTmp > 1 && fabs(etaTmp) < 2.5 && chTmp !=0)
        {
            if(nTracks_> nTracksMax) break;
            LorentzPoint vtxLab = (*it)->labVertex();
            float dx = vtxLab.x()/micrometer ;
            float dy = vtxLab.y()/micrometer ;
            float dz = vtxLab.z()/micrometer ;
            float dr = sqrt(dx*dx + dy*dy);
            float dp = vtxLab.phi();
            
            if(dr < 100) // zero all coordinates for vertices that are not easy
            {
		dr = 0; dp = 0; dz = 0;
	    }
  
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
           

	    pt_[nTracks_]      = ptTmp;
	    eta_[nTracks_]     = etaTmp;
	    phi_[nTracks_]     = phiTmp;
	    dr_[nTracks_]      = dr;
	    dp_[nTracks_]      = dp;
	    dz_[nTracks_]      = dz;
            ID_[nTracks_]      = (*it)->id(); 
           // MID_[nTracks_]     = (*it)->parents().size() > 0 ? ((*it)->parents()[0])->id():0; 
            MID_[nTracks_]     = MID; 
            GMID_[nTracks_]    = GMID; 
            charge_[nTracks_]  = chTmp; 
            nTracks_++;

        }
    }
 

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
  myfile_ ->Close();
}

void lambo::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
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

