// -*- C++ -*-
#ifndef ThePEG_lambo_H
#define ThePEG_lambo_H
//
// This is the declaration of the lambo class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector> 
#include <algorithm>

#define nTracksMax 4000



namespace ThePEG {

/**
 * Here is the documentation of the lambo class.
 *
 * @see \ref lamboInterfaces "The interfaces"
 * defined for lambo.
 */
class lambo: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  lambo();

  /**
   * The destructor.
   */
  virtual ~lambo();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Return a LorentzTransform which would put the event in the
   * desired Lorentz frame.
   * @param event a pointer to the Event to be considered.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tcEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the weight of the current event.
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight the weight of the current event.
   */
  virtual void analyze(tPPtr particle, double weight);
  //@}

protected:

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  void findDecayProducts(tcPPtr myPart ,std::vector<tcPPtr> &products);

  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  lambo & operator=(const lambo &);


  private:
  TTree * events_;
  TFile * myfile_;
  int       loop_;
  float     eveW_;
  float       mW_;
  int      nTracks_;
  float    pt_[nTracksMax];
  float    eta_[nTracksMax];
  float    phi_[nTracksMax];
  float    dr_[nTracksMax]; // rho of the lab vertex
  float    dp_[nTracksMax]; // phi of the lab vertex
  float    dz_[nTracksMax]; // rho of the lab vertex
  int      ID_[nTracksMax];
  int      MID_[nTracksMax];
  int      GMID_[nTracksMax];
  int      charge_[nTracksMax];
  bool     isFromW_[nTracksMax];
  float    invM_;
  float    ptW_;
  float    etaW_;
  int      nCh_;
  int      nNu_;

};

}

#endif /* ThePEG_lambo_H */
