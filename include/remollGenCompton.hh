#ifndef __REMOLLGENCOMPTON_HH
#define __REMOLLGENCOMPTON_HH

#include "remollVEventGen.hh"

class remollGenCompton : public remollVEventGen {
public:
  remollGenCompton();
  ~remollGenCompton();

private:
  void SamplePhysics(remollVertex *, remollEvent *);
};

#endif
