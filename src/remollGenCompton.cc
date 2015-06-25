#include "remollGenCompton.hh" //Need to define class structure in a .hh file

#include "CLHEP/Random/RandFlat.h"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "remolltypes.hh"

remollGenCompton::remollGenCompton(){
  fThCoM_min = 0.0*deg;
  fThCoM_max = 90.0*deg;

  fApplyMultScatt = false;
}

remollGenCompton::~remollGenCompton(){
}

void remollGenCompton::SamplePhysics(remollVertex *vert, remollEvent *evt){
  double beamE = vert->GetBeamE();
  double me = electron_mass_c2;

  double beta_com  = sqrt( (beamE - me)/(beamE + me) );
    double gamma_com = 1.0/sqrt(1.0 - beta_com*beta_com);

    double e_com = me*gamma_com;
    double thcom = acos(CLHEP::RandFlat::shoot(cos(fThCoM_max), cos(fThCoM_min)));
    double phcom = CLHEP::RandFlat::shoot(0.0, 2.0*pi);

    double alpha = 1.348;
    double k_photon = 2.3305*eV;

    //Replace with cross section for compton interaction
    double sigma = alpha*alpha*pow(3.0+cos(thcom)*cos(thcom),2.0)*hbarc*hbarc/pow(sin(thcom),4.0)/(2.0*me*beamE); // units of area

    double V = 2.0*pi*(cos(fThCoM_min) - cos(fThCoM_max));

    //This will need to be replaced with actual cross section of compton interactoin
    evt->SetEffCrossSection(sigma*V*vert->GetMaterial()->GetZ()/2.0);

    //electronPosition = G4ThreeVector(0.0, 0.0, 0.0);
    //opticPhotonPosition = G4ThreeVector(0.0, 0.0, 0.0); //Need to be in same position in the centre of mass frame
    //electronMomentum = G4ThreeVector(0.0,0.0,e_com); //Need to figure out the coordinate transformation for this and the optical photon momentum in com frame
    //opticPhotonMomentum = G4ThreeVector(k_photon*sin(alpha),0.0,k_photon*cos(alpha));

    evt->ProduceNewParticle(G4ThreeVector(0.0, 0.0, 0.0), G4ThreeVector(0.0,0.0,e_com), "e-");
    evt->ProduceNewParticle(G4ThreeVector(0.0, 0.0, 0.0), G4ThreeVector(k_photon*sin(alpha),0.0,k_photon*cos(alpha)), "opticalphoton"); 
    return;
}
