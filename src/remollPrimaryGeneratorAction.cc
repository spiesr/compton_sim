#include "remollPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "remollIO.hh"
#include "remollVEventGen.hh"
#include "remollEvent.hh"
#include "remollRun.hh"
#include "remollRunData.hh"
#include "remolltypes.hh"
#include "globals.hh"

#include "remollGenMoller.hh"
#include "remollGenpElastic.hh"
#include "remollGenpInelastic.hh"
#include "remollGenPion.hh"
#include "remollGenBeam.hh"
#include "remollGenFlat.hh"
#include "remollGenCompton.hh"

remollPrimaryGeneratorAction::remollPrimaryGeneratorAction() {
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);


  fDefaultEvent = new remollEvent();
  fDefaultEvent->ProduceNewParticle( 
	  G4ThreeVector(0.*cm,0.*cm,-100.*cm),
	  G4ThreeVector(0.0,0.0, gDefaultBeamE),
	  "e-" );

  double kinE = sqrt(fDefaultEvent->fPartMom[0].mag()*fDefaultEvent->fPartMom[0].mag() 
	  + fDefaultEvent->fPartType[0]->GetPDGMass()*fDefaultEvent->fPartType[0]->GetPDGMass() )
      -  fDefaultEvent->fPartType[0]->GetPDGMass();

  // Default generator data
  fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
  fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
  fParticleGun->SetParticleEnergy( kinE  );
  fParticleGun->SetParticlePosition( fDefaultEvent->fPartPos[0] );

  fEventGen = NULL;
}

remollPrimaryGeneratorAction::~remollPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fDefaultEvent;
}

void remollPrimaryGeneratorAction::SetGenerator(G4String genname) {

    fEventGen = NULL;

    if( genname == "moller" ){
	fEventGen = new remollGenMoller();
    }

    if( genname == "elastic" ){
	fEventGen = new remollGenpElastic();
    }

    if( genname == "inelastic" ){
	fEventGen = new remollGenpInelastic();
    }

    if( genname == "pion" ){
	fEventGen = new remollGenPion();
    }

    if( genname == "beam" ){
	fEventGen = new remollGenBeam();
    }

    if( genname == "flat" ){
	fEventGen = new remollGenFlat();
    }

    if( genname == "compton" ){
      fEventGen = new remollGenCompton();
    }

    if( !fEventGen ){
	G4cerr << __FILE__ << " line " << __LINE__ << " - ERROR generator " << genname << " invalid" << G4endl;
	exit(1);
    } else {
	G4cout << "Setting generator to " << genname << G4endl;
    }

    remollRun::GetRun()->GetData()->SetGenName(genname.data());

    return;
}

void remollPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  /*  Generate event, set IO data */

  remollEvent *thisev = NULL;
  if( fEventGen ){  // Specified our own generator
      thisev = fEventGen->GenerateEvent();
      for( unsigned int pidx = 0; pidx < thisev->fPartType.size(); pidx++ ){

	  double kinE = sqrt(thisev->fPartMom[pidx].mag()*thisev->fPartMom[pidx].mag() + 
		             thisev->fPartType[pidx]->GetPDGMass()*thisev->fPartType[pidx]->GetPDGMass())
	              -  thisev->fPartType[pidx]->GetPDGMass();

	  fParticleGun->SetParticleDefinition(thisev->fPartType[pidx]);
	  fParticleGun->SetParticleMomentumDirection(thisev->fPartMom[pidx].unit());
	  fParticleGun->SetParticleEnergy( kinE  );
	  fParticleGun->SetParticlePosition( thisev->fPartPos[pidx] );

	  fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      if( thisev->fPartType.size() > 0 ){
	  fIO->SetEventData(thisev);
      }
  } else { // Use default, static single generator
      // Update this just in case things changed
      // from the command user interface
      fDefaultEvent->Reset();
      fDefaultEvent->ProduceNewParticle( 
	      fParticleGun->GetParticlePosition(),
	      fParticleGun->GetParticleMomentumDirection()*
	      fParticleGun->GetParticleMomentum(),
	      fParticleGun->GetParticleDefinition()->GetParticleName() );
      fIO->SetEventData(fDefaultEvent);

      fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if( thisev != NULL ){
      delete thisev;
  }
}

G4ParticleGun* remollPrimaryGeneratorAction::GetParticleGun() {
  return fParticleGun;
} 

