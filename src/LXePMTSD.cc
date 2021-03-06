//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file optical/LXe/src/LXePMTSD.cc
/// \brief Implementation of the LXePMTSD class
//
//
#include "LXePMTSD.hh"
#include "LXePMTHit.hh"
#include "LXeDetectorConstruction.hh"
#include "LXeUserTrackInformation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePMTSD::LXePMTSD(G4String name)
  : G4VSensitiveDetector(name),SiPMHitCollection(0),SiPMPositionsX(0)
  ,SiPMPositionsY(0),SiPMPositionsZ(0)
{
  collectionName.insert("sipmHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePMTSD::~LXePMTSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePMTSD::Initialize(G4HCofThisEvent* hitsCE){
  SiPMHitCollection = new LXePMTHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //Store collection with event and keep ID
  static G4int hitCID = -1;
  if(hitCID<0){
    hitCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitCID, SiPMHitCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LXePMTSD::ProcessHits(G4Step* ,G4TouchableHistory* ){
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode

G4bool LXePMTSD::ProcessHits_constStep(const G4Step* aStep,
                                       G4TouchableHistory*, G4TrackVector* trackVec){

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition()
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
 
  //User replica number 1 since photocathode is a daughter volume
  //to the pmt which was replicated
  //!!double check this
  G4int SiPM_ID=  
    4*aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1)
    +aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(0);
  //printf("SiPM_ID = %d volum1 = %d volume2 = %d\n",SiPM_ID,aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(0),aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1));
  G4VPhysicalVolume* physVol=
    aStep->GetPostStepPoint()->GetTouchable()->GetVolume(0);

  //Find the correct hit collection
  G4int n=SiPMHitCollection->entries();
  LXePMTHit* hit=NULL;
  for(G4int i=0;i<n;i++){
    if((*SiPMHitCollection)[i]->GetID()==SiPM_ID){
      hit=(*SiPMHitCollection)[i];
      break;
    }
  }
 
  if(hit==NULL){//this pmt wasnt previously hit in this event
    hit = new LXePMTHit(); //so create new hit
    hit->SetPMTNumber(SiPM_ID);
    hit->SetPMTPhysVol(physVol);
    SiPMHitCollection->insert(hit);
    hit->SetSiPMPos((*SiPMPositionsX)[SiPM_ID],(*SiPMPositionsY)[SiPM_ID],
                   (*SiPMPositionsZ)[SiPM_ID]);
  }

  hit->IncPhotonCount(); //increment hit for the selected pmt
  hit->AddHitTime(aStep->GetTrack()->GetGlobalTime()/ns); // or Local????
  hit->AddHitPosition(aStep->GetTrack()->GetPosition());
    
  if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="Scintillation") hit->AddHitProcess(1);
  else if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="Cerenkov") hit->AddHitProcess(0);
  else hit->AddHitProcess(-1);

  hit->SetDrawit(true);
    
  if(false){
        
        //Add optical cross talk
        G4double gain = 1.0e4;
        G4double MeanNumberOfPhotons = 3.0E-5*gain;
        G4int NumPhotons;
        
        if (MeanNumberOfPhotons > 10.){
            G4double sigma = std::sqrt(MeanNumberOfPhotons);
            NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);
        }else{
            NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
        }
        
        for (G4int i = 0; i < NumPhotons; i++) {
            
            // Determine photon energy
            G4double sampledEnergy = 1.15*eV+0.25*eV*G4UniformRand();
            
            // Generate random photon direction
            G4double cost = 1.-2.*G4UniformRand();
            G4double sint = std::sqrt((1.-cost)*(1.+cost));
            G4double phi = twopi*G4UniformRand();
            G4double sinp = std::sin(phi);
            G4double cosp = std::cos(phi);
            
            // Create photon momentum direction vector
            G4ParticleMomentum photonMomentum(sint*cosp, sint*sinp, cost);
            
            // Determine polarization of new photon
            G4ThreeVector photonPolarization(cost*cosp,cost*sinp,-sint);
            G4ThreeVector perp = photonMomentum.cross(photonPolarization);
            
            phi = twopi*G4UniformRand();
            sinp = std::sin(phi);
            cosp = std::cos(phi);
            
            photonPolarization = cosp * photonPolarization + sinp * perp;
            photonPolarization = photonPolarization.unit();
            
            // Generate a new photon
            G4DynamicParticle* aScintillationPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
            aScintillationPhoton->SetPolarization(photonPolarization.x(),photonPolarization.y(),photonPolarization.z());
            aScintillationPhoton->SetKineticEnergy(sampledEnergy);
            
            G4Track* aSecondaryTrack = new G4Track(aScintillationPhoton,
                                                   aStep->GetTrack()->GetGlobalTime(),
                                                   aStep->GetTrack()->GetPosition ());
            aSecondaryTrack->SetTouchableHandle(aStep->GetPreStepPoint()->GetTouchableHandle());
            aSecondaryTrack->SetParentID(aStep->GetTrack()->GetTrackID());
            
            trackVec->push_back(aSecondaryTrack);
        }
    }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePMTSD::EndOfEvent(G4HCofThisEvent* ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePMTSD::clear() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePMTSD::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePMTSD::PrintAll() {}
