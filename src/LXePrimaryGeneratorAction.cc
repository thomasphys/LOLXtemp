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
/// \file optical/LXe/src/LXePrimaryGeneratorAction.cc
/// \brief Implementation of the LXePrimaryGeneratorAction class
//
//
#include "LXePrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::LXePrimaryGeneratorAction(){
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e+");
//  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="alpha");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(0*MeV);
  eventmode = 0;
  loadbetadecay = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction(){
    if(loadbetadecay){ 
	    	delete fDBetaParticleGun1;
    		delete fDBetaParticleGun2;
    }else{
	    delete fParticleGun;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::GetEnergyFraction(double r){
	int bin = r*10000;
	return dbetaenergysplit[bin];
}

void LXePrimaryGeneratorAction::GenerateDoubleBetaDecay(G4Event* anEvent){

	if(loadbetadecay){
		delete fParticleGun;
		G4int n_particle = 2;
		fDBetaParticleGun1 = new G4ParticleGun(1);
		fDBetaParticleGun2 = new G4ParticleGun(1);
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
		G4ParticleDefinition* particle = particleTable->FindParticle("e-");
		fDBetaParticleGun1->SetParticleDefinition(particle);
		fDBetaParticleGun2->SetParticleDefinition(particle);

		//load spectrum from file
		loadbetadecay = true;
	}

	//create random vertex position
	G4ThreeVector vertex = G4ThreeVector(0.0,0.0,0.0);

	//Split energy between neutrinos and electrons
	G4double Q = 2.479;
	G4double M_e = 0.511;

	G4double NeutrinoEn = Q*(1.-GetEnergyFraction(rand->Uniform()));
        G4double ElectronEn = Q-NeutrinoEn;
        G4double NeutrinoEn1 = NeutrinoEn*rand->Uniform();
	G4double NeutrinoEn2 = NeutrinoEn-NeutrinoEn1;	

	G4double maximum_neutrino_net_momenum = 2.0*TMath::Sqrt((0.5*ElectronEn+M_e)*(0.5*ElectronEn+M_e)-M_e*M_e);
	G4double net_momentum = maximum_neutrino_net_momenum*rand->Uniform();

	double theta = TMath::ACos(-1.+2.*rand->Uniform());
	double phi = TMath::Pi()*2.*rand->Uniform();	
	G4ThreeVector netmomentumdir = G4ThreeVector(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta));

	theta = TMath::ACos(-1.+2.*rand->Uniform());
        phi = TMath::Pi()*2.*rand->Uniform();
        G4ThreeVector dir1 = G4ThreeVector(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta));
	G4ThreeVector dir2 = G4ThreeVector(-TMath::Sin(theta)*TMath::Cos(phi),-TMath::Sin(theta)*TMath::Sin(phi),-TMath::Cos(theta));

	fDBetaParticleGun1->SetParticlePosition(vertex);
  	fDBetaParticleGun1->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  	fDBetaParticleGun1->SetParticleEnergy(0*MeV);

	fDBetaParticleGun2->SetParticlePosition(vertex);
        fDBetaParticleGun2->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        fDBetaParticleGun2->SetParticleEnergy(0*MeV);
}

void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

	if(eventmode == 1){
	}else{
		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
}
