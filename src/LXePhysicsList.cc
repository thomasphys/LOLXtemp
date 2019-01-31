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
//
//
//
//



#include "LXePhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4PhysicsListHelper.hh"
//#include "G4NuclideTable.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
//#include "G4NuclearLevelData.hh"
//#include "G4DeexPrecoParameters.hh"
#include "G4ParticleTypes.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4OpticalPhysics.hh"

#include "G4SystemOfUnits.hh"


#include "G4Cerenkov.hh"
#include "G4OpticalProcessIndex.hh"
#include "G4ParticleDefinition.hh"

#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4IonConstructor.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePhysicsList::LXePhysicsList()
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  defaultCutValue=0.1*mm;//0.5*eV;
  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // EM physics
  RegisterPhysics(new G4EmLowEPPhysics());
  //RegisterPhysics(new G4EmStandardPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // Optical Physics
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

   // Not required parameters are commented out
  //opticalPhysics->SetWLSTimeProfile("exponential");
  opticalPhysics->SetFiniteRiseTime(false);
  opticalPhysics->SetScintillationYieldFactor(1.0);// Percentage of scintillation light we actually produce
  opticalPhysics->SetScintillationExcitationRatio(0.0);
  opticalPhysics->SetMaxNumPhotonsPerStep(10);
  opticalPhysics->SetMaxBetaChangePerStep(1.0);
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);

  RegisterPhysics( opticalPhysics );

  G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  G4int MaxNumPhotons = 30;
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePhysicsList::~LXePhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCutsWithDefault();
}
