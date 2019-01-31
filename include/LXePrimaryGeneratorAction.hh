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
/// \file optical/LXe/include/LXePrimaryGeneratorAction.hh
/// \brief Definition of the LXePrimaryGeneratorAction class
//
//
#ifndef LXePrimaryGeneratorAction_h
#define LXePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include "LXeDetectorConstruction.hh"
#include "LXePrimaryGeneratorMessenger.hh"
#include "G4Navigator.hh"

class G4ParticleGun;
class G4Event;

class LXePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    LXePrimaryGeneratorAction();
    virtual ~LXePrimaryGeneratorAction();
 
  public:

    virtual void GeneratePrimaries(G4Event* anEvent);
    void Set_generator(G4String generator);
    void Set_nCaptureXe_Method(G4String method);
    void Set_Xe_Component(G4String component);
    void Set_bb2n_CutOffMin(G4double frac);
    void Set_bb2n_CutOffMax(G4double frac);
    void Set_GenInCenter(G4bool b);
    void SetRandDimensions();
    void GetIsotropicDirection(G4ThreeVector& pos);
    void GetUnifRandPosInLXe(G4ThreeVector& pos);
    void Generate_bb0n(G4Event* anEvent);
    G4double D_bb0n_spectral_max(G4double k);
    G4double D_bb0n_spectrum(G4double k, G4double d);
    G4double Fermi_function(G4int z, G4double ke);
    void Generate_bb2n(G4Event* anEvent);
    G4double BB2n_sum_spectrum(G4double k);
    G4double D_spectrum(G4double k, G4double d);
    G4double SimpsonsRule(G4double x0, G4double xn,G4int n, G4double f[]);
    void Set_norm();
    void Generate_nCaptureXe136(G4Event* anEvent);
    void SetGeometry(G4int _geometry){fGeometry = _geometry;}
    void Generate_singleelec(G4Event* anEvent);
    void Generate_doubleelec(G4Event* anEvent);
    void Generate_tripleelec(G4Event* anEvent);
    void SetEnergy(G4double energy){fQ_value = energy*keV;}

  private:

    int eventmode;
    int seed;
    bool loadbetadecay;
    G4String fnCaptureXe_Method;
    G4String fXeComponent;
    void GenerateDoubleBetaDecay(G4Event* anEvent);
    G4double GetEnergyFraction(double rand);
    std::vector<double> dbetaenergysplit;
    TRandom3 *rand;
    G4ParticleGun* fParticleGun1;
    G4ParticleGun* fParticleGun3;
    G4GeneralParticleSource* fParticleGun2;
    G4double fB8NeutrinoMaxEnergy;
    G4double fB8NeutrinoMaxY;
    G4double fB8NeutrinoBinWidth;
    G4int fCuIsotope;
    G4double fBb2nCutOffMinFraction;
    G4double fBb2nCutOffMaxFraction;
    G4bool fGenInCenter;
    G4double fRandPos;
    G4double fRandHalfLength;
    G4double fRandRadius;
    G4double electron_mass_c2;
    G4double fQ_value;
    G4double fD_spectral_max;
    G4double fFF_factor;
    G4double fNormalization;
    G4double fK_spectral_max;
    G4String fGenerator;
    LXePrimaryGeneratorMessenger *fMessenger;
    G4Navigator* theNavigator;
    G4double fGeometry;
    
};

#endif
