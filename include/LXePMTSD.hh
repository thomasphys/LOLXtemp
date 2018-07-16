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
/// \file optical/LXe/include/LXePMTSD.hh
/// \brief Definition of the LXePMTSD class
//
//
#ifndef LXePMTSD_h
#define LXePMTSD_h 1

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "LXePMTHit.hh"

class G4Step;
class G4HCofThisEvent;

class LXePMTSD : public G4VSensitiveDetector
{

  public:

    LXePMTSD(G4String name);
    virtual ~LXePMTSD();
 
    virtual void Initialize(G4HCofThisEvent* );
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* );
 
    //A version of processHits that keeps aStep constant
    G4bool ProcessHits_constStep(const G4Step* ,
                                 G4TouchableHistory* );
    virtual void EndOfEvent(G4HCofThisEvent* );
    virtual void clear();
    void DrawAll();
    void PrintAll();
 
    //Initialize the arrays to store pmt possitions
    inline void InitPMTs(G4int nSiPMs){
      if(SiPMPositionsX)delete SiPMPositionsX;
      if(SiPMPositionsY)delete SiPMPositionsY;
      if(SiPMPositionsZ)delete SiPMPositionsZ;
      SiPMPositionsX=new G4DataVector(nSiPMs);
      SiPMPositionsY=new G4DataVector(nSiPMs);
      SiPMPositionsZ=new G4DataVector(nSiPMs);
    }

    //Store a pmt position
    inline void SetPMTPos(G4int n,G4double x,G4double y,G4double z){
      if(SiPMPositionsX)SiPMPositionsX->insertAt(n,x);
      if(SiPMPositionsY)SiPMPositionsY->insertAt(n,y);
      if(SiPMPositionsZ)SiPMPositionsZ->insertAt(n,z);
    }

  private:

    LXePMTHitsCollection* SiPMHitCollection;

    G4DataVector* SiPMPositionsX;
    G4DataVector* SiPMPositionsY;
    G4DataVector* SiPMPositionsZ;
};

#endif
