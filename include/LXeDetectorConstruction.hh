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
/// \file optical/LXe/include/LXeDetectorConstruction.hh
/// \brief Definition of the LXeDetectorConstruction class
//
//
#ifndef LXeDetectorConstruction_H
#define LXeDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class LXePMTSD;
class LXeScintSD;

#include "G4Material.hh"
#include "LXeDetectorMessenger.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "LXeScintSD.hh"
#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume_Sphere.hh"
#include "LXeMainVolume_Cylindrical.hh"


#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "TMath.h"
#include "TH1.h"
#include "TFile.h"

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    LXeDetectorConstruction();
    virtual ~LXeDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    //Functions to modify the geometry
    void SetDimensions(G4ThreeVector );
    void SetHousingThickness(G4double );
    void SetDefaults();

    //Get values
    G4double GetScintX(){return fScint_x;}
    G4double GetScintY(){return fScint_y;}
    G4double GetScintZ(){return fScint_z;}
    G4double GetHousingThickness(){return fD_mtl;}
 
    //rebuild the geometry based on changes. must be called
    void UpdateGeometry();
    G4bool GetUpdated(){return fUpdated;}

    void SetHousingReflectivity(G4double r){fRefl=r; fUpdated=true;}
    G4double GetHousingReflectivity(){return fRefl;}

    void SetMainVolumeOn(G4bool b){fMainVolume=b; fUpdated=true;}
    G4bool GetMainVolumeOn(){return fMainVolume;}

    void SetMainScintYield(G4double );
    void SetGeometry(G4int geo){fGeometry = geo;}

  private:

    G4VPhysicalVolume* ConstructDetector();
    
    LXeDetectorMessenger* fDetectorMessenger;

    G4bool fUpdated;
 
    G4Box* fExperimentalHall_box;
    G4LogicalVolume* fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;

    //Geometry
    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;
    G4double fD_mtl;
    G4double fRefl;
    G4bool fMainVolume;
    G4int fGeometry;

    G4MaterialPropertiesTable* fLXe_mt;
    static LXeScintSD* fScint_SD;

};

#endif
