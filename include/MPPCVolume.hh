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
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef MPPCVolume_H
#define MPPCVolume_H 1

#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "LXeMaterials.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SystemOfUnits.hh"

class MPPCVolume : public G4LogicalVolume
{
  public:

    MPPCVolume();
    
    G4LogicalVolume* GetChipLogicVolume(){return MPPClogic;}
    G4ThreeVector GetChipPosition(G4int chipnum,G4RotationMatrix* rotation,G4ThreeVector translation);
    G4double GetWindowFaceZ(){return WindowFaceZ;}
    G4double GetTopFace(){return Package_sizeZ/2.;}

  private:

    void VisAttributes();
    void SurfaceProperties();

    void CopyValues();

    G4bool fUpdated;
    
    G4LogicalVolume* MPPClogic;
    G4ThreeVector MPPC_pos[4];
    G4LogicalVolume* Ceramiclogic;
    
    static G4double Package_sizeXY;  // size of the  MPPC package initially made from Xenon.
    static G4double Package_sizeZ;  // everything is housed within this volume
    static G4double Package_border;
    G4double WindowFaceZ;


};

#endif
