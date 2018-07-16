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
#ifndef FilterVolume_H
#define FilterVolume_H 1

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SubtractionSolid.hh"
#include "LXeDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SystemOfUnits.hh"
#include "LXeMaterials.hh"
#include "G4PVPlacement.hh"

class FilterVolume : public G4LogicalVolume
{
  public:

    FilterVolume(G4double Window_sizeXY,G4double Window_sizeZ,bool boarder=false);

  private:

    void SurfaceProperties();

    G4bool fUpdated;

};

#endif
