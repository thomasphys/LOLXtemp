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
/// \file optical/LXe/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeScintSD* LXeDetectorConstruction::fScint_SD=NULL;

LXeDetectorConstruction::LXeDetectorConstruction()
: fLXe_mt(NULL)
{

  fUpdated = false;
  
  fGeometry = 0;

  SetDefaults();

  fDetectorMessenger = new LXeDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct(){
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  G4double world_sizeXY = 88*mm;
  G4double world_sizeZ  = 88*mm;
  bool checkOverlaps = true;
  
  G4Box* solidWorld =
    new G4Box("World",              //its name
              0.5*world_sizeXY,
              0.5*world_sizeXY,
              0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        DMaterials::Get_fe_mat(),//DMaterials::Get_world_mat(),           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
    
  G4double xenon_sizeXY = 78*mm; //was 59 // was 19
  G4double xenon_sizeZ  = 78*mm; // was 38
  
  G4Box* solidXenon = new G4Box("Xenon",0.5*xenon_sizeXY,0.5*xenon_sizeXY,0.5*xenon_sizeZ);
  G4LogicalVolume* logicXenon = new G4LogicalVolume(solidXenon,DMaterials::Get_xenon_mat(),"Xenon");

  new G4PVPlacement(0,G4ThreeVector(),logicXenon,"Xenon",logicWorld,false,0,checkOverlaps);

  
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  logicXenon->SetVisAttributes(G4VisAttributes::Invisible);

//put in two geometries
  //Place the main volume
  if(fGeometry ==  1){
    new LXeMainVolume_Sphere(0,G4ThreeVector(),logicXenon);
  }else{
    new LXeMainVolume_Cylindrical(0,G4ThreeVector(),logicXenon);
  }
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if(!fScint_SD){//determine if it has already been created
      fScint_SD = new LXeScintSD("/LXeDet/scintSD");
      SDman->AddNewDetector(fScint_SD);
  }
  logicXenon->SetSensitiveDetector(fScint_SD);
  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims){
  this->fScint_x=dims[0];
  this->fScint_y=dims[1];
  this->fScint_z=dims[2];
  fUpdated=true;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl){
  this->fD_mtl=d_mtl;
  fUpdated=true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDefaults(){
  //Resets to default values
  fD_mtl=0.0635*cm;

  fScint_x = 10.0*cm;
  fScint_y = 10.0*cm;
  fScint_z = 10.0*cm;

  fRefl=1.0;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/LXe/detector/scintYieldFactor 1.");

  if(fLXe_mt)fLXe_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV);

  fUpdated=true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::UpdateGeometry(){

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();
  G4SurfaceProperty::CleanSurfacePropertyTable();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  fUpdated=false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y){
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD",y/MeV);
}
