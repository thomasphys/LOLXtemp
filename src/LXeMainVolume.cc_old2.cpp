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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "LXeMainVolume.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "G4SystemOfUnits.hh"
#include "TMath.h"
#include "TVector3.h"

LXePMTSD* LXeMainVolume::SiPM_SD=NULL;

G4LogicalVolume* LXeMainVolume::fHousing_log=NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix* BuildRotation(G4double eta,G4double theta,G4double phi){

    G4RotationMatrix *RotMat = new G4RotationMatrix();
    //RotMat->rotateX(eta*rad);
    RotMat->rotateZ(-phi*rad);
    //RotMat->rotateZ((-eta)*rad);
    //RotMat->rotateY(theta*rad);
    RotMat->rotateY(-theta*rad);
    //RotMat->rotateZ(phi*rad);
    RotMat->rotateZ(eta*rad);
    //RotMat->rotate(eta*rad,G4ThreeVector(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)));
    return RotMat;
}

LXeMainVolume::LXeMainVolume(G4RotationMatrix *pRot,
                             const G4ThreeVector &tlate,
                             G4LogicalVolume *pMotherLogical,
                             G4bool pMany,
                             G4int pCopyNo,
                             LXeDetectorConstruction* c)
  //Pass info to the G4PVPlacement constructor
  :G4PVPlacement(pRot,tlate,
                 //Temp logical volume must be created here
                 //new G4LogicalVolume(new G4Box("temp",1,1,1),
                 //                    G4Material::GetMaterial("Vacuum"),
                 //                    "temp",0,0,0),
                // "housing",pMotherLogical,pMany,pCopyNo),fConstructor(c)
                 new G4LogicalVolume(new G4Box("temp",1,1,1),
                                     c->fVacuum,
                                     "temp",0,0,0),
                 "housing",pMotherLogical,pMany,pCopyNo),fConstructor(c)
{
  CopyValues();
    
    SiPMholder_w = 15.0*mm;
    SiPMholder_h = 2.5*mm;
    SiPMholder_r = 28.355044*mm+SiPMholder_h/2.+0.02;
    SiPMQuartz_w = 12.0*mm;
    SiPMQuartz_h = 0.5*mm;
    SiPMQuartz_r = SiPMholder_r-SiPMholder_h/2.+SiPMQuartz_h/2.;
    SiPMChip_w = 12.0*mm;
    SiPMChip_h = 1.3*mm;
    SiPMChip_r = SiPMholder_r-SiPMholder_h/2.+SiPMQuartz_h+SiPMChip_h/2.;
    filter2_w = 13.4*mm;
    filter2_h = 1.0-0.03*mm;
    filter1_w = 13.4*mm;
    filter1_h = 0.03*mm;
    filter1_r = SiPMholder_r-SiPMholder_h/2.-filter2_h-filter1_h/2.;
    filter2_r = SiPMholder_r-SiPMholder_h/2.-filter2_h/2.;
    
    
    
    std::vector<G4double> rot;
    std::vector<G4double> theta;
    std::vector<G4double> phi;
    std::vector<G4ThreeVector> pos_sipmholder;
    std::vector<G4ThreeVector> pos_sipmquartz;
    std::vector<G4ThreeVector> pos_sipmchip;
    std::vector<G4ThreeVector> pos_filter1;
    std::vector<G4ThreeVector> pos_filter2;
    std::vector<bool> hasfilter;
    
    //Put this in lookup file
    rot.push_back(pi/4.); theta.push_back(pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);

    rot.push_back(pi/4.); theta.push_back(pi*55./180.); phi.push_back(pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*55./180.); phi.push_back(3.*pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*55./180.); phi.push_back(5.*pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*55./180.); phi.push_back(7.*pi/4.); hasfilter.push_back(true);

    rot.push_back(pi/4.); theta.push_back(3.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(3.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(3.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(3.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);

    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(3.*pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(5.*pi/8.); hasfilter.push_back(false);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(7.*pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(9.*pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(11.*pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(13.*pi/8.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi/2.); phi.push_back(15.*pi/8.); hasfilter.push_back(true);

    rot.push_back(pi/4.); theta.push_back(5.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(5.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(5.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(5.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);


    rot.push_back(pi/4.); theta.push_back(pi*125./180.); phi.push_back(pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*125./180.); phi.push_back(3.*pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*125./180.); phi.push_back(5.*pi/4.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(pi*125./180.); phi.push_back(7.*pi/4.); hasfilter.push_back(true);

    rot.push_back(pi/4.); theta.push_back(7.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(7.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(7.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    rot.push_back(pi/4.); theta.push_back(7.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);
    
    n_sipm = rot.size();
    
    for(int i=0; i<n_sipm; ++i){
        pos_sipmholder.push_back(G4ThreeVector(SiPMholder_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                         SiPMholder_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                         SiPMholder_r*TMath::Cos(theta[i])));
        pos_sipmquartz.push_back(G4ThreeVector(SiPMQuartz_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                         SiPMQuartz_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                         SiPMQuartz_r*TMath::Cos(theta[i])));
        pos_sipmchip.push_back(G4ThreeVector(SiPMChip_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                         SiPMChip_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                         SiPMChip_r*TMath::Cos(theta[i])));
        pos_filter1.push_back(G4ThreeVector(filter1_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                           filter1_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                           filter1_r*TMath::Cos(theta[i])));
        pos_filter2.push_back(G4ThreeVector(filter2_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                           filter2_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                           filter2_r*TMath::Cos(theta[i])));
    }

    G4double housing_x=fScint_x+fD_mtl;
    G4double housing_y=fScint_y+fD_mtl;
    G4double housing_z=fScint_z+fD_mtl;
 
    //*************************** housing and scintillator
    fScint_box = new G4Box("scint_box",fScint_x/2.,fScint_y/2.,fScint_z/2.);
    fHousing_box = new G4Box("housing_box",housing_x/2.,housing_y/2.,
                            housing_z/2.);
 
    //fScint_log = new G4LogicalVolume(fScint_box,G4Material::GetMaterial("LXe"),"scint_log",0,0,0);
    fScint_log = new G4LogicalVolume(fScint_box,c->fLXe,"scint_log",0,0,0);
    //fHousing_log = new G4LogicalVolume(fHousing_box,G4Material::GetMaterial("Al"),"housing_log",0,0,0);
    fHousing_log = new G4LogicalVolume(fHousing_box,c->fAl, "housing_log",0,0,0);
 
    //*************** Miscellaneous sphere to demonstrate skin surfaces
 
    CADMesh * mesh = new CADMesh("/Users/thomasmcelroy/lolxsim/data/LoLXSphere.stl",mm,G4ThreeVector(-32.35,-32.35,-32.35),false);
    LoLXSphere = mesh->TessellatedMesh();
    //LoLXSphere_log = new G4LogicalVolume(LoLXSphere,G4Material::GetMaterial("PMMA"),"LoLXSphere_log");
    LoLXSphere_log = new G4LogicalVolume(LoLXSphere,c->fPMMA,"LoLXSphere_log");
   // LoLXSphere_log = new G4LogicalVolume(LoLXSphere,c->fLXe,"LoLXSphere_log");
    new G4PVPlacement(0,G4ThreeVector(),LoLXSphere_log,"LoLXSphere_Physics",fScint_log,false,0);
 
    //****************** Build SiPMs and filters
    SiPMChip = new G4Box("sipmChip",SiPMChip_w/2.,SiPMChip_w/2.,SiPMChip_h/2.);
    SiPMQuartz = new G4Box("sipmQuartz",SiPMQuartz_w/2.,SiPMQuartz_w/2.,SiPMQuartz_h/2.);
    
    G4Box* SiPM_sub = new G4Box("sipm_sub",SiPMChip_w/2.,SiPMChip_w/2.,(SiPMChip_h+SiPMQuartz_h)/2.);
    G4Box* SiPM_main = new G4Box("sipm_min",SiPMholder_w/2.,SiPMholder_w/2.,SiPMholder_h/2.);
    
    SiPMholder = new G4SubtractionSolid("SiPMholder",SiPM_main,SiPM_sub,0,G4ThreeVector(0.,0.,-SiPMholder_h/2.+(SiPMChip_h+SiPMQuartz_h)/2.));
    
    filter1 = new G4Box("sipm",filter1_w/2.,filter1_w/2.,filter1_h/2.);
    filter2 = new G4Box("sipm",filter2_w/2.,filter2_w/2.,filter2_h/2.);
 
    SiPMChip_log = new G4LogicalVolume(SiPMChip,c->fSi,"SiPMChip_log");
    SiPMQuartz_log = new G4LogicalVolume(SiPMQuartz,c->fQuartz,"SiPMQuartz_log");
    SiPMholder_log = new G4LogicalVolume(SiPMholder,c->fCeramic,"SiPMholder_log");
    
    filter1_log = new G4LogicalVolume(filter1,c->fGlass,"filter1_log");
    filter2_log = new G4LogicalVolume(filter2,c->fGlass,"filter2_log");
 
    //***********Arrange pmts around the outside of housing**********
    //---pmt sensitive detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if(!SiPM_SD){
      SiPM_SD = new LXePMTSD("/LXeDet/SiPM_SD");
      SDman->AddNewDetector(SiPM_SD);
      //Created here so it exists as pmts are being placed
    }
    SiPM_SD->InitPMTs(32); //let pmtSD know # of pmts
    //-------
    Printf("Setting SiPM locations");
    
   // G4OpticalSurface* OpSurface = new G4OpticalSurface("name");
   // OpSurface -> SetType(dielectric_dichroic);
   // OpSurface -> SetModel(unified);
   // OpSurface -> SetFinish(polished);
    
    int k=0;
    for(int i=0; i<n_sipm; ++i){
        new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_sipmholder[i],SiPMholder_log,"sipmtholder",fScint_log,false,i);
        new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_sipmquartz[i],SiPMQuartz_log,"sipmtquartz",fScint_log,false,i);
        new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_sipmchip[i],SiPMChip_log,"sipmtchip",fScint_log,false,i);
    
        SiPM_SD->SetPMTPos(i,pos_sipmchip[i].x(),pos_sipmchip[i].y(),pos_sipmchip[i].z());
        
        if(hasfilter[i]){
            G4PVPlacement* filter1pv = new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_filter1[i],filter1_log,"uvfilter1",fScint_log,false,k++);
            G4PVPlacement* filter2pv = new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_filter2[i],filter2_log,"uvfilter2",fScint_log,false,k++);

            //G4LogicalBorderSurface* Surface = new G4LogicalBorderSurface("filteroptics",filter1pv,filter2pv,OpSurface);
  
        }
    }
    
    new G4PVPlacement(0,G4ThreeVector(),fScint_log,"scintillator",
                                   fHousing_log,false,0);
    
    //**********Setup Sensitive Detectors***************
    if(!fScint_SD){//determine if it has already been created
       Printf("Setting Sensitive Detector");
      fScint_SD = new LXeScintSD("/LXeDet/scintSD");
      SDman->AddNewDetector(fScint_SD);
    }
    Printf("Setting LXe sensitive volume");
    fScint_log->SetSensitiveDetector(fScint_SD);
 
    SiPMChip_log->SetSensitiveDetector(SiPM_SD);
    //SiPMQuartz_log->SetSensitiveDetector(SiPM_SD);
    VisAttributes();
    SurfaceProperties();
    Printf("Done");
   
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  housing_va->SetVisibility(false);
  fHousing_log->SetVisAttributes(housing_va);

  SetLogicalVolume(fHousing_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues(){
  fUpdated=fConstructor->GetUpdated();

  fScint_x=fConstructor->GetScintX();
  fScint_y=fConstructor->GetScintY();
  fScint_z=fConstructor->GetScintZ();
  fD_mtl=fConstructor->GetHousingThickness();
  fRefl=fConstructor->GetHousingReflectivity();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes(){
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  fHousing_log->SetVisAttributes(housing_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::SurfaceProperties(){
  //put in proper tables.
  const G4int num = 2;
  G4double ephoton[num] = {7.0*eV, 7.14*eV};

  //**Scintillator housing properties
  G4double reflectivity[num] = {0.0, 0.0};
  G4double efficiency[num] = {0.0, 0.0};
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency, num);
  G4OpticalSurface* OpScintHousingSurface = new G4OpticalSurface("HousingSurface",unified,polished,dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);
 
  G4double nm2ev = 1240.;
  
  G4double sipm_EFF_wavelength[40] = {060.0,121.6,124.6,127.6,129.4,131.0,132.4,134.8,137.1,138.3,139.6,141.0,143.0,144.5,146.8,148.8,152.0,153.2,154.8,157.0,
                                      158.6,159.5,160.7,161.8,162.7,163.7,164.3,165.3,166.4,167.8,169.4,171.0,173.2,176.4,180.3,185.2,191.1,197.2,200.0,700.0};

  G4double sipm_EFF_Energy[40];
  for(int i=0; i<40; ++i) sipm_EFF_Energy[39-i] = nm2ev*eV/sipm_EFF_wavelength[i];
    
  G4double sipm_EFF[40]={0.000,0.107,0.123,0.141,0.156,0.166,0.176,0.192,0.207,0.214,0.219,0.226,0.236,0.239,0.239,0.237,0.233,0.230,0.225,0.209,
                         0.194,0.187,0.185,0.188,0.199,0.219,0.227,0.233,0.236,0.236,0.235,0.236,0.238,0.238,0.238,0.238,0.238,0.238,0.238,0.238};
  G4double sipm_EFF_EnergyOrder[40];
  for(int i=0; i<40; ++i) sipm_EFF_EnergyOrder[39-i] = 0.96;//sipm_EFF[i];
                         
  G4double sipm_ReR[40];
  G4double sipm_ImR[40];
  
  for(int i=0; i<40; ++i){
    sipm_ReR[i]=0.0;//1.92;
    sipm_ImR[i]=1.69;
  }
  
  G4MaterialPropertiesTable* sipm_mt = new G4MaterialPropertiesTable();
  sipm_mt->AddProperty("EFFICIENCY",sipm_EFF_Energy,sipm_EFF_EnergyOrder,40);
  sipm_mt->AddProperty("REFLECTIVITY",sipm_EFF_Energy,sipm_ReR,40);
  //sipm_mt->AddProperty("REALRINDEX",sipm_EFF_Energy,sipm_ReR,40);
  //sipm_mt->AddProperty("IMAGINARYRINDEX",sipm_EFF_Energy,sipm_ImR,40);
  G4OpticalSurface* sipm_opsurf= new G4OpticalSurface("sipm_opsurf",glisur,polished,dielectric_metal);
  sipm_opsurf->SetMaterialPropertiesTable(sipm_mt);
  
  G4OpticalSurface* pmma_opsurf= new G4OpticalSurface("sipm_opsurf",glisur,ground,dielectric_dielectric);
  G4OpticalSurface* quartz_opsurf= new G4OpticalSurface("quartz_opsurf",glisur,polished,dielectric_dielectric);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("housing_surf",fHousing_log,OpScintHousingSurface);
  //new G4LogicalSkinSurface("sipm_surf",SiPMChip_log,sipm_opsurf);
  new G4LogicalSkinSurface("sipm_surf",SiPMQuartz_log,sipm_opsurf);
  new G4LogicalSkinSurface("pmma_surf",LoLXSphere_log,pmma_opsurf);
  //new G4LogicalSkinSurface("Quartz_surf",SiPMQuartz_log,quartz_opsurf);
  
}
