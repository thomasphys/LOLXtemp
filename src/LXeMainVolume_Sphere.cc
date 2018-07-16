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
#include "LXeMainVolume_Sphere.hh"

LXePMTSD* LXeMainVolume_Sphere::SiPM_SD=NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix* BuildRotation(G4double eta,G4double theta,G4double phi){

    G4RotationMatrix *RotMat = new G4RotationMatrix();
    RotMat->rotateZ(pi*rad-phi*rad);
    //RotMat->rotateY(-theta*rad);
    RotMat->rotateY(pi*rad+theta*rad);
    RotMat->rotateZ(eta*rad);
    return RotMat;
}

LXeMainVolume_Sphere::LXeMainVolume_Sphere(G4RotationMatrix *pRot,const G4ThreeVector &tlate,G4LogicalVolume *pMotherLogical)
  //:G4PVPlacement(pRot,tlate,new G4LogicalVolume(new G4Box("temp",1,1,1),DMaterials::Get_xenon_mat(),"temp",0,0,0), "housing",pMotherLogical,0,0)
{
    
    std::vector<G4double> rot;
    std::vector<G4double> theta;
    std::vector<G4double> phi;
    std::vector<G4ThreeVector> pos_sipm;
    std::vector<G4ThreeVector> pos_filter;
    std::vector<bool> hasfilter;
    
    MPPCVolume *SiPM_volume = new MPPCVolume();
    FilterVolume *filter_volume = new FilterVolume(13.4*mm,0.97*mm);
    
    double SiPMholder_r = 28.355044*mm + SiPM_volume->GetTopFace()+0.1;
    double filter_r = SiPMholder_r - SiPM_volume->GetTopFace() - 1.0*mm/2.0;
    
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
        pos_sipm.push_back(G4ThreeVector(SiPMholder_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                         SiPMholder_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                         SiPMholder_r*TMath::Cos(theta[i])));
        pos_filter.push_back(G4ThreeVector(filter_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                           filter_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                           filter_r*TMath::Cos(theta[i])));
    }
 
    CADMesh * mesh = new CADMesh("/home/tmcelroy/lolxsim/data/LoLXSphere.stl",mm,G4ThreeVector(-32.35,-32.35,-32.35),false);
    LoLXSphere = mesh->TessellatedMesh();
    LoLXSphere_log = new G4LogicalVolume(LoLXSphere,DMaterials::Get_fPMMA(),"LoLXSphere_log");
    new G4PVPlacement(0,G4ThreeVector(),LoLXSphere_log,"LoLXSphere_Physics",pMotherLogical,false,0);
 
    //***********Arrange pmts around the outside of housing**********
    //---pmt sensitive detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if(!SiPM_SD){
      SiPM_SD = new LXePMTSD("/LXeDet/SiPM_SD");
      SDman->AddNewDetector(SiPM_SD);
      //Created here so it exists as pmts are being placed
    }
    SiPM_SD->InitPMTs(n_sipm*4); //let pmtSD know # of pmts
    //-------
    Printf("Setting SiPM locations");
    
    int k=0;
    for(int i=0; i<n_sipm; ++i){
    //for(int i=0; i<1; ++i){
       // if(i<16 or i>20) continue;
        new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_sipm[i],SiPM_volume,"sipmtholder",pMotherLogical,false,i,true);

        G4ThreeVector chippos = SiPM_volume->GetChipPosition(0,BuildRotation(rot[i],theta[i],phi[i]),pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(1,BuildRotation(rot[i],theta[i],phi[i]),pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+1,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(2,BuildRotation(rot[i],theta[i],phi[i]),pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+2,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(3,BuildRotation(rot[i],theta[i],phi[i]),pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+3,chippos.x(),chippos.y(),chippos.z());

        if(hasfilter[i]){
            new G4PVPlacement(BuildRotation(rot[i],theta[i],phi[i]),pos_filter[i],filter_volume,"uvfilter1",pMotherLogical,false,k++,true);
        }
    }
 
    SiPM_volume->GetChipLogicVolume()->SetSensitiveDetector(SiPM_SD);
    SurfaceProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume_Sphere::SurfaceProperties(){
  //put in proper tables.
  G4OpticalSurface* pmma_opsurf= new G4OpticalSurface("sipm_opsurf",glisur,ground,dielectric_dielectric);
  new G4LogicalSkinSurface("pmma_surf",LoLXSphere_log,pmma_opsurf);
  
}
