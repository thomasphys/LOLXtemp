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
#include "LXeMainVolume_Cylindrical.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "G4SystemOfUnits.hh"
#include "TMath.h"
#include "TVector3.h"
#include "LXeMaterials.hh"

LXePMTSD* LXeMainVolume_Cylindrical::SiPM_SD=NULL;

G4RotationMatrix* BuildRotation(G4double theta,G4double phi){

    G4RotationMatrix *RotMat = new G4RotationMatrix();
    RotMat->rotateZ(-phi*rad);
    RotMat->rotateY(pi*rad-theta*rad);
    return RotMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume_Cylindrical::LXeMainVolume_Cylindrical(G4RotationMatrix *pRot,const G4ThreeVector &tlate,G4LogicalVolume *pMotherLogical)
  //Pass info to the G4PVPlacement constructor
  //:G4PVPlacement(pRot,tlate,new G4LogicalVolume(new G4Box("temp",19*mm,19*mm,38*mm),DMaterials::Get_xenon_mat(),"temp",0,0,0),"housing",pMotherLogical,0,0)
{
    
    bool checkOverlaps = true;
    
    G4RotationMatrix *RotMat = new G4RotationMatrix();
    RotMat->rotateX(-pi/2.);
    RotMat->rotateY(pi*(22.5/180.));
    G4RotationMatrix *RotMat2 = new G4RotationMatrix();
    RotMat2->rotateX(-pi/2.);
    RotMat2->rotateY(pi*(22.5/180.));
    RotMat2->rotateZ(-pi);
    //(,z,) 65.11/2.
 /*   CADMesh * mesh = new CADMesh(Form("%s/LoLXCylinder%d.stl",std::getenv("LOLXSTLDIR"),1),mm,G4ThreeVector(-30.08,-20.5,-30.08),false);
    LoLXCylinder = mesh->TessellatedMesh();
    LoLXCylinder_log = new G4LogicalVolume(LoLXCylinder,DMaterials::Get_fPMMA(),"LoLXCylinder_log");
    new G4PVPlacement(RotMat,G4ThreeVector(0.,0.,0.),LoLXCylinder_log,"LoLXCylinder_Physics",pMotherLogical,false,0);
*/
    
    std::ifstream in(Form("%s/%s",std::getenv("LOLXSTLDIR"),"STLOffsets.txt"));
    
    if(!in){
        std::cerr << "Cannot open the File : "<<Form("%s/%s",std::getenv("LOLXSTLDIR"),"STLOffsets.txt")<<std::endl;
        return;
    }
    
    std::string str;
    std::getline(in, str);
    double MiddleX = std::atof(str.data());
    std::getline(in, str);
    double MiddleY = std::atof(str.data());
    std::getline(in, str);
    double MiddleZ = std::atof(str.data());
    std::getline(in, str);
    double TopX = std::atof(str.data());
    std::getline(in, str);
    double TopY = std::atof(str.data());
    std::getline(in, str);
    double TopZ = std::atof(str.data());
    std::getline(in, str);
    double BottomX = std::atof(str.data());
    std::getline(in, str);
    double BottomY = std::atof(str.data());
    std::getline(in, str);
    double BottomZ = std::atof(str.data());
    in.close();
    
    printf("MiddleX = %f MiddleY = %f MiddleZ = %f\n",MiddleX,MiddleY,MiddleZ);
    printf("TopX = %f TopY = %f TopZ = %f\n",MiddleX,MiddleY,MiddleZ);
    printf("BottomX = %f BottomY = %f BottomZ = %f\n",MiddleX,MiddleY,MiddleZ);
    
    CADMesh * meshside = new CADMesh(Form("%s/MiddleOctagon_New.STL",std::getenv("LOLXSTLDIR"),1),mm,G4ThreeVector(MiddleX,MiddleY,MiddleZ),false);
    LoLXCylinderSide = meshside->TessellatedMesh();
    LoLXCylinder_logSide = new G4LogicalVolume(LoLXCylinderSide,DMaterials::Get_fPMMA(),"LoLXCylinder_logSide");
    new G4PVPlacement(RotMat,G4ThreeVector(0.,0.,0.),LoLXCylinder_logSide,"LoLXCylinder_PhysicsSide",pMotherLogical,false,0,checkOverlaps);
    
    CADMesh * meshtop = new CADMesh(Form("%s/TopOctagon_New.STL",std::getenv("LOLXSTLDIR"),1),mm,G4ThreeVector(TopX,TopY,TopZ),false);
    LoLXCylinderTop = meshtop->TessellatedMesh();
    LoLXCylinder_logTop = new G4LogicalVolume(LoLXCylinderTop,DMaterials::Get_fPMMA(),"LoLXCylinder_logTop");
    new G4PVPlacement(RotMat2,G4ThreeVector(0.,0.,20.0),LoLXCylinder_logTop,"LoLXCylinder_PhysicsTop",pMotherLogical,false,0,checkOverlaps);
    
    CADMesh * meshbottom = new CADMesh(Form("%s/BottomOctagon_New.STL",std::getenv("LOLXSTLDIR"),1),mm,G4ThreeVector(BottomX,BottomY,BottomZ),false);
    LoLXCylinderBot = meshbottom->TessellatedMesh();
    LoLXCylinder_logBot = new G4LogicalVolume(LoLXCylinderBot,DMaterials::Get_fPMMA(),"LoLXCylinder_logBot");
    new G4PVPlacement(RotMat2,G4ThreeVector(0.,0.,-20.0),LoLXCylinder_logBot,"LoLXCylinder_PhysicsBottom",pMotherLogical,false,0,checkOverlaps);
    
    G4double xenon_sizeXY = 19*mm; // was 19  and then 14
    G4double xenon_sizeZ  = 38*mm; // was 38  and then 31
    G4double Package_sizeXY = 15*mm;  // size of the  MPPC package initially made from Xenon.
    G4double Package_sizeZ = 2.5*mm;  // everything is housed within this volume
    G4double Package_border = 1.5*mm; // this is the size of the border that is left by the indent
                                      // in the ceramic holding structure.
    G4double Window_gap_scaler = 0.5; // normall 1.5 how much larger the distance from the edge of the Packagelogic will be compared to the Package_border
    
    G4double Value_XY = xenon_sizeXY/2. + Package_sizeZ/2.;
    G4double Value_Z = xenon_sizeZ/2.;
    
    G4double offset = Package_sizeXY*sqrt(2)-1.5*mm;
    G4double zoffset = Package_sizeXY/2.+0.1*mm;
    G4double sqroot2 = sqrt(2);
    G4double octagon_offset = zoffset*2+0.5*mm+Package_sizeZ/2.;
    G4double gap = 0.1*mm;
    G4double Window_sizeZ = 0.97*mm;

    MPPCVolume *SiPM_volume = new MPPCVolume();
    double filterthickness = 0.97*mm;
    double filterwidth = 13.4*mm;
    FilterVolume *filter_volume = new FilterVolume(filterwidth,filterthickness,false);

    std::vector<G4double> theta;
    std::vector<G4double> phi;
    std::vector<G4ThreeVector> pos_sipm;
    std::vector<G4ThreeVector> pos_filter;
    std::vector<G4double> rad;
    std::vector<bool> hasfilter;
    
    G4ThreeVector vert[2];
    vert[0] = G4ThreeVector(0.0,0.0,8.35);
    vert[1] = G4ThreeVector(0.0,0.0,-8.35);
    
    double toprot = -pi*(22.5/180.);
    
    G4ThreeVector quad[4];
    //quad[0] = G4ThreeVector( 7.85, 7.85,0.0);
    //quad[1] = G4ThreeVector( 7.85,-7.85,0.0);
    //quad[2] = G4ThreeVector(-7.85, 7.85,0.0);
    //quad[3] = G4ThreeVector(-7.85, -7.85,0.0);
    quad[0] = G4ThreeVector( 8.4, 8.4,0.0);
    quad[1] = G4ThreeVector( 8.4,-8.4,0.0);
    quad[2] = G4ThreeVector(-8.4, 8.4,0.0);
    quad[3] = G4ThreeVector(-8.4, -8.4,0.0);
    
    //Put this in lookup file
    G4double rad_temp = 20.6*mm + SiPM_volume->GetTopFace();
    rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi/2.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi/2.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(3.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(3.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(pi); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(5.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(5.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(3.*pi/2.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(3.*pi/2.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(7.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(7.*pi/4.); theta.push_back(pi/2.); hasfilter.push_back(true);
    //rad_temp = zoffset*2+0.5*mm+Package_sizeZ/2.;
    rad_temp = 20.2*mm + SiPM_volume->GetTopFace();
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(0.0); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(0.0); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(0.0); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(0.0); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(pi); hasfilter.push_back(false);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(pi); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(pi); hasfilter.push_back(true);
    rad.push_back(rad_temp); phi.push_back(toprot); theta.push_back(pi); hasfilter.push_back(true);

    for(int i=0; i<24; ++i){
        double x = TMath::Sin(theta[i])*TMath::Cos(phi[i]);
        double y = TMath::Sin(theta[i])*TMath::Sin(phi[i]);
        double z = TMath::Cos(theta[i]);
        double sink = -0.0*mm;
        if(i<16){
            double radius = rad[i];
            pos_sipm.push_back(G4ThreeVector(x*radius+vert[i%2].x(),y*radius+vert[i%2].y(),z*radius+vert[i%2].z()));
            radius -=  SiPM_volume->GetTopFace() + Window_sizeZ/2.0;
            pos_filter.push_back(G4ThreeVector(x*radius+vert[i%2].x(),y*radius+vert[i%2].y(),z*radius+vert[i%2].z()));
        }else{
            double radius = rad[i];
            TVector3 pos = TVector3(x*radius+quad[i%4].x(),y*radius+quad[i%4].y(),z*radius+quad[i%4].z());
            pos.Rotate(toprot,TVector3(0.0,0.0,1.0));
            pos_sipm.push_back(G4ThreeVector(pos.x(),pos.y(),pos.z()));
            radius -= SiPM_volume->GetTopFace() + Window_sizeZ/2.0+sink;
            pos = TVector3(x*radius+quad[i%4].x(),y*radius+quad[i%4].y(),z*radius+quad[i%4].z());
            pos.Rotate(toprot,TVector3(0.0,0.0,1.0));
            pos_filter.push_back(G4ThreeVector(pos.x(),pos.y(),pos.z()));
        }
    }
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if(!SiPM_SD){
      SiPM_SD = new LXePMTSD("/LXeDet/SiPM_SD");
      SDman->AddNewDetector(SiPM_SD);
      //Created here so it exists as pmts are being placed
    }
    
    SiPM_SD->InitPMTs(pos_sipm.size()*4);

    for (int i = 0; i < pos_sipm.size(); i++){
        new G4PVPlacement(BuildRotation(theta[i],phi[i]), pos_sipm[i],SiPM_volume,"Package",pMotherLogical,false,i,checkOverlaps);
        G4ThreeVector chippos = SiPM_volume->GetChipPosition(0,BuildRotation(theta[i],phi[i]), pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(1,BuildRotation(theta[i],phi[i]), pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+1,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(2,BuildRotation(theta[i],phi[i]), pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+2,chippos.x(),chippos.y(),chippos.z());
        chippos = SiPM_volume->GetChipPosition(3,BuildRotation(theta[i],phi[i]), pos_sipm[i]);
        SiPM_SD->SetPMTPos(4*i+3,chippos.x(),chippos.y(),chippos.z());
        if(hasfilter[i])new G4PVPlacement(BuildRotation(theta[i],phi[i]), pos_filter[i],filter_volume,"uvfilter1",pMotherLogical,false,i,checkOverlaps);
    }

    SiPM_volume->GetChipLogicVolume()->SetSensitiveDetector(SiPM_SD);
    SurfaceProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume_Cylindrical::SurfaceProperties(){
    //put in proper tables.
    G4OpticalSurface* pmma_opsurf= new G4OpticalSurface("pmma_opsurf",glisur,ground,dielectric_dielectric);
    new G4LogicalSkinSurface("pmma_surf",LoLXCylinder_logSide,pmma_opsurf);
    new G4LogicalSkinSurface("pmma_surf",LoLXCylinder_logTop,pmma_opsurf);
}
