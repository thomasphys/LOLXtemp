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
#include "MPPCVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MPPCVolume::Package_sizeXY = 15*mm;  // size of the  MPPC package initially made from Xenon.
G4double MPPCVolume::Package_sizeZ = 2.5*mm;  // everything is housed within this volume
G4double MPPCVolume::Package_border = 1.1*mm;


MPPCVolume::MPPCVolume()
  //Pass info to the G4PVPlacement constructor
  :G4LogicalVolume(new G4Box("Package", Package_sizeXY/2, Package_sizeXY/2, Package_sizeZ/2),DMaterials::Get_xenon_mat(), "Package")
{

    bool checkOverlaps = true;
    
    // create a logical volume, where will put everything

    G4double MPPC_sizeXY = 5.9*mm; // size of individual MPPC, which is 5.85 x 5.95 in actuallity (approximated by a square of 5.9x5.9 here)
    G4double MPPC_sizeZ = 0.0025*mm; // thickness of MPPC
    // position of MPPC arrays inside the package (Packagelogic)

    MPPC_pos[0] = G4ThreeVector(    (MPPC_sizeXY+0.5*mm)/2.,     (MPPC_sizeXY+0.5*mm)/2., MPPC_sizeZ/2.);
    MPPC_pos[1] = G4ThreeVector(    (MPPC_sizeXY+0.5*mm)/2., -1.*(MPPC_sizeXY+0.5*mm)/2., MPPC_sizeZ/2.);
    MPPC_pos[2] = G4ThreeVector(-1.*(MPPC_sizeXY+0.5*mm)/2., -1.*(MPPC_sizeXY+0.5*mm)/2., MPPC_sizeZ/2.);
    MPPC_pos[3] = G4ThreeVector(-1.*(MPPC_sizeXY+0.5*mm)/2.,     (MPPC_sizeXY+0.5*mm)/2., MPPC_sizeZ/2.);

    //MPPC
    G4Box* SolidMPPC = new G4Box("MPPCsold", MPPC_sizeXY/2., MPPC_sizeXY/2., MPPC_sizeZ/2.);
    MPPClogic = new G4LogicalVolume(SolidMPPC, DMaterials::Get_silicon_mat(), "MPPClog"); // MPPC, made out of silicon

    // To form a ceramic holding structure need to do a boolean subtraction of volumes
    // first creating two boxes, one the size of the Packagelogic, the other is half depth (z) and is smaller by the Package_border

    G4double delta = 0.0*mm;

    G4Box* box1 = new G4Box("Box1",Package_sizeXY/2.,Package_sizeXY/2.,Package_sizeZ/2.); // size of the Packagelogic
    G4Box* box2 = new G4Box("Box2",(Package_sizeXY-Package_border)/2.,(Package_sizeXY-Package_border)/2.,Package_sizeZ/4.+delta); // half depth
                                                                                               // and smaller in size by Package_border
    G4ThreeVector  translation_ceramic(0.0,0.0,(Package_sizeZ/4.)+delta/2.); // this vectro is needed for volume subtraction
    G4SubtractionSolid* b1minusb2 = new G4SubtractionSolid("box1-box2",box1,box2,0,translation_ceramic); // creating a solid for ceramic holding structure
    // by subtracting one G4Box from another G4Box.
    //Ceramiclogic = new G4LogicalVolume(b1minusb2, DMaterials::Get_silicon_mat(), "ceramic"); // logical volume for the ceramic holding structure
    //Ceramiclogic = new G4LogicalVolume(box1, DMaterials::Get_xenon_mat(), "ceramic");
    Ceramiclogic = new G4LogicalVolume(b1minusb2, DMaterials::Get_fCeramic(), "ceramic");
   // new G4PVPlacement(0, G4ThreeVector(), Ceramiclogic, "Ceramic", this,0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), Ceramiclogic, "Ceramic", this,0, checkOverlaps);

    // creating a quarts window
    G4double Window_sizeZ = 0.5*mm; // window thickness  nominally was .5mm. I changed it to .25mm
    G4double Window_gap_scaler = 1.; // normall 1.5 how much larger the distance from the edge of the Packagelogic will be compared to the Package_border
    G4double Window_sizeXY = Package_sizeXY-Window_gap_scaler*Package_border;

    G4Box* window_solid = new G4Box("QuartzWindow",Window_sizeXY/2.,Window_sizeXY/2.,Window_sizeZ/2.);
    G4LogicalVolume* Windowlogic = new G4LogicalVolume(window_solid, DMaterials::Get_quartz_mat(), "Window");

    new G4PVPlacement(0, MPPC_pos[0], MPPClogic, "MPPC", this,false,0, checkOverlaps); // individual MPPCs
    new G4PVPlacement(0, MPPC_pos[1], MPPClogic, "MPPC", this,false,1, checkOverlaps);
    new G4PVPlacement(0, MPPC_pos[2], MPPClogic, "MPPC", this,false,2, checkOverlaps);
    new G4PVPlacement(0, MPPC_pos[3], MPPClogic, "MPPC", this,false,3, checkOverlaps);
 
    WindowFaceZ = Window_sizeZ/2.+MPPC_sizeZ;
    
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,WindowFaceZ), Windowlogic, "Window", this,0, checkOverlaps);
    SurfaceProperties();
}

G4ThreeVector MPPCVolume::GetChipPosition(G4int num,G4RotationMatrix* rot,G4ThreeVector trans){

    G4double x = MPPC_pos[num].x()*rot->xx()+MPPC_pos[num].y()*rot->xy()+MPPC_pos[num].z()*rot->xz()+trans.x();
    G4double y = MPPC_pos[num].x()*rot->yx()+MPPC_pos[num].y()*rot->yy()+MPPC_pos[num].z()*rot->yz()+trans.y();
    G4double z = MPPC_pos[num].x()*rot->zx()+MPPC_pos[num].y()*rot->zy()+MPPC_pos[num].z()*rot->zz()+trans.z();

//	printf("%f %f %f\n%f %f %f\n%f %f %f\n",rot->xx(),rot->xy(),rot->xz(),rot->yx(),rot->yy(),rot->yz(),rot->zx(),rot->zy(),rot->zz());

    //G4double x = MPPC_pos[num].x()*rot->xx()+MPPC_pos[num].y()*rot->yx()+MPPC_pos[num].z()*rot->zx()+trans.x();
    //G4double y = MPPC_pos[num].x()*rot->xy()+MPPC_pos[num].y()*rot->yy()+MPPC_pos[num].z()*rot->zy()+trans.y();
    //G4double z = MPPC_pos[num].x()*rot->xz()+MPPC_pos[num].y()*rot->yz()+MPPC_pos[num].z()*rot->zz()+trans.z();

    return G4ThreeVector(x,y,z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MPPCVolume::SurfaceProperties(){

    const G4int NUMENTRIES1 = 850;  //from 150 to 1000 
    G4double LXe_PP1[NUMENTRIES1];  //energies
    G4double LXe_SCINT1[NUMENTRIES1]; //emission probability
    G4double LXe_RIND1[NUMENTRIES1]; // Refraction index
    G4double LXe_ABSL1[NUMENTRIES1]; // Absorption length
    G4double LXe_Rayleigh1[NUMENTRIES1]; // Rayleigh scattering length
    
    for(int iE=0; iE<NUMENTRIES1; iE++){
      G4double WL=1000-(iE); // starting from 1000nm WL, ending with 150nm  


      LXe_PP1[iE]=1240/(WL)*eV; // tabulate energy
      LXe_SCINT1[iE] = 1./sqrt(2*3.14)/5.*exp(-pow((WL-178.)/5.,2.)/2.);//approximating with a gaus (178,5)
      LXe_RIND1[iE] = sqrt(1.5+
			0.38*WL*WL/(WL*WL-146.9*146.9)+
			0.009*WL*WL/(WL*WL-827*827));//from https://arxiv.org/pdf/1502.04213.pdf  AS Feb 21, 2018
      if(WL>700) LXe_RIND1[iE]=1.365;// avoid the fall of ref index in the formula
  //    if(WL<150) LXe_RIND1[iE]=1.7;
      LXe_ABSL1[iE] = 200.*cm;
      LXe_Rayleigh1[iE] = 35.*cm; // need to code proper formula but should not be critical for LoLX
      G4cout << " WL= " << WL<< " and energy="<<LXe_PP1[iE]<< " iE="<<iE<<G4endl; //debug output
    }


    G4OpticalSurface* OpMPPCSurface = new G4OpticalSurface("OpMPPCSurface");
    OpMPPCSurface->SetModel(glisur); //glisur
    OpMPPCSurface->SetType(dielectric_metal); //_metal
    OpMPPCSurface->SetFinish(polished);

    G4double reflectivity[NUMENTRIES1];for(double &r: reflectivity) r=0.; // non Reflective
    G4double efficiency[NUMENTRIES1];for(double &r: efficiency) r=1.0; // Perfect efficiency
    G4MaterialPropertiesTable* MPPCTable = new G4MaterialPropertiesTable();
    MPPCTable->AddProperty("REFLECTIVITY",LXe_PP1,reflectivity,NUMENTRIES1);
    MPPCTable->AddProperty("EFFICIENCY",LXe_PP1,efficiency,NUMENTRIES1);
    OpMPPCSurface->SetMaterialPropertiesTable(MPPCTable);
    new G4LogicalSkinSurface("OpMPPCSurface",MPPClogic,OpMPPCSurface);

    G4OpticalSurface* OpCeramicSurface = new G4OpticalSurface("OpCeramicSurface");
    OpCeramicSurface->SetModel(glisur);
    OpCeramicSurface->SetType(dielectric_metal);
    OpCeramicSurface->SetFinish(ground);

    G4double reflectivity_ceramic[NUMENTRIES1];for(double &r: reflectivity_ceramic) r=0.8; // 80% Reflective ceramic MPPC holding structure
    G4double efficiency_ceramic[NUMENTRIES1];for(double &r: efficiency_ceramic) r=0; // zero efficiency
    G4MaterialPropertiesTable* CeramicTable = new G4MaterialPropertiesTable();
    CeramicTable->AddProperty("REFLECTIVITY",LXe_PP1,reflectivity_ceramic,NUMENTRIES1);
    CeramicTable->AddProperty("EFFICIENCY",LXe_PP1,efficiency_ceramic,NUMENTRIES1);
    OpCeramicSurface->SetMaterialPropertiesTable(CeramicTable);
    new G4LogicalSkinSurface("OpCeramicSurface",Ceramiclogic,OpCeramicSurface); //This surface is for the Ceramic, hence for Ceramiclogic
  
}
