#include "FilterVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FilterVolume::FilterVolume(G4double Window_sizeXY,G4double Window_sizeZ,bool boarder)
  //Pass info to the G4PVPlacement constructor
  :G4LogicalVolume(new G4Box("temp",Window_sizeXY/2.,Window_sizeXY/2.,Window_sizeZ/2.),DMaterials::Get_xenon_mat(),"temp",0,0,0)
{
    bool checkOverlaps = true;
    
    G4double Package_sizeZ = 2.5*mm;  // everything is housed within this volume
    G4double Package_border = 1.5*mm;  // this is the size of the border that is left by the indent
                                    // in the ceramic holding str
    G4double film_thickness = 0.1*mm;
    if(boarder) Window_sizeXY -= Package_border;
    
    G4Box* film_solid = new G4Box("",Window_sizeXY/2.,Window_sizeXY/2.,film_thickness/2.);
    G4Box* window_solid = new G4Box("QuartzWindow",Window_sizeXY/2.,Window_sizeXY/2.,(Window_sizeZ-film_thickness)/2.);

    G4LogicalVolume* FilterUSlogic = new G4LogicalVolume(window_solid, DMaterials::Get_xenon_mat(), "Window");
    G4LogicalVolume* FilterDSlogic = new G4LogicalVolume(film_solid, DMaterials::Get_xenon_mat(), "Window");
    
    G4VPhysicalVolume* volume1 =  new G4PVPlacement(0,G4ThreeVector(0,0,-film_thickness/2.),
                                                    FilterUSlogic,"WindowUS",this,0,checkOverlaps);
    
    G4VPhysicalVolume* volume2 =  new G4PVPlacement(0,G4ThreeVector(0,0,(Window_sizeZ-film_thickness)/2.),
                                                    FilterDSlogic,"WindowDS",this,0,checkOverlaps);

   if(boarder){
	printf("adding boarder\n");
	G4Box* box1 = new G4Box("Box1",Window_sizeXY/2.+0.05*mm,Window_sizeXY/2.+0.05*mm,Window_sizeZ/2.);
    	G4Box* box2 = new G4Box("Box2",(Window_sizeXY+Package_border)/2.,(Window_sizeXY+Package_border)/2.,Window_sizeZ/2.);
    	G4ThreeVector  translation(0.0,0.0,0.0);
    	G4SubtractionSolid* b1minusb2 = new G4SubtractionSolid("box1-box2",box2,box1,0,translation);
	G4LogicalVolume* Boarderlogic = new G4LogicalVolume(b1minusb2, DMaterials::Get_fPMMA(), "boarderlog");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Boarderlogic,"boarder",this,0,checkOverlaps);

	G4OpticalSurface* pmma_opsurf= new G4OpticalSurface("sipm_opsurf",glisur,ground,dielectric_dielectric);
  	new G4LogicalSkinSurface("pmma_surf",Boarderlogic,pmma_opsurf);
    }

    G4OpticalSurface* OpSurface = new G4OpticalSurface("name");

    G4LogicalBorderSurface* Surface1 = new G4LogicalBorderSurface("name1",volume1,volume2,OpSurface);
    G4LogicalBorderSurface* Surface2 = new G4LogicalBorderSurface("name2",volume2,volume1,OpSurface);

    G4double sigma_alpha = 0.1;
    OpSurface -> SetType(dielectric_dichroic);
    OpSurface -> SetModel(unified);
    OpSurface -> SetFinish(polished);
}
