#ifndef LXeMaterials_H
#define LXeMaterials_H 1

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "LXeDetectorMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"

class LXeMaterials
{
  public:

    LXeMaterials();
    virtual ~LXeMaterials();
    
    G4Material* fLXe;
    G4Material* xenon_mat;
    G4Material* world_mat;
    G4Material* fAl;
    G4Material* fSi;
    G4Element* fN;
    G4Element* fO;
    G4Material* fAir;
    G4Material* fVacuum;
    G4Element* fC;
    G4Element* fCl;
    G4Element* fH;
    G4Element* fSi_e;
    G4Material* fCeramic;
    G4Material* fGlass;
    G4Material* fPMMA;
    G4Material* fQuartz;
    G4Material* fe_mat;
    G4Material* silicon_mat;
    G4Material* quartz_mat;

  private:

    G4MaterialPropertiesTable* fLXe_mt;

};

namespace DMaterials {
extern LXeMaterials* DetectorMaterials;

void BuildMaterials();

G4Material* Get_LXe();
G4Material* Get_xenon_mat();
G4Material* Get_world_mat();
G4Material* Get_fAl();
G4Material* Get_fSi();
G4Element* Get_fN();
G4Element* Get_fO();
G4Material* Get_fAir();
G4Material* Get_fVacuum();
G4Element* Get_fC();
G4Element* Get_fCl();
G4Element* Get_fH();
G4Element* Get_fSi_e();
G4Material* Get_fCeramic();
G4Material* Get_fGlass();
G4Material* Get_fPMMA();
G4Material* Get_fQuartz();
G4Material* Get_fe_mat();
G4Material* Get_silicon_mat();
G4Material* Get_quartz_mat();
}

#endif
