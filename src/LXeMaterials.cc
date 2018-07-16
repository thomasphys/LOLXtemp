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

#include "LXeMaterials.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"

#include "globals.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "TMath.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMaterials::LXeMaterials()
{
    G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA = 3+2*polyPMMA;
  G4int nH_PMMA = 6+2*polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth = 2*polyeth;
  G4int nH_eth = 4*polyeth;
  
  G4NistManager* nist = G4NistManager::Instance();

  //***Elements
  fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  fO = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);
  fSi_e = new G4Element("Si"  , "Si", z=14., a= 28.09*g/mole);
  fCl = new G4Element("Cl", "Cl", z=17., a=35.45*g/mole);

  //***Materials
  //Liquid Xenon
  
  G4double nm2ev = 1240.;
  
  G4Material* fLXe = nist->FindOrBuildMaterial("G4_lXe");
  const G4int NUMENTRIES1 = 850;  //from 150 to 1000
  G4double LXe_PP1[NUMENTRIES1];  //energies
  G4double LXe_SCINT1[NUMENTRIES1]; //emission probability
  G4double LXe_RIND1[NUMENTRIES1]; // Refraction index
  G4double LXe_ABSL1[NUMENTRIES1]; // Absorption length
  G4double LXe_Rayleigh1[NUMENTRIES1]; // Rayleigh scattering length
  // defining properties of LiXe
  for(int iE=0; iE<NUMENTRIES1; iE++){
    G4double WL=1000-(iE); // starting from 1000nm WL, ending with 150nm

    LXe_PP1[iE]=nm2ev/(WL)*eV; // tabulate energy
    LXe_SCINT1[iE] = 1./sqrt(2*3.14)/5.*exp(-pow((WL-178.)/5.,2.)/2.);//approximating with a gaus (178,5)
    //from https://arxiv.org/pdf/1502.04213.pdf  AS Feb 21, 2018
    LXe_RIND1[iE] = sqrt(1.5+ 0.38*WL*WL/(WL*WL-146.9*146.9)+ 0.009*WL*WL/(WL*WL-827*827));
    if(WL>700) LXe_RIND1[iE]=1.365;// avoid the fall of ref index in the formula
    if(WL<147.) LXe_RIND1[iE]= sqrt(1.5+ 0.38*147.*147./(147.*147.-146.9*146.9)+ 0.009*147.*147./(147.*147.-827*827));
    LXe_ABSL1[iE] = 200.*cm;
    LXe_Rayleigh1[iE] = 35.*cm; // need to code proper formula but should not be critical for LoLX
    //G4cout << " WL= " << WL<< " and energy="<<LXe_PP1[iE]<< " iE="<<iE<<G4endl; //debug output
  }
    
  G4MaterialPropertiesTable* LXe_MPT = new G4MaterialPropertiesTable();
  LXe_MPT -> AddProperty("FASTCOMPONENT",LXe_PP1, LXe_SCINT1, NUMENTRIES1);
  LXe_MPT -> AddProperty("SLOWCOMPONENT",LXe_PP1, LXe_SCINT1, NUMENTRIES1);
  LXe_MPT -> AddProperty("RINDEX", LXe_PP1, LXe_RIND1, NUMENTRIES1);
  LXe_MPT -> AddProperty("ABSLENGTH",LXe_PP1, LXe_ABSL1, NUMENTRIES1);
  LXe_MPT -> AddProperty("RAYLEIGH",LXe_PP1, LXe_Rayleigh1, NUMENTRIES1);
  LXe_MPT -> AddConstProperty("RESOLUTIONSCALE", 1.0);
  LXe_MPT -> AddConstProperty ("SCINTILLATIONYIELD",46300/MeV);// 46300 was 68000 was 100 before should be 46296, based on 21.6 eV for Beta or 17.9 eV for alpha
  //LXe_MPT -> AddConstProperty ("SCINTILLATIONYIELD",0.0/MeV);
  LXe_MPT -> AddConstProperty("FASTTIMECONSTANT",2.2*ns); // from http://www.pd.infn.it/~conti/images/LXe/tabellatau.pdf
  LXe_MPT -> AddConstProperty("SLOWTIMECONSTANT",34.*ns);
  LXe_MPT -> AddConstProperty("YIELDRATIO",.5); // fraction of fast to slow components

  fLXe -> SetMaterialPropertiesTable(LXe_MPT);
  //LXe->GetIonisation()->SetBirksConstant(0.126*mm/MeV); //no saturation 

  xenon_mat = fLXe;
  
  //Silicone
  fSi = new G4Material("Si",z=14.,a=28.084*g/mole,density=2.3290*g/cm3);//put in better numbers
  //Aluminum
  fAl = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);
  //Air
  fAir = new G4Material("Air", density= 1.29*mg/cm3, 2);
  fAir->AddElement(fN, 70*perCent);
  fAir->AddElement(fO, 30*perCent);

  world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4MaterialPropertiesTable *world_mt = new G4MaterialPropertiesTable();
  G4double world_rindex[NUMENTRIES1];for(double &r: world_rindex) r=1.;
  world_mt->AddProperty("RINDEX",LXe_PP1,world_rindex,NUMENTRIES1);
  world_mat->SetMaterialPropertiesTable(world_mt);
  
  fe_mat = nist->FindOrBuildMaterial("G4_Fe");
  G4MaterialPropertiesTable *fe_mt = new G4MaterialPropertiesTable();
  G4double fe_rindex[NUMENTRIES1];for(double &r: fe_rindex) r=2.9;
  fe_mt->AddProperty("RINDEX",LXe_PP1,fe_rindex,NUMENTRIES1);
  fe_mat->SetMaterialPropertiesTable(fe_mt);
  
  silicon_mat = nist->FindOrBuildMaterial("G4_Si");
  G4MaterialPropertiesTable *silicon_mt = new G4MaterialPropertiesTable();
  G4double silicon_rindex[NUMENTRIES1];for(double &r: silicon_rindex) r=3.4;
  silicon_mt->AddProperty("RINDEX",LXe_PP1,silicon_rindex,NUMENTRIES1);
  silicon_mat->SetMaterialPropertiesTable(silicon_mt);
  
  quartz_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();
  G4double quartz_rindex[NUMENTRIES1];for(double &r: quartz_rindex) r=1.42;
  quartz_mt->AddProperty("RINDEX",LXe_PP1,quartz_rindex,NUMENTRIES1);
  quartz_mat->SetMaterialPropertiesTable(quartz_mt);

  //Glass
  fGlass = new G4Material("Glass", density=1.032*g/cm3,2);
  fGlass->AddElement(fC,91.533*perCent);
  fGlass->AddElement(fH,8.467*perCent);
  //Fiber(PMMA)
  fPMMA = new G4Material("PMMA", density=1190*kg/m3,3);
  fPMMA->AddElement(fH,nH_PMMA);
  fPMMA->AddElement(fC,nC_PMMA);
  fPMMA->AddElement(fO,2);
  
  fCeramic = new G4Material("Ceramic", density=2.0*g/cm3,2);
  fCeramic->AddElement(fSi_e,50*perCent);
  fCeramic->AddElement(fC,50*perCent);
  
  fQuartz = new G4Material("Quartz", density=2.2*g/cm3,2);
  fQuartz->AddElement(fSi_e,20*perCent);
  fQuartz->AddElement(fCl,80*perCent);
 
  //***Material properties tables
  
  
  G4double glass_AbsLength_wavelength[57]={202.968,208.880,214.793,231.362,246.094,254.982,269.816,284.655,302.454,323.232, 346.970, 361.805, 385.547, 406.317, 435.998, 456.776, 474.575, 501.288, 522.066, 551.740,
                                       575.474, 599.212, 628.897, 667.482, 700.123, 720.890, 759.466, 792.115, 827.732, 860.377, 904.897, 949.413, 996.889,1032.494,1074.039,1115.583,1171.970,1219.458,1278.816,1329.268,
                                      1388.627,1436.110,1501.400,1566.687,1634.937,1694.284,1744.736,1818.925,1890.151,1958.409,2014.796,2091.961,2172.094,2243.323,2299.710,2367.976,2468.890};
  G4double glass_AbsLength_Energy[57];
  for(int i=0; i<57; ++i) glass_AbsLength_Energy[56-i] = nm2ev*eV/glass_AbsLength_wavelength[i];
  

  G4double glass_Trasmittance[57]={0.000006,0.01505,0.030094,0.839635,0.909826,0.919861,0.922384,0.9224,0.927432,0.924948,0.92748,0.930003,0.930029,0.932557,0.930083,0.9276,0.932632,0.930155,
                                       0.927671,0.93021,0.935248,0.93778,0.9328,0.930336,0.932878,0.937913,0.940461,0.93799,0.935523,0.935558,0.933101,0.933149,0.938213,0.943265,0.945816,0.948368,
                                       0.948429,0.945975,0.943533,0.943588,0.941147,0.941198,0.941269,0.943847,0.948934,0.954011,0.954066,0.956653,0.956731,0.956805,0.956866,0.954444,0.952025,
                                       0.949597,0.949658,0.94472,0.937311};
  G4double maxtrans = 0.0;
  for(int i=0; i<57; ++i) if(glass_Trasmittance[i]>maxtrans)maxtrans=glass_Trasmittance[i];
  maxtrans *= 1.01;
  for(int i=0; i<57; ++i) glass_Trasmittance[i] /= maxtrans;

  G4double glass_AbsLength[57];
  for(int i=0; i<57; ++i) glass_AbsLength[56-i] = 1000*mm;
  for(int i=0; i<57; ++i)if(glass_AbsLength_wavelength[i]<260.) glass_AbsLength[56-i] = 0.001*mm;//-0.97*mm/TMath::Log(glass_Trasmittance[i]);
 
  G4double glass_RIND[57];
  for(int i=0; i<57; ++i) glass_RIND[i] = 2.242;
  
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",glass_AbsLength_Energy,glass_AbsLength,57);
  glass_mt->AddProperty("RINDEX",glass_AbsLength_Energy,glass_RIND,57);
  fGlass->SetMaterialPropertiesTable(glass_mt);
  
  G4double PMMA_RIND_Wavelength[62] =  { 60.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,
                                        390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,
                                        590.,600.,610.,620.,630.,640.,650.,660.,670.,680.,690.,700.,710.,720.,730.,740.,750.,760.,770.,780.,
                                        790.,800.};
 G4double PMMA_RIND_Energy[62];
  for(int i=0; i<62; ++i) PMMA_RIND_Energy[61-i] =  nm2ev*eV/PMMA_RIND_Wavelength[i];
  
  G4double PMMA_RIND[62] =  {1.597,1.597,1.584,1.573,1.564,1.556,1.550,1.544,1.539,1.534,1.531,1.527,1.524,1.521,1.519,1.516,1.514,1.512,1.510,1.509,
                             1.507,1.506,1.505,1.503,1.502,1.501,1.500,1.499,1.499,1.498,1.497,1.496,1.496,1.495,1.494, 1.494,1.493,1.493,1.492,1.492,
                             1.491,1.491,1.490,1.490,1.490,1.489,1.489,1.488,1.488,1.488,1.488,1.487,1.487,1.487,1.486,1.486,1.486,1.486,1.485,1.485,
                             1.485,1.485};
    
  G4double PMMA_RINDEnergyOrdered[62];
  for(int i=0; i<62; ++i) PMMA_RINDEnergyOrdered[61-i] = PMMA_RIND[i];

  G4double PMMA_ABS_Wavelength[423] =  { 65.,200.,280.,281.,282.,283.,284.,285.,286.,287.,288.,289.,290.,291.,292.,293.,294.,295.,296.,297.,
                                     298.,299.,300.,301.,302.,303.,304.,305.,306.,307.,308.,309.,310.,311.,312.,313.,314.,315.,316.,317.,
                                     318.,319.,320.,321.,322.,323.,324.,325.,326.,327.,328.,329.,330.,331.,332.,333.,334.,335.,336.,337.,
                                     338.,339.,340.,341.,342.,343.,344.,345.,346.,347.,348.,349.,350.,351.,352.,353.,354.,355.,356.,357.,
                                     358.,359.,360.,361.,362.,363.,364.,365.,366.,367.,368.,369.,370.,371.,372.,373.,374.,375.,376.,377.,
                                     378.,379.,380.,381.,382.,383.,384.,385.,386.,387.,388.,389.,390.,391.,392.,393.,394.,395.,396.,397.,
                                     398.,399.,400.,401.,402.,403.,404.,405.,406.,407.,408.,409.,410.,411.,412.,413.,414.,415.,416.,417.,
                                     418.,419.,420.,421.,422.,423.,424.,425.,426.,427.,428.,429.,430.,431.,432.,433.,434.,435.,436.,437.,
                                     438.,439.,440.,441.,442.,443.,444.,445.,446.,447.,448.,449.,450.,451.,452.,453.,454.,455.,456.,457.,
                                     458.,459.,460.,461.,462.,463.,464.,465.,466.,467.,468.,469.,470.,471.,472.,473.,474.,475.,476.,477.,
                                     478.,479.,480.,481.,482.,483.,484.,485.,486.,487.,488.,489.,490.,491.,492.,493.,494.,495.,496.,497.,
                                     498.,499.,500.,501.,502.,503.,504.,505.,506.,507.,508.,509.,510.,511.,512.,513.,514.,515.,516.,517.,
                                     518.,519.,520.,521.,522.,523.,524.,525.,526.,527.,528.,529.,530.,531.,532.,533.,534.,535.,536.,537.,
                                     538.,539.,540.,541.,542.,543.,544.,545.,546.,547.,548.,549.,550.,551.,552.,553.,554.,555.,556.,557.,
                                     558.,559.,560.,561.,562.,563.,564.,565.,566.,567.,568.,569.,570.,571.,572.,573.,574.,575.,576.,577.,
                                     578.,579.,580.,581.,582.,583.,584.,585.,586.,587.,588.,589.,590.,591.,592.,593.,594.,595.,596.,597.,
                                     598.,599.,600.,601.,602.,603.,604.,605.,606.,607.,608.,609.,610.,611.,612.,613.,614.,615.,616.,617.,
                                     618.,619.,620.,621.,622.,623.,624.,625.,626.,627.,628.,629.,630.,631.,632.,633.,634.,635.,636.,637.,
                                     638.,639.,640.,641.,642.,643.,644.,645.,646.,647.,648.,649.,650.,651.,652.,653.,654.,655.,656.,657.,
                                     658.,659.,660.,661.,662.,663.,664.,665.,666.,667.,668.,669.,670.,671.,672.,673.,674.,675.,676.,677.,
                                     678.,679.,680.,681.,682.,683.,684.,685.,686.,687.,688.,689.,690.,691.,692.,693.,694.,695.,696.,697.,
                                     698.,699.,700.};
  G4double PMMA_ABS_Energy[423];
  for(int i=0; i<423; ++i) PMMA_ABS_Energy[422-i] = nm2ev*eV/PMMA_ABS_Wavelength[i];
    
  G4double PMMA_ABS_Length[423] =  {0.001,000.001,0067.58,0042.92,0037.46,0034.36,0075.65,0014.24,0004.84,0002.19,0001.47,0001.38,0001.38,0001.50,0001.85,0002.18,0002.43,0002.73,0003.12,0003.54,
                                    03.94,0004.44,0005.05,0005.87,0006.81,0007.88,0009.05,0010.53,0012.07,0013.64,0015.29,0017.04,0018.86,0020.56,0022.25,0023.91,0025.55,0027.02,0028.61,0030.28,
                                    32.11,0033.89,0035.91,0038.03,0040.34,0042.81,0045.36,0048.11,0050.93,0053.87,0057.30,0060.15,0063.11,0065.91,0068.55,0071.27,0073.78,0076.08,0078.37,0080.59,
                                    82.79,0084.96,0087.26,0089.71,0092.21,0094.80,0097.55,0100.50,0103.53,0106.60,0110.01,0113.77,0117.46,0121.36,0125.36,0129.90,0134.40,0139.22,0144.55,0149.91,
                                   155.81,0162.14,0168.81,0176.11,0184.06,0192.85,0202.31,0213.17,0225.05,0237.41,0250.85,0266.65,0283.61,0302.01,0322.04,0343.89,0367.01,0391.49,0417.07,0445.36,
                                   475.02,0508.50,0538.37,0572.88,0607.17,0648.79, 694.59,0746.68,0791.99,0835.87,0882.24,0928.49,0979.52,1030.52,1081.12,1127.95,1169.07,1212.70,1257.31,1302.66,
                                  1339.17,1370.10,1407.10,1440.62,1472.00,1495.92,1528.44,1553.69,1580.66,1606.40,1625.76,1658.86,1672.81,1698.49,1719.23,1742.58,1757.11,1774.78,1812.56,1823.23,
                                  1836.98,1851.56,1879.05,1894.44,1914.24,1933.27,1963.48,1969.48,1996.37,2024.76,2033.99,2042.78,2057.53,2086.12,2125.80,2141.81,2146.43,2162.58,2183.73,2202.56,
                                  2221.08,2240.20,2277.44,2290.00,2321.74,2340.21,2373.72,2389.63,2389.13,2404.46,2450.94,2475.93,2488.57,2498.67,2541.45,2580.17,2616.36,2662.25,2670.56,2675.02,
                                  2699.55,2728.09,2790.08,2759.08,2777.12,2824.90,2843.81,2854.17,2871.54,2889.72,2910.04,2931.55,2926.42,2953.43,2989.92,3023.97,3034.70,3035.75,3062.08,3047.45,
                                  3082.25,3104.92,3102.34,3131.03,3175.08,3173.87,3184.51,3207.41,3234.83,3249.92,3274.21,3285.67,3317.77,3353.93,3350.44,3359.06,3407.99,3416.40,3423.90,3426.44,
                                  3432.33,3460.35,3495.39,3497.56,3565.50,3582.85,3567.63,3582.17,3655.17,3674.32,3668.88,3715.38,3729.30,3745.51,3770.15,3793.91,3855.13,3916.19,3928.38,3938.47,
                                  3959.27,4006.08,4054.41,4076.31,4090.77,4109.07,4143.65,4178.79,4209.61,4246.41,4270.02,4287.30,4316.69,4319.79,4311.64,4292.95,4326.17,4351.55,4377.48,4383.87,
                                  4366.35,4354.38,4361.28,4402.13,4407.61,4404.43,4407.98,4476.82,4525.52,4541.79,4596.28,4589.10,4556.81,4619.94,4656.46,4671.57,4765.87,4797.87,4866.13,4921.54,
                                  4933.44,4959.98,5064.61,5150.38,5106.72,5146.09,5170.84,5181.52,5243.49,5295.43,5344.13,5391.89,5396.31,5388.61,5430.96,5499.15,5466.23,5412.68,5430.21,5469.04,
                                  5463.16,5529.56,5527.87,5586.64,5603.18,5527.96,5620.59,5652.23,5644.46,5694.86,5704.38,5727.56,5748.94,5729.26,5695.56,5624.84,5730.34,5697.05,5660.15,5697.30,
                                  5698.42,5650.42,5622.80,5680.14,5579.96,5510.08,5486.90,5404.21,5363.16,5288.98,5223.99,5145.40,5039.10,4898.50,4818.82,4726.98,4621.64,4555.02,4508.72,4429.16,
                                  4339.17,4278.00,4214.02,4201.53,4212.28,4221.29,4261.50,4314.62,4380.30,4436.89,4480.61,4599.42,4687.90,4782.77,4889.66,5084.54,5219.08,5341.56,5482.05,5667.66,
                                  5806.92,5913.54,6174.61,6341.41,6414.37,6512.86,6635.74,6717.40,6762.14,6790.97,6990.80,7129.37,7098.97,7039.25,6998.95,7074.30,7026.00,6956.08,6965.11,6908.33,
                                  6784.48,6698.72,6638.20,6573.63,6520.23,6528.01,6372.79,6330.25,6351.60,6226.50,6145.43,6112.99,6111.87,6097.08,6148.07,6161.22,6228.73,6135.23,6080.73,6216.11,
                                  6100.26,6129.49,6012.49,6065.03,6091.34,6096.55,6035.70,6057.97,6073.05,6075.99,6032.31,6004.06,5990.53,5960.05,5954.05,5886.62,5813.18,5744.58,5572.41,5428.94,
                                  5272.57,5128.36,4966.45};
    
  G4double PMMA_ABS_LengthEnergyOrdered[423];
  for(int i=0; i<423; ++i) PMMA_ABS_LengthEnergyOrdered[422-i] = 0.001*mm;//PMMA_ABS_Length[i]*mm;
    
  G4MaterialPropertiesTable *PMMA_mt = new G4MaterialPropertiesTable();
  PMMA_mt->AddProperty("ABSLENGTH",PMMA_ABS_Energy,PMMA_ABS_LengthEnergyOrdered,423);
  PMMA_mt->AddProperty("RINDEX",PMMA_RIND_Energy,PMMA_RINDEnergyOrdered,62);
  fPMMA->SetMaterialPropertiesTable(PMMA_mt);
  
  G4double Quartz_Wavelength[23] = {142.946,158.016,160.396,161.983,163.966,165.155,171.104,172.691,174.277,175.467,
                                    178.640,183.399,187.762,192.521,196.090,204.419,210.368,217.507,236.940,266.288,
                                    295.240,400.000,991.517};
    
  G4double Quartz_Transmitance[23] = {0.006,0.006,0.018,0.034,0.058,0.095,0.718,0.760,0.788,0.820,
                                      0.841,0.855,0.865,0.893,0.900,0.902,0.906,0.916,0.918,0.923,
                                      0.930,0.939,0.955};
    
  G4double Quartz_Energy[23];
  G4double Quartz_Index[23];
  G4double Quartz_Abs[23];
  
  for(int i=0; i<23; ++i){
      Quartz_Energy[22-i] = nm2ev*eV/Quartz_Wavelength[i];
      Quartz_Index[i] = 1.46;
      Quartz_Abs[22-i] = -10./TMath::Log(Quartz_Transmitance[i]);
  }
    
  G4MaterialPropertiesTable *Quartz_mt = new G4MaterialPropertiesTable();
  //Quartz_mt->AddProperty("ABSLENGTH",Quartz_Energy,Quartz_Abs,13);
  Quartz_mt->AddProperty("RINDEX",Quartz_Energy,Quartz_Index,13);
  fQuartz->SetMaterialPropertiesTable(Quartz_mt);
  
  G4double Si_Energy[3] = {0.1*eV,1.0*eV,10.0*eV};
  G4double Si_Abs[3] = {0.1*mm,0.1*mm,0.1*mm};
  
  G4MaterialPropertiesTable *Si_mt = new G4MaterialPropertiesTable();
  Si_mt->AddProperty("ABSLENGTH",Si_Energy,Si_Abs,3);
  fSi->SetMaterialPropertiesTable(Si_mt);
  
  G4double Ceramic_Energy[3] = {0.1*eV,1.0*eV,10.0*eV};
  G4double Ceramic_Abs[3] = {0.1*mm,0.1*mm,0.1*mm};
  
  G4MaterialPropertiesTable *ceramic_mt = new G4MaterialPropertiesTable();
  ceramic_mt->AddProperty("ABSLENGTH",Ceramic_Energy,Ceramic_Abs,3);
  fCeramic->SetMaterialPropertiesTable(ceramic_mt);

  G4double vacuum_Energy[3]={2.0*eV,7.0*eV,7.14*eV};
  G4double vacuum_RIND[3]={1.,1.,1.};
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND,3);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);//Give air the same rindex
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMaterials::~LXeMaterials() {}

namespace DMaterials {

LXeMaterials* DetectorMaterials = new LXeMaterials();

void BuildMaterials(){
//return;    if(!DetectorMaterials) DetectorMaterials = new LXeMaterials();
}

G4Material* Get_LXe(){BuildMaterials(); return DetectorMaterials->fLXe;}
G4Material* Get_xenon_mat(){BuildMaterials(); return DetectorMaterials->xenon_mat;}
G4Material* Get_world_mat(){BuildMaterials(); return DetectorMaterials->world_mat;}
G4Material* Get_fAl(){BuildMaterials(); return DetectorMaterials->fAl;}
G4Material* Get_fSi(){BuildMaterials(); return DetectorMaterials->fSi;}
G4Element* Get_fN(){BuildMaterials(); return DetectorMaterials->fN;}
G4Element* Get_fO(){BuildMaterials(); return DetectorMaterials->fO;}
G4Material* Get_fAir(){BuildMaterials(); return DetectorMaterials->fAir;}
G4Material* Get_fVacuum(){BuildMaterials(); return DetectorMaterials->fAir;}
G4Element* Get_fC(){BuildMaterials(); return DetectorMaterials->fC;}
G4Element* Get_fCl(){BuildMaterials(); return DetectorMaterials->fCl;}
G4Element* Get_fH(){BuildMaterials(); return DetectorMaterials->fH;}
G4Element* Get_fSi_e(){BuildMaterials(); return DetectorMaterials->fSi_e;}
G4Material* Get_fCeramic(){BuildMaterials(); return DetectorMaterials->fCeramic;}
G4Material* Get_fGlass(){BuildMaterials(); return DetectorMaterials->fGlass;}
G4Material* Get_fPMMA(){BuildMaterials(); return DetectorMaterials->fPMMA;}
G4Material* Get_fQuartz(){BuildMaterials(); return DetectorMaterials->fQuartz;}
G4Material* Get_fe_mat(){BuildMaterials(); return DetectorMaterials->fe_mat;}
G4Material* Get_silicon_mat(){BuildMaterials(); return DetectorMaterials->silicon_mat;}
G4Material* Get_quartz_mat(){BuildMaterials(); return DetectorMaterials->quartz_mat;}
}
