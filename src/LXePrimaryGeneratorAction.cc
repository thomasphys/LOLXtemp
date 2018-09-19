#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4OpticalPhoton.hh"

#include "Randomize.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "LXePrimaryGeneratorAction.hh"
#include "LXePrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector randomDirection() {
    double phi = 2.*M_PI*G4UniformRand();
    double costheta = 2.*G4UniformRand()-1.;
    double sintheta = sqrt(1.-costheta*costheta);
    return G4ThreeVector(cos(phi)*sintheta, sin(phi)*sintheta, costheta);
}


LXePrimaryGeneratorAction::LXePrimaryGeneratorAction()
{
  fMessenger   = new LXePrimaryGeneratorMessenger(this);

  fQ_value = 2457.8*keV;
  fBb2nCutOffMinFraction = 0.;
  fBb2nCutOffMaxFraction = 1.;

  fFF_factor = 0;
  fK_spectral_max = fD_spectral_max = fNormalization = 0;

  fGenerator = G4String(" ");

  fXeComponent = G4String(" ");

  fParticleGun1 = new G4ParticleGun();
  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");

  fParticleGun1->SetParticleDefinition(particle);
  fParticleGun1->SetParticleTime(0.0*ns);
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(0,0,0));
  fParticleGun1->SetParticleEnergy(6.966292135*eV);
    
    fParticleGun3 = new G4ParticleGun(1);
    G4ParticleDefinition* particle2 = particleTable->FindParticle("e+");
    
    fParticleGun3->SetParticleDefinition(particle2);
    fParticleGun3->SetParticleTime(0.0*ns);
    fParticleGun3->SetParticleMomentumDirection(G4ThreeVector(0,0,0));
    fParticleGun3->SetParticleEnergy(0*MeV);

  fParticleGun2 = new G4GeneralParticleSource();

  // Solar Boron-8 Neutrino Electron Recoil
  fB8NeutrinoMaxEnergy = 16.56; // MeV
  fB8NeutrinoMaxY = 0.133; // arb.
  fB8NeutrinoBinWidth = 0.02; // MeV
  electron_mass_c2 = 0.511;
    
    theNavigator = NULL;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction()
{
  delete fMessenger;
  delete fParticleGun1;
  delete fParticleGun2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if ( fGenerator == "gps" || fGenerator == "GeneralParticleSource" ){

     fParticleGun2->GeneratePrimaryVertex(anEvent);
     for(int i=0; i<anEvent->GetNumberOfPrimaryVertex();i++){
         auto vertex = anEvent->GetPrimaryVertex(i);
         for(int j=0; j<vertex->GetNumberOfParticle(); j++){
             auto particle = vertex->GetPrimary(j);
             if(particle->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()
                     && particle->GetPolarization().mag2()==0){
                 G4ThreeVector polar = particle->GetMomentumDirection().cross(randomDirection());
                 particle->SetPolarization(polar);
             }
         }
     }
  }
  else if ( fGenerator == "gun" ) fParticleGun3->GeneratePrimaryVertex(anEvent);
  else if ( fGenerator == "bb0n" ) Generate_bb0n(anEvent);
  else if ( fGenerator == "single" ) Generate_singleelec(anEvent);
  else if ( fGenerator == "bb2n" ) Generate_bb2n(anEvent);
  else if ( fGenerator == "nCaptureXe136" )
              Generate_nCaptureXe136(anEvent);
  else {
    G4double cosTheta = 2*G4UniformRand() - 1., phi = 2*3.14159265358979323846*G4UniformRand();
    G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    G4double ux = cosTheta,
             uy = sinTheta*std::sin(phi),
             uz = sinTheta*std::cos(phi);

    fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));

    phi = 2*3.14159265358979323846*G4UniformRand();
    G4double r = std::sqrt(G4UniformRand())*650*mm;
    G4double z = r*std::cos(phi),
             y = r*std::sin(phi),
             x = (G4UniformRand()-0.5)*2*650*mm;

    fParticleGun1->SetParticlePosition(G4ThreeVector(x,y,z));

    fParticleGun1->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Set_generator(G4String generator)
{
  if ( generator == "gps"                      ||
       generator == "gun"                      ||
       generator == "GeneralParticleSource"    ||
       generator == "gun"                      ||
       generator == "bb0n"                     ||
       generator == "single"                   ||
       generator == "bb2n"                     ||
       generator == "nCaptureXe136")
  {
      G4cout << "Using " << generator << " Generator" << G4endl;
      fGenerator = generator;
      Set_norm();

  } else {
      G4cout << "Generator "<< generator <<" not recognized" << G4endl;
      fGenerator = G4String("gun");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Set_nCaptureXe_Method(G4String method)
{
  G4cout << "***************************************************" << G4endl;

  if (method != "InternalConversions" &&
      method != "RandomGammas" &&
      method != "ImbalancedCascade")
  {
    G4cout << "Error: must be InternalConversions, RandomGammas, or ImbalancedCascade" << G4endl;
  } else {
    G4cout << "setting the method to " << method << G4endl;
    fnCaptureXe_Method = method;
  }

  G4cout << "***************************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Set_Xe_Component(G4String component)
{
  G4cout << "***************************************************" << G4endl;

  if (component != "ActiveLXe" && component != "InactiveLXe") {
    G4cout << "Error: must be ActiveLXe or InactiveLXe" << G4endl;
  } else {
    G4cout << "setting the component to " << component << G4endl;
    fXeComponent = component;
  }

  G4cout << "***************************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Set_bb2n_CutOffMin(G4double frac)
{
  fBb2nCutOffMinFraction = frac;
  G4cout << "Setting the fraction of bb2n min cut off to " << fBb2nCutOffMinFraction << " = " << fBb2nCutOffMinFraction*fQ_value << " keV\n";
  Set_norm();
}

void LXePrimaryGeneratorAction::Set_bb2n_CutOffMax(G4double frac)
{
  fBb2nCutOffMaxFraction = frac;
  G4cout << "Setting the fraction of bb2n max cut off to " << fBb2nCutOffMaxFraction << " = " << fBb2nCutOffMaxFraction*fQ_value << " keV\n";
  Set_norm();
}

void LXePrimaryGeneratorAction::Set_GenInCenter(G4bool b)
{
    fGenInCenter = b;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::GetIsotropicDirection(G4ThreeVector& pos){
  G4double cost = 2*G4UniformRand() - 1.;
  G4double sint = sqrt(1. - cost*cost);
  G4double phi = 2*pi*G4UniformRand();
  G4double vdx = sint*cos(phi);
  G4double vdy = sint*sin(phi);
  G4double vdz = cost;
  pos.set(vdx,vdy,vdz);
}

void LXePrimaryGeneratorAction::GetUnifRandPosInLXe(G4ThreeVector& pos)
{
  // Return a random position in LXe volume following a uniform distribution through this volume
  // Assumes LXe volume is contained in TPC vessel which is a conical (G4Cons) along z
    
    if(!theNavigator){
        theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    }

   if(fGenInCenter){
       pos = G4ThreeVector(0.,0.,0.);
       return;
  }
    
    G4VPhysicalVolume* volume = theNavigator->LocateGlobalPointAndSetup(pos);
    
    G4String LXevolume = volume->GetName();

  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    
    G4bool isoutsidevolume = true;
    
  while ( isoutsidevolume )
  {
      if(fGeometry == 0){
          G4double z0 = ((2.0*G4UniformRand() - 1.0)*38*mm/2.);
          G4double phi0 = 2*pi*G4UniformRand();
          G4double r0 = 19*mm*sqrt(G4UniformRand());
          pos.set(r0*cos(phi0),r0*sin(phi0),z0);
      }else{
          G4double costheta = 2.0*G4UniformRand() - 1.0;
          G4double sintheta = TMath::Sqrt(1.-costheta*costheta);
          G4double phi0 = 2*pi*G4UniformRand();
          G4double r0 = 28.355044*mm*TMath::Power(G4UniformRand(),1./3.);
          pos.set(r0*cos(phi0)*sintheta,r0*sin(phi0)*sintheta,r0*costheta);
      }

      volume = theNavigator->LocateGlobalPointAndSetup(pos);

      if(volume->GetName() == LXevolume) isoutsidevolume = false;
  }

  return;
}

void LXePrimaryGeneratorAction::Generate_bb0n(G4Event* anEvent)
{
  //
  // K=T1+T2 (sum of kinetic energies T1 & T2).
  // D=T1-T2 (their difference; -K < D < +K; D=0 => T1=T2=T0/2)
  //
  // For bb0n events, K is the Q value divided by electron mass.
  // The sum electron spectrum of the bb0n mode is very simple; it is just a
  // delta-function peak at the end point energy E0 (or T0).
  //
  // Algorithm to generate the energy of the two electrons for the bb0n decay:
  // 1. Randomly generate d that satisfies the phase space constraint
  //
  //    P1*E1*FermiFunction(Z,E1)*P2*E2*FermiFunction(Z,E2)
  //
  // 2. With this set of K and D, calculate T1 and T2
  //    T1 = (K+D)/2
  //    t2 = (K-D)/2
  //
  // To generate an arbitrary probabilistic distribution function f(x),
  // use the Acceptance-Rejection Monte Carlo Method.
  // (cf. Frank Porter's lecture on the Acceptance-Rejection Method.)
  //
  // http://www.hep.caltech.edu/~fcp/statistics/sluo00/sluolec5.pdf
  //
  // Max-------------------
  //    |                 |
  //    |    *** f(x)     |
  //    |   *   **        |
  //    |  *      **      |
  //    | *         **    |
  //    |*            **  |
  //  0 *---------------***
  //    a                 b
  //
  // Generate randomly a point within the box bounded by (a->b, 0->Max).
  // If the point falls under f(x), accept it, otherwise, reject.
  // If it is rejected, repeat the process until a point falls under f(x).
  //

  G4int max_iterations = 100; //The max no. of iterations to satisfy
                              //the Acceptance-Rejection condition

  fParticleGun1->SetParticleDefinition(G4Electron::Electron());

  // Choose electron energies

  G4double t0 = fQ_value/electron_mass_c2;

  G4double k = t0;
  G4double t1 = 0.;
  G4double t2 = 0.;

  // Find the peak of the D spectrum first

  fD_spectral_max = D_bb0n_spectral_max(k);

  // Acceptance-Rejection iteration loop

  G4int n = 0;
  while ( n < max_iterations ) {
    G4double d = k*(2*G4UniformRand()-1.0);   // -k < d < +k
    G4double d_spectral = fD_spectral_max * G4UniformRand();

    // Test whether or not a randomly generated point is
    // below the D spectrum
    if ( d_spectral <= D_bb0n_spectrum(k,d) ) {
       t1 = 0.5*(k+d);
       t2 = 0.5*(k-d);
       break; // Condition is met; stop the loop.
    }
    n++;
  }

  t1 = t1*electron_mass_c2; // scale by the electron mass
  t2 = t2*electron_mass_c2;

  //*****Want to generate a primary vertex within the LXe boundary,
  //*****test whether or not it is actually in the xenon.
  //***** Remember that the axis of the TPC lies along the Z coordinate!
  
  G4ThreeVector decayVertex;
  GetUnifRandPosInLXe(decayVertex);
  //  G4cout << "Final point = " << decayVertex << G4endl;
  fParticleGun1->SetParticlePosition(decayVertex);

  // GENERATE FIRST ELECTRON ------------------------------------

  G4double theta = acos(2*(G4UniformRand()-0.5));
  G4double phi   = 2*pi*G4UniformRand();
  G4double ux1   = sin(theta)*cos(phi);
  G4double uy1   = sin(theta)*sin(phi);
  G4double uz1   = cos(theta);

  fParticleGun1->SetParticleEnergy(t1);
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux1,uy1,uz1));
  fParticleGun1->GeneratePrimaryVertex(anEvent);

  // GENERATE SECOND ELECTRON ----------------------------------

  // For the bb0n decay, the angular distribution of the electrons
  // is of the form
  // f(cos)= 1 - beta1*beta2*cos(theta12) with beta = P/E,
  // for the 0+ --> 0+ transition (cf. Boehm & Vogel p.146).
  //
  // This distribution favors a negative cos(theta12), i.e., when the
  // electrons are pointing at opposite hemispheres.
  //
  // The range of f(cos) is:
  // 1 - beta1*beta2 < f(cos) < 1 + beta1*beta2
  //
  // When the direction of the first electron is chosen, the direction
  // of the second electron has to be generated according to the above
  // angular distribution.
  //

  G4double e1 = t1 + electron_mass_c2;
  G4double e2 = t2 + electron_mass_c2;
  G4double p1 = sqrt(e1*e1-electron_mass_c2*electron_mass_c2);
  G4double p2 = sqrt(e2*e2-electron_mass_c2*electron_mass_c2);
  G4double beta1 = p1/e1;
  G4double beta2 = p2/e2;

  G4double theta12 = 0., phi2 = 0.;

  n = 0;
  while ( n < max_iterations) {
    theta12 = acos(2*(G4UniformRand()-0.5));
    G4double angular_distribution = 1.0 + beta1*beta2*(2*G4UniformRand()-1.0);
    // 1 - beta1*beta2 < f(cos) < 1 + beta1*beta2
    if ( angular_distribution <= ( 1.0 - beta1*beta2*cos(theta12) ) ) {
       // Acceptance-Rejection condition for the second electron to satisfy
       // the angular distribution

       phi2 = 2*pi*G4UniformRand();
       break;  // Condition is met; stop the loop.
    }
    n++;
  }

  G4double ux2 = sin(theta12)*cos(phi2);
  G4double uy2 = sin(theta12)*sin(phi2);
  G4double uz2 = cos(theta12);

  // Rotate second electron to detector coordinates

  G4double ct = cos(theta);
  G4double st = sin(theta);
  G4double cp = cos(phi);
  G4double sp = sin(phi);

  G4double uxp = ct*ux2 + st*uz2;
  G4double uyp = uy2;
  G4double uzp = -1.0*st*ux2 + ct*uz2;
  ux2 = uxp;
  uy2 = uyp;
  uz2 = uzp;

  uxp = cp*ux2 - sp*uy2;
  uyp = sp*ux2 + cp*uy2;
  uzp = uz2;
  ux2 = uxp;
  uy2 = uyp;
  uz2 = uzp;

  fParticleGun1->SetParticleEnergy(t2);
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux2,uy2,uz2));
  fParticleGun1->GeneratePrimaryVertex(anEvent);
}

void LXePrimaryGeneratorAction::Generate_singleelec(G4Event* anEvent)
{
    
    fParticleGun1->SetParticleDefinition(G4Electron::Electron());
    G4ThreeVector decayVertex;
    GetUnifRandPosInLXe(decayVertex);
    fParticleGun1->SetParticlePosition(decayVertex);
    
    G4double theta = acos(2*(G4UniformRand()-0.5));
    G4double phi   = 2*pi*G4UniformRand();
    G4double ux1   = sin(theta)*cos(phi);
    G4double uy1   = sin(theta)*sin(phi);
    G4double uz1   = cos(theta);
    
    fParticleGun1->SetParticleEnergy(fQ_value);
    fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux1,uy1,uz1));
    fParticleGun1->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::D_bb0n_spectral_max(G4double k)
{
  // This function calculates the maximum value of the bb0n D spectrum.
  // fFF_factor = 1 for xe136

  fFF_factor = 1;

  //
  // Energies (E1, E2),  momenta (P1, P2), and kinetic energies (T1, T2)
  // are all divided by the electron mass.
  //
  // K = T0 = fQ_value/electron_mass_c2
  // D = T1-T2 (their difference; Range:  -K < D < +K )
  //
  // T1 = (K+D)/2
  // T2 = (k-D)/2
  //

  G4int z = 54;

  G4int nbins = 100;    // nbins is the # of bins.

  G4double a = -k;
  G4double b =  k;
  G4double dD = (b-a)/nbins;

  G4double phase_space[101] = {0.};

  for (G4int n = 0; n < nbins+1; n++) {

      G4double d;

      if (n != nbins) {
        d = a + n*dD;
      } else {
        d = b;
      }

      G4double e1 = 0.5*(k+d)+1.0;
      G4double e2 = 0.5*(k-d)+1.0;

      G4double p1 = sqrt(e1*e1-1.0);
      G4double p2 = sqrt(e2*e2-1.0);

      G4double t1 = e1 - 1.0;
      G4double t2 = e2 - 1.0;

      phase_space[n] = p1*e1*Fermi_function(z,t1*electron_mass_c2) *
                       p2*e2*Fermi_function(z,t2*electron_mass_c2);
  }

  //search for the peak value of phase space[i]

  G4double d_spectral_max = 0.;

  for (G4int i=0; i<nbins+1; i++) {
      if (phase_space[i] > d_spectral_max) d_spectral_max = phase_space[i];
  }

  return d_spectral_max;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::D_bb0n_spectrum(G4double k, G4double d)
{
    G4int z = 54;

    G4double e1 = 0.5*(k+d)+1.0;
    G4double e2 = 0.5*(k-d)+1.0;

    G4double p1 = sqrt(e1*e1-1.0);
    G4double p2 = sqrt(e2*e2-1.0);

    G4double t1 = e1-1.0;
    G4double t2 = e2-1.0;

    return p1*e1*Fermi_function(z,t1*electron_mass_c2) *
           p2*e2*Fermi_function(z,t2*electron_mass_c2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::Fermi_function(G4int z, G4double ke)
{
  // From Schenter+Vogel in Nucl.Sci.Eng,83,393(1983)
  // W is full energy in electron mass units

  G4double alpha = 7.2974e-3; // fine structure constant

  G4double total_energy = ke + electron_mass_c2;
  G4double w = total_energy/electron_mass_c2;

  G4double z0 = z + fFF_factor*2;

  if ( w <= 1 ) w = 1 + 1e-4;

  G4double a = -0.811+4.46e-2*z0+1.08e-4*z0*z0;
  G4double b = 0.673-1.82e-2*z0+6.38e-5*z0*z0;

  if ( w > 2.2 ) {
     a = -8.46e-2+2.48e-2*z0+2.37e-4*z0*z0;
     b = 1.15e-2+3.58e-4*z0-6.17e-5*z0*z0;
  }

  G4double x = sqrt(w-1);
  G4double p = sqrt(w*w-1);

  G4double result = exp(a+b*x)*w/p;

  if (p<=0) result = 1; //just to be consistent with the old Fermi Function code

  if (fFF_factor == -1) { // for double positron decays
     G4double v = p/w;
     G4double y = 2*pi*z0*alpha/v;
     G4double yy = 1./exp(y);
     result = result*yy;
  }

  return result;
}

void LXePrimaryGeneratorAction::Generate_bb2n(G4Event* anEvent)
{
  //
  // K=T1+T2 (sum of kinetic energies T1 and T2; 0 < K < T0)
  // D=T1-T2 (their difference; -K < D < +K; D=0 => T1=T2=T0/2)
  //
  // Algorithm to generate the energy of the two electrons for the bb2n decay:
  //
  //    P1*E1*FermiFunction(Z,E1)*P2*E2*FermiFunction(Z,E2)*pow(T0-K,5)
  //
  // 3. With this set of K and D, calculate T1 and T2
  //    T1 = (K+D)/2
  //    T2 = (K-D)/2
  //
  // This set of K & D (or T1, T2) guarantees to satisfy both
  // the sum electron spectrum dn/dk and the single electron spectra dn/dt0.
  // (cf. Physics of Massive Neutrinosy Felix Boehm & Petr Vogel,pp.145-151.)
  //
  // To generate an arbitrary probabilistic distribution function f(x),
  // use the Acceptance-Rejection Monte Carlo Method.
  // (cf. Frank Porter's lecture on the Acceptance-Rejection Method.)
  //
  // http://www.hep.caltech.edu/~fcp/statistics/sluo00/sluolec5.pdf
  //
  // Max-------------------
  //    |                 |
  //    |    *** f(x)     |
  //    |   *   **        |
  //    |  *      **      |
  //    | *         **    |
  //    |*            **  |
  //  0 *---------------***
  //    a                 b
  //
  // Generate randomly a point within the box bounded by (a->b, 0->Max).
  // If the point falls under f(x), accept it, otherwise, reject.
  // If it is rejected, repeat the process until a point falls under f(x).
  //

  G4int max_iterations = 100; //The max no. of iterations to satisfy
                              //the Acceptance-Rejection condition

  fParticleGun1->SetParticleDefinition(G4Electron::Electron());

  // Choose electron energies

  G4double t0 = fQ_value/electron_mass_c2;

  G4double t1 = 0.;
  G4double t2 = 0.;

  //G4cout << "K_spectral_max = " << fK_spectral_max << G4endl;

  // G4double fK_spectral_max=0.0233;
  // k_spectral_max is now calculated once at initialization.
  // Value of k_spectral_max assumes that the sum electron spectrum dn/dk
  // is normalized to 1, assuming the # of bins=100 (if the # of bins
  // is changed, change also this spectral max value.) With this normalization
  // the max value is slightly less than 0.0233.

  G4double range_width = (fBb2nCutOffMaxFraction - fBb2nCutOffMinFraction);

  G4int i = 0;
  while ( i < max_iterations ) {

    G4double k = t0*(range_width*G4UniformRand() + fBb2nCutOffMinFraction);

    G4double k_spectrum = fK_spectral_max*G4UniformRand();

    if ( k_spectrum <= BB2n_sum_spectrum(k)) {
       // Acceptance-Rejection condition for the sum electron spectrum dn/dk

       G4int n = 0;
       while ( n < max_iterations ) {

         G4double d = k*(2*G4UniformRand()-1.0);   // -k < d < +k

         G4double d_spectral = fD_spectral_max * G4UniformRand();

         //Acceptance-Rejection condition for the single electron spectrum dn/t0

         if ( d_spectral <= D_spectrum(k,d)) {
            t1 = 0.5*(k+d);
            t2 = 0.5*(k-d);
            break; // Condition is met; stop the loop.

         }
         n++;
       }
       break;
    }
    i++;
  }

  t1 = t1*electron_mass_c2; // scale by the electron mass
  t2 = t2*electron_mass_c2;

  // Beta source vertex----------------------------------------

  //*****Want to generate a primary vertex within the LXe boundary,
  //*****test whether or not it is actually in the xenon.
  //***** Remember that the axis of the TPC lies along the Z coordinate!

  G4ThreeVector decayVertex;
  GetUnifRandPosInLXe(decayVertex);
  //  G4cout << "Final point = " << decayVertex << G4endl;
  fParticleGun1->SetParticlePosition(decayVertex);

  // GENERATE FIRST ELECTRON ------------------------------------

  G4double theta = acos(2*(G4UniformRand()-0.5));
  G4double phi   = 2*pi*G4UniformRand();
  G4double ux1   = sin(theta)*cos(phi);
  G4double uy1   = sin(theta)*sin(phi);
  G4double uz1   = cos(theta);

  fParticleGun1->SetParticleEnergy(t1);
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux1,uy1,uz1));
  fParticleGun1->GeneratePrimaryVertex(anEvent);

  // GENERATE SECOND ELECTRON ----------------------------------

  // For the bb2n decay, the angular distribution of the electrons
  // is of the form
  // f(cos)= 1 - beta1*beta2*cos(theta12) with beta = P/E,
  // for the 0+ --> 0+ transition (cf. Boehm & Vogel p.146).
  //
  // This distribution favors a negative cos(theta12), i.e., when the
  // electrons are pointing at opposite hemispheres.
  //
  // The range of f(cos) is:
  // 1 - beta1*beta2 < f(cos) < 1 + beta1*beta2
  //
  // When the direction of the first electron is chosen, the direction
  // of the second electron has to be generated according to the above
  // angular distribution.
  //

  G4double e1 = t1 + electron_mass_c2;
  G4double e2 = t2 + electron_mass_c2;
  G4double p1 = sqrt(e1*e1-electron_mass_c2*electron_mass_c2);
  G4double p2 = sqrt(e2*e2-electron_mass_c2*electron_mass_c2);
  G4double beta1 = p1/e1;
  G4double beta2 = p2/e2;
  G4double theta12 = 0., phi2 = 0.;

  G4int n = 0;
  while ( n < max_iterations) {
    theta12 = acos(2*(G4UniformRand()-0.5));
    G4double angular_distribution = 1.0 + beta1*beta2*(2*G4UniformRand()-1.0);
    // 1 - beta1*beta2 < f(cos) < 1 + beta1*beta2
    if ( angular_distribution <= ( 1.0 - beta1*beta2*cos(theta12) ) ) {
       // Acceptance-Rejection condition for the second electron to satisfy
       // the angular distribution

      phi2 = 2*pi*G4UniformRand();
      break;   // Condition is met; stop the loop.
    }
    n++;
  }

  G4double ux2 = sin(theta12)*cos(phi2);
  G4double uy2 = sin(theta12)*sin(phi2);
  G4double uz2 = cos(theta12);

  // Rotate second electron to detector coordinates

  G4double ct = cos(theta);
  G4double st = sin(theta);
  G4double cp = cos(phi);
  G4double sp = sin(phi);

  G4double uxp = ct*ux2 + st*uz2;
  G4double uyp = uy2;
  G4double uzp = -1.0*st*ux2 + ct*uz2;
  ux2 = uxp;
  uy2 = uyp;
  uz2 = uzp;

  uxp = cp*ux2 - sp*uy2;
  uyp = sp*ux2 + cp*uy2;
  uzp = uz2;
  ux2 = uxp;
  uy2 = uyp;
  uz2 = uzp;

  fParticleGun1->SetParticleEnergy(t2);
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(ux2,uy2,uz2));
  fParticleGun1->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::BB2n_sum_spectrum(G4double k)
{
// This function calculates the sum electron spectrum dn/dk for the bb2n decay.

  //
  // Energies (E1, E2),  momenta (P1, P2), and kinetic energies (T1, T2)
  // are all divided by the electron mass.
  //
  // K = T1+T2 (Sum of kinetic energies; Range: 0 < K < T0 )
  // D = T1-T2 (their difference; Range:  -K < D < +K )
  //
  // T1 = (K+D)/2
  // t2 = (K-D)/2
  //

  G4int z = 54;

  G4int nbins = 100;    // nbins  is the # of bins.

  G4double t0 = fQ_value/electron_mass_c2;

  G4double a = -k;
  G4double b =  k;
  G4double dD = (b-a)/nbins;

  G4double phase_space[101]={0.};

  for (G4int n = 0; n < nbins+1; n++) {

      G4double d;

      if (n != nbins) {
        d = a + n*dD;
      } else {
        d = b;
      }

      G4double e1 = 0.5*(k+d)+1.0;
      G4double e2 = 0.5*(k-d)+1.0;

      G4double p1 = sqrt(e1*e1-1.0);
      G4double p2 = sqrt(e2*e2-1.0);

      G4double t1 = e1 - 1.0;
      G4double t2 = e2 - 1.0;

      phase_space[n] = p1*e1*Fermi_function(z,t1*electron_mass_c2) *
                       p2*e2*Fermi_function(z,t2*electron_mass_c2)*pow(t0-k,5);
  }

    //
    // d_spectral_max is the peak value for the phase_space[n] spectrum.
    // This value is needed in generating the energy of the second electron.
    //

    G4double d_spectral_max = 0.;

    //search for the peak value of phase space[i]
    for (G4int i=0; i<nbins+1; i++) {
        if (phase_space[i] > d_spectral_max) d_spectral_max = phase_space[i];
    }

    // G4double normalization = 134958.7414*nbins;
    // Normalization is now calculated once at initialization.
    // This factor normalizes the integral dn/dk (sum electron spectrum) to 1.
    // nbins = # of bins. If nbins is changed, normalization should also be
    // adjusted.

    G4double normalized_spectrum = SimpsonsRule(a, b, nbins, phase_space) /
                                   fNormalization;

    fD_spectral_max = d_spectral_max;

    return normalized_spectrum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::D_spectrum(G4double k, G4double d)

// This is for the bb2n spectrum
{
    G4int z = 54;

    G4double t0 = fQ_value/electron_mass_c2;

    G4double e1 = 0.5*(k+d)+1.0;
    G4double e2 = 0.5*(k-d)+1.0;

    G4double p1 = sqrt(e1*e1-1.0);
    G4double p2 = sqrt(e2*e2-1.0);

    G4double t1 = e1 - 1.0;
    G4double t2 = e2 - 1.0;

    return p1*e1*Fermi_function(z,t1*electron_mass_c2) *
           p2*e2*Fermi_function(z,t2*electron_mass_c2)*pow(t0-k,5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LXePrimaryGeneratorAction::SimpsonsRule(G4double x0, G4double xn,
                                                  G4int n, G4double f[])
{
  //
  // Simpson's Rule:
  // Partition [a, b] into intervals all of the same width.
  // We must use an even number of intervals, so n will be even.
  // xk = a + kx = a + k (b-a)/n
  //
  // integral(x0, xn)
  // = (xn-x0)/3n *[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) +
  //                         ... + 2f(xn-2) + 4f(xn-1) + f(xn)]
  //
  // Sum all the odd terms and then multiply by 4.
  // Sum all the even terms and then multiply by 2.
  //

    if ( n%2 != 0 ) G4cout << "SimpsonsRule: N is not even";

    G4double sum_odd = 0.;
    for (G4int i = 0; i < n/2; i++) sum_odd += f[2*i+1];
    sum_odd = 4.0*sum_odd;

    G4double sum_even = 0.;
    for (G4int j = 1; j < n/2; j++) sum_even += f[2*j];
    sum_even = 2.0*sum_even;

    return (xn-x0)*(f[0]+sum_odd+sum_even+f[n])/(3.0*n);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Set_norm()
{
  // Calculate the normalization factor (normalization) as well
  // as the maximum value of the k distribution (k_spectral_max) which
  // varies for different isotopes.
  // This function should only run once in the beginning.

  // This function should only run once in the beginning.

   G4double t0 = fQ_value/electron_mass_c2;
   G4double ts = fBb2nCutOffMinFraction*t0;

   fNormalization = 1; //to make BB2n_sum_spectrum run the first time

   G4int nbins = 100;//00;

   G4double bb2n_sum_max = 0;
   G4double* bb2n_sum_spec = new G4double[nbins+1];// = {0.};

   G4cout << "Evaluation of bb2n spectrum...\n";
   for (G4int n = 0; n <= nbins; n++) {
     G4double sum_ene_me = ts + (t0 - ts)*n/nbins;
       bb2n_sum_spec[n] = BB2n_sum_spectrum(sum_ene_me);
       if (bb2n_sum_spec[n] > bb2n_sum_max)
         bb2n_sum_max = bb2n_sum_spec[n];
       //std::cout << electron_mass_c2*sum_ene_me << " " << bb2n_sum_spec[n] << std::endl;
   }
   G4cout << "Max: " << bb2n_sum_max << G4endl;

   fNormalization = SimpsonsRule(0, t0, nbins, bb2n_sum_spec);

   G4cout << "normalization " << fNormalization << G4endl;

   fK_spectral_max = bb2n_sum_max*1.01/fNormalization;

   G4cout << "K_spectral_max " << fK_spectral_max << G4endl;

   delete bb2n_sum_spec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::Generate_nCaptureXe136(G4Event* anEvent)
{
  struct nCapture_Xe136Info
  {
    G4int jump;
    G4double level, gamma;
    G4double ratio;
  };

  fParticleGun1->SetParticleDefinition(G4Gamma::Gamma());

  //******** Set Vertex of Capture + Excited State Decay *************//

  /*
  G4double inactive_Xenon_R = 65.*cm;
  G4double inactive_Xenon_Z = 65.*cm;

  while ( true ) {
    G4double phi0 = 2*pi*G4UniformRand();
    G4double r0 = inactive_Xenon_R*sqrt(G4UniformRand());
    G4double x0 = inactive_Xenon_Z*(2*(G4UniformRand() - 0.5));
    G4double y0 = r0*cos(phi0);
    G4double z0 = r0*sin(phi0);

    G4ThreeVector point = G4ThreeVector(x0, y0, z0);
    G4Navigator* theNavigator =
                 G4TransportationManager::GetTransportationManager()->
                                                     GetNavigatorForTracking();
    G4VPhysicalVolume* volume = theNavigator->LocateGlobalPointAndSetup(point);

    if (volume->GetName()== fXeComponent) {
      //G4cout << XeComponent << " chosen" << G4endl;
      fParticleGun1->SetParticlePosition(G4ThreeVector(x0, y0, z0));
      break;
    }
  } */

  // Draw a random generation position in LXe (Active and Inactive)
  
  G4ThreeVector decayVertex;
  GetUnifRandPosInLXe(decayVertex);
  //G4cout << "Final point = " << decayVertex << G4endl;
  fParticleGun1->SetParticlePosition(decayVertex);
  //G4cout << "fnCaptureXe_Method = " << fnCaptureXe_Method << G4endl;

  if (fnCaptureXe_Method == "InternalConversions" ||
      fnCaptureXe_Method == "RandomGammas") {

     nCapture_Xe136Info decay[] = {
       {  3, 4025,    0.00, 100/135. }, // skip to seen lines
       { 32, 4025,-3424.47,  23/135. }, // special to 601
       { 30, 4025,-3039.33,   1. }, // special to 986
       //
       { 30, 4025, 3424.47, 46/87. }, // 601
       { 28, 4025, 3039.33, 15/87. }, // 986
       { 17, 4025, 2184.04, 12/87. }, // 1842
       { 11, 4025, 2088.93,  5/87. }, // 1936
       {  9, 4025, 1829.38,  1/87. }, // 2196
       {  4, 4025, 1535.15,  7/87. }, // 2490
       {  1, 4025, 1416.68,  1. }, // 2609
       { 23, 2609, 2007.80, 10/110. }, // 601
       { 13, 2609,  893.30, 1. }, // 1716
       { 22, 2490, 2490.38, 84/267. }, // 0
       { 20, 2490, 1889.21, 25/267. }, // 601
       { 18, 2490, 1504.30, 100/267. }, // 986
       { 14, 2490, 1187.55, 1. }, // 1303
       { 13, 2196,  893.42, 1. }, // 1303
       { 17, 1936, 1936.05, 63/340. }, // 0
       { 15, 1936, 1335.00, 71/340. }, // 601
       { 13, 1936,  949.85, 100/340. }, // 986
       {  9, 1936,  633.32, 46/340. }, // 1303
       {  6, 1936,  267.92, 1. }, // 1668
       { 12, 1842, 1841.49, 30/130. }, // 0
       {  6, 1842,  538.76, 1. }, // 1303
       { 10, 1716, 1715.55, 100/170. }, // 0
       {  8, 1716, 1114.50, 46/170. }, // 601
       {  3, 1716,  412.82, 1. }, // 1303
       {  6, 1668, 1067.08, 100/167. }, // 601
       {  4, 1668,  681.93, 1. }, // 986
       {  5, 1303, 1302.73, 100/112. }, // 0
       {  3, 1303,  701.68, 10/112. }, // 601
       {  1, 1303,  316.53, 1. }, // 986
       {  1,  986,  385.15, 1. }, // 601
       {  1,  601,  601.05, 1. }, // 0
       {  0,    0,    0.00, 1. },
     };

     G4int decay_num = sizeof(decay)/sizeof(decay[0]);
     //G4int decay_num = GetNumArrayElements();

     //G4cout << "At beginning of loop" << G4endl;
     G4double p = G4UniformRand();
     //G4cout << "p = " << p << G4endl;
     for (G4int j = 0; 0<= j && j < decay_num;) {
        if (p < decay[j].ratio) {
           G4double energy = decay[j].gamma;
           //G4cout << "Energy = " << energy << G4endl;
           if (energy > 0.) {
              //G4cout << "In Energy > 0 " << G4endl;
              fParticleGun1->SetParticleEnergy(energy*keV);
              G4ThreeVector momentumDirection;
              GetIsotropicDirection(momentumDirection);
              fParticleGun1->SetParticleMomentumDirection(momentumDirection);
              fParticleGun1->GeneratePrimaryVertex(anEvent);
           }
           if ( energy < 0.) {
              if (fnCaptureXe_Method == "InternalConversions") {
                 //G4cout << "electron being generated" << G4endl;
                 fParticleGun1->SetParticleDefinition(G4Electron::Electron());
                 energy = -decay[j].gamma;
                 G4ThreeVector momentumDirection;
                 GetIsotropicDirection(momentumDirection);
                 fParticleGun1->SetParticleMomentumDirection(momentumDirection);
                 fParticleGun1->SetParticleEnergy(energy*keV);
                 fParticleGun1->GeneratePrimaryVertex(anEvent);
                 fParticleGun1->SetParticleDefinition(G4Gamma::Gamma());
              }
              else if (fnCaptureXe_Method == "RandomGammas") {
                 G4double level = (G4UniformRand() + G4UniformRand() +1.)*1000.;
                 G4double target = decay[j].level+decay[j].gamma;
                 if (target > 610.) level += 500.;
                 energy = decay[j].level - level;
                 //G4cout << "Energy in 1st 'else' statement" << G4endl;
                 //G4cout << "Energy = " << energy << G4endl;
                 fParticleGun1->SetParticleEnergy(energy*keV);
                 G4ThreeVector momentumDirection;
                 GetIsotropicDirection(momentumDirection);
                 fParticleGun1->SetParticleMomentumDirection(momentumDirection);
                 fParticleGun1->GeneratePrimaryVertex(anEvent);
                 energy = level - target;
                 //G4cout << "Energy in 2nd 'else' statement" << G4endl;
                 //G4cout << "Energy = " << energy << G4endl;
                 fParticleGun1->SetParticleEnergy(energy*keV);
                 GetIsotropicDirection(momentumDirection);
                 fParticleGun1->SetParticleMomentumDirection(momentumDirection);
                 fParticleGun1->GeneratePrimaryVertex(anEvent);
              }
           }
           if (decay[j].jump == 0) break;
           else j += decay[j].jump;
           //G4cout << "Creating new Random Number " << G4endl;
           p = G4UniformRand();
        } else {
           p -= decay[j].ratio;
           j++;
        }
     }

  } else if (fnCaptureXe_Method == "ImbalancedCascade") {
  //**** Start choosing energies of gammas from 4025keV cascade ****//
     //G4cout << "imbalance" << G4endl;
     nCapture_Xe136Info decay[] = {
       {4,    4025,      0.0,  85./111.},
       {12,   4025,      0.0,   3./111.},
       {20,   4025,      0.0,   8./111.},
       {20,   4025,      0.0,   1.     },
       {19,   4025,   3.4245,   46./85.},
       {17,   4025,   3.0394,   15./85.},
       {10,   4025,   2.1839,   12./85.},
       {6,    4025,   2.0889,   5.3/85.},
       {1,    4025,   1.5351,   1.     },
       {15,   2490,   2.4905,   1.6/5.1},
       {13,   2490,   1.8892,   0.5/5.1},
       {11,   2490,   1.5043,   1.9/5.1},
       {8,    2490,   1.1876,   1.     },
       {10,   1937,   1.3355,   2.5/8.1},
       {8,    1937,   0.9502,    3.5/8.1},
       {3,    1937,   0.2684,    1.     },
       {8,    1841,   1.8415,   3.0/12.9},
       {3,    1841,   0.5389,    1.     },
       {5,    1668,   1.0671,   1.9/2.9},
       {3,    1668,   0.6820,    1.     },
       {4,    1303,   1.3027,   12./14.3},
       {2,    1303,   0.7012,    1.     },
       {1,    986,    0.3852,    1.     },
       {1,    601,    0.6010,    1.     },
       {0,    0,      0.,        1.     },
     };

     G4int decay_num = sizeof(decay)/sizeof(decay[0]);
     //G4int decay_num = GetNumArrayElements();

     G4double p = G4UniformRand();
     for (G4int j = 0; 0<= j && j < decay_num;){
         if (p < decay[j].ratio){
            G4double energy = decay[j].gamma;
            if (energy > 0.){
               fParticleGun1->SetParticleEnergy(energy*MeV);
               G4ThreeVector momentumDirection;
               GetIsotropicDirection(momentumDirection);
               fParticleGun1->SetParticleMomentumDirection(momentumDirection);
               // If energy > 0, a gamma is generated, but that doesn't mean
               // the cascade is finished!
               fParticleGun1->GeneratePrimaryVertex(anEvent);
            }
            // Only break out of loop when 'jump' == 0
                  if (decay[j].jump == 0) break;
                   else j += decay[j].jump;
            p = G4UniformRand();
         } else {
                  p -= decay[j].ratio;
            j++;
         }
     }
  } else {
    G4cout << "WARNING: Unknown fnCaptureXe_Method: \"" << fnCaptureXe_Method << "\"" <<  G4endl;
  }
}

