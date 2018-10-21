#include "TMath.h"
#include "TRandom3.h"
#include "LOLXReadData.hh"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

double m_e = 511000.;
double nm2cm = 1.0e-7;
double density = 2.942;//g/cm^3

double d2N_dxdl(double wl,double beta){
	double n = LOLXReadData::GetXenon_n_stitch(1240./wl);
	//if(beta<(1./n)) return 0.;
	double value = 2.*TMath::Pi()*(1./137.);
	value *= 1.-1./(beta*beta*n*n);
	value *= 1./wl*wl;
    return TMath::Max(0.0,value);
}

double Getbeta(double en){
	double beta = 1. + en/m_e;
	beta = 1./(beta*beta);
	beta = 1.-beta;
	return TMath::Sqrt(beta);
}

int main(){

    TFile *fout = new TFile("MCElectronCherenkov.root","RECREATE");
    
    TH2F* CherenkovSpectrum = new TH2F("CherenkovSpectrum","",1000.,0,1000.,1000,0.,1000.);
    TH1F* RadiativeStoppingPower = new TH1F("RadiativeStoppingPower","",5000,0.,5);
    TH1F* RadiativeStoppingPower_MC = new TH1F("RadiativeStoppingPower_MC","",5000,0.,5);
    TH1F* XeIndex_formula = new TH1F("XeIndex_formula","",1100,0.,11.);
    TH1F* XeIndex_data = new TH1F("XeIndex_data","",1100,0.,11.);
    
    double electronenergy = 1000.;
    double dx = 1.0;
    double totalcherenkov = 0.0;
    double bin[1000];
    for(int i=0; i<1000; ++i)bin[i]=0.0;
    
    for(int i=1; i<=XeIndex_formula->GetNbinsX(); ++i){
        double eV = XeIndex_formula->GetBinCenter(i);
        XeIndex_formula->SetBinContent(i,LOLXReadData::GetXenon_n(eV));
        if(eV>8.0) XeIndex_data->SetBinContent(i,LOLXReadData::GetXenon_IndexofRefraction(eV));
    }
 
    while(electronenergy<5.0e6){
        printf("electron energy = %f\n",electronenergy);
        double newEnergy = electronenergy+density*dx*nm2cm*LOLXReadData::GeteStop_Total(electronenergy/1.0e6)*1.0e6;
        double beta = Getbeta(0.5*(electronenergy+newEnergy));
        double radiativeloss = 0.0;
        for(int i = 20; i<1000; ++i){
            double wl = (double)i;
            double dN = d2N_dxdl(wl,beta)*dx;
            radiativeloss += dN*(1240./wl);
            bin[i] = dN;
            totalcherenkov += dN;
        }
        for(int i = 148; i<1000; ++i){
            CherenkovSpectrum->Fill(totalcherenkov,(double)i,bin[i]);
        }
        
        RadiativeStoppingPower->Fill(electronenergy,LOLXReadData::GeteStop_Radiation(electronenergy/1.0e6));
        RadiativeStoppingPower_MC->Fill(electronenergy,radiativeloss/(1.0e6*density*dx*nm2cm));
        
        electronenergy = newEnergy;
    }
    fout->cd();
    CherenkovSpectrum->Write();
    RadiativeStoppingPower->Write();
    RadiativeStoppingPower_MC->Write();
    XeIndex_formula->Write();
    XeIndex_data->Write();
    fout->Close();
    
    return 0;
}

