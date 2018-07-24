#include "TFile.h"
#include "TTree.h"
#include "Event.h"
#include <vector>
#include <utility>
#include "TMath.h"
#include "MPPC.h"
#include "TChain.h"
#include "TH2.h"
#include "TVector3.h"

bool Comparison_pair_int_double_second(const std::pair<int,double> &a, const std::pair<int,double> &b)
{
	return a.second>b.second;
}

int main(int argc, char** argv)
{

    double pi = 3.141592653589793;
    
    double SiPM_w = 15.0;
    double SiPM_h = 2.5;
    double SiPM_r = 28.355044+SiPM_h/2.;
    double filter_w = 13.5;
    double filter_h = 1.0;
    double filter_r = 28.355044-filter_h/2.;
    
    std::vector<double> theta;
    std::vector<double> phi;
    std::vector<TVector3> pos_sipm;
    std::vector<bool> hasfilter;
    std::vector< std::vector<std::pair<int,double> > > neighbour;

    if(argc >3 && std::atoi(argv[3]) == 1){
    	theta.push_back(pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    	theta.push_back(pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    	theta.push_back(pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    	theta.push_back(pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);

    	theta.push_back(pi*55./180.); phi.push_back(pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*55./180.); phi.push_back(3.*pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*55./180.); phi.push_back(5.*pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*55./180.); phi.push_back(7.*pi/4.); hasfilter.push_back(true);

    	theta.push_back(3.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    	theta.push_back(3.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    	theta.push_back(3.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    	theta.push_back(3.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);

    	theta.push_back(pi/2.); phi.push_back(pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(3.*pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(5.*pi/8.); hasfilter.push_back(false);
    	theta.push_back(pi/2.); phi.push_back(7.*pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(9.*pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(11.*pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(13.*pi/8.); hasfilter.push_back(true);
    	theta.push_back(pi/2.); phi.push_back(15.*pi/8.); hasfilter.push_back(true);

    	theta.push_back(5.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    	theta.push_back(5.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    	theta.push_back(5.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    	theta.push_back(5.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);


    	theta.push_back(pi*125./180.); phi.push_back(pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*125./180.); phi.push_back(3.*pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*125./180.); phi.push_back(5.*pi/4.); hasfilter.push_back(true);
    	theta.push_back(pi*125./180.); phi.push_back(7.*pi/4.); hasfilter.push_back(true);

    	theta.push_back(7.*pi/8.); phi.push_back(0.0); hasfilter.push_back(true);
    	theta.push_back(7.*pi/8.); phi.push_back(pi/2.); hasfilter.push_back(true);
    	theta.push_back(7.*pi/8.); phi.push_back(pi); hasfilter.push_back(true);
    	theta.push_back(7.*pi/8.); phi.push_back(3.*pi/2.); hasfilter.push_back(true);

    	for(int i=0; i<32; ++i){
        	pos_sipm.push_back(TVector3(SiPM_r*TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                            SiPM_r*TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                            SiPM_r*TMath::Cos(theta[i])));
    	}

    }else{
	std::vector<double> rad;
	double Package_sizeXY = 15.0;
    	double Package_sizeZ = 4.0;
    	double gap = 0.1;
	TVector3 vert[2];
    	vert[0] = TVector3(0.0,0.0,(Package_sizeXY+gap)/2.);
    	vert[1] = TVector3(0.0,0.0,-(Package_sizeXY+gap)/2.);

    	TVector3 quad[4];
    	quad[0] = TVector3( (Package_sizeXY+gap)/2., (Package_sizeXY+gap)/2.,0.0);
    	quad[1] = TVector3( (Package_sizeXY+gap)/2.,-(Package_sizeXY+gap)/2.,0.0);
    	quad[2] = TVector3(-(Package_sizeXY+gap)/2., (Package_sizeXY+gap)/2.,0.0);
    	quad[3] = TVector3(-(Package_sizeXY+gap)/2.,-(Package_sizeXY+gap)/2.,0.0);

    	double rad_temp = Package_sizeXY*sqrt(2)-1.5;
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
    	rad_temp = Package_sizeXY+1.0+0.5+Package_sizeZ/2.;
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(0.0); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(0.0); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(0.0); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(0.0); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi); hasfilter.push_back(false);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi); hasfilter.push_back(true);
    	rad.push_back(rad_temp); phi.push_back(0.0); theta.push_back(pi); hasfilter.push_back(true);

	for(int i=0; i<24; ++i){
        	double x = TMath::Sin(theta[i])*TMath::Cos(phi[i]);
        	double y = TMath::Sin(theta[i])*TMath::Sin(phi[i]);
        	double z = TMath::Cos(theta[i]);
        	if(i<16){
            		double radius = rad[i];
            		pos_sipm.push_back(TVector3(x*radius+vert[i%2].x(),y*radius+vert[i%2].y(),z*radius+vert[i%2].z()));
        	}else{
            		double radius = rad[i];
            		pos_sipm.push_back(TVector3(x*radius+quad[i%4].x(),y*radius+quad[i%4].y(),z*radius+quad[i%4].z()));
        	}
    	}	

    }
    int n_sipm = theta.size();

    for(int i=0; i<n_sipm; ++i){
	   neighbour.push_back(std::vector< std::pair<int,double> >(0,std::pair<int,double>(0,0.)));
	   for(int j=0; j<n_sipm; ++j){
		neighbour[i].push_back(std::pair<int,double>(j,pos_sipm[i].Unit()*pos_sipm[j].Unit()));	
	   }
	   std::sort(neighbour[i].begin(),neighbour[i].end(), Comparison_pair_int_double_second);
	   if(i==0){
		for(int j=0; j<n_sipm; ++j)
		Printf(Form("%d: %d %f",i,j,neighbour[i][j].second));
	   }
    }

    TFile *fout = new TFile(argv[2],"RECREATE");
    
    TH2F* Filtered_NonFiltered = new TH2F("Filtered_NonFiltered","",10000,0.,10000.,10000,0,10000);
    TH2F *fMaxHist = new TH2F("fMax","",10000,0.,10000.,100,0.,1.0);
    TH2F *GroupfMaxHist = new TH2F("GroupfMax","",10000,0.,10000.,100,0.,1.0);
    TH2F *MCScintillation_Vs_MCCherenkov = new TH2F("MCScintillation_Vs_MCCherenkov","",10000,0.,10000.,10000,0,10000.);
    TH2F *NMCCherenkov_Vs_CherSpectrum = new TH2F("NMCCherenkov_Vs_CherSpectrum","",10000,0.,10000.,700,0.,700.);
    TH2F *NonFiltered_vs_CherSpectrum = new TH2F("NonFiltered_vs_CherSpectrum","",10000,0.,10000.,700,0.,700.);
    TH1F *SiPM_TotalHits = new TH1F("SiPM_TotalHits","",n_sipm,-0.5,(double)n_sipm -0.5);
    TH1F *SiPM_ScintillationHits = new TH1F("SiPM_ScintillationHits","",n_sipm,-0.5,(double)n_sipm -0.5);
    TH1F *SiPM_CherenkovHits = new TH1F("SiPM_CherenkovHits","",n_sipm,-0.5,(double)n_sipm -0.5);

    TChain* t = new TChain("T");
    t->AddFile(argv[1]);
    Event* ds = NULL;
    MPPC* mppc = NULL;
    
    t->SetBranchAddress("EV",&ds);
    
    int nentries = t->GetEntries();
    
    for(int i=0; i<nentries; ++i){
        t->GetEntry(i);
        double integral_nonfiltered = 0.;
        double integral_filtered = 0.;
       
	std::vector<double> groupfMax;
	std::vector<double> fMax;

	for(int j=0; j<n_sipm; ++j){
		groupfMax.push_back(0.0);
		fMax.push_back(0.0);
	}

        int n_mppc = ds->GetMPPCCount();
        for(int i_mppc=0; i_mppc<n_mppc; ++i_mppc){
            mppc = ds->GetMPPC(i_mppc);
            int id = mppc->GetID();
            int siPMid = id/4;
            //Printf(Form("%d: %d has %d photons",i,siPMid,mppc->GetHitCount()));
	    if(hasfilter[siPMid]){
                integral_filtered += mppc->GetHitCount();
		fMax[siPMid] += mppc->GetHitCount();
		for(int jj=0; jj<n_sipm; ++jj){
            		for(int j=0; j<5; ++j){
				if(siPMid == neighbour[jj][j].first){
					groupfMax[jj] += mppc->GetHitCount();
				}	
            		}
		}
            }else{
                integral_nonfiltered += mppc->GetHitCount();
            }
	    SiPM_TotalHits->Fill(siPMid,mppc->GetHitCount());
    	    SiPM_ScintillationHits->Fill(siPMid,mppc->GetScintillationHitCount());
    	    SiPM_CherenkovHits->Fill(siPMid,mppc->GetCherenkovHitCount());
        }
	double max_groupfmax = 0.0;
	double max_fmax = 0.0;
	for(int jj=0; jj<n_sipm; ++jj){
		 groupfMax[jj] /= integral_filtered;
		 fMax[jj] /= integral_filtered;
		 if(groupfMax[jj]>max_groupfmax) max_groupfmax = groupfMax[jj];
		 if(fMax[jj]>max_fmax) max_fmax = fMax[jj];
	}
        //Printf(Form("%d: Filtered = %f unfiltered %f",i,integral_filtered,integral_nonfiltered)); 
        Filtered_NonFiltered->Fill(integral_nonfiltered, integral_filtered);
	fMaxHist->Fill(integral_filtered,max_fmax);
	GroupfMaxHist->Fill(integral_filtered,max_groupfmax);

	double fillfactor = 1.0;
	//if(argc >3 && std::atoi(argv[3]) == 1) fillfactor = 0.53;

	MCScintillation_Vs_MCCherenkov->Fill(fillfactor*(ds->GetMCSintillation()*(1.0/(n_sipm))+ds->GetMCCherenkov()*(1.0/(n_sipm))),
					     fillfactor*(ds->GetMCCherenkov()*(((double)(n_sipm-1))/(n_sipm))));
	int nsamples = ds->GetCherenkovSpectrumCount();
	for(int jj=0; jj<nsamples; ++jj){
		NMCCherenkov_Vs_CherSpectrum->Fill(ds->GetMCCherenkov(),ds->GetCherenkovWavelength(jj));
    		NonFiltered_vs_CherSpectrum->Fill(integral_nonfiltered,ds->GetCherenkovWavelength(jj));
	}

    }

    Filtered_NonFiltered->Write();
    fMaxHist->Write();
    GroupfMaxHist->Write();
    MCScintillation_Vs_MCCherenkov->Write();
    NMCCherenkov_Vs_CherSpectrum->Write();
    NonFiltered_vs_CherSpectrum->Write();
    SiPM_TotalHits->Write();
    SiPM_ScintillationHits->Write();
    SiPM_CherenkovHits->Write();
    fout->Close();

    return 0;

}
