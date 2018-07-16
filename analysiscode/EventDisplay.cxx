#include "TFile.h"
#include "TTree.h"
#include "Event.h"
#include <vector>
#include "TMath.h"
#include "MPPC.h"
#include "TChain.h"
#include "TH2.h"
#include "TVector3.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TStyle.h"

TH2F *Hist;
std::vector<TVector3> pos_sipm;
void Fill(int iSiPM,double val){
	int xmax = Hist->GetXaxis()->GetNbins();
	int ymax = Hist->GetYaxis()->GetNbins();
	for(int x=1; x<=xmax; ++x){
		double phi = Hist->GetXaxis()->GetBinCenter(x);
		for(int y=1; y<=ymax; ++y){
			double theta = TMath::ACos(Hist->GetYaxis()->GetBinCenter(y));
			TVector3 pos = TVector3(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta));
			if(pos*pos_sipm[iSiPM] > TMath::Cos(TMath::ATan(7.5/30.0))) Hist->SetBinContent(x,y,val);
		}
	}
}

int main(int argc, char** argv)
{

    double pi = 3.141592653589793;

    double r = 30.0;

    std::vector<double> theta;
    std::vector<double> phi;
    std::vector<bool> hasfilter;

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
    
    int n_sipm = theta.size();
    
    for(int i=0; i<n_sipm; ++i){
        pos_sipm.push_back(TVector3(TMath::Sin(theta[i])*TMath::Cos(phi[i]),
                                         TMath::Sin(theta[i])*TMath::Sin(phi[i]),
                                         TMath::Cos(theta[i])));
    }

    
    Hist = new TH2F("EventDisplay","",1000,-pi,pi,1000,-1.0,1.0);
    
    TChain* t = new TChain("T");
    t->AddFile(argv[1]);
    Event* ds = NULL;
    MPPC* mppc = NULL;
    
    t->SetBranchAddress("EV",&ds);
    
    int nentries = t->GetEntries();
    t->GetEntry(std::atoi(argv[2]));
    
    std::vector<double> SiPM(phi.size(),0.0);
    std::vector<TVector3> plotpos;

    int n_mppc = ds->GetMPPCCount();
    for(int i_mppc=0; i_mppc<n_mppc; ++i_mppc){
        mppc = ds->GetMPPC(i_mppc);
    	int id = mppc->GetID();
	int siPMid = id/4;
        Printf(Form("%d,%d has %d photons",siPMid,id-siPMid*4,mppc->GetHitCount()));
        if(hasfilter[siPMid]){
                SiPM[siPMid] += mppc->GetHitCount();
    	}
    }

    for(int i=0; i<SiPM.size(); ++i){
	    Fill(i,SiPM[i]);
    }

    gStyle->SetPalette(79);//kFall

    TImage *img = TImage::Create();
    TCanvas *c1 = new TCanvas("c1","c1",500,500);
    c1->cd();
    Hist->SetStats(false);
    Hist->Draw("colz");
    img->FromPad(c1);
    img->WriteImage(Form("Event_%d.png",std::atoi(argv[2])),TImage::kPng);

    return 0;

}
