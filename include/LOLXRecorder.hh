#ifndef _LOLXRecorder_
#define _LOLXRecorder_
#include "LXeRecorderBase.hh"
#include "LXeUserEventInformation.hh"
#include "G4SDManager.hh"
#include "LXePMTHit.hh"
#include "TFile.h"
#include "TTree.h"
#include "Event.h"
#include "globals.hh"
#include "TMath.h"
#include "G4SystemOfUnits.hh"

class LOLXRecorder : public LXeRecorderBase{

  public:

  LOLXRecorder();
  ~LOLXRecorder(){};
  
  virtual void RecordBeginOfRun(const G4Run*) override {
        
		fout = new TFile(LOLXOutputFileName,"RECREATE");
		tree = new TTree("T","T");
		event = new Event();
		tree->Branch("EV",&event);
	}

    virtual void RecordEndOfRun(const G4Run*) override {
        
		fout->cd();
		tree->Write();
		fout->Close();
	}
    
    virtual void RecordBeginOfEvent(const G4Event*){
        
		event->Clear();
	}

    virtual void RecordEndOfEvent(const G4Event* anEvent){

         LXeUserEventInformation* eventInformation =(LXeUserEventInformation*)anEvent->GetUserInformation();	    

        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        LXePMTHitsCollection* pmtHC = 0;
        G4HCofThisEvent* hitsCE = anEvent->GetHCofThisEvent();
        pmtHC = (LXePMTHitsCollection*)(hitsCE->GetHC(SDman->GetCollectionID("sipmHitCollection")));
        
		event->SetEventID((long int)anEvent->GetEventID());
		if(pmtHC){
            G4int n_sipm=pmtHC->entries();
            for(G4int i_sipm=0;i_sipm<n_sipm;i_sipm++){
                MPPC* mppc = new MPPC((*pmtHC)[i_sipm]->GetID());
                G4int n_photons = (*pmtHC)[i_sipm]->GetPhotonCount();
                for(G4int i_photon=0;i_photon<n_photons;i_photon++){
                    mppc->AddSensorHit((*pmtHC)[i_sipm]->GetHitTime(i_photon));
		    if((*pmtHC)[i_sipm]->GetHitProcess(i_photon) == 1) 
			mppc->AddScintillationHit((*pmtHC)[i_sipm]->GetHitTime(i_photon));
		    else if((*pmtHC)[i_sipm]->GetHitProcess(i_photon) == 0)
			mppc->AddCherenkovHit((*pmtHC)[i_sipm]->GetHitTime(i_photon));

                }
                event->AddMPPC(mppc);
            }
        }

   	event->SetMCSintillation(eventInformation->GetPhotonCount_Scint());
        event->SetMCCherenkov(eventInformation->GetPhotonCount_Ceren());

	 int nsamples = eventInformation->GetCherenkovSamples();
	 for(int i=0; i<nsamples; ++i){
    		double photonenergy = eventInformation->GetCherenkovPhotonEnergy(i);
    		double electronenergy = eventInformation->GetElectronMomentum(i);
		photonenergy = 1240./(photonenergy*1.0e6);
		electronenergy/=.511;
		electronenergy += 1.0;
		electronenergy = 1./electronenergy;
		electronenergy = TMath::Sqrt(electronenergy);
		electronenergy = 1.0 - electronenergy;
		event->AddCherenkovWavelength(photonenergy);
        	event->AddElectronBeta(TMath::Sqrt(electronenergy));
	}

        tree->Fill();
        counter++;
        Printf(Form("%d Events saved",counter));
    }

    virtual void RecordTrack(const G4Track*) {return;}
    virtual void RecordStep(const G4Step*) {return;}
    
    void SetFileName(char* name){LOLXOutputFileName = name;}

  private:

  TFile *fout;
  TTree *tree;
  Event *event;
  int counter;
  char* LOLXOutputFileName;

};

#endif
