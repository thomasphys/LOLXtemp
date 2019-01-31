LOLXDAQSimulation::LOLXDAQSimulation(){

    for(){
        SiPM.push_back(new SiPMPixelMap());
    }
    
    for(){
        SiPMGrouping.push_back(std::vector<int>(0,0));
        for(){
            SiPMGrouping[i].push_back(i);
        }
            
    }
}
LOLXDAQSimulation::~LOLXDAQSimulation();

void LOLXDAQSimulation::ResetForNewEvent(){
    for(){
        SiPM[i]->ClearSiPM();
    }
}

void LOLXDAQSimulation::GenerateWaveforms(){
    for(int channelid=0; channelid<SiPMGroup.size(); ++channelid){
        for(int i=0; i<SiPMGroup[channelid].size(); ++i){
            SiPM[SiPMGroup[channelid][i]]->TimeOrderHits();
            SiPM[SiPMGroup[channelid][i]]->FillOverlaps();
        }
        nsamples = (timemax - timemin+timeoffset*2.)/16.0;
        for(int j=0; j<nsamples; ++j){
            double time = j*16.+timemax-timeoffset;
            double value = noise->GetWhiteNoise();
            for(int i=0; i<SiPMGroup[channelid].size(); ++i){
                value += SiPM[SiPMGroup[channelid][i]]->GetSiPMCurrent(time);
            }
            fullwave[channelid]->AddPoint(value);
        }
    }
}

bool LOLXDAQSimulation::FindTrigger(){
    //simple threshold
    for(int j=0; j<nsamples; ++j){
        for(int channelid=0; channelid<fullwave.size(); ++channelid){
            if(fullwave[channelid]->GetSample(j)<threshold){
                triggerbin = j;
                return triggered;
            }
        }
    }
    return false;
}

Waveform* GetChannelTriggeredWaveform(G4Int channelid){
    fullwave[channelid]->GetWaveformSegment(triggerbin-presamples,triggerbin+postsamples,noisegen);
}

void LOLXDAQSimulation::AddPhotonHit(G4Int siPMid, G4double time, G4ThreeVector pos, G4double wavelength){
    timemin = time<timemin ? timemin:time;
    timemax = time>timemax ? timemax:time;
    SiPM[siPMid]->AddPhotonHit(time,wavelength,pos);
}

