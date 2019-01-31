void SiPM_pixelSim(){

	//rate_recovery_pitch_size

	TRandom3 *r = new TRandom3();
	TGraph *measuredrate = new TGraph();

	int npix = 6*1000/50;
	double SiPM[npix][npix];
	double recovertime= 200.;

	TFile* fout = new TFile(Form("SiPM_Steady_Rate_%f.root",recovertime),"RECREATE");

	for(int px=0; px<npix; ++px)
		for(int py=0; py<npix; ++py)
			SiPM[px][py]=-recovertime;

	int i=0;
	for(double rate = 0.1; rate<1000; rate+=0.1){
		 int counts = 0;
		for(double time=0; time<1e9;){
			int px = r->Uniform(0,npix);
			int py = r->Uniform(0,npix);
			time += r->Exp(1./rate);
		
			if(time - SiPM[px][py] > recovertime){
				++counts;
				SiPM[px][py] = time;
			}
		}
		double fluxrate = (rate/recovertime)/(npix*npix);
		measuredrate->SetPoint(i++,fluxrate,(((double)counts)/1.0e9)/rate);
	}

	fout->cd();
	measuredrate->Write("measuredrate");
	fout->Close();
}
