#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"

TRandom3 *randomGen;
double pixelpitch = 50.;
double pixeldepth = 50.;
TH1F* AveSpecHist;
TGraph *AveSpecGraph;
TGraph *AbsLengthGraph;
TGraph *SiPMPDE;

double GetPixelX(double x){
	return x/pixelpitch;
}

double GetPixelY(double y){
	return y/pixelpitch;
}

TVector3 GetRandomPosition(){

	double x = pixelpitch*0.5*(1.-2.*randomGen->Uniform());
	double y = pixelpitch*0.5*(1.-2.*randomGen->Uniform());
	double z = pixeldepth*TMath::Power(randomGen->Uniform(),7.0);
	return TVector3(x,y,z) + TVector3(3000.,3000.,0.);
}

TVector3 GetRandomDirection(){

	double theta = TMath::ACos(1.0-2.0*randomGen->Uniform());
	double phi = TMath::Pi()*2.*randomGen->Uniform();
	return TVector3(TMath::Cos(phi)*TMath::Sin(theta),TMath::Sin(phi)*TMath::Sin(theta),TMath::Cos(theta));

}
 
double SiAbsWave[63]={251.241,280.135,289.164,298.194,305.417,314.446,327.088,339.729,352.370,359.593,368.623,372.234,377.652,383.069,386.681,397.516,408.352,419.187,431.828,444.469,460.722,473.363,489.616,505.869,525.733,543.792,563.656,585.327,605.191,626.862,650.338,673.814,695.485,720.767,747.855,776.749,807.449,838.148,865.237,886.907,913.995,933.860,955.530,977.200,998.871,1016.930,1033.182,1045.823,1058.465,1071.106,1087.358,1107.223,1125.282,1141.534,1154.176,1163.205,1170.428,1174.040,1177.652,1181.264,1188.487,1193.905,1700.0};
double SiAbsLength[63]={0.0052953,0.00423106,0.00439229,0.00529538,0.00638416,0.00714213,0.00829456,0.00893873,0.00963293,0.00963293,0.0129923,0.0181911,0.0254702,0.0370207,0.0558598,0.0908315,0.142276,0.191894,0.268679,0.349078,0.470818,0.635015,0.765580,0.958161,1.19918,1.44574,1.80942,2.10138,2.53344,2.94223,3.54718,4.27651,4.96655,5.98772,7.49393,9.37902,12.6499,16.4352,21.3533,27.7431,37.4184,50.4680,70.6623,98.9372,154.972,233.834,352.828,512.833,833.899,1306.198,1898.551,3326.915,5615.907,9841.004,17244.830,29109.659,40757.620,66274.385,103810.503,156637.445,219314.469,318772.081,318772.081};
double AvalancheWave[72]={400.000,462.069,477.155,494.397,507.328,524.569,537.500,554.741,567.672,580.603,595.690,604.310,615.086,625.862,636.638,649.569,660.345,668.966,675.431,684.052,692.672,703.448,712.069,720.690,729.310,740.086,748.707,761.638,776.724,793.966,811.207,826.293,841.379,852.155,865.086,878.017,888.793,906.034,921.121,942.672,977.155,994.397,1013.793,1031.034,1041.810,1052.586,1076.293,1097.845,1108.621,1128.017,1151.724,1175.431,1203.448,1237.931,1276.724,1315.517,1362.931,1390.948,1425.431,1442.672,1459.914,1490.086,1524.569,1554.741,1574.138,1580.603,1587.069,1595.690,1604.310,1617.241,1625.862,1656.034};
double AvalancheAmp[72]={0.000,0.805,0.805,0.604,1.208,1.611,2.215,2.819,3.624,4.832,5.839,6.644,7.651,9.262,10.470,12.081,13.893,15.101,16.711,18.121,19.732,21.141,23.356,24.966,26.376,27.785,28.389,28.993,29.396,29.799,30.201,30.805,32.013,33.624,35.034,36.443,38.255,38.859,39.664,40.067,38.658,40.268,41.074,42.282,46.107,48.926,51.141,52.550,54.564,56.779,57.383,56.376,54.765,55.168,53.758,52.752,48.926,46.510,43.289,42.282,45.906,46.309,43.691,40.671,35.839,30.201,25.369,18.523,10.671,5.034,0.000,0.000};

double SiPMEffwave[63]={60.000,121.538,127.132,130.209,132.727,137.202,140.139,142.937,144.475,146.713,148.251,150.629,153.146,154.965,156.223,158.181,160.000,161.118,161.958,163.916,165.034,166.293,168.111,170.069,172.727,175.524,178.041,180.979,184.335,188.251,192.307,195.804,199.580,326.627,330.177,333.727,338.461,345.562,357.396,365.680,376.331,385.798,392.899,408.284,426.035,442.603,460.355,485.207,505.325,521.893,550.295,572.781,605.917,636.686,666.272,700.591,739.644,775.147,814.201,862.721,899.408,1000.000,1700.00};
double SiPMEff[63]={0.0000,0.0533,0.0688,0.0799,0.0890,0.1035,0.1108,0.1180,0.1200,0.1200,0.1190,0.1176,0.1153,0.1114,0.1071,0.0990,0.0932,0.0923,0.0961,0.1120,0.1168,0.1183,0.1183,0.1183,0.1189,0.1198,0.1199,0.1199,0.1195,0.1195,0.1196,0.1196,0.1197,0.1208,0.1554,0.1682,0.1820,0.1929,0.1999,0.2069,0.2207,0.2346,0.2425,0.2515,0.2604,0.2664,0.2675,0.2647,0.2579,0.2471,0.2255,0.2059,0.1794,0.1559,0.1363,0.1128,0.0903,0.0717,0.0541,0.0356,0.0259,0.0000,0.000};

void SetAvalancheSpec(){

	AveSpecGraph = new TGraph();
	for(int i=0; i<72; ++i){
		AveSpecGraph->SetPoint(i,AvalancheWave[i],AvalancheAmp[i]);
	}
	AbsLengthGraph = new TGraph();
	for(int i=0; i<63; ++i){
                AbsLengthGraph->SetPoint(i,SiAbsWave[i],SiAbsLength[i]);
        }
	AveSpecHist = new TH1F("AveSpecHist","AveSpecHist",1100,500.,1600.);
	for(int i=1; i<1101; ++i){
		AveSpecHist->SetBinContent(i,AveSpecGraph->Eval(AveSpecHist->GetXaxis()->GetBinCenter(i)));
	}
	SiPMPDE = new TGraph();
	for(int i=0; i<63; ++i){
		SiPMPDE->SetPoint(i,SiPMEffwave[i],SiPMEff[i]);
	}
}

double GetRandomWavelength(){
	return AveSpecHist->GetRandom();
}

double GetRandomDistancetoAbs(double wl){

	return -AbsLengthGraph->Eval(wl)*TMath::Log(randomGen->Uniform());
}

int GetHistBin(double wl){
	return (wl-500)/50;
}

int main()
{
	randomGen = new TRandom3();
	SetAvalancheSpec();
	TFile *fout = new TFile("OpticalCrossTalkParameters.root","RECREATE");
	//500 to 1600 in steps of 50 nm
	TH2F* AngularDistEmission[22];
        for(int i=0; i<22; ++i) AngularDistEmission[i] = new TH2F(Form("PositionExcape_%d_%d",500+50*i,550+50*i),Form("PositionExcape_%d_%d",500+50*i,550+50*i),120,0,119,120,0,119);
	TH2F* PositionAbsorbed[22];
	for(int i=0; i<22; ++i) PositionAbsorbed[i] = new TH2F(Form("PositionAbsorbed_%d_%d",500+50*i,550+50*i),Form("PositionAbsorbed_%d_%d",500+50*i,550+50*i),120,0,119,120,0,119);
	TH1F* AbsorptionSpec = new TH1F("AbsorptionSpec","AbsorptionSpec",22,500,1600);
	TH1F* EmissionSpec = new TH1F("EmissionSpec","EmissionSpec",22,500,1600);
	TH1F* RedetectionSpectrum = new TH1F("RedetectionSpectrum","RedetectionSpectrum",22,500,1600);

	int nsim = 1e7;
	int nemit = 0;
	int nabs = 0;

	for(int i=0; i<nsim; ++i){

		if(nsim%1000==0) printf("%d / %d\n",i,nsim);
		TVector3 pos = GetRandomPosition();
		TVector3 dir = GetRandomDirection();

		double wl = GetRandomWavelength();
		double DistToAbs = GetRandomDistancetoAbs(wl);
		double DistToEdge = 0.0;

		if(dir.z() == 0.0){
			DistToEdge = 10.*DistToAbs;
		}else if(dir.z()>0.0){
			DistToEdge = (pixeldepth-pos.z())/dir.z();
		}else{
			DistToEdge = -pos.z()/dir.z();
		}

		if(DistToAbs<DistToEdge){
			pos = pos + DistToAbs*dir;
			if(pos.x()>0. && pos.x()<6000. && pos.y()>0. && pos.y()<6000.) ++nabs;
			AbsorptionSpec->Fill(wl);
			PositionAbsorbed[GetHistBin(wl)]->Fill(GetPixelX(pos.x()),GetPixelY(pos.y()));
		}else if(dir.z()>0.0){
			pos = pos + DistToEdge*dir;
                        if(pos.x()>0. && pos.x()<6000. && pos.y()>0. && pos.y()<6000.) ++nemit;
                        EmissionSpec->Fill(wl);
                        AngularDistEmission[GetHistBin(wl)]->Fill(GetPixelX(pos.x()),GetPixelY(pos.y()));
		}

	}

	for(int i=1; i<23; ++i){
		RedetectionSpectrum->SetBinContent(i,SiPMPDE->Eval(RedetectionSpectrum->GetXaxis()->GetBinCenter(i))*EmissionSpec->GetBinContent(i));
	}

	printf("Absorbed %d / %d \n",nabs,nsim);
	printf("Emmit %d / %d \n",nemit,nsim);
	fout->cd();
	for(int i=0; i<22; ++i) AngularDistEmission[i]->Write();
	for(int i=0; i<22; ++i) PositionAbsorbed[i]->Write();
	AbsorptionSpec->Write();
	EmissionSpec->Write();
	RedetectionSpectrum->Write();
	AveSpecHist->Write();
	AbsLengthGraph->Write("AbsLength");
	SiPMPDE->Write("SiPMPDE");
	fout->Close();

	return 0;
}
