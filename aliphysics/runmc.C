
#include <stdio.h>
#include <signal.h>

#include <TFile.h>
#include <TRandom3.h>
//#include <TDatabasePDG.h>

#include <AliBasicParticle.h>
#include <AliCFContainer.h>
#include <AliUEHistograms.h>

#include <Pythia8/Pythia.h>

typedef unsigned int uint;

using namespace Pythia8;

volatile sig_atomic_t termFlag = 0;
static void SignalHandler(int sig){
	printf("Caught SIGTERM. Exiting...\n");
	termFlag = 1;
}

int main() {
	uint seed = 1000;
	const TString &output = "AnalysisResults.root";
	setenv("PYTHIA8DATA","/root/alice/sw/ubuntu2204_x86-64/AliRoot/latest/PYTHIA8/pythia8243/xmldoc",1);
	struct sigaction sa;
	memset(&sa,0,sizeof(sa));
	sa.sa_handler = SignalHandler;
	sigaction(SIGTERM,&sa,0);

	TList *phistosOut = new TList();
	phistosOut->SetOwner(kFALSE);

	Pythia pythia;
	Event &event = pythia.event;
	pythia.readFile("PythiaTuneMonash.cmnd");
	pythia.readString("Random:setSeed = on");
	pythia.readString(Form("Random:seed = %u",seed));
	
	pythia.init();
	pythia.settings.listChanged();

	//TRandom3 *prng = new TRandom3(151221);
	//TFormula tf("detaAcceptance","1-(x < 0.0)*0.5");

	const char *pcustomBinning = "p_t_assoc: 1.0, 2.0, 3.0, 4.0\np_t_leading_course: 1.0, 2.0, 3.0, 4.0\neta: -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2\nmultiplicity:-0.5, 2.750, 5.250, 7.750, 12.750, 17.750, 22.750, 27.750, 32.750, 37.750, 42.750, 47.750, 52.750, 57.750, 62.750, 67.750, 72.750, 77.750, 82.750, 87.750, 92.750, 97.750, 250.000\n";
	AliUEHistograms *pueSame = new AliUEHistograms("AliUEHistogramsSame","5RD",pcustomBinning);
	AliUEHistograms *pueMixed = new AliUEHistograms("AliUEHistogramsMixed","5RD",pcustomBinning);

	pueSame->SetTwoTrackCutMinRadius(0.8);
	pueMixed->SetCheckEventNumberInCorrelation(kFALSE);
	
	AliUEHistograms *pueAll[] = {pueSame,pueMixed};
	for(auto &pue : pueAll){
		pue->SetSelectCharge(0);
		pue->SetSelectTriggerCharge(0);
		pue->SetSelectAssociatedCharge(0);
		pue->SetTriggerRestrictEta(-1);
		pue->SetOnlyOneEtaSide(0);
		pue->SetOnlyOneAssocEtaSide(0);
		pue->SetEtaOrdering(kFALSE);
		pue->SetPairCuts(-1,-1);
		pue->SetCustomCut(-1,-1,-1,-1);
		//pue->SetCutOnPhi(-1);
		//pue->SetCutOnRho(-1);
		pue->SetCutOnK0s(-1);
		pue->SetCutOnLambda(-1);
		pue->SetRejectResonanceDaughters(-1);
		//pue->SetTrackEtaCut(1.1);
		//pue->SetTrackEtaCut(0.9);
		pue->SetTrackEtaCut(4.0);
		pue->SetWeightPerEvent(kFALSE);
		pue->SetPtOrder(kTRUE);
		phistosOut->Add(pue);
	}

	Int_t nCentralityBins = pueSame->GetUEHist(2)->GetEventHist()->GetNBins(1);
	
#define MAX_TRACKS 2048
	//previous-event tracks for mixed correlations (for each mult bin)
	AliBasicParticle (*pmultParticles)[MAX_TRACKS] = new AliBasicParticle[nCentralityBins][MAX_TRACKS];
	TObjArray *ptracksPrevEvent = new TObjArray[nCentralityBins];

	//current event tracks
	AliBasicParticle particles[MAX_TRACKS];
	TObjArray tracks;

	//Generate events until SIGTERM
	for(uint i = 0; termFlag == 0; ++i){
		if(!pythia.next())
			break;

		//leading particle pT event selection
		/*double lppt = 0.0;
		for(uint j = 0, l = event.size(); j < l; ++j){
			if(!event[j].isFinal() || !event[j].isCharged() || TMath::Abs(event[j].eta()) > 0.9)
				continue;
			lppt = TMath::Max(lppt,event[j].pT());
		}
		if(lppt < 0.0)
			continue;*/

		uint mult = 0;

		tracks.Clear();
		for(uint j = 0, k = 0, l = event.size(); j < l; ++j){
			if(!event[j].isFinal() || !event[j].isCharged())// || !event[j].isHadron())
				continue;

			//include pT > 1.0 GeV particles only in the correlations
			if(event[j].pT() > 1.0){ //
				double phi = event[j].phi() < 0.0?2.0*TMath::Pi()+event[j].phi():event[j].phi(); //Add 2pi if negative
				new(&particles[k]) AliBasicParticle(event[j].eta(),phi,event[j].pT(),1);
				tracks.Add(&particles[k]);
				particles[k].SetUniqueID(i*MAX_TRACKS+j);
				++k;
			}

			//
			//if(prng->Uniform(0,1) > tf.Eval(event[j].eta()))
			//	continue;

			if(TMath::Abs(event[j].eta()) < 1.0 && event[j].pT() > 0.2)
				++mult; //count multiplicity in a simple manner
		}

		//get multiplicity bin
		double multD = (double)mult;
		Int_t bin = pueSame->GetUEHist(2)->GetEventHist()->GetAxis(1,0)->FindBin(multD)-1;
		if(bin >= nCentralityBins-1)
			continue;

		//calculate and fill the correlations
		double zvertex = 0.0;
		pueSame->FillEvent(multD,AliUEHist::kCFStepReconstructed); //Step does not matter in this case. Step 0 (MC truth) would be better label

		float weight = 1.0f;
		pueSame->FillCorrelations(multD,zvertex,AliUEHist::kCFStepReconstructed,&tracks,0,weight);
		if(ptracksPrevEvent[bin].GetEntriesFast() > 0)
			pueMixed->FillCorrelations(multD,zvertex,AliUEHist::kCFStepReconstructed,&tracks,&ptracksPrevEvent[bin],weight);

		//copy the current event tracks to previous event of the corresponding bin
		ptracksPrevEvent[bin].Clear();
		for(uint j = 0, l = tracks.GetEntriesFast(); j < l; ++j){
			pmultParticles[bin][j] = particles[j];
			ptracksPrevEvent[bin].Add(&pmultParticles[bin][j]);
		}
	}

	TFile *pfo = new TFile(output.Data(),"recreate");
	TDirectory *pmcDir = pfo->mkdir("PWG4_PhiCorrelations");
	pfo->cd("PWG4_PhiCorrelations");

	phistosOut->Write("histosPhiCorrelations",TObject::kSingleKey); 
	pfo->Close();
	delete pfo;

	delete pueSame;
	delete pueMixed;

	//delete prng;

	delete []ptracksPrevEvent;
	delete []pmultParticles;

	delete phistosOut;
}

