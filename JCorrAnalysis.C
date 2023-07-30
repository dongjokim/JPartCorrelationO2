#include "TRandom.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
#include "src/JBaseTrack.h"
#include "src/JTreeDataManager.h"
#include "src/JBaseEventHeader.h"
#include "src/JHistos.h"
// AliPhysics includes 
#include <AliBasicParticle.h>
#include <AliCFContainer.h>
#include <AliUEHistograms.h>

typedef unsigned int uint;
using namespace std;

int main(int argc, char **argv){

	if ( argc < 2 ) {
		cout << "Usuage : " << endl;
		cout<<"+  "<<argv[0]<<" [outputFileList] [outputfile]"<<endl;
		exit(1);
	}

	char *inputfile = argv[1];
	char *outputfile = argv[2];

	TFile *pfo = new TFile(outputfile,"recreate","Final analysis");
	TClonesArray *trackList = new TClonesArray("JBaseTrack",1000);
	TClonesArray *inputList = NULL;
	JTreeDataManager* dmg = new JTreeDataManager();

	JHistos *histos = new JHistos();
	histos->CreateQAHistos();

	JBaseEventHeader *eventHeader;

	//==== Read the Data files =====
	dmg->ChainInputStream(inputfile);
	AliJBaseEventHeader *eventHeader;
	//------------------------------------------------------------------------------
	// For AliPhysics code from Jasper
	//------------------------------------------------------------------------------
	TList *phistosOut = new TList();
	phistosOut->SetOwner(kFALSE);
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
	//------------------------------------------------------------------------------
	// End For AliPhysics code from Jasper
	//------------------------------------------------------------------------------


	int numberEvents = dmg->GetNEvents();
	int ieout = numberEvents/20;
	int EventCounter = 0;

	if (ieout<1) ieout=1;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Loop over all events
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	TStopwatch timer;
	timer.Start();

	for(Int_t evt=0; evt<numberEvents; evt++) {
		if(evt % ieout == 0) cout << evt << "\t" << int(float(evt)/numberEvents*100) << "%" << endl ;
		dmg->LoadEvent(evt);
		uint mult = 0; // UE
		tracks.Clear();

		eventHeader  = dmg->GetEventHeader();
		//cout << eventHeader->GetEventID() << endl;
		if(!dmg->IsGoodEvent()) continue;  // Vertex cut applied in IsGoodEvent and histo saved there too

		trackList->Clear();
		dmg->RegisterList(trackList, NULL);
		inputList = trackList;
		//---- assign input list ----
		int noAllTracks = inputList->GetEntries();
		// calculate multiplicity
		for(int itrack=0; itrack<noAllTracks; itrack++){
			JBaseTrack *trk = (JBaseTrack*)inputList->At(itrack);
			histos->fhPt->Fill(trk->Pt());
			histos->fheta->Fill(trk->Eta());
			//include pT > 1.0 GeV particles only in the correlations
			if(trk->Pt() > 1.0){ //
				double phi = trk->Phi() < 0.0?2.0*TMath::Pi()+trk->Phi():trk->Phi(); //Add 2pi if negative
				new(&particles[k]) AliBasicParticle(trk->Eta(),phi,trk->Pt(),1);
				tracks.Add(&particles[k]);
				particles[k].SetUniqueID(i*MAX_TRACKS+j);
				++k;
			}
			if(TMath::Abs(trk->Eta()) < 1.0 && trk->Pt() > 0.2)
				++mult; //count multiplicity in a simple manner
			//trk->Print("dmg");
		}
		cout << "mult = " << mult << endl;
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
		EventCounter++;
	}
	
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
	
	timer.Print();

	cout << "Total Number of Event scanned from input = " << numberEvents << endl;
	cout << "Total Number of Event used for analysis  = " << EventCounter << endl;
	cout<<"All files properly closed. Good Bye!"<<endl;

}

