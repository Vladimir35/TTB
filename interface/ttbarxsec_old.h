#ifndef TTBARAN_H
#define TTBARAN_H
#include <iostream>
#include <list>
#include <set>
#include "AnalyzerBase.h"
#include "URStreamer.h"
#include "URDriver.h"
#include "IDMuon.h"
#include "IDElectron.h"
#include "IDJet.h"
#include "IDMet.h"
#include "GenObject.h"
#include "TTBarGenPlots.h"
#include "TTBarPlots.h"
#include "TTBarSolver.h"
#include "TTBarResponse.h"
#include "TTBarResponse2D.h"
#include "TTBarResponse2Dvar.h"
#include "Permutation.h"
#include "BtagEff.h"
#include "JetScaler.h"
#include "BTagWeight.h"
#include "BHadronDecayWeights.h"

#include <TGraphErrors.h>

using namespace std;
class PDFuncertainty;

class ttbar : public AnalyzerBase
{
	friend class TTBarGenPlots;
    friend class TTBarPlots;
    friend class TTBarResponse;
    friend class BTagWeight;
    friend class JetScaler;

	private:
		int isDA = 0;
		map<int, set<int> >  runinfo;
		PDFuncertainty* pdfunc;
		//Collections
		//Gen:
		bool FULLHAD;
		bool SEMILEP;
		bool FULLLEP;
		bool SEMILEPACC;
		GenObject gps[8];
		GenObject rivetobjs[8];
		list<GenObject> sgenparticles;
	  	vector<GenObject*> genalljets;
	  	vector<GenObject*> genbjets;
	  	vector<GenObject*> gencjets;
	  	vector<GenObject*> genljets;
	  	vector<GenObject*> plleptons;
	  	vector<GenObject*> plphotons;
		TLorentzVector gentq;
		TLorentzVector gentqbar;
		TLorentzVector gentqlep;
		TLorentzVector gentqhad;
		TLorentzVector gent;
		TLorentzVector gentbar;
		TLorentzVector gentlep;
		TLorentzVector genthad;

		Permutation genallper;
		Permutation psper;
		Permutation* genper;
		Permutation rightper;
		Permutation testper;
		Permutation bestper;

		//reco
		list<IDJet> sjets;
		vector<IDJet*> cleanedjets;
		list<IDMuon> smuons;
		vector<IDMuon*> loosemuons;
		vector<IDMuon*> tightmuons;
		list<IDElectron> selectrons;
		vector<IDElectron*> looseelectrons;
		vector<IDElectron*> mediumelectrons;
		IDMet met;

		TLorentzVector* lep = nullptr;
		int leppdgid = 0;
		//hists
		TH1DCollection gen1d;
		TH2DCollection gen2d;
		TH1DCollection reco1d;
		TH2DCollection reco2d;
		TH1DCollection truth1d;
		TH2DCollection truth2d;

        TTBarGenPlots ttp_genall;
        TTBarGenPlots ttp_genacc;
        TTBarGenPlots ttp_genjet;

		TTBarPlots ttp_truth;
		TTBarPlots ttp_right;
		TTBarPlots ttp_wrong;
		TTBarPlots ttp_semi;
		TTBarPlots ttp_other;
		TTBarPlots ttp_all;

	//	TTBarPlots ttp_jetspos_right;
	//	TTBarPlots ttp_jetspos_wrong;
	//	TTBarPlots ttp_hadjets_right;
	//	TTBarPlots ttp_hadjets_wrong;
	//	TTBarPlots ttp_jets_right;
	//	TTBarPlots ttp_jets_wrong;
	//	TTBarPlots ttp_blep_right;
	//	TTBarPlots ttp_blep_wrong;
		//TTBarPlots ttp_whad_right;
		//TTBarPlots ttp_whad_wrong;

		TTBarPlots ttp_tlepthad_right;
		TTBarPlots ttp_tlep_right;
		TTBarPlots ttp_thad_right;
		TTBarPlots ttp_nn_right;
		TTBarPlots ttp_nsemi_right;
		TTBarPlots ttp_alljets_right;
		TTBarPlots ttp_alljets_wrong;

		TTBarResponse response;
		TTBarResponse2D response2d;
		TTBarResponse2Dvar response2dvar;
		TTBarResponse response_ps;

		BtagEff btageff;
		BTagWeight btagweight;

		//ttbar solver
		TTBarSolver ttsolver;

		JetScaler jetscaler;

		//configuration
		bool PDFTEST;
		bool PSEUDOTOP;
		bool LHCPS;
		bool BTAGMODE;
		bool JETSCALEMODE;
		bool MUONS;
		bool ELECTRONS;
		double B_TIGHT;
		double B_MEDIUM;
		double B_LOOSE;
		string cnbtag;
		size_t cnusedjets;
		double clikelihoodcut;
		double cwjetptsoft;
		double cwjetpthard;
		double cbjetptsoft;
		double cbjetpthard;
		double cjetetamax;
		double clptmin;
		double cletamax;
		double cpwjetptsoft;
		double cpwjetpthard;
		double cpbjetptsoft;
		double cpbjetpthard;
		double cpjetetamax;
		double cplptmin;
		double cpletamax;
		double cpjetsep;
		double csigmajet;
		double csigmajetwj;
		double cscalejetwj;
		string cjecuncertainty;
		string cBTaggingSF;
		string cBTaggingEff;
		int cjetres;
		double csigmamet;
		double csigmalep;
		double ctopptweight = 0.;
		double ctoprapweight = 0.;
		double cttptweight = 0.;
		double cjetptweight = 0.;
		double cbdecay = 0.;
		double cbfrag = 0.;
		double cbsplitting = 1.;
		double ccsplitting = 1.;
		int cfacscale;
		int crenscale;
		int chdamp;
		int cbtagunc;
		int cltagunc;
		int cpileup;
		string cLeptonScaleFactor;
		string cJetEnergyUncertainty;
		string cJetResolution;
		string cJetResolutionSF;
		vector<double> cel27eff;
		bool TTMC;
		bool HERWIGPP;
		bool PYTHIA6;
		bool ISRUP;
		bool ISRDOWN;
		bool FSRUP;
		bool FSRDOWN;
		bool TUNEUP;
		bool TUNEDOWN;
		//
		double jetptmin;
	
		int lastlumi = -1;
		double weight;
		double mcweight;
		double puweight;

		BHadronDecayWeights bdecayweights;
		BFragWeights bfragweights;

		//binning vectors
		vector<double> topptbins;
		vector<double> topybins;
		vector<double> ttmbins;
		vector<double> ttybins;
		vector<double> ttptbins;
		vector<double> metbins;
		vector<double> htbins;
		vector<double> evtmassbins;
		vector<double> jetbins;
		vector<double> nobins;
		vector<double> dybins;
		vector<double> btagpt;
		vector<double> alljetbins;
		vector<vector<double> > jetptbins;
		vector<vector<double> > jetetabins;
		vector<vector<double> > jetdrbins;
		vector<vector<double> > jetdrtopbins;
		vector<string> jetnames = {"b_{l}", "b_{h}", "j_{W1}" , "j_{W2}", "j_{1}", "j_{2}","j_{3}", "j_{4}"};

		TH1D* puhist;
		TH2D* musfhist;
		TH2D* elsfhist;
		TH2D* mutrgsfhist;
		TH2D* eltrgsfhist;

		bool STUDENT = false;
		TFile* stud_tf = nullptr;
		TTree* stud_tr = nullptr;
		Float_t gen_px[2];
		Float_t gen_py[2];
		Float_t gen_pz[2];
		Float_t gen_e[2];
		Int_t gen_type[2];
		Float_t det_px[2];
		Float_t det_py[2];
		Float_t det_pz[2];
		Float_t det_e[2];
		Int_t det_type[2];
		UInt_t num_gen = 2;
		UInt_t num_det = 2;

	public:

		ttbar(const std::string output_filename);
		~ttbar();

		//This method is called once per job at the beginning of the analysis
		//book here your histograms/tree and run every initialization needed
		virtual void begin();
		//virtual void end();
		virtual void analyze();

		void SelectGenParticles(URStreamer& event);
		void SelectRivetPS(URStreamer& event);
		void SelectPseudoTopLHC(URStreamer& event);
		void SelectPseudoTop(URStreamer& event);
		void SelectRecoParticles(URStreamer& event);
		void AddGenJetSelection(URStreamer& event);
		void ttanalysis(URStreamer& event);
		void reconstruction();
		double lepeffweight(TLorentzVector* lep, URStreamer& event);

		static void setOptions() {}
};

#endif
