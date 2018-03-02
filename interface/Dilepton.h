#ifndef DILEPTON_H
#define DILEPTON_H
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <TLorentzVector.h>
#include <NeutrinoSolver.h>
#include <DiNeutrinoSolver.h>
#include <helper.h>
#include <IDJet.h>
using namespace std;

class ttbar;
class DiNeutrinoSolver;

class Dilepton
{

	private:
		TLorentzVector* lt_ = 0;
		TLorentzVector* ltbar_ = 0;
		TLorentzVector top_;
		TLorentzVector tbar_ ;
		TLorentzVector* bt_ = 0;
		TLorentzVector* btbar_ = 0;
		TLorentzVector* nut_ = 0;
		TLorentzVector* nutbar_ = 0;
		TLorentzVector wt_ ;
		TLorentzVector wtbar_;
		TLorentzVector tt_ ;
		int pdgid_lt_ = 0;
		int pdgid_ltbar_ = 0;
		TLorentzVector* met_;
		TLorentzVector genmet_;
		TLorentzVector nunu_;
		vector<IDJet*> addjets_;
		bool calculated ;
		bool isgen_ = false;
		double jet_catch_deltaR = 0.4;

		//class formerly known as DilepNu

		bool solved_ =false;
		int nbins_;
		vector<vector<double>> minlist_;
		int nmins_ =-2;
		int minbins_size_ =10;
		vector<double> min_deltar_;
		vector<double> min_netdeltar_;
		vector<double> min_deltaphi_;
		vector<double> min_netdeltaphi;
		vector<double> min_deltap_;
		vector<double> min_deltapxpy_;
		double tstep_ ;
		int min_deltar_ind_;
		int min_netdeltar_ind_;
		int min_deltaphi__ind_;
		int min_netdeltaphi_ind_;
		int min_deltap_ind_;
		double MW = 80.4;
		double Mt = 172.5;
		double twopi = 6.2832;
		vector<TLorentzVector> nut_solns_ ;
		vector<TLorentzVector> nutbar_solns_;
		DiNeutrinoSolver DNS_;


	public:

//inits, resets, key utils
		Dilepton(){};
		void Init(TLorentzVector* lt, TLorentzVector* ltbar, TLorentzVector* bt, TLorentzVector* btbar, int pdgid_lt, int pdgid_ltbar, TLorentzVector* met);
		void InitGen(TLorentzVector* top, TLorentzVector* tbar, TLorentzVector* lt, TLorentzVector* ltbar, TLorentzVector* bt, TLorentzVector* btbar, int pdgid_lt, int pdgid_ltbar, TLorentzVector* nut, TLorentzVector* nutbar);
		void Reset();
		void Solve();
		void Copy(Dilepton* di);
		bool Error() {return DNS_.Error();}


//get basic elements

		TLorentzVector* Lt() const {return(lt_);}
		TLorentzVector* Ltbar() const {return(ltbar_);}
		int IdLt() {return pdgid_lt_;}
		int IdLtbar() {return pdgid_ltbar_;}
		TLorentzVector* Bt() const {return(bt_);}
		TLorentzVector* Btbar() const {return(btbar_);}
		TLorentzVector* Nut() {return(nut_);}
		TLorentzVector* Nutbar() {return(nutbar_);}
		TLorentzVector& Wt() {return wt_;}
		TLorentzVector& Wtbar() {return wtbar_;}
		TLorentzVector Top() {return top_;}
		TLorentzVector Tbar() { return tbar_;}
		TLorentzVector& TT() { return tt_;}
		TLorentzVector *MET() {
			if (isgen_) {
				return(&genmet_);
			}
		  return (met_);
		}
		vector<IDJet*> AddJets() {return addjets_;}
		TLorentzVector& NuNu() {return(nunu_);}
		bool IsGen() {return isgen_;}
		TLorentzVector* NuMET() {return &genmet_;}


//Get top physics stuff
		void Calculate();
		double Mtt();
		double DY();
		double IM();


//set elements

		void SetMET(TLorentzVector* met)
		{
			met_ =met;
			return;
		}
		void SetNeutrinos(TLorentzVector nut, TLorentzVector nutbar);
		void SetBs(TLorentzVector* bt, TLorentzVector* btbar)
		{
			bt_ = bt;
			btbar_ = btbar;
			nmins_ = -2;
			solved_ = false;
		}
		void SetAddJets(vector<IDJet*> addjets) {addjets_ = addjets;}

//Neutrino Utils


		DiNeutrinoSolver DNS() {return DNS_;}
		vector<TLorentzVector> NutSolns() {return nut_solns_;}
		vector<TLorentzVector> NutbarSolns() {return nutbar_solns_;}
		bool Solved() {return solved_;}
		void SetNeutrinos(int i);
		void AverageNeutrinos();
		vector<int> SolnsSorted(Dilepton gendi);
		void SortSolutions(Dilepton gendi);
		void SortSolutionsPt();
		void NuPlot(TH2DCollection nugraphs, string s, int i, int steps);
		double MinDeltaP(Dilepton gendi);
		double SumDeltaP(int ind, Dilepton gendi) {return 0.;};
		int LowPIndex();
		int LowPtIndex();
		int LowPzIndex();
		int LowPsIndex();
		int LowPsIndex2();
		int LowPmaxIndex();
		int LowMttIndex();
		void NumericalNuSolve(int nbins) ;
		vector<vector<double>> MinList()
		{
			if (!solved_){return {{0}};}
			return minlist_;
		}
		int MinNum()
		{
			if (!solved_){return -1;}
			return nut_solns_.size();
		}





//b-jet matching tools

		int BJetsCorrect(const Dilepton& other) const
		{
			int bs=0;
			for(int j = 0 ; j < NumBJets() ; j++)
			{
				if(GetJet(j)->DeltaR(*other.GetJet(j)) < jet_catch_deltaR) {bs++;}
			}
			return(bs);
		}
		int BJetsSwapped(const Dilepton& other) const
		{
			int bs=0;
			for(int j = 0 ; j < NumBJets() ; j++)
			{
				if(GetJet(1-j)->DeltaR(*other.GetJet(j)) < jet_catch_deltaR) {bs++;}
			}
			return(bs);
		}
		bool BtCorrect(const Dilepton& other) const
		{
			if(GetJet(0)->DeltaR(*other.GetJet(0)) < jet_catch_deltaR) {return true;}
			return(false);
		}
		bool BtbarCorrect(const Dilepton& other) const
		{
			if(GetJet(1)->DeltaR(*other.GetJet(1)) < jet_catch_deltaR) {return true;}
			return(false);
		}
		bool BtSwapped(const Dilepton& other) const
		{
			if(GetJet(0)->DeltaR(*other.GetJet(1)) < jet_catch_deltaR) {return true;}
			return(false);
		}
		bool BtbarSwapped(const Dilepton& other) const
		{
			if(GetJet(1)->DeltaR(*other.GetJet(0)) < jet_catch_deltaR) {return true;}
			return(false);
		}

//leftovers from Permutation class

		bool IsComplete() const {return(lt_ != 0 && ltbar_ != 0 && bt_ != 0 && btbar_ != 0 && (met_->E() != 0 || nut_->E() + nutbar_->E() > -0.5)&& bt_ != btbar_);}
		int NumBJets() const {return((bt_ != 0 ? 1 : 0) + (btbar_ != 0 ? 1 : 0));}
		int IsJetIn(const TLorentzVector* jet) const
		{
			if(jet->DeltaR(*Bt()) < 0.2) return 0;
			if(jet->DeltaR(*Btbar()) < 0.2) return 1;
			for(size_t i = 0 ; i < addjets_.size() ; ++i)
			{
				if(jet->DeltaR(*addjets_[i]) < 0.2) return i+2;
			}
			return -1;
		}
		TLorentzVector* GetJet(size_t n) const
		{
			if(n == 0) return Bt();
			if(n == 1) return Btbar();
			if(n > 1 && n-2 < addjets_.size()) return addjets_[n-2];

			return nullptr;
		}



};

#endif
