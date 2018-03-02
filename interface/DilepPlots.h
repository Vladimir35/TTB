#ifndef DIPLOTS
#define DIPLOTS
#include <helper.h>
#include <string>
#include <iostream>
#include <TLorentzVector.h>
#include <Dilepton.h>

class ttbar;
class Dilepton;

class DilepPlots
{
	protected:
		string prefix_;
		TH1DCollection plot1d;
		TH2DCollection plot2d;
		ttbar* an;

	public:
		DilepPlots(string prefix);
		~DilepPlots();
		void Init(ttbar* analysis);
		void Fill(Dilepton& di, double weight);
		void FillNuInfo(Dilepton& di, Dilepton& gendi, double weight);
		void FillErrorInfo(Dilepton& di,Dilepton& dilo, Dilepton& dihi, Dilepton& gendi, double weight);
		double AveSolnDist(vector<TLorentzVector> nutsolnsA, vector<TLorentzVector> nutbarsolnsA, vector<TLorentzVector> nutsolnsB, vector<TLorentzVector> nutbarsolnsB);

};



#endif
