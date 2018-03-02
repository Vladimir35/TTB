#ifndef HELPER_H
#define HELPER_H
#include <string>
#include <sstream>
#include <map>
#include <vector>

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMatrixD.h"

class Dilepton;
using namespace std;

//Path and Filenames
string FNbody(string filename);
string FNfilename(string filename);
string FNpath(string filename);
string FNnoext(string filename);

vector<string> string_split(const string& in, const vector<string> splits = {" "});

template<typename T> T stringtotype(string s)
{
	T i;
	istringstream(s) >> i;
	return(i);
}

// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
				std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(),
				std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

//Stat
double ProbToSigma(double prob);
double SigmaToProb(double sigma);

//Histgramm Collections

class TH1DCollection
{
	private:
		string prefix;
		map<string, TH1D*> hmap;
	public:
		TH1DCollection(string _prefix);
		TH1D* AddHist(string name, Int_t bins, Double_t min, Double_t max, string xlabel, string ylabel);
		TH1D* AddHist(string name, const vector<Double_t>& bins, string xlabel, string ylabel);

		TH1D* Get(string name);
		TH1D* operator[](string name);
};

class TH2DCollection
{
	private:
		string prefix;
		map<string, TH2D*> hmap;
	public:
		TH2DCollection(string _prefix);
		TH2D* AddHist(string name, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax, string xlabel, string ylabel);
		TH2D* AddHist(string name, Int_t xbins, Double_t xmin, Double_t xmax, const vector<Double_t>& ybins, string xlabel, string ylabel);
		TH2D* AddHist(string name, const vector<Double_t>& xbins, Int_t ybins, Double_t ymin, Double_t ymax, string xlabel, string ylabel);
		TH2D* AddHist(string name, const vector<Double_t>& xbins, const vector<Double_t>& ybins, string xlabel, string ylabel);
		TH2D* AddNeutrinoHist(string name, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax, string xlabel, string ylabel);

		TH2D* Get(string name);
		TH2D* operator[](string name);
};

//Sort
double getThetaStar(TLorentzVector Z1, TLorentzVector Z2);



template<typename T> bool sortpt(const T* A, const T* B) {return A->Pt() >  B->Pt();};
//bool tlvsortpt(const TLorentzVector& JA, const TLorentzVector& JB);
//bool tlvsortetaabs(const TLorentzVector& JA, const TLorentzVector& JB);

//bool ptlvsortpt(const TLorentzVector* JA, const TLorentzVector* JB);
//bool ptlvsortetaabs(const TLorentzVector* JA, const TLorentzVector* JB);
//bool ptlvsorteta(const TLorentzVector* JA, const TLorentzVector* JB);


//BINNER
class Bin
{
	private:
		double min_;
		double max_;
	public:
		Bin(double min, double max) : min_(min), max_(max) {}
		Bin(double val) : min_(val), max_(val) {}

		double min() const {return min_;}
		double max() const {return max_;}

};

bool operator<(const Bin& A, const Bin& B);

//Evan's stuff

double AddInQ(double d1, double d2);
int DYbin(double dy);
void MMPrint(TMatrixD M);
void MMPrint(TVectorD V);
void MMPrint(TVector3 V);
void MMPrint(TLorentzVector V);
void NuDetailPrint(Dilepton di, string s);
TLorentzVector Smear(TLorentzVector b, double sig);
TLorentzVector ESmear(TLorentzVector b, double sig);
TLorentzVector METsmear(TLorentzVector b, double sig);
double deltaP(TLorentzVector Z1, TLorentzVector Z2);
double deltaPxy(TLorentzVector Z1, TLorentzVector Z2);
double deltaPz(TLorentzVector Z1, TLorentzVector Z2);

#endif
