#ifndef HNSOLVER
#define HNSOLVER
#include <TMatrixD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;
using namespace TMath;

class NeutrinoSolver
{
	private:
		double Mt;
		double Mw;
		double Ml;
		double Mb;
		double Mn;

		bool ERROR = false;
		TMatrixD H_ =TMatrixD(3,3);
		TMatrixD Hp_ =TMatrixD(3,3);
		TMatrixD T =TMatrixD(3,1);
		TMatrixD MET =TMatrixD(2,1);
		TMatrixD VM =TMatrixD(2,2);
		TMatrixD RotationX(double a);
		TMatrixD RotationY(double a);
		TMatrixD RotationZ(double a);

		double IM_;

		void Solve(double t);
		TMatrixD GetPtSolution(double t);
		double Chi2(double t);
		pair<double, double> Extrem(double t, bool MIN = true);

	public:
		TMatrixD H(){return H_;}
		TMatrixD Hp(){return Hp_;}
		double IM(){return IM_;}
		TLorentzVector GetSolution(double t);
		double Chi2(TLorentzVector L, double t);
		void SetMET(double metx, double mety, double metxerr, double metyerr, double metxyrho);
		void Build(const TLorentzVector* lep, const TLorentzVector* bjet, double MW, double MT);
		NeutrinoSolver(){}
		TLorentzVector GetBest(double metx, double mety, double metxerr, double metyerr, double metxyrho, double& test, bool INFO = false);
		bool Error() {return ERROR;}

};


#endif
