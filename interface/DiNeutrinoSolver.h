#ifndef DNSOLVER
#define DNSOLVER

#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <helper.h>
#include <TMatrixDSymEigen.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <NeutrinoSolver.h>
#include <DiNeutrinoSolver.h>
#include <iostream>
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Functor.h"

using namespace std;
using namespace TMath;

class Dilepton;

class DiNeutrinoSolver
{
	private:
		double MW = 80.4;
		double Mt = 172.5;
		vector<TLorentzVector> nut_solns_ ;
		vector<TLorentzVector> nutbar_solns_;
		bool error_ =false;
		bool solved_;
		int nmins_;
		double minEllipseDist_;

		TMatrixD U = TMatrixD(3,3);
		TMatrixD Gamma = TMatrixD(3,3);
		TMatrixD D = TMatrixD(3,3);

		double metx;
		double mety;

		bool problem = false;

		double Cofactor(TMatrixD M, int row, int col);
		TMatrixD MatCross(TVector3 L, TMatrixD M);
		int swapI(int i);
		TMatrixD SwapMat(TMatrixD G);
		TVector3 SwapLine(TVector3 L);
		vector<TVector3> Lfactor(TMatrixD M);
		vector<TVector3> IntersectLineEllipse(TVector3 line, TMatrixD ellipse);
		vector<TVector3> IntersectEllipses(TMatrixD ellipse1, TMatrixD ellipse2);

	public:

		static ROOT::Math::Minimizer* mini;

		double MinEllipseDist() {return minEllipseDist_;}
		int eps[3][3][3];
		void NumericalBestSoln(double guess1, double guess2);
		double DistFromExact(const double *tt);
		void buildeps();
		vector<TVector3> BurtsApprox(TMatrixD Xp, TMatrixD XpBar);
		vector<TLorentzVector> NutSolns() {return nut_solns_;}
		vector<TLorentzVector> NutbarSolns() {return nutbar_solns_;}
		void Init(Dilepton* di);
		DiNeutrinoSolver(){};
		NeutrinoSolver NuSolve;
		NeutrinoSolver NuBarSolve;

		TMatrixD GammaMat()
		{
			return Gamma;
		}
		TMatrixD Hpnu()
		{
			return NuSolve.Hp();
		}
		TMatrixD Hpnub()
		{
			return NuBarSolve.Hp();
		}


		void Solve();

		int MinNum()
		{
			if (!solved_)
			{
				return 0;
			}
			return nut_solns_.size();
		}
		void Reset()
		{
			nut_solns_.clear();
			nutbar_solns_.clear();
			solved_ = false;

		}
		TLorentzVector Nu(int i)
		{
			return nut_solns_[i];
		}
		bool Error()
		{
			return error_;
		}



};


#endif
