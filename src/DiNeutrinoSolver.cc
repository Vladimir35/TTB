#include "DiNeutrinoSolver.h"
#include "Dilepton.h"
#include "Math/GSLMinimizer.h"

//define a matrix cofactor

double DiNeutrinoSolver::Cofactor(TMatrixD M, int row, int col)
{
return(M((row + 1) % 3 , (col + 1) % 3) * M((row + 2) % 3,(col + 2) % 3) - M((row + 1) % 3,(col + 2) % 3) * M((row + 2) % 3,(col + 1) % 3));
}

//define the levi-civita tensor

void DiNeutrinoSolver::buildeps()
{
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			for (size_t k = 0; k < 3; k++) {
				if (i==j || j== k || k==i)
				{
					eps[i][j][k]= 0;
				}
				else if (i==0)
				{
					if (j==1) {eps[i][j][k]= 1;}
					else {eps[i][j][k]= -1;}
				}
				else if (i==1)
				{
					if (j==2) {eps[i][j][k]= 1;}
					else {eps[i][j][k]= -1;}
				}
				else if (i==2)
				{
					if (j==0) {eps[i][j][k]= 1;}
					else {eps[i][j][k]= -1;}
				}
			}
		}
	}
}

//Define the "cross product" of a line with a matrix, used to find intersections of line with ellipse

TMatrixD DiNeutrinoSolver::MatCross(TVector3 L, TMatrixD M)
{
	TMatrixD D(3, 3);
	for (int i = 0; i < 3; ++i)
	{
		for (int m = 0; m < 3; ++m)
		{
			D(i,m)=0;
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					D(i,m)+=eps[i][j][k]*L(j)*M(k,m);
				}
			}
		}
	}
	return D;
}

//tool to swap first and second indices, which sometimes you need to do
int DiNeutrinoSolver::swapI(int i)
	{if (i==0) return 1;
	 if (i==1) return 0;
	 if (i==2) return 2;
	 else return 3;
	}

//returns matrix with swapped first and second indices
TMatrixD DiNeutrinoSolver::SwapMat(TMatrixD G)
{
		TMatrixD Gs(3,3);
		for (int x = 0; x<3; ++x)
		{
			for (int y = 0; y < 3; ++y)
			{
			Gs[x][y]=G[swapI(x)][swapI(y)];
			}
		}
		return Gs;
}

//returns line with swapped first and second indices

TVector3 DiNeutrinoSolver::SwapLine(TVector3 L)
{
		TVector3 Ls;
		for (int x = 0; x<3; ++x)
		{
			Ls(x)=L(swapI(x));
		}
		return Ls;
}

//factors degenrate conic into outer product of two lines. These lines can then be intersected with an ellipse to get the intersection of the conic and ellipse

vector<TVector3> DiNeutrinoSolver::Lfactor(TMatrixD M)
{
	TVector3 Lp;
	TVector3 Lm;
	vector<TVector3> Ls;
	double small=10e-10;
	bool swapQ = abs(M[0][0])>abs(M[1][1]);
	if (swapQ) {M = SwapMat(M);}
	M *= (1/M(1,1));
	if (abs(M[0][0])<small && abs(M[1][1])<small)
	{
		Lp=TVector3(M(0,1),0,M(1,2));
		Lm=TVector3(0,M(0,1),M(0,2)-M(1,2));
	}
	else if (abs(Cofactor(M,2,2))<small && abs(M[1][1])>small)
	{
		if (Cofactor(M,0,0)>0) {return Ls;}
		Lp=TVector3(M(0,1),M(1,1),M(1,2)+sqrt(-Cofactor(M,0,0)));
		Lm=TVector3(M(0,1),M(1,1),M(1,2)-sqrt(-Cofactor(M,0,0)));
	}
	else if (abs(M[1][1])>small)
	{
		if (Cofactor(M,2,2)>0) {return Ls;}
		Lp=TVector3(M(0,1)+sqrt(-Cofactor(M,2,2)),M(1,1),
		-1*(Cofactor(M,1,2)/Cofactor(M,2,2)+(Cofactor(M,0,2)/Cofactor(M,2,2))*(M(0,1)+sqrt(-Cofactor(M,2,2)))));
		Lm=TVector3(M(0,1)-sqrt(-Cofactor(M,2,2)),M(1,1),
		-1*(Cofactor(M,1,2)/Cofactor(M,2,2)+(Cofactor(M,0,2)/Cofactor(M,2,2))*(M(0,1)-sqrt(-Cofactor(M,2,2)))));
	}
	if (swapQ) {Ls.push_back(SwapLine(Lp));Ls.push_back(SwapLine(Lm));}
	else {Ls.push_back(Lp); Ls.push_back(Lm);}

	return Ls;
}

//the intersections of a line with ellipse is found by looking at eigenvectors of the 'matrix cross product'
//the parameter k has been somewhat problematic here. Burt says it should be zero for valid solutions, but sonetimes it is of order 1e-3 for valid solutions

vector<TVector3> DiNeutrinoSolver::IntersectLineEllipse(TVector3 L, TMatrixD M)
{
	if (problem) {
		cout<<endl<<"EL running: ";
	}
	vector<TVector3> list;
	vector<TVector3> emptylist;
	TMatrixD MT = M;
	MT.T();
	M = 0.5 *(M+MT);

	TMatrixDEigen Mx = MatCross(L, M);
	TVectorD evRe = Mx.GetEigenValuesRe();


	TVectorD evIm = Mx.GetEigenValuesIm();

	if (problem) {
		cout<<"###"<<endl;
		evRe.Print();
	}

	TMatrixD evecs = Mx.GetEigenVectors();
	for (int i = 0; i < 3; ++i)
	{
		if(abs(evIm(i))>1e-6 || abs(evRe(i))<1e-12){continue;}
		TVector3 v = TVector3(evecs(0,i),evecs(1,i),evecs(2,i));
		TMatrixDColumn vcol = TMatrixDColumn(evecs,i);
		double k = pow(L.Dot(v),2)+v.Dot(M*v);
		if (problem) {
			cout<<endl<<"k= "<<k;
		}
		if(abs(k)>1e-2){/*cout<<" (>>>k="<<k<<") ";*/continue;}
		else {/*cout<<"k="<<k <<endl;*/ list.push_back((1/v(2))*v);}
	}

	if (problem) {
		cout<<endl<<"N_EL: "<<list.size()<<" ";
	}
	return list;
}

//Mysterious approximation from Burt's paper
//Some fixes to implementation --the paper's sugguestions do not always work reliably

vector<TVector3> DiNeutrinoSolver::BurtsApprox(TMatrixD Xp, TMatrixD XpBar)
{

		TMatrixD Mp = TMatrixD(TMatrixD::kTransposed, Xp*D) + Xp*D;
		TMatrixD MpBar = TMatrixD(TMatrixD::kTransposed, XpBar*D) + XpBar*D;

		vector<TVector3>  Ts = IntersectEllipses(Mp,U);
		vector<TVector3>  TBars = IntersectEllipses(MpBar,U);

		if (Ts.size()==0 || TBars.size()==0 ) {
			problem = true;
			cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
			vector<TVector3>  probTs = IntersectEllipses(Mp,U);
			vector<TVector3>  probTBars = IntersectEllipses(MpBar,U);
			cout<<endl<<"Mp=";
			MMPrint(Mp);
			cout<<endl<<"Mpb=";
			MMPrint(MpBar);
			cout<<endl;
			cout<<probTs.size()<<", "<<probTBars.size();
			cout<<endl;
			cout<<Ts.size()<<", "<<TBars.size();
		}


		double bestD = 1e10;
		vector<size_t> bestIJ = {0,0};

		for (size_t i = 0; i < Ts.size(); i++) {
			for (size_t j = 0; j < TBars.size(); j++) {
				double currentD =
					(NuSolve.Hp()*Ts[i] + NuBarSolve.Hp()*TBars[j] -TVector3(metx,mety,0) ).Mag2();
				if (currentD < bestD) {
					bestD = currentD;
					bestIJ = {i,j};
				}
			}
		}
		vector<TVector3> approxMinTs;

		if (Ts.size()==0 || TBars.size()==0 ) {
			approxMinTs = {TVector3(1,0,0),TVector3(0,1,0)};
		}
		else{
			approxMinTs = {Ts[bestIJ[0]],TBars[bestIJ[1]]};
		}

		return approxMinTs;

}

//Distance (squared) from exact solution. to be minimized in NumericalBestSoln()

double DiNeutrinoSolver::DistFromExact(const double *tt)
{
	const double_t t1 = tt[0];
	const double_t t2 = tt[1];
	const TVector3 nu1 = NuSolve.Hp()*TVector3(cos(t1),sin(t1),1);
	const TVector3 nu2 = NuBarSolve.Hp()*TVector3(cos(t2),sin(t2),1)-TVector3(metx,mety,0);
	const TVector3 distvec = nu1+nu2;
	return distvec.Mag2();
}

//Numerical solver part of the approximate solution method

void DiNeutrinoSolver::NumericalBestSoln(double guess1, double guess2)
{


   ROOT::Math::Functor distFunc(this,&DiNeutrinoSolver::DistFromExact,2);


   mini->SetMaxFunctionCalls(1000000);
   mini->SetMaxIterations(1000000);
   mini->SetTolerance(0.01);

   double step[2] = {0.01,0.01};
   double variable[2] = {guess1,guess2};

   mini->SetFunction(distFunc);

   mini->SetVariable(0,"t1",variable[0], step[0]);
   mini->SetVariable(1,"t2",variable[1], step[1]);
	mini->SetVariableLimits(0,-6.4,12.8);
	mini->SetVariableLimits(1,-6.4,12.8);

   mini->Minimize();


	if(mini->Edm()<10)
	{
		minEllipseDist_ = mini->MinValue();
		const double *ts = mini->X();
		double t1f = ts[0];
		double t2f = ts[1];

		TVector3 nu = NuSolve.H()*TVector3(cos(t1f),sin(t1f),1);
		TVector3 nubar = NuBarSolve.H()*TVector3(cos(t2f),sin(t2f),1);
		nut_solns_.push_back(TLorentzVector(nu,nu.Mag()));
		nutbar_solns_.push_back(TLorentzVector(nubar,nubar.Mag()));
	}
	else
	{
		TVector3 gv1 = NuSolve.H()*TVector3(cos(guess1),sin(guess1),1);
		TVector3 gv2 = NuBarSolve.H()*TVector3(cos(guess2),sin(guess2),1);
		nut_solns_.push_back(TLorentzVector(gv1,gv1.Mag()));
		nutbar_solns_.push_back(TLorentzVector(gv2,gv2.Mag()));
	}
	return;
}

//Form the Pencil, pick a value of lambda for which it is degenerate, factor, intersect line with ellipse.

vector<TVector3> DiNeutrinoSolver::IntersectEllipses(TMatrixD A, TMatrixD B)
{
	TMatrixD AT = A;
	TMatrixD BT = B;
	AT.T();
	BT.T();
	A=0.5*(A+AT);
	B=0.5*(B+BT);
	if (abs(A.Determinant())>abs(B.Determinant()))
	{
		TMatrixD M = A;
		A = B;
		B = M;
	}
	TMatrixD Ai = A;
	Ai.Invert();

	TMatrixDEigen Mat = Ai*B;
	TVectorD evRe = Mat.GetEigenValuesRe();
	TVectorD evIm = Mat.GetEigenValuesIm();

	vector<TVector3> pointlist;

	for(int evindex=0 ; evindex<3; evindex++)
	{
		double minsize	= .01;

		if(abs(evIm(evindex))>10e-10 ){continue;}
		TMatrixD DegMat = B - evRe(evindex)*A;
		vector<TVector3> lines = Lfactor(DegMat);
		for(size_t i=0;i<lines.size();i++)
		{
			vector<TVector3> points = IntersectLineEllipse(lines[i],A);
			for(size_t j=0; j<points.size(); j++)
			{
				bool newmin =true;
				for(size_t k=0; k<pointlist.size(); k++)
				{
					if ((points[j]-pointlist[k]).Mag()<minsize){newmin=false; continue; }
				}
				if (newmin) {pointlist.push_back(TVector3(points[j]));}
			}
		}
	}
	if (problem) {
		cout<<endl<<"N_EE= "<<pointlist.size()<<" ";
	}
	return pointlist;
}

//Given a Dilepton object, forms two NeutrinoSolver objects and gets the solution

	void DiNeutrinoSolver::Solve()
	{
		//there has got to be a better way to set TMatrixD elements, but when I fed it an array of the elements in the TMatrix constructor, it didn't want to fix the size of the matrix correctly and kept giving matrix multiplication errors!
		nut_solns_.clear();
		nutbar_solns_.clear();

		D(0, 0) = 0.;
		D(1, 0) = -1.;
		D(2, 0) = 0.;
		D(0, 1) = 1.;
		D(1, 1) = 0.;
		D(2, 1) = 0.;
		D(0, 2) = 0.;
		D(1, 2) = 0.;
		D(2, 2) = 0.;

		U(0, 0) = 1.;
		U(1, 0) = 0.;
		U(2, 0) = 0.;
		U(0, 1) = 0.;
		U(1, 1) = 1.;
		U(2, 1) = 0.;
		U(0, 2) = 0.;
		U(1, 2) = 0.;
		U(2, 2) = -1.;


		Gamma(0, 0) = -1.;
		Gamma(1, 0) = 0.;
		Gamma(2, 0) = 0;
		Gamma(0, 1) = 0.;
		Gamma(1, 1) = -1.;
		Gamma(2, 1) = 0;
		Gamma(0, 2) = metx;
		Gamma(1, 2) = mety;
		Gamma(2, 2) = 1.;

		TMatrixD GammaT = Gamma;
		GammaT.T();

		if (NuSolve.Error() || NuBarSolve.Error() ) {
			error_ = true;
			solved_= true;
			return;
		}


		TMatrixD HPerp = NuSolve.Hp();

		TMatrixD HPerpT(TMatrixD::kTransposed, HPerp);
		TMatrixD HPerpI(TMatrixD::kInverted, HPerp);
		TMatrixD HPerpIT(TMatrixD::kTransposed, HPerpI);

		TMatrixD HBarPerp = NuBarSolve.Hp();

		TMatrixD HBarPerpT(TMatrixD::kTransposed, HBarPerp);
		TMatrixD HBarPerpI(TMatrixD::kInverted, HBarPerp);
		TMatrixD HBarPerpIT(TMatrixD::kTransposed, HBarPerpI);

		TMatrixD NuPerp = (HPerpIT*U)*HPerpI;
		TMatrixD NuBarPerp = HBarPerpIT*U*HBarPerpI;
		TMatrixD NuBarPerpPrime = GammaT*NuBarPerp*Gamma;


		vector<TVector3>  solns2d  = IntersectEllipses(NuPerp,NuBarPerpPrime);

		if (solns2d.size()==0) {
			TMatrixD Xp1 = HPerpT*NuBarPerpPrime*HPerp;
			TMatrixD NuPerpPrime = GammaT*NuPerp*Gamma;
			TMatrixD Xp2 = HBarPerpT*NuPerpPrime*HBarPerp;
			vector<TVector3> guesses = BurtsApprox(Xp1,Xp2);
			double t1guess = atan2(guesses[0](1),guesses[0](0));
			double t2guess = atan2(guesses[1](1),guesses[1](0));
			NumericalBestSoln(t1guess,t2guess);
		}
		else for (size_t i = 0; i < solns2d.size(); ++i)
		{
			TVector3 nu_p = solns2d[i];
			TVector3 nubar_p = Gamma*nu_p;
			TVector3 nu = TVector3((NuSolve.H()*HPerpI)*nu_p);
			TVector3 nubar = NuBarSolve.H()*HBarPerpI*nubar_p;
			nut_solns_.push_back(TLorentzVector(nu,nu.Mag()));
			nutbar_solns_.push_back(TLorentzVector(nubar,nubar.Mag()));
			minEllipseDist_ = 0.;
		}

		solved_ = true;
	}

//init.

void DiNeutrinoSolver::Init(Dilepton* di)
{
	Reset();
	NuSolve.Build(di->Lt(),di->Bt(),MW,Mt);
	NuBarSolve.Build(di->Ltbar(),di->Btbar(),MW,Mt);
	error_=false;
	solved_=false;
	buildeps();
	TLorentzVector* met = di->MET();
	metx = met->X();
	mety = met->Y();
}

ROOT::Math::Minimizer* DiNeutrinoSolver::mini = 0;
