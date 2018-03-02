#include "Dilepton.h"
#include <helper.h>
#include <ttbarxsec.h>
#include "TMath.h"

using namespace std;

void Dilepton::InitGen(TLorentzVector* top, TLorentzVector* tbar, TLorentzVector* lt, TLorentzVector* ltbar, TLorentzVector* bt, TLorentzVector* btbar, int pdgid_lt, int pdgid_ltbar, TLorentzVector* nut, TLorentzVector* nutbar)
{
	Reset();
	top_ = *top;
	tbar_ = *tbar;
	lt_ = lt;
	ltbar_ = ltbar;
	bt_ = bt;
	btbar_ = btbar;
	nut_ = nut;
	nutbar_ = nutbar;
	nunu_ = *nut_ + *nutbar_;
	pdgid_lt_ = pdgid_lt;
	pdgid_ltbar_ = pdgid_ltbar;
	this->Calculate();
	isgen_ = true;
	solved_ = true;
}

void Dilepton::Init(TLorentzVector* lt, TLorentzVector* ltbar, TLorentzVector* bt, TLorentzVector* btbar, int pdgid_lt, int pdgid_ltbar, TLorentzVector* met)
{
	Reset();
	lt_ = lt;
	ltbar_ = ltbar;
	bt_ = bt;
	btbar_ = btbar;
	pdgid_lt_ = pdgid_lt;
	pdgid_ltbar_ = pdgid_ltbar;
	met_ = met;
}

void Dilepton::Copy(Dilepton* di)
{
	Reset();
	lt_ = di->Lt();
	ltbar_ = di->Ltbar();
	pdgid_lt_ = di->IdLt();
	pdgid_ltbar_ = di->IdLtbar();
	bt_ = di->Bt();
	btbar_ = di->Btbar();
	met_ = di->MET();
}


void Dilepton::Calculate()
{
	nunu_=*nut_+*nutbar_;
	wt_=*lt_+*nut_;
	wtbar_=*ltbar_+*nutbar_;

	if (!isgen_) {
		top_ = wt_ + *bt_;
		tbar_ = wtbar_ + *btbar_;
	}

	tt_=top_+tbar_;

	if (isgen_) {
		genmet_ = nunu_;
		genmet_.SetPz(0);
		genmet_.SetE(genmet_.Pt());
		//genmet_.SetE(sqrt(pow(genmet_.E(),2)-pow(genmet_.Z(),2)));
	}

	return;
}

void Dilepton::Reset()
{
	lt_ = nullptr;
	ltbar_ = nullptr;
	bt_ = nullptr;
	btbar_ = nullptr;
	nut_ = nullptr;
	nutbar_ = nullptr;
	top_ = TLorentzVector(0,0,0,0);
	tbar_ = TLorentzVector(0,0,0,0);
	addjets_.clear();
	pdgid_lt_ = 0;
	pdgid_ltbar_ = 0;
	calculated=false;


	minlist_.clear();
	nut_solns_.clear();
	nutbar_solns_.clear();
	solved_ = false;
}

void Dilepton::Solve()
{
	nut_solns_.clear();
	nutbar_solns_.clear();
	DNS_.Init(this);

	DNS_.Solve();

	nut_solns_ = DNS_.NutSolns();
	nutbar_solns_ = DNS_.NutbarSolns();
	solved_ =true;
	nmins_=DNS_.MinNum();
}

double Dilepton::Mtt()
{
	if (isgen_)
	{
		return (this->TT()).M();
	}
	TLorentzVector total = *nut_+*nutbar_+*bt_+*btbar_+*lt_+*ltbar_;
	for (int i = 0; i < addjets_.size(); ++i)
	{
		total += *(addjets_[i]);
	}
	return total.M();
}

double Dilepton::DY()
{
	TLorentzVector total1 = *nut_+*bt_+*lt_;
	TLorentzVector total2 = *nutbar_+*btbar_+*ltbar_;

	return total1.Rapidity() - total2.Rapidity();
}



double Dilepton::MinDeltaP(Dilepton gendi)
{

		vector<double> mindists;
		if (this->MinNum()<=0)
		{
			return 100000.;
		}
		for (int i=0; i<nmins_; ++i)
		{
			mindists.push_back(deltaP(*gendi.Nut(),nut_solns_[i]) + deltaP(*gendi.Nutbar(),nutbar_solns_[i]));
		}

		double bestmindist =mindists[0];
		int bestminindex = 0;
		for(int i=1; i<nmins_; ++i)
		{
			if (mindists[i]<bestmindist)
			{
				bestmindist = mindists[i];
				bestminindex = i;
			}
		}
		return mindists[bestminindex];


}



void Dilepton::NumericalNuSolve(int nbins)
		{
			if (solved_ == false)
			{
				NeutrinoSolver ntguess;
				ntguess.Build(lt_,bt_,MW,Mt);
				NeutrinoSolver ntbarguess;
				ntbarguess.Build(ltbar_,btbar_,MW,Mt);
				nbins_ = nbins;
				tstep_ = twopi/nbins_;
				TLorentzVector n1;
				TLorentzVector n2;
				TLorentzVector pn;
				TLorentzVector* met = this->MET();
				double dx;
				double dy;
				double dist;
				vector<vector<double>> minbins={10,{100000,20,20}};
				vector<vector<double>> minlist;

				for (double t1 = 0.; t1 < twopi; t1=t1+tstep_)
				{
					for (double t2 = 0.; t2 < twopi; t2=t2+tstep_)
					{
						n1 = ntguess.GetSolution(t1);
						n2 = ntbarguess.GetSolution(t2);
						pn = n1 + n2 - *met;
						TLorentzVector diff = pn;
						dx = pow(diff[0],2);
						dy = pow(diff[1],2);


						dist = dx+dy;

						if (dist<minbins[minbins_size_-1][0])
						{
							int k=minbins_size_-1;
							while(dist<minbins[k-1][0])
							{
								if(k==1){k=0;break;}
								else {k=k-1;}
							}
							for (int i = minbins_size_-1; i>k; i=i-1)
							{
								minbins[i]=minbins[i-1];
							}
							minbins[k]={dist,t1,t2};
						}
					}
				}
				//check for error
				if (DNS_.Error())
				{
					solved_ = true;
					nmins_ = 0;
					return;
				}

				for (size_t i = 0; i < minbins.size(); ++i)
				{
					bool newmin=true;
					vector<double> currentmin=minbins[i];
					for (vector<double> min : minlist)
					{
						if ((abs(min[1]-currentmin[1]) < minbins_size_*tstep_ || abs(abs(min[1]-currentmin[1]) -twopi) < minbins_size_*tstep_)
							&& (abs(min[2]-currentmin[2]) < minbins_size_*tstep_|| abs(abs(min[2]-currentmin[2]) -twopi) < minbins_size_*tstep_))
						{
							newmin = false;
						}
					}
					if(newmin)
					{
						minlist.push_back(currentmin);
						nut_solns_.push_back(ntguess.GetSolution(currentmin[1]));
						nutbar_solns_.push_back(ntbarguess.GetSolution(currentmin[2]));
					}
				}
				minlist_ = minlist;
				nmins_ = minlist.size();
				solved_ = true;


			}
		}

void Dilepton::NuPlot(TH2DCollection nugraphs, string s, int i, int steps)
	{
		int neutrino_bins=steps;
		double twopi = 6.2832;
		double neutrino_step = twopi/neutrino_bins;
		//nugraphs.AddNeutrinoHist("neutrino_solver"+to_string(i)+s,neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,"t1","t2");

		NeutrinoSolver ntguess;
		ntguess.Build(this->Lt(),this->Bt(),MW,Mt);
		NeutrinoSolver ntbarguess;
		ntbarguess.Build(this->Ltbar(),this->Btbar(),MW,Mt);
		TLorentzVector met = *this->MET();
		TLorentzVector n1;
		TLorentzVector n2;
		TLorentzVector pn;
		double dx;
		double dy;
		double pdist;


		for (double t1 = 0.; t1 < 6.283; t1=t1+neutrino_step)
		{
			for (double t2 = 0.; t2 < 6.283; t2=t2+neutrino_step)
			{
				n1 = ntguess.GetSolution(t1);
				n2 = ntbarguess.GetSolution(t2);
				pn = n1 + n2;
				dx = pow((pn[0]-met.Px()),2) ;
				dy = pow((pn[1]-met.Py()),2) ;
				pdist = sqrt(dx+dy);

				nugraphs["neutrino_solver"+to_string(i)+"_"+s]->Fill(t1,t2,pdist);

			}
		}
		//print some info
				cout<<"Neutrino Graph "<<i<<"_"<<s<<endl;
				for (int i = 0; i < nmins_; ++i)
				{
					cout<<minlist_[i][1]<<", "<<minlist_[i][2]<<endl;
				}
				cout<<endl;
	}


vector<int> Dilepton::SolnsSorted(Dilepton gendi)
{
	std::vector<int> inds;
	for (int i = 0; i < nmins_; ++i)
	{
		inds.push_back(i);
	}

	//have to do this for some reason
	const TLorentzVector* correct_nut = gendi.Nut();
	const TLorentzVector* correct_nutbar = gendi.Nutbar();

	std::sort(inds.begin(), inds.end(),
	[=](int& i1,  int& i2) {
	return deltaP(nut_solns_[i1] , *correct_nut) + deltaP(nutbar_solns_[i1] , *correct_nutbar)
		< deltaP(nut_solns_[i2] , *correct_nut) + deltaP(nutbar_solns_[i2] , *correct_nutbar);});



	return inds;
}

void Dilepton::SetNeutrinos(int i)
{
	nut_= &nut_solns_[i];
	nutbar_= &nutbar_solns_[i];
	Calculate();
}

void Dilepton::SetNeutrinos(TLorentzVector nut, TLorentzVector nutbar)
{
	nut_solns_.clear();
	nutbar_solns_.clear();
	nut_solns_.push_back(nut);
	nutbar_solns_.push_back(nutbar);
	nut_= &nut_solns_[0];
	nutbar_= &nutbar_solns_[0];

	Calculate();
	solved_ = true;
	nmins_ = 1;

}

void Dilepton::AverageNeutrinos()
{
	size_t s = nut_solns_.size();
	TLorentzVector avenut = TLorentzVector(0,0,0,0);
	TLorentzVector avenutbar = TLorentzVector(0,0,0,0);
	if (s<1) {
		return;
	}
	for (size_t i = 0; i < s; i++) {
		avenut = avenut + nut_solns_[i];
		avenutbar = avenutbar + nutbar_solns_[i];
	}
	avenut = 1/(double)s * avenut;
	avenutbar = 1/(double)s * avenutbar;
	this->SetNeutrinos(avenut,avenutbar);
}

void Dilepton::SortSolutions(Dilepton gendi)
{
	if (nut_solns_.size()<2) {
		return;
	}
	TLorentzVector* tru_nut = gendi.Nut();
	TLorentzVector* tru_nutbar = gendi.Nutbar();



	std::vector<size_t> inds;
	for (size_t i = 0; i < nut_solns_.size(); ++i)
	{
		inds.push_back(i);
	}

	//have to do this for some reason
	const TLorentzVector* correct_nut = gendi.Nut();
	const TLorentzVector* correct_nutbar = gendi.Nutbar();

	std::sort(inds.begin(), inds.end(),
	[=](size_t& i1,  size_t& i2) {
	return deltaP(nut_solns_[i1],  *correct_nut) + deltaP(nutbar_solns_[i1] , *correct_nutbar)
		< deltaP(nut_solns_[i2],  *correct_nut) + deltaP(nutbar_solns_[i2] , *correct_nutbar);});

	vector<TLorentzVector> sorted_nuts;
	vector<TLorentzVector> sorted_nutbars;
	for (size_t i =0; i<nut_solns_.size() ; i++) {
		sorted_nuts.push_back(nut_solns_[inds[i]]);
		sorted_nutbars.push_back(nutbar_solns_[inds[i]]);
		nut_solns_ = sorted_nuts;
		nutbar_solns_ = sorted_nutbars;
	}
}

void Dilepton::SortSolutionsPt()
{
	std::vector<size_t> inds;
	for (size_t i = 0; i < nmins_; ++i)
	{
		inds.push_back(i);
	}


	std::sort(inds.begin(), inds.end(),
	[=](size_t& i1,  size_t& i2) {
	return nut_solns_[i1].Pt()  +nutbar_solns_[i1].Pt()
		< nut_solns_[i2].Pt()  +nutbar_solns_[i2].Pt();});

	vector<TLorentzVector> sorted_nuts;
	vector<TLorentzVector> sorted_nutbars;
	for (size_t i : inds) {
		sorted_nuts.push_back(nut_solns_[i]);
		sorted_nutbars.push_back(nutbar_solns_[i]);
		nut_solns_ = sorted_nuts;
		nutbar_solns_ = sorted_nutbars;
	}
}

int Dilepton::LowPIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (nut_solns_[i].P()<lowp) {
			lowp = nut_solns_[i].P();
			lowpi = i;
		}
	}
	return(lowpi);
}

int Dilepton::LowPsIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (nut_solns_[i].P()+ nutbar_solns_[i].P()<lowp) {
			lowp = nut_solns_[i].P() +nutbar_solns_[i].P();
			lowpi = i;
		}
	}
	return(lowpi);
}

int Dilepton::LowPsIndex2()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (pow(nut_solns_[i].P(),2)+ pow(nutbar_solns_[i].P(),2)<lowp) {
			lowp = pow(nut_solns_[i].P(),2)+ pow(nutbar_solns_[i].P(),2);
			lowpi = i;
		}
	}
	return(lowpi);
}

int Dilepton::LowPmaxIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (max(nut_solns_[i].P(),nutbar_solns_[i].P()) < lowp) {
			lowp = max(nut_solns_[i].P(),nutbar_solns_[i].P());
			lowpi = i;
		}
	}
	return(lowpi);
}



int Dilepton::LowPtIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (nut_solns_[i].Pt()<lowp) {
			lowp = nut_solns_[i].Pt();
			lowpi = i;
		}
	}
	return(lowpi);
}

int Dilepton::LowPzIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		if (nut_solns_[i].Pz()<lowp) {
			lowp = nut_solns_[i].Pz();
			lowpi = i;
		}
	}
	return(lowpi);
}

int Dilepton::LowMttIndex()
{
	int lowpi = -1;
	double lowp = 1000000000000;
	for (size_t i = 0; i < nut_solns_.size(); i++) {
		this->SetNeutrinos(i);
		if (this->Mtt()<lowp) {
			lowp = this->Mtt();
			lowpi = i;
		}
	}
	return(lowpi);
}
