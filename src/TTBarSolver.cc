#include "TTBarSolver.h"
#include "NeutrinoSolver.h"
#include "IDMet.h"
#include "IDJet.h"


TTBarSolver* TTBarSolver::TTBS = 0;

void myfuncln(Int_t& npar, Double_t* deriv, Double_t& val, Double_t* par, Int_t flag)
{
	val = TTBarSolver::TTBS->Test(par);
}

TTBarSolver::TTBarSolver() : minuit(9), probfile(0)
{
}

TTBarSolver::~TTBarSolver()
{
	if(probfile != 0) {probfile->Close();}
}

void TTBarSolver::Init(bool pseudo, string filename, bool usebtag, bool usens, bool usemass)
{
	USEBTAG = usebtag;
	USENS = usens;
	USEMASS = usemass;
	TDirectory* dir = gDirectory;
	probfile = new TFile(filename.c_str(), "READ");
	N_right = dynamic_cast<TH1D*>(probfile->Get("D_right"));
	WTmass_right = dynamic_cast<TH2D*>(probfile->Get("MtMW_right"));
	masscutoff = WTmass_right->GetMaximum()/1000.;
//	WTmass_right = dynamic_cast<TH2D*>(probfile->Get("TRUTH/truth_Wmasshad_tmasshad_right"));
//	WTmass_right->RebinX(5);
//	WTmass_right->RebinY(5);
//	WTmass_right->Scale(1./WTmass_right->Integral("width"));
	//BTag_right = dynamic_cast<TH1D*>(probfile->Get("TRUTH/truth_btag_true"));
	//BTag_right->Scale(1./BTag_right->Integral("width"));
	//BTag_wrong = dynamic_cast<TH1D*>(probfile->Get("TRUTH/truth_btag_wrong"));
	//BTag_wrong->Scale(1./BTag_wrong->Integral("width"));
//	N_right = dynamic_cast<TH1D*>(probfile->Get("TRUTH/truth_nschi_right"));
	//N_right->Scale(1./N_right->Integral("width"));
	//N_wrong = dynamic_cast<TH1D*>(probfile->Get("TRUTH/truth_nschi_wrong"));
	//N_wrong->Scale(1./N_wrong->Integral("width"));
	if(pseudo)
	{
		c_mt = 171.5;
		c_mw = 83.1;
		c_rt = 1.2376E-3;
		c_rw = 2.6318E-3;
		c_rwt = -9.761E-4;
	}
	else
	{
		c_mt = 168.218;
		c_mw = 81.0519;
		c_rt = 2.105E-3;
		c_rw = 5.444E-3;
		c_rwt = -2.1583E-3;
	}
	norm = -1.* Log(Sqrt(c_rw*c_rt - c_rwt*c_rwt)/Pi());
	masstestmax = norm + c_mw*c_mw*c_rw + 2.*c_mw*c_mt*c_rwt + c_mt*c_mt*c_rt;
	dir->cd();
}

void TTBarSolver::Solve(TLorentzVector* bhad, TLorentzVector* j1had, TLorentzVector* j2had, TLorentzVector* blep, TLorentzVector* llep, IDMet* met, bool kinfit)
{
	kinfit_ = kinfit;
	bhad_ = dynamic_cast<IDJet*>(bhad);
	j1had_ = dynamic_cast<IDJet*>(j1had);
	j2had_ = dynamic_cast<IDJet*>(j2had);
	blep_ = dynamic_cast<IDJet*>(blep);
	llep_ = llep;
	met_ = met;
	ubhad_ = 0.1;
	uj1had_ = 0.1;
	uj2had_ = 0.1;
	ublep_ = 0.1;
	ullep_ = 0.01;
	umetx_ = 0.2;//met->pxunctot();
	umety_ = 0.2;//met->pyunctot();
	//umety_ = 0.05*met->Py();//Sqrt(met_->pyUnc());
	//umetx_ = 25;//Sqrt(met_->pxUnc());
	//umety_ = 25;//Sqrt(met_->pyUnc());
	rhomet_ = 0.;//met_->pxpyUnc()/(umetx_*umety_);

	nschi = -1;
	res = 1.E10;
	nstest = 1.E10;
	masstest = 1.E10;
	btagtest = 0;
	//cout << bhad->csvIncl() << " " << blep->csvIncl() << " " << j1had->csvIncl() << " " << j2had->csvIncl() << endl;
	//btagtest = -1.*Log(BTag_right->Interpolate(bhad_->csvIncl())/BTag_wrong->Interpolate(bhad_->csvIncl()));
	//btagtest -= Log(BTag_right->Interpolate(blep_->csvIncl())/BTag_wrong->Interpolate(blep_->csvIncl()));
	//btagtest -= Log(BTag_wrong->Interpolate(j1had_->csvIncl())/BTag_right->Interpolate(j1had_->csvIncl()));
	//btagtest -= Log(BTag_wrong->Interpolate(j2had_->csvIncl())/BTag_right->Interpolate(j2had_->csvIncl()));

	//minuit.mnparm(9, "metz", neutrino.Pz(), 1., 0., 10000., flags);
	//minuit.mnparm(10, "wn mass", 80., 1., 10, 1000, flags);
	//if(!LEPMASS)
	//{
	//	minuit.FixParameter(0);
	//	minuit.FixParameter(1);
	//}
	//if(!kinfit_)
	//{
	//	minuit.FixParameter(2);
	//	minuit.FixParameter(3);
	//	minuit.FixParameter(4);
	//	minuit.FixParameter(5);
	//	minuit.FixParameter(6);
	//	minuit.FixParameter(7);
	//	minuit.FixParameter(8);
	//}
	//else
	//{
	//	Int_t err;
	//	Double_t tmp[1];
	//	tmp[0] =  2+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  3+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  4+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  5+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  6+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  7+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//	tmp[0] =  8+1; minuit.mnexcm( "REL", tmp,  1,  err );
	//}
	if(!kinfit_)
	{
		double par[9];
		par[0] = 172.5;
		par[1] = 80.4;
		par[2] = 1.;
		par[3] = 1.;
		par[4] = 1.;
		par[5] = 1.;
		par[6] = 1.;
		par[7] = 1.;
		par[8] = 1.;
		Test(par);
	}
	else
	{
		TTBS = this;
		minuit.SetFCN(myfuncln);
		minuit.SetPrintLevel(-1);
		Int_t flags = 0;
		minuit.mnparm(0, "top mass", 172.5, 1., 10, 1000, flags);
		minuit.mnparm(1, "w mass", 80.4, 1., 10, 1000, flags);
		minuit.mnparm(2, "ubhad", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(3, "uj1had", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(4, "uj2had", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(5, "ublep", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(6, "ullep", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(7, "umetx", 1, 0.01, 0.5, 1.5, flags);
		minuit.mnparm(8, "umety", 1, 0.01, 0.5, 1.5, flags);
		minuit.FixParameter(0);
		minuit.FixParameter(1);
		minuit.SetMaxIterations(500);
		minuit.Migrad();
	}
}

double TTBarSolver::Test(double* par)
{
	bhadT_ = TLorentzVector(bhad_->Px()*par[2], bhad_->Py()*par[2], bhad_->Pz()*par[2], bhad_->E()*par[2]);
	j1hadT_ = TLorentzVector(j1had_->Px()*par[3], j1had_->Py()*par[3], j1had_->Pz()*par[3], j1had_->E()*par[3]);
	j2hadT_ = TLorentzVector(j2had_->Px()*par[4], j2had_->Py()*par[4], j2had_->Pz()*par[4], j2had_->E()*par[4]);
	blepT_ = TLorentzVector(blep_->Px()*par[5], blep_->Py()*par[5], blep_->Pz()*par[5], blep_->E()*par[5]);
	llepT_ = TLorentzVector(llep_->Px()*par[6], llep_->Py()*par[6], llep_->Pz()*par[6], llep_->E()*par[6]);
	NeutrinoSolver NS;
	NS.Build(&llepT_, &blepT_, par[1], par[0]);
	//metT_ = TLorentzVector(NS.GetBest(met_->Px()*par[7], met_->Py()*par[8], umetx_, umety_, rhomet_, nschi));
	metT_ = TLorentzVector(NS.GetBest(met_->Px()*par[7], met_->Py()*par[8], 1., 1., 0., nschi));
	// //Evan trying to do something
	// if (metT_[0]==0)
	// {
	// 	cout<<"NSFAIL IM "<<NS.IM()<<endl;
	// }
	//cout << nschi << " NS " << metT_.E() << " " << metT_.Pt() << endl;
	if(nschi > 0. && nschi < 150.*150.)
	{
		nstest = -1.*Log(N_right->Interpolate(Sqrt(nschi)));
	}

	double mwhad = (j1hadT_ + j2hadT_).M();
	double mthad = (j1hadT_ + j2hadT_ + bhadT_).M();
	//cout << mwhad << " M " << mthad << endl;


	if(false)
	{
		masstest = norm + ((mwhad-c_mw)*(mwhad-c_mw)*c_rw + 2.*(mwhad-c_mw)*(mthad-c_mt)*c_rwt + (mthad-c_mt)*(mthad-c_mt)*c_rt);
		if(masstest - masstestmax > Log(0.001)) masstest = 1.E10;
	}
	else
	{
		if(mthad < 500. && mwhad < 500.)
		{
			double massdisval = WTmass_right->Interpolate(mthad, mwhad);
			if(massdisval > masscutoff) {masstest = -1.*Log(massdisval);}
		}
	}

	res = 0.;
	if(kinfit_)
	{
		//res += Power((par[0]-173)/(20), 2);
		//res += Power((par[1]-80)/(20), 2);
		res += Power((par[2]-1.)/ubhad_ , 2)/2.;
		res += Power((par[3]-1.)/uj1had_ , 2)/2.;
		res += Power((par[4]-1.)/uj2had_ , 2)/2.;
		res += Power((par[5]-1.)/ublep_ , 2)/2.;
		res += Power((par[6]-1.)/ullep_ , 2)/2.;
		res += Power((par[7]-1.)/umetx_ , 2)/2.;
		res += Power((par[8]-1.)/umety_ , 2)/2.;
	}
	if(USEMASS) {res += masstest;}
	if(USENS) {res += nstest;}
	if(USEBTAG) {res += btagtest;}
	//cout << res << " " << par[2]-1 << " " << par[3]-1 <<endl;
	return(res);
}
