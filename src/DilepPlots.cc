#include <DilepPlots.h>
#include <ttbarxsec.h>
#include <helper.h>
#include <Permutation.h>


using namespace std;
using namespace TMath;

DilepPlots::DilepPlots(string prefix) : prefix_(prefix), plot1d(prefix), plot2d(prefix)
{

}

DilepPlots::~DilepPlots()
{

}

void DilepPlots::Init(ttbar* analysis)
{
	an = analysis;
	plot1d.AddHist("b&lep_inv_mass",100,0,500.,"b + lep invariant mass","Events");
	plot1d.AddHist("b&lep_sum_deltar",25,0,10.,"bt.deltaR(lept) + btbar.deltaR(leptbar)","Events");
//	plot1d.AddHist("bestmin_pxy",50,0,100.,"best min delta_pxy","Events");
	plot1d.AddHist("allmin_Pt",50,0,300,"nu Pt","Events");
	plot1d.AddHist("allmin_Pz",50,0,300,"nu Pt","Events");
	plot1d.AddHist("allmin_P",50,0,500,"nu Pt","Events");
	plot1d.AddHist("dns_min_num",6,-0.5,5.5,"number of mins","Events");

	plot1d.AddHist("dns_bestmin_deltaP_sum",30,0,500,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaP_each",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_bestmin_diffP_each",100,-300,300,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaP_each_percent",50,0,3.,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaP_nunu",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaPxy_sum",30,0,300,"delta Pxy","Events");
	plot1d.AddHist("dns_bestmin_deltaPxy_each",30,0,200,"delta Pxy","Events");
	plot1d.AddHist("dns_bestmin_deltaPxy_nunu",30,0,200,"delta Pxy","Events");
	plot1d.AddHist("dns_bestmin_deltaPz_sum",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaPz_each",30,0,200,"delta P","Events");
	plot1d.AddHist("dns_bestmin_deltaPz_nunu",30,0,200,"delta P","Events");
	plot1d.AddHist("dns_allmin_diffP_each",100,-300,300,"delta P","Events");

	plot1d.AddHist("dns_allmin_deltaP_sum",30,0,500,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaP_each",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaP_nunu",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPxy_sum",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPxy_each",30,0,200,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPxy_nunu",30,0,200,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPz_sum",30,0,300,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPz_each",30,0,200,"delta P","Events");
	plot1d.AddHist("dns_allmin_deltaPz_nunu",30,0,200,"delta P","Events");

	// plot1d.AddHist("dns_avemin_deltaP_sum",30,0,500,"delta P","Events");
	// plot1d.AddHist("dns_avemin_deltaP_each",30,0,300,"delta P","Events");
	// plot1d.AddHist("dns_avemin_deltaP_nunu",30,0,300,"delta P","Events");
	// plot1d.AddHist("dns_avemin_deltaPxy_sum",30,0,300,"delta Pxy","Events");
	// plot1d.AddHist("dns_avemin_deltaPxy_each",30,0,200,"delta Pxy","Events");
	// plot1d.AddHist("dns_avemin_deltaPxy_nunu",30,0,200,"delta Pxy","Events");
	// plot1d.AddHist("dns_avemin_deltaPz_sum",30,0,300,"delta Pz","Events");
	// plot1d.AddHist("dns_avemin_deltaPz_each",30,0,200,"delta Pz","Events");
	// plot1d.AddHist("dns_avemin_deltaPz_nunu",30,0,200,"delta Pz","Events");


	plot1d.AddHist("dns_nu_P_error_best",100,-500,500.,"IM","Events");
	plot1d.AddHist("dns_TT_IM_error_best",100,-500,500.,"IM","Events");
	plot1d.AddHist("dns_TT_DeltaY_error_best",50,0.,2.5,"delta Y","Events");

	plot1d.AddHist("dns_TT_IM_percenterror_best",50,-1.5,1.5,"IM percenterror","Events");
	plot1d.AddHist("dns_TT_DeltaY_percenterror_best",50,0.,1.,"delta Y percenterror","Events");

	plot1d.AddHist("dns_bestmin_momentum_percenterror",50,0,3.,"momentum percenterror","Events");

	plot2d.AddHist("dns_nu_vs_nubar_deltaP",100,0,300.,100,0,300.,"nu delta P","nu bar delta P");
	plot2d.AddHist("dns_deltaP_vs_Pt",100,0,300.,100,0,300.,"soln pt","nu soln delta P");
	plot2d.AddHist("dns_deltaP_vs_Pz",100,0,500.,100,0,300.,"soln pz","nu soln delta P");
	plot2d.AddHist("dns_deltaP_vs_nummins",100,0,500.,15,-.5,4.5,"nu soln delta P (each)","min num");

	plot1d.AddHist("W_mass_approx",50,60.,100.,"mass","events");
	plot1d.AddHist("Top_mass_approx",50,140.,210.,"mass","events");
	plot1d.AddHist("W_mass_exact",50,70.,90.,"mass","events");
	plot1d.AddHist("Top_mass_exact",50,150.,200.,"mass","events");
	plot1d.AddHist("W_mass_vsgen_approx",50,-50,50.,"mass","events");
	plot1d.AddHist("Top_mass_vsgen_approx",50,-50,50.,"mass","events");
	plot1d.AddHist("W_mass_vsgen_exact",50,-50,50,"mass","events");
	plot1d.AddHist("Top_mass_vsgen_exact",50,-50.,50.,"mass","events");

	plot1d.AddHist("bjerr_change_in_best_sumdeltaP",50,0.,500.,"change in delta P (GeV)","events");
	plot1d.AddHist("bjerr_change_in_best_sumdeltaP_hilo",50,0.,500.,"change in delta P (GeV)","events");
	plot1d.AddHist("bjerr_avedist",50,0.,500.,"change in delta P (GeV)","events");
	plot1d.AddHist("bjerr_avedist_lo",50,0.,500.,"change in delta P (GeV)","events");

	plot1d.AddHist("bjerr_avedist_hi",50,0.,500.,"change in delta P (GeV)","events");

	plot1d.AddHist("bjerr_avedist_hilo",50,0.,500.,"change in delta P (GeV)","events");
	plot1d.AddHist("bjerr_btest_hi",50,0.,100.,"change in P of bjet (GeV)","events");
	plot1d.AddHist("bjerr_btest_lo",50,0.,100.,"change in P of bjet (GeV)","events");
	plot1d.AddHist("bjerr_btest_hilo",50,0.,100.,"change in P of bjet (GeV)","events");



	plot2d.AddHist("bestmin_deltaP_vs_avedist_hilo",100,0,500.,100,0,500,"hilo avedist","nu soln delta P (each)");
	plot2d.AddHist("avedist_hilo_vs_numminsA",100,0,500.,15,-.5,4.5,"nu soln delta P (each)","min num");
	plot2d.AddHist("avedist_hilo_vs_numminsB",100,0,500.,15,-.5,4.5,"nu soln delta P (each)","min num");

	// plot1d.AddHist("dns_worsemins_deltap_each",25,0,400,"delta P","Events");
	// plot1d.AddHist("dns_bestmin_deltaPxy",25,0,200,"number of mins","Events");
	// plot1d.AddHist("dns_bestmin_deltaR",25,0,10,"delta R","Events");
	// plot1d.AddHist("dns_deltaPxy_test",25,0,20,"deltaPxy (Net, sanity test)","Events");
	// plot2d.AddHist("deltaP_min_comparison_each",25,0,500,10,0.5,2.5,"deltaP","1=best, 2=other solns");
	// plot2d.AddHist("deltaRlep_vs_correct_min",25,0,9,10,0.5,2.5,"deltaR lep+Nu sum","1=best, 2=other solns");
	// plot2d.AddHist("deltaRbl_vs_correct_min",25,0,16,10,0.5,2.5,"deltaR lep+Nu sum","1=best, 2=other solns");
	// plot2d.AddHist("sum_Pt_vs_sum_deltaP",25,0,500,25,0,500,"sum of nu Pt","sum soln delta P");
	// plot2d.AddHist("deltaZ_vs_deltaP",25,0,250,25,0,500,"nu delta z","sum soln delta P");
	// plot1d.AddHist("dns_sum_Pt",25,0,500,"sum nu Pt","events");
	// plot1d.AddHist("dns_sum_Pz",25,0,500,"sum nu Pt","events");
	// plot2d.AddHist("allbestsoln_MinDeltaPxy_vs_DeltaP",25,0,500,25,0,600,"min delta pxy","sum soln delta P");
   //
   //
	// plot1d.AddHist("bestsoln_sumdeltaP",25,0.,400,"best minimum sum delta P", "Events");
	// plot1d.AddHist("othersolns_sumdeltaP",25,0.,400,"other minima sum delta P", "Events");
	// plot1d.AddHist("bestsoln_Pt",25,0.,300,"best minimum sum Pt (GeV)", "Events");
	// plot1d.AddHist("othersolns_Pt",25,0.,300,"other minima sum Pt (GeV)", "Events");
	// plot1d.AddHist("bestsoln_deltaRlep",25,0.,7.5,"best minimum delta R", "Events");
	// plot1d.AddHist("othersolns_deltaRlep",25,0.,7.5,"other minima  delta R", "Events");
	// plot1d.AddHist("bestsoln_netPz",25,0.,400,"best minimum net Pz", "Events");
	// plot1d.AddHist("othersolns_netPz",25,0.,400,"other minima net Pz", "Events");
	// plot1d.AddHist("bestsoln_Pz",25,0.,300,"best minimum pz", "Events");
	// plot1d.AddHist("othersolns_Pz",25,0.,300,"other minima pz", "Events");
	// plot2d.AddHist("bestsoln_sumPt_vs_sum_deltaP",25,0,500,25,0,500,"sum of nu Pt (GeV)","sum soln delta P");
	// plot2d.AddHist("othersolns_sumPt_vs_sum_deltaP",25,0,500,25,0,500,"sum of nu Pt (GeV)","sum soln delta P");
	// plot1d.AddHist("bestsoln_angle_lep",25,0.,6.3,"angle with lep", "Events");
	// plot1d.AddHist("othersolns_angle_lep",25,0.,6.3,"angle with lep", "Events");
	// plot2d.AddHist("bestsoln_MinDeltaPxy_vs_DeltaP",25,0,500,25,0,600,"min delta pxy","sum soln delta P");
	// plot2d.AddHist("othersolns_MinDeltaPxy_vs_DeltaP",25,0,500,25,0,600,"min delta pxy","sum soln delta P");
	// plot2d.AddHist("othersolns_DeltaPbest_vs_DeltaPgen",25,0,500,25,0,600,"min delta pxy (GeV)","sum soln delta P");
	// plot1d.AddHist("othersolns_diff_Pt",25,-250,250,"diff in Pt (GeV)", "Events");
	// plot1d.AddHist("othersolns_diff_sum_Pt",25,-250,250,"diff in Pt (GeV)", "Events");
	// plot1d.AddHist("othersolns_diff_max_Pt",25,-250,250,"diff in Pt (GeV)", "Events");
   //
	// plot1d.AddHist("othersolns_diff_absPz",25,-250,250,"diff in abs(Pz)", "Events");


}

void DilepPlots::Fill(Dilepton& di, double weight)
{
	plot1d["b&lep_inv_mass"]->Fill((*di.Bt()+*di.Lt()).M(),weight);
	plot1d["b&lep_inv_mass"]->Fill((*di.Btbar()+*di.Ltbar()).M(),weight);
	plot1d["b&lep_sum_deltar"]->Fill(di.Bt()->DeltaR(*di.Lt())+di.Btbar()->DeltaR(*di.Ltbar()), weight);
	plot1d["dns_min_num"]->Fill(di.MinNum(),weight);
}

void DilepPlots::FillNuInfo(Dilepton& di, Dilepton& gendi, double weight)
{
	if ( di.Error() || !di.Solved() ) {
		return;
	}
	di.SortSolutions(gendi);

	di.SetNeutrinos(0);

	vector<TLorentzVector> nuts = di.NutSolns();
	vector<TLorentzVector> nutbars = di.NutbarSolns();
	vector<TLorentzVector> nunus;

	vector<TLorentzVector> better_nuts;
	vector<TLorentzVector> better_nutbars;
	vector<TLorentzVector> worse_nuts;
	vector<TLorentzVector> worse_nutbars;

	vector<TLorentzVector> good_nuts;
	vector<TLorentzVector> good_nutbars;
	vector<TLorentzVector> bad_nuts;
	vector<TLorentzVector> bad_nutbars;


	double ave = 1/(nuts.size());
	TLorentzVector avenut(0,0,0,0);
	TLorentzVector avenutbar(0,0,0,0);
	for (size_t i = 0; i < nuts.size(); i++) {
		avenut += ave*nuts[i];
		avenutbar += ave*nutbars[i];
		nunus.push_back(nuts[i]+nutbars[i]);

	}

	if (nuts.size()==1) {
		TLorentzVector Wt = *di.Lt() + nuts[0];
		TLorentzVector Wtbar = *di.Ltbar() + nutbars[0];
		TLorentzVector Top = *di.Bt() + Wt;
		TLorentzVector Tbar = *di.Btbar() + Wtbar;
		plot1d["W_mass_approx"]->Fill(Wt.M(),weight);
		plot1d["W_mass_approx"]->Fill(Wtbar.M(),weight);
		plot1d["Top_mass_approx"]->Fill(Top.M(),weight);
		plot1d["Top_mass_approx"]->Fill(Tbar.M(),weight);
		plot1d["W_mass_vsgen_approx"]->Fill(Wt.M()-gendi.Wt().M(),weight);
		plot1d["W_mass_vsgen_approx"]->Fill(Wtbar.M()-gendi.Wtbar().M(),weight);
		plot1d["Top_mass_vsgen_approx"]->Fill(Top.M()-gendi.Top().M(),weight);
		plot1d["Top_mass_vsgen_approx"]->Fill(Tbar.M()-gendi.Tbar().M(),weight);
	}

	if (nuts.size()>1) {
		for (size_t i = 0; i < nuts.size(); i++) {
			TLorentzVector Wt = *di.Lt() + nuts[i];
			TLorentzVector Wtbar = *di.Ltbar() + nutbars[i];
			TLorentzVector Top = *di.Bt() + Wt;
			TLorentzVector Tbar = *di.Btbar() + Wtbar;
			plot1d["W_mass_exact"]->Fill(Wt.M(),weight);
			plot1d["W_mass_exact"]->Fill(Wtbar.M(),weight);
			plot1d["Top_mass_exact"]->Fill(Top.M(),weight);
			plot1d["Top_mass_exact"]->Fill(Tbar.M(),weight);
			plot1d["W_mass_vsgen_exact"]->Fill(Wt.M()-gendi.Wt().M(),weight);
			plot1d["W_mass_vsgen_exact"]->Fill(Wtbar.M()-gendi.Wtbar().M(),weight);
			plot1d["Top_mass_vsgen_exact"]->Fill(Top.M()-gendi.Top().M(),weight);
			plot1d["Top_mass_vsgen_exact"]->Fill(Tbar.M()-gendi.Tbar().M(),weight);

		}

	}



	TLorentzVector avenunu = avenut+avenutbar;
	TLorentzVector trunut = *gendi.Nut();
	TLorentzVector trunutbar = *gendi.Nutbar();
	TLorentzVector trununu = gendi.NuNu();


	plot1d["dns_bestmin_deltaP_sum"]->Fill(deltaP(nuts[0],trunut)+deltaP(nutbars[0],trunutbar),weight);
	plot1d["dns_bestmin_deltaP_nunu"]->Fill(deltaP(nunus[0],trununu),weight);
	plot1d["dns_bestmin_deltaP_each"]->Fill(deltaP(nuts[0],trunut),weight);
	plot1d["dns_bestmin_deltaP_each"]->Fill(deltaP(nutbars[0],trunutbar),weight);
	plot1d["dns_bestmin_diffP_each"]->Fill(nuts[0].P()-trunut.P(),weight);
	plot1d["dns_bestmin_diffP_each"]->Fill(nutbars[0].P()-trunutbar.P(),weight);

	plot1d["dns_bestmin_deltaP_each_percent"]->Fill(deltaP(nuts[0],trunut)/trunut.P(),weight);
	plot1d["dns_bestmin_deltaP_each_percent"]->Fill(deltaP(nutbars[0],trunutbar)/trunutbar.P(),weight);

	plot1d["dns_bestmin_deltaPxy_sum"]->Fill(deltaPxy(nuts[0],trunut)+deltaPxy(nutbars[0],trunutbar),weight);
	plot1d["dns_bestmin_deltaPxy_nunu"]->Fill(deltaPxy(nunus[0],trununu),weight);
	plot1d["dns_bestmin_deltaPxy_each"]->Fill(deltaPxy(nuts[0],trunut),weight);
	plot1d["dns_bestmin_deltaPxy_each"]->Fill(deltaPxy(nutbars[0],trunutbar),weight);
	plot1d["dns_bestmin_deltaPz_sum"]->Fill(deltaPz(nuts[0],trunut)+deltaPz(nutbars[0],trunutbar),weight);
	plot1d["dns_bestmin_deltaPz_nunu"]->Fill(deltaPz(nunus[0],trununu),weight);
	plot1d["dns_bestmin_deltaPz_each"]->Fill(deltaPz(nuts[0],trunut),weight);
	plot1d["dns_bestmin_deltaPz_each"]->Fill(deltaPz(nutbars[0],trunutbar),weight);

	for (size_t i=0; i< nuts.size(); i++) {
		plot1d["allmin_Pt"]->Fill(nuts[i].Pt(),weight);
		plot1d["allmin_Pt"]->Fill(nutbars[i].Pt(),weight);
		plot1d["allmin_Pz"]->Fill(nuts[i].Pz(),weight);
		plot1d["allmin_Pz"]->Fill(nutbars[i].Pz(),weight);
		plot1d["allmin_P"]->Fill(nuts[i].P(),weight);
		plot1d["allmin_P"]->Fill(nutbars[i].P(),weight);
		plot1d["dns_allmin_diffP_each"]->Fill(nuts[i].P()-trunut.P(),weight);
		plot1d["dns_allmin_diffP_each"]->Fill(nutbars[i].P()-trunutbar.P(),weight);
		plot1d["dns_allmin_deltaP_sum"]->Fill(deltaP(nuts[i],trunut)+deltaP(nutbars[i],trunutbar),weight);
		plot1d["dns_allmin_deltaP_nunu"]->Fill(deltaP(nunus[i],trununu),weight);
		plot1d["dns_allmin_deltaP_each"]->Fill(deltaP(nuts[i],trunut),weight);
		plot1d["dns_allmin_deltaP_each"]->Fill(deltaP(nutbars[i],trunutbar),weight);
		plot1d["dns_allmin_deltaPxy_sum"]->Fill(deltaPxy(nuts[i],trunut)+deltaPxy(nutbars[i],trunutbar),weight);
		plot1d["dns_allmin_deltaPxy_nunu"]->Fill(deltaPxy(nunus[i],trununu),weight);

		plot1d["dns_allmin_deltaPxy_each"]->Fill(deltaPxy(nuts[i],trunut),weight);
		plot1d["dns_allmin_deltaPxy_each"]->Fill(deltaPxy(nutbars[i],trunutbar),weight);
		plot1d["dns_allmin_deltaPz_sum"]->Fill(deltaPz(nuts[i],trunut)+deltaPz(nutbars[i],trunutbar),weight);
		plot1d["dns_allmin_deltaPz_nunu"]->Fill(deltaPz(nunus[i],trununu),weight);
		plot1d["dns_allmin_deltaPz_each"]->Fill(deltaPz(nuts[i],trunut),weight);
		plot1d["dns_allmin_deltaPz_each"]->Fill(deltaPz(nutbars[i],trunutbar),weight);
		plot2d["dns_nu_vs_nubar_deltaP"]->Fill(deltaP(nuts[i],trunut),deltaP(nutbars[i],trunutbar),weight);
		plot2d["dns_deltaP_vs_Pt"]->Fill(nuts[i].Pt(),deltaP(nuts[i],trunut),weight);
		plot2d["dns_deltaP_vs_Pt"]->Fill(nutbars[i].Pt(),deltaP(nutbars[i],trunutbar),weight);

		plot2d["dns_deltaP_vs_Pz"]->Fill(nuts[i].Pz(),deltaP(nuts[i],trunut),weight);
		plot2d["dns_deltaP_vs_Pz"]->Fill(nutbars[i].Pz(),deltaP(nutbars[i],trunutbar),weight);

		plot2d["dns_deltaP_vs_nummins"]->Fill(deltaP(nuts[i],trunut),nuts.size(),weight);
		plot2d["dns_deltaP_vs_nummins"]->Fill(deltaP(nutbars[i],trunutbar),nuts.size(),weight);



	}

	// plot1d["dns_avemin_deltaP_sum"]->Fill(deltaP(avenut,trunut)+deltaP(avenutbar,trunutbar),weight);
	// plot1d["dns_avemin_deltaP_nunu"]->Fill(deltaP(avenunu,trununu),weight);
	// plot1d["dns_avemin_deltaP_each"]->Fill(deltaP(avenut,trunut),weight);
	// plot1d["dns_avemin_deltaP_each"]->Fill(deltaP(avenutbar,trunutbar),weight);
	// plot1d["dns_avemin_deltaPxy_sum"]->Fill(deltaPxy(avenut,trunut)+deltaPxy(avenutbar,trunutbar),weight);
	// plot1d["dns_avemin_deltaPxy_nunu"]->Fill(deltaPxy(avenunu,trununu),weight);
	// plot1d["dns_avemin_deltaPxy_each"]->Fill(deltaPxy(avenut,trunut),weight);
	// plot1d["dns_avemin_deltaPxy_each"]->Fill(deltaPxy(avenutbar,trunutbar),weight);
	// plot1d["dns_avemin_deltaPz_sum"]->Fill(deltaPz(avenut,trunut)+deltaPz(avenutbar,trunutbar),weight);
	// plot1d["dns_avemin_deltaPz_nunu"]->Fill(deltaPz(avenunu,trununu),weight);
	// plot1d["dns_avemin_deltaPz_each"]->Fill(deltaPz(avenut,trunut),weight);
	// plot1d["dns_avemin_deltaPz_each"]->Fill(deltaPz(avenutbar,trunutbar),weight);



	plot1d["dns_TT_IM_error_best"]->Fill( ( di.TT().M()-gendi.TT().M() ),weight);
	plot1d["dns_TT_IM_percenterror_best"]->Fill( ( di.TT().M()-gendi.TT().M() )/gendi.TT().M(),weight);
	plot1d["dns_TT_DeltaY_error_best"]->Fill( abs( di.Top().Rapidity()-di.Tbar().Rapidity() ) - abs( gendi.Top().Rapidity()-gendi.Tbar().Rapidity() ),weight);
	plot1d["dns_TT_DeltaY_percenterror_best"]->Fill( (abs( di.Top().Rapidity()-di.Tbar().Rapidity() ) - abs( gendi.Top().Rapidity()-gendi.Tbar().Rapidity() ))
	 			/ abs( gendi.Top().Rapidity()-gendi.Tbar().Rapidity() ),weight);

	plot1d["dns_bestmin_momentum_percenterror"]->Fill( ( deltaP(*di.Nut(),*gendi.Nut() )/ gendi.Nut()->P() ),weight);
	plot1d["dns_bestmin_momentum_percenterror"]->Fill( ( deltaP(*di.Nutbar(),*gendi.Nutbar() )/ gendi.Nutbar()->P() ),weight);


	// if (.9< deltaP(*di.Nut(),*gendi.Nut() )/ gendi.Nut()->P() && deltaP(*di.Nut(),*gendi.Nut() )/ gendi.Nut()->P() < 1.1) {
	// 	cout<<"###############################"<<endl;
	// 	gendi.Nut()->Print();
	// 	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	// 	di.Nut()->Print();
	// }

}

void DilepPlots::FillErrorInfo(Dilepton& di, Dilepton& dilo, Dilepton& dihi, Dilepton& gendi, double weight)
{
	di.SortSolutions(gendi);
	dilo.SortSolutions(gendi);
	dihi.SortSolutions(gendi);

	vector<TLorentzVector> nuts = di.NutSolns();
	vector<TLorentzVector> nutbars = di.NutbarSolns();

	TLorentzVector trunut = *gendi.Nut();
	TLorentzVector trunutbar = *gendi.Nutbar();

	if (nuts.size()==0) {
		return;
	}
	vector<TLorentzVector> nuts_lo = dilo.NutSolns();
	vector<TLorentzVector> nutbars_lo = dilo.NutbarSolns();
	vector<TLorentzVector> nuts_hi = dihi.NutSolns();
	vector<TLorentzVector> nutbars_hi = dihi.NutbarSolns();

	if (nuts_lo.size()>0) {
		plot1d["bjerr_change_in_best_sumdeltaP"]->Fill(deltaP(nuts_lo[0],nuts[0])+deltaP(nutbars_lo[0],nutbars[0]),weight);
	}
	if (nuts_hi.size()>0) {
		plot1d["bjerr_change_in_best_sumdeltaP"]->Fill(deltaP(nuts_hi[0],nuts[0])+deltaP(nutbars_hi[0],nutbars[0]),weight);
	}
	if (nuts_hi.size()>0 && nuts_lo.size()>0) {
		plot1d["bjerr_change_in_best_sumdeltaP_hilo"]->Fill(deltaP(nuts_hi[0],nuts_lo[0])+deltaP(nutbars_hi[0],nutbars_lo[0]),weight);
	}


	if (nuts.size()>0 && nuts_lo.size()>0) {
		plot1d["bjerr_btest_lo"]->Fill(deltaP(*dilo.Bt(),*di.Bt()) ,weight);
		plot1d["bjerr_btest_lo"]->Fill(deltaP(*dilo.Btbar(),*di.Btbar()) ,weight);
		plot1d["bjerr_avedist_lo"]->Fill(AveSolnDist(nuts,nutbars,nuts_lo,nutbars_lo) , weight);
	}
	if (nuts_hi.size()>0 && nuts.size()>0) {

		// plot1d["bjerr_avedist_hi"]->Fill(AveSolnDist(nuts,nutbars,nuts_hi,nutbars_hi) , weight);
		plot1d["bjerr_btest_hi"]->Fill(deltaP(*dihi.Bt(),*di.Bt()) ,weight);
		plot1d["bjerr_btest_hi"]->Fill(deltaP(*dihi.Btbar(),*di.Btbar()) ,weight);
		plot1d["bjerr_avedist_hi"]->Fill(AveSolnDist(nuts_hi,nutbars_hi,nuts,nutbars) , weight);
	}
	if (nuts_hi.size()>0 && nuts_lo.size()>0) {
		// plot1d["bjerr_avedist_hilo"]->Fill(AveSolnDist(nuts_lo,nutbars_lo,nuts_hi,nutbars_hi) , weight);
		plot1d["bjerr_avedist_hilo"]->Fill( AveSolnDist(nuts_lo,nutbars_lo,nuts_hi,nutbars_hi) , weight);
		plot1d["bjerr_btest_hilo"]->Fill(deltaP(*dihi.Bt(),*dilo.Bt()) ,weight);
		plot1d["bjerr_btest_hilo"]->Fill(deltaP(*dihi.Btbar(),*dilo.Btbar()) ,weight);
		plot2d["bestmin_deltaP_vs_avedist_hilo"]->Fill(AveSolnDist(nuts_lo,nutbars_lo,nuts_hi,nutbars_hi), deltaP(nuts[0],trunut) + deltaP(nutbars[0],trunutbar) ,weight);
		plot2d["avedist_hilo_vs_numminsA"]->Fill(AveSolnDist(nuts_lo,nutbars_lo,nuts_hi,nutbars_hi),nuts_lo.size(),weight);
		plot2d["avedist_hilo_vs_numminsB"]->Fill(AveSolnDist(nuts_lo,nutbars_lo,nuts_hi,nutbars_hi),nuts_hi.size(),weight);
	}

}

double DilepPlots::AveSolnDist(vector<TLorentzVector> nutsolnsA, vector<TLorentzVector> nutbarsolnsA, vector<TLorentzVector> nutsolnsB, vector<TLorentzVector> nutbarsolnsB)
{
	vector<TLorentzVector> less_nuts;
	vector<TLorentzVector> less_nutbars;
	vector<TLorentzVector> more_nuts;
	vector<TLorentzVector> more_nutbars;

	if (nutsolnsA.size() <= nutsolnsB.size()) {
		less_nuts = nutsolnsA;
		less_nutbars = nutbarsolnsA;
		more_nuts = nutsolnsB;
		more_nutbars = nutbarsolnsB;
	}
	else
	{
		less_nuts = nutsolnsB;
		less_nutbars = nutbarsolnsB;
		more_nuts = nutsolnsA;
		more_nutbars = nutbarsolnsA;
	}


	double avedist = 0;

	for (size_t n=0; n< less_nuts.size(); n++) {
		double this_soln_mindist = 10000;
		size_t besti = 0;
		for (size_t i = 0; i < more_nuts.size() ; i++) {
			double thisdist = deltaP(less_nuts[n], more_nuts[i]) + deltaP(less_nutbars[n], more_nutbars[i]);
			if (thisdist < this_soln_mindist) {
				this_soln_mindist = thisdist;
				besti = i;
			}
		}
		more_nuts.erase(more_nuts.begin() + besti);
		more_nutbars.erase(more_nutbars.begin() + besti);
		avedist += this_soln_mindist;
	}
	avedist /= less_nuts.size();

	return avedist;
}
