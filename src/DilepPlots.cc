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
	plot1d.AddHist("dns_min_num",6,-0.5,5.5,"number of mins","Events");
	plot1d.AddHist("nu_Pt",100,0,500,"nu Pt","Events");
	plot1d.AddHist("nu_Pz",100,0,500,"nu Pz","Events");
	plot1d.AddHist("nu_P",100,0,600,"nu P","Events");
	plot1d.AddHist("M_tt",200,0,2000,"nu Pt","Events");

	plot1d.AddHist("nu_Pt_diff",100,-250,250,"nu Pt: reco-gen (GeV)","Events");
	plot1d.AddHist("nu_Pz_diff",100,-250,250,"nu Pz: reco-gen (GeV)","Events");
	plot1d.AddHist("nu_P_diff",100,-300,300,"nu P: reco-gen (GeV)","Events");
	plot1d.AddHist("M_tt_diff",200,-1000,1000,"M_{tt} (reco) - M_{tt} (gen)","Events");

	plot1d.AddHist("nu_Pt_err",100,-3.,3.,"nu Pt error (fractional)","Events");
	plot1d.AddHist("nu_Pz_err",100,-3.,3,"nu Pz error (fractional)","Events");
	plot1d.AddHist("nu_P_err",100,-3,3,"nu P error (fractional)","Events");
	plot1d.AddHist("M_tt_err",200,-2.5,2.5,"M_tt error (fractional)","Events");

}

void DilepPlots::Fill(Dilepton& di, double weight)
{
	if (!di.IsComplete()) {
		return;
	}

	plot1d["b&lep_inv_mass"]->Fill((*di.Bt()+*di.Lt()).M(),weight);
	plot1d["b&lep_inv_mass"]->Fill((*di.Btbar()+*di.Ltbar()).M(),weight);
	plot1d["b&lep_sum_deltar"]->Fill(di.Bt()->DeltaR(*di.Lt())+di.Btbar()->DeltaR(*di.Ltbar()), weight);
	plot1d["nu_Pt"]->Fill(di.Nut()->Pt(),weight);
	plot1d["nu_Pz"]->Fill(di.Nut()->Pz(),weight);
	plot1d["nu_P"]->Fill(di.Nut()->P(),weight);
	plot1d["nu_Pt"]->Fill(di.Nutbar()->Pt(),weight);
	plot1d["nu_Pz"]->Fill(di.Nutbar()->Pz(),weight);
	plot1d["nu_P"]->Fill(di.Nutbar()->P(),weight);
	plot1d["M_tt"]->Fill(di.Mtt(),weight);
	plot1d["dns_min_num"]->Fill(di.MinNum(),weight);
}

void DilepPlots::FillNuInfo(Dilepton& di, Dilepton& gendi, double weight)
{
	if ( di.Error() || !di.Solved() || !di.IsComplete()) {
		return;
	}

	plot1d["nu_Pt_diff"]->Fill(di.Nut()->Pt()-gendi.Nut()->Pt(),weight);
	plot1d["nu_Pz_diff"]->Fill(di.Nut()->Pz()-gendi.Nut()->Pz(),weight);
	plot1d["nu_P_diff"]->Fill(di.Nut()->P()-gendi.Nut()->P(),weight);
	plot1d["M_tt_diff"]->Fill(di.Mtt()-gendi.Mtt(),weight);

	plot1d["nu_Pt_err"]->Fill((di.Nut()->Pt()-gendi.Nut()->Pt())/abs(gendi.Nut()->Pt()),weight);
	plot1d["nu_Pz_err"]->Fill((di.Nut()->Pz()-gendi.Nut()->Pz())/abs(gendi.Nut()->Pz()),weight);
	plot1d["nu_P_err"]->Fill((di.Nut()->P()-gendi.Nut()->P())/abs(gendi.Nut()->P()),weight);
	plot1d["M_tt_err"]->Fill((di.Mtt()-gendi.Mtt())/abs(gendi.Mtt()),weight);


}

void DilepPlots::FillErrorInfo(Dilepton& di, Dilepton& dilo, Dilepton& dihi, Dilepton& gendi, double weight)
{

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
