#include "ttbarxsec.h"
#include "TRandom3.h"

#include "Permutation.h"
#include "PDFuncertainty.h"
#include "NeutrinoSolver.h"
#include "ConfigParser.h"

using namespace std;

ttbar::ttbar(const std::string output_filename):
	AnalyzerBase("ttbar", output_filename),
	FULLHAD(false),
	SEMILEP(false),
	FULLLEP(false),
	SEMILEPACC(false),
	gen1d("gen"),
	gen2d("gen"),
	ttp_genall("genall"),
	ttp_genacc("genacc"),
	ttp_genjet("genjet"),
	reco1d("reco"),
	reco2d("reco"),
	truth1d("truth"),
	truth2d("truth"),
	ttp_truth("truth"),
	ttp_right("right"),
	ttp_wrong("wrong"),
	ttp_semi("semi"),
	ttp_other("other"),
	ttp_all("all"),
	//ttp_jetspos_right("jetspos_right"),
	//ttp_jetspos_wrong("jetspos_wrong"),
	//	ttp_hadjets_right("hadjets_right"),
	//	ttp_hadjets_wrong("hadjets_wrong"),
	//	ttp_jets_right("jets_right"),
	//	ttp_jets_wrong("jets_wrong"),
	//	ttp_blep_right("blep_right"),
	//	ttp_blep_wrong("blep_wrong"),
	//ttp_whad_right("whad_right"),
	//ttp_whad_wrong("whad_wrong"),
	ttp_tlepthad_right("ttp_tlepthad_right"),
	ttp_tlep_right("ttp_tlep_right"),
	ttp_thad_right("ttp_thad_right"),
	ttp_nn_right("ttp_nn_right"),
	ttp_nsemi_right("ttp_nsemi_right"),
	ttp_alljets_right("ttp_alljets_right"),
	ttp_alljets_wrong("ttp_alljets_wrong"),
	response("response", this),
	response2d("response"),
	response2dvar("response"),
	response_ps("response_ps", this),
	PDFTEST(false),
	PSEUDOTOP(true),
	LHCPS(true),
	BTAGMODE(false), //set true for the b-tag efficiency measurement
	JETSCALEMODE(false), //set true for jet scale measurement
	ELECTRONS(true),
	MUONS(true),
	//B_TIGHT(0.935),
	//B_MEDIUM(0.800),
	//B_LOOSE(0.460),
	B_TIGHT(0.9535),
	B_MEDIUM(0.8484),
	B_LOOSE(0.5426),
	cnbtag("sigml"), //1: one thight b-jet, 2: two medium b-jets
	cnusedjets(1000), //only nused jets, ordered by pT are used for the permutations
	cwjetptsoft(25.), //min pT of softer W-jet
	cwjetpthard(35.), //min pT of harder W-jet 
	cbjetptsoft(25.), //min pT of softer b-jets
	cbjetpthard(35.), //min pT of harder b-jets
	cjetetamax(2.4),//max |eta| for jets
	clptmin(30.), //min pT of lepton (el or mu)
	cletamax(2.1),//max |eta| of leptons (max allowed value is 2.4) 
	cpwjetptsoft(25.), //min pT of softer W-jet
	cpwjetpthard(35.), //min pT of harder W-jet 
	cpbjetptsoft(25.), //min pT of softer b-jets
	cpbjetpthard(35.), //min pT of harder b-jets
	cpjetetamax(2.4),//max |eta| for jets
	cplptmin(30.), //min pT of lepton (el or mu)
	cpletamax(2.1),//max |eta| of leptons (max allowed value is 2.4) 
	cpjetsep(0.),
	csigmajet(0.),
	cjetres(0),
	csigmamet(0.),
	crenscale(0),
	cfacscale(0),
	chdamp(0),
	cbtagunc(0),
	cltagunc(0),
	cpileup(0),
	TTMC(false),
	HERWIGPP(false),
	PYTHIA6(false),
	ISRUP(false),
	ISRDOWN(false),
	FSRUP(false),
	FSRDOWN(false),
	TUNEUP(false),
	TUNEDOWN(false)
{
	ConfigParser CP("ttbarxsec.cfg");
	cLeptonScaleFactor = CP.Get<string>("LeptonScaleFactor");
	cJetEnergyUncertainty = CP.Get<string>("JetEnergyUncertainty");
	cBTaggingSF = CP.Get<string>("BTaggingSF");
	cBTaggingEff = CP.Get<string>("BTaggingEff");
	cJetResolution = CP.Get<string>("JetResolution");
	cJetResolutionSF = CP.Get<string>("JetResolutionSF");
	cel27eff = CP.GetVector<double>("el27eff");
	PSEUDOTOP = CP.Get<bool>("PSEUDOTOP");
	LHCPS = CP.Get<bool>("LHCPS");
	BTAGMODE = CP.Get<bool>("BTAGMODE");
	ELECTRONS = CP.Get<bool>("ELECTRONS");
	MUONS = CP.Get<bool>("MUONS");
	IDMuon::USEISO = CP.Get<bool>("LEPTONISO");
	IDElectron::USEISO = CP.Get<bool>("LEPTONISO");
	cnbtag = CP.Get<string>("nbtag");
	clikelihoodcut = CP.Get<double>("likelihoodcut");

	cwjetptsoft = CP.Get<double>("wjetptsoft"); 
	cwjetpthard = CP.Get<double>("wjetpthard"); 
	cbjetptsoft = CP.Get<double>("bjetptsoft"); 
	cbjetpthard = CP.Get<double>("bjetpthard"); 
	cjetetamax = CP.Get<double>("jetetamax");
	clptmin = CP.Get<double>("lptmin"); 
	cletamax = CP.Get<double>("letamax");

	cpwjetptsoft = CP.Get<double>("Pwjetptsoft"); 
	cpwjetpthard = CP.Get<double>("Pwjetpthard"); 
	cpbjetptsoft = CP.Get<double>("Pbjetptsoft"); 
	cpbjetpthard = CP.Get<double>("Pbjetpthard"); 
	cpjetetamax = CP.Get<double>("Pjetetamax");
	cplptmin = CP.Get<double>("Plptmin"); 
	cpletamax = CP.Get<double>("Pletamax");
	cpjetsep = CP.Get<double>("Pjetsep");

	csigmajet = CP.Get<double>("sigmajet");
	csigmajetwj = CP.Get<double>("sigmajetwj");
	cscalejetwj = CP.Get<double>("scalejetwj");
	cjecuncertainty = CP.Get<string>("jecuncertainty");
	cjetres = CP.Get<int>("jetres");
	csigmamet = CP.Get<double>("sigmamet");
	csigmalep = CP.Get<double>("sigmalep");
	ctopptweight = CP.Get<double>("topptweight");
	ctoprapweight = CP.Get<double>("toprapweight");
	cttptweight = CP.Get<double>("ttptweight");
	cjetptweight = CP.Get<double>("jetptweight");
	if(output_filename.find("tt_PowhegP8") != string::npos)
	{
		cbdecay = CP.Get<double>("bdecay");
		cbfrag = CP.Get<double>("bfrag");
		cbsplitting = CP.Get<double>("bsplitting");
		ccsplitting = CP.Get<double>("csplitting");
		cfacscale = CP.Get<int>("facscale");
		crenscale = CP.Get<int>("renscale");
		chdamp = CP.Get<int>("hdamp");
		PDFTEST = CP.Get<bool>("PDFTEST");
	}
	cbtagunc = CP.Get<int>("btagunc");
	cltagunc = CP.Get<int>("ltagunc");
	cpileup = CP.Get<int>("pileupunc");

	cout << output_filename << endl;
	if(output_filename.find("tt_") == 0){TTMC = true;}
	if(CP.Get<bool>("mcspecificcorrections"))
	{
		if(output_filename.find("Hpp") != string::npos){HERWIGPP = true;}
		if(output_filename.find("P6") != string::npos){PYTHIA6 = true;}
		if(output_filename.find("isrup") != string::npos){ISRUP = true;}
		if(output_filename.find("isrdown") != string::npos){ISRDOWN = true;}
		if(output_filename.find("fsrup") != string::npos){FSRUP = true;}
		if(output_filename.find("fsrdown") != string::npos){FSRDOWN = true;}
		if(output_filename.find("tuneup") != string::npos){TUNEUP = true;}
		if(output_filename.find("tunedown") != string::npos){TUNEDOWN = true;}
	}
	isDA = 0;
	if(output_filename.find("DATAEL") == 0){isDA = 11;}
	if(output_filename.find("DATAMU") == 0){isDA = 13;}

	if(STUDENT)
	{
		stud_tf = new TFile(("tree_"+output_filename).c_str(), "recreate");
		stud_tr = new TTree("events", "ttbar events");
		stud_tr->Branch("num_gen", &num_gen, "num_gen/i");
		stud_tr->Branch("gen_px", gen_px, "gen_px[num_gen]");
		stud_tr->Branch("gen_py", gen_py, "gen_py[num_gen]");
		stud_tr->Branch("gen_pz", gen_pz, "gen_pz[num_gen]");
		stud_tr->Branch("gen_e", gen_e, "gen_e[num_gen]");
		stud_tr->Branch("gen_type", gen_type, "gen_type[num_gen]/I");
		stud_tr->Branch("num_det", &num_det, "num_det/i");
		stud_tr->Branch("det_px", det_px, "det_px[num_det]");
		stud_tr->Branch("det_py", det_py, "det_py[num_det]");
		stud_tr->Branch("det_pz", det_pz, "det_pz[num_det]");
		stud_tr->Branch("det_e", det_e, "det_e[num_det]");
		stud_tr->Branch("det_type", det_type, "det_type[num_det]/I");
	}
	jetptmin = min(cwjetptsoft, cbjetptsoft);
// binning TOP-16-008
//	topptbins = {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 400.0, 800.0};
//	topybins = {0.0, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6, 2.5};
//	ttmbins = {300.0, 375.0, 450.0, 530.0, 625.0, 740.0, 850.0, 1100.0, 2000.0};
//	ttybins = {0.0, 0.2, 0.4, 0.6, 0.9, 1.3, 2.3};
//	ttptbins = {0.0, 35.0, 80.0, 140.0, 200.0, 500.0};
//	metbins = {0.0, 30.0, 45.0, 60.0, 80.0, 120.0, 580.0};
//	jetbins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5};
//	dybins = {-2.0, -1.5, -1., -0.5, 0., 0.5, 1.0, 1.5, 2.0};
//	nobins = {0., 13000.};

	//topptbins = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 290.0, 335.0, 380.0, 430.0, 500.0, 800.0};
	topptbins = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 430.0, 500.0, 800.0};
	topybins = {0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5};
	//ttmbins = {300.0, 350.0, 420.0, 500.0, 580.0, 680.0, 780.0, 900.0, 1100.0, 1200.0, 1400.0, 2000.0},
	ttmbins = {300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 1500.0, 2500.0},
	ttybins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4};
	//ttptbins = {0.0, 35.0, 70.0, 110.0, 165.0, 230.0, 300.0, 370.0, 450.0, 550.0};
	ttptbins = {0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 1000.0};
	metbins = {0.0, 30.0, 45.0, 60.0, 80.0, 120.0, 580.0};
	htbins = {120.0, 180.0, 240.0, 300.0, 360.0, 420.0, 480.0, 540.0, 600.0, 700.0, 800.0, 900.0, 1100.0, 1500.0, 2500.0};
	evtmassbins = {300., 400., 500., 600., 700., 800., 900., 1000., 1200.0, 1500.0, 2000., 2500., 5500.};
	jetbins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
	dybins = {-2.0, -1.5, -1., -0.5, 0., 0.5, 1.0, 1.5, 2.0};
	nobins = {0., 13000.};

	//btagpt = {30., 50., 70., 100., 140., 200., 300., 670.};
	btagpt = {30., 40., 50., 60., 80., 100., 120., 140., 170., 200., 250., 300., 400., 670.};
	//vector<string> testpdf = {"CT10", "CT10as", "NNPDF30_nnlo_as_0118"};
	//vector<string> testpdf = {"CT10nlo", "NNPDF30_nlo_as_0118", "MMHT2014nlo68clas118"};
	vector<string> testpdf;
	if(PDFTEST) {pdfunc = new PDFuncertainty("NNPDF30_nlo_as_0118", 0, testpdf);}


}

void ttbar::begin()
{
	outFile_.cd();
	TDirectory* dir_gen = outFile_.mkdir("GEN");
	dir_gen->cd();
	gen1d.AddHist("TYP", 5, 0., 5., "Decay TYP", "Events");
	gen1d.AddHist("tpt", 500, 0, 1000, "p_{T}(t) [GeV]", "Events");
	gen1d.AddHist("ty", 500, 0, 10, "y(t)", "Events");
	gen1d.AddHist("ttpt", 500, 0, 1000, "p_{T}(t#bar{t}) [GeV]", "Events");
	gen1d.AddHist("tty", 500, 0, 10, "y(t#bar{t})", "Events");
	gen2d.AddHist("plmuiso", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plmuisoprompt", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plmuisott", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("pleliso", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plelisoprompt", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plelisott", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plphiso", 100, 0, 500,200, 0, 2, "phiso", "Events");
	gen2d.AddHist("plphisoprompt", 100, 0, 500,200, 0, 2, "phiso", "Events");
	gen2d.AddHist("plphisoreco", 100, 0, 500,200, 0, 2, "phiso", "Events");
	gen2d.AddHist("plmuisoreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plmuisopromptreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plmuisottreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plelisoreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plelisopromptreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("plelisottreco", 100, 0, 500, 200, 0, 2, "lepiso", "Events");
	gen2d.AddHist("met_mt", 500, 0, 500, 500, 0, 500, "met", "mt");
    ttp_genall.Init(this);
    ttp_genacc.Init(this);
    ttp_genjet.Init(this);

	TDirectory* dir_truth = outFile_.mkdir("TRUTH");
	dir_truth->cd();
	truth1d.AddHist("counter", 20, 0., 20., "counter", "Events");
	truth1d.AddHist("npu", 100, 0., 100., "npu", "Events");
	truth1d.AddHist("npuorig", 100, 0., 100., "npuorig", "Events");
	truth1d.AddHist("found", 10, 0., 10., "found", "Events");
	truth1d.AddHist("foundgen", 10, 0., 10., "found", "Events");
	truth2d.AddHist("btag2d_true", 100, 0., 1., 100, 0., 1., "btag", "btag");
	truth1d.AddHist("btag_true", 100, 0., 1., "btag", "Events");
	truth1d.AddHist("btag_wrong", 100, 0., 1., "btag", "Events");
	truth1d.AddHist("wjet_sep", 100, 0., 100, "sep", "Events");	
	truth1d.AddHist("bjet_sep", 100, 0., 100, "sep", "Events");
	truth2d.AddHist("tt_jets", 3, 0., 3., 3, 0., 3., "num b-jets", "num w-jets");
	truth2d.AddHist("Jetstt_JetsAll", 5, 0., 5., 20, 0., 20., "jet tt-system", "all jets");
	truth2d.AddHist("Nmu_Ne", 5, 0., 5., 5, 0., 5., "# #mu", "# e");
	truth2d.AddHist("Ne_Nmu", 5, 0., 5., 5, 0., 5., "# e", "# #mu");
	truth2d.AddHist("jetscale_lB", 800, 0., 800., 200, -1., 1., "p_{T} [GeV]", "#Deltap_{T}/p_{T}");
	truth2d.AddHist("jetscale_lE", 800, 0., 800., 200, -1., 1., "p_{T} [GeV]", "#Deltap_{T}/p_{T}");
	truth2d.AddHist("jetscale_bB", 800, 0., 800., 200, -1., 1., "p_{T} [GeV]", "#Deltap_{T}/p_{T}");
	truth2d.AddHist("jetscale_bE", 800, 0., 800., 200, -1., 1., "p_{T} [GeV]", "#Deltap_{T}/p_{T}");
	truth1d.AddHist("dRNu_right", 200, 0., 10., "#DeltaR(#nu_{gen}, #nu_{rec})", "Events");
	truth1d.AddHist("dPtNu_right", 200, -2., 2., "#Deltap_{T}/p_{T}", "Events");
	truth1d.AddHist("dPzNu_right", 200, -2., 2., "#Deltap_{z}/p_{z}", "Events");
	truth2d.AddHist("dPzNu_dPhi_right", 200, -2., 2., 100, 0., 3.15, "#Deltap_{z}/p_{z}", "#Delta#Phi(#nu, met)");
	truth2d.AddHist("dPzNu_chi2_right", 200, -2., 2., 100, 0., 100, "#Deltap_{z}/p_{z}", "#Delta#Phi(#nu, met)");
	truth2d.AddHist("RES_Mtt_right", 200, -1., 1., 50, 300, 1300 , "#Delta M_{tt}/M_{tt}", "M_{tt} [GeV]");
	truth2d.AddHist("RES_Mtt_all", 200, -1., 1., 50, 300, 1300 , "#Delta M_{tt}/M_{tt}", "M_{tt} [GeV]");
	truth2d.AddHist("RES_dbeta_right", 200, -1., 1., 50, 0, 2 , "#Delta (#delta #beta)/#delta #beta", "#delta #beta");
	truth2d.AddHist("RES_dbeta_all", 200, -1., 1., 50, 0, 2 , "#Delta (#delta #beta)/#delta #beta", "#delta #beta");
	truth1d.AddHist("dRNuMet_right", 200, 0., 10., "#DeltaR(#nu_{gen}, #nu_{rec})", "Events");
	truth1d.AddHist("dPtNuMet_right", 200, -2., 2., "#Deltap_{T}/p_{T}", "Events");
	truth2d.AddHist("Wmasshad_tmasshad_right", 500, 0., 500., 500, 0., 500, "M(t_{h}) [GeV]", "M(W_{h}) [GeV]");
	truth2d.AddHist("Wmtlep_tmtlep_right", 500, 0., 500., 500, 0., 500, "M_{T}(t_{l}) [GeV]", "M(W_{l}) [GeV]");
	truth2d.AddHist("Wmasshad_tmasshad_wrong", 500, 0., 500., 500, 0., 500, "M(t) [GeV]", "M(W) [GeV]");
	truth2d.AddHist("Wmtlep_tmtlep_wrong", 500, 0., 500., 500, 0., 500, "M_{T}(t_{l}) [GeV]", "M(W_{l}) [GeV]");
	truth2d.AddHist("Wmasshad_tmasshad_rightw", 500, 0., 500., 500, 0., 500, "M(W) [GeV]", "M(t) [GeV]");
	truth2d.AddHist("Wmasshad_tmasshad_wrongw", 500, 0., 500., 500, 0., 500, "M(W) [GeV]", "M(t) [GeV]");
	truth1d.AddHist("btagtest_wrong", 1000, -100, 100., "-Log(p) btag-test", "Events");
	truth1d.AddHist("masstest_wrong", 1000, -100, 100., "-Log(p) mass-test", "Events");
	truth1d.AddHist("nstest_wrong", 200, 0, 20., "neutrino-test", "Events");
	truth1d.AddHist("nschi_wrong", 75, 0., 150., "#chi neutrino-test", "Events");
	truth1d.AddHist("comtest_wrong", 1000, -100, 100., "-Log(p)", "Events");
	truth1d.AddHist("btagtest_right", 1000, -100, 100., "-Log(p) btag-test", "Events");
	truth1d.AddHist("masstest_right", 1000, -100, 100., "-Log(p) mass-test", "Events");
	truth1d.AddHist("nstest_right", 200, 0, 20., "neutrino-test", "Events");
	truth1d.AddHist("nschi_right", 75, 0., 150., "#chi neutrino-test", "Events");
	truth1d.AddHist("comtest_right", 1000, -100, 100., "-Log(p)", "Events");
	truth2d.AddHist("Elrho_iso_1", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Elrho_iso_2", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Elrho_iso_3", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Elrho_iso_4", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Murho_iso_1", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Murho_iso_2", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Murho_iso_3", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth2d.AddHist("Murho_iso_4", 10, 0., 50., 100, 0., 150, "rho", "iso");
	truth1d.AddHist("TTRECO", 20, 0, 20, "ttreco", "Events");
	truth1d.AddHist("Eff_BpassingM", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_BpassingL", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_Ball", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_CpassingM", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_CpassingL", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_Call", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_LpassingM", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_LpassingL", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Eff_Lall", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("Mu", 100, 0, 100, "#mu", "Events");
	truth1d.AddHist("MuWeighted", 100, 0, 100, "#mu", "Events");
	truth1d.AddHist("EffL_BpassingM", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("EffL_BpassingL", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("EffL_Ball", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("EffH_BpassingM", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("EffH_BpassingL", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("EffH_Ball", btagpt, "p_{T} [GeV]", "Events");
	truth1d.AddHist("xb", 150, 0, 1.5, "x_{b}", "Events");
	truth1d.AddHist("xc", 150, 0, 1.5, "x_{c}", "Events");
	truth1d.AddHist("bids_had", 10000, 0, 10000, "pdgids b", "Events");
	truth1d.AddHist("bids_lep", 10000, 0, 10000, "pdgids b", "Events");

	response.AddMatrix("thardpt", topptbins, topptbins, "high p_{T}(t) [GeV]");
	response.AddMatrix("tsoftpt", topptbins, topptbins, "soft p_{T}(t) [GeV]");
	response.AddMatrix("thadpt", topptbins, topptbins, "p_{T}(t_{h}) [GeV]");
	response.AddMatrix("thady", topybins, topybins, "|y(t_{h})|");
	response.AddMatrix("tleppt", topptbins, topptbins, "p_{T}(t_{l}) [GeV]");
	response.AddMatrix("tlepy", topybins, topybins, "|y(t_{l})|");
	response.AddMatrix("ttm", ttmbins, ttmbins, "m(t#bar{t}) [GeV]");
	response.AddMatrix("ttpt", ttptbins, ttptbins, "p_{T}(t#bar{t}) [GeV]");
	response.AddMatrix("tty", ttybins, ttybins, "|y(t#bar{t})|");
	response.AddMatrix("njet", jetbins, jetbins, "n-jets", true);
	response.AddMatrix("nobin", nobins, nobins, "total", true);
	response.AddMatrix("dymp", dybins, dybins, "y(t)-y(#bar{t})");
	response.AddMatrix("dy", dybins, dybins, "|y(t)|-|y(#bar{t})|");
	response.AddMatrix("ht", htbins, htbins, "H_{T} [GeV]");
	response.AddMatrix("evtmass", evtmassbins, evtmassbins, "M_{evt} [GeV]");

	vector<double> jetbins_large = {-0.5, 0.5, 1.5, 2.5, 3.5};
	vector<double> topptbins_large = {0.0, 90.0, 180.0, 270.0, 800.0};
	vector<double> topybins_large = {0.0, 0.5, 1.0, 1.5, 2.5};
	vector<double> ttmbins_large = {300.0, 450.0, 625.0, 850.0, 2000.0};
	vector<double> ttybins_large = {0.0, 0.4, 0.8, 1.2, 2.3};
	vector<double> ttptbins_large = {0.0, 35.0, 80.0, 140.0, 500.0};
	//vector<double> topptbins_large = {0.0, 60.0, 100.0, 140.0, 180.0, 250.0, 500.0};
	response2d.AddMatrix("njets_thadpt", topptbins, jetbins_large, topptbins, jetbins_large, false, true);
	//vector<double> ttptbins_large = {0.0, 25.0, 45.0, 70.0, 110.0, 200., 500.0};
	response2d.AddMatrix("njets_ttpt", ttptbins, jetbins_large, ttptbins, jetbins_large, false, true);
	response2d.AddMatrix("njets_ttm", ttmbins, jetbins_large, ttmbins, jetbins_large, false, true);
	response2d.AddMatrix("thady_thadpt", topptbins, topybins_large, topptbins, topybins_large);
	response2d.AddMatrix("ttpt_thadpt", topptbins, ttptbins_large, topptbins, ttptbins_large);
	response2d.AddMatrix("thadpt_ttm", ttmbins, topptbins_large, ttmbins, topptbins_large);
	response2d.AddMatrix("ttpt_ttm", ttmbins, ttptbins_large, ttmbins, ttptbins_large);
	response2d.AddMatrix("ttm_tty", ttybins, ttmbins_large, ttybins, ttmbins_large);
	response2d.AddMatrix("ttm_dy", dybins, ttmbins_large, dybins, ttmbins_large);

	vector< vector<double> > bins_njet_thadpt = {
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0}
	};
	response2dvar.AddMatrix("njet+thadpt", jetbins_large, bins_njet_thadpt, true, false);

	vector< vector<double> > bins_njet_ttpt = {
	{0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 1000.0},
	{0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 1000.0},
	{0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 1000.0},
	{0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 1000.0}
	};
	response2dvar.AddMatrix("njet+ttpt", jetbins_large, bins_njet_ttpt, true, false);

	vector< vector<double> > bins_njet_ttm = {
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 2000.0},
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 2000.0},
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 2000.0},
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 2000.0}
	};
	response2dvar.AddMatrix("njet+ttm", jetbins_large, bins_njet_ttm, true, false);

	vector< vector<double> > bins_thady_thadpt = {
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0}
	};
	response2dvar.AddMatrix("thady+thadpt", topybins_large, bins_thady_thadpt, false, false);

	vector< vector<double> > bins_thadpt_ttm = {
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 2000.0},
	{300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 2000.0},
	{300.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200., 2000.0},
	{300.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200., 2000.0}
	};
	response2dvar.AddMatrix("thadpt+ttm", topptbins_large, bins_thadpt_ttm, false, false);

	vector< vector<double> > bins_ttpt_ttm = {
	{300.0, 350.0, 420.0, 500.0, 580.0, 680.0, 780.0, 900.0, 1100.0, 2000.0},
	{300.0, 350.0, 420.0, 500.0, 580.0, 680.0, 780.0, 900.0, 1100.0, 2000.0},
	{300.0, 350.0, 420.0, 500.0, 580.0, 680.0, 780.0, 900.0, 1100.0, 2000.0},
	{300.0, 350.0, 420.0, 500.0, 580.0, 680.0, 780.0, 900.0, 1100.0, 2000.0}
	};
	response2dvar.AddMatrix("ttpt+ttm", ttptbins_large, bins_ttpt_ttm, false, false);

	vector< vector<double> > bins_ttpt_thadpt = {
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
	{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0}
	};
	response2dvar.AddMatrix("ttpt+thadpt", ttptbins_large, bins_ttpt_thadpt, false, false);

	vector< vector<double> > bins_ttm_tty = {
	{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.4},
	{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.4},
	{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.4},
	{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.4}
	};
	response2dvar.AddMatrix("ttm+tty", ttmbins_large, bins_ttm_tty, false, false);

	vector< vector<double> > bins_ttm_dy = {
	{-2.0, -1.4, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.4, 2.0},
	{-2.0, -1.4, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.4, 2.0},
	{-2.0, -1.4, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.4, 2.0},
	{-2.0, -1.4, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.4, 2.0}
	};
	response2dvar.AddMatrix("ttm+dy", ttmbins_large, bins_ttm_dy, false, false);


	alljetbins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	jetptbins = {
	{30., 50., 75., 100., 150., 200., 350.},
	{30., 50., 75., 100., 150., 200., 350.},
	{30., 50., 75., 100., 150., 200., 350.},
	{30., 50., 75., 100., 250.},
	{30., 50., 75., 100., 125., 150., 175., 200., 250., 320., 500.},
	{30., 50., 75., 100., 125., 150., 180., 350.},
	{30., 50., 75., 100., 250.},
	{30., 50., 75., 100., 200.}};
	response2dvar.AddMatrix("jet+jetpt", alljetbins, jetptbins, false, false);
	jetptbins = {
	{30., 50., 75., 100., 150., 200., 350.,351},
	{30., 50., 75., 100., 150., 200., 350.,351},
	{30., 50., 75., 100., 150., 200., 350.,351},
	{30., 50., 75., 100., 150.,151},
	{30., 50., 75., 100., 125., 150., 175., 200., 250., 320., 400.,401},
	{30., 50., 75., 100., 125., 150., 180., 350.,351},
	{30., 50., 75., 100., 200.,201},
	{30., 50., 75., 100., 120.,121}};
	response2dvar.AddMatrix("jet+jetptOF", alljetbins, jetptbins, false, true);
	jetetabins = {
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5},
	{0., 0.5, 1.0, 1.5, 2.0, 2.5},
	{0., 0.5, 1.0, 1.5, 2.0, 2.5}};
	response2dvar.AddMatrix("jet+jeteta", alljetbins, jetetabins, false, false);
	jetdrbins = {
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2., 2.5, 4.5},
	{0.4, 0.8, 1.2, 1.6, 2., 2.5, 4.5},
	{0.4, 0.8, 1.2, 1.6, 2., 2.5, 4.5}};
	response2dvar.AddMatrix("jet+jetdr", alljetbins, jetdrbins, false, false);
	jetdrtopbins = {
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.3, 0.6, 0.9, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.4, 0.8, 1.2, 1.5, 2., 2.5, 4.5},
	{0., 0.4, 0.8, 1.2, 1.5, 2., 2.5, 4.5}};
	response2dvar.AddMatrix("jet+jetdrtop", alljetbins, jetdrtopbins, false, false);

	response_ps.AddMatrix("thardpt", topptbins, topptbins, "hard p_{T}(t_{h}) [GeV]");
	response_ps.AddMatrix("tsoftpt", topptbins, topptbins, "soft p_{T}(t_{h}) [GeV]");
	response_ps.AddMatrix("thadpt", topptbins, topptbins, "p_{T}(t_{h}) [GeV]");
	response_ps.AddMatrix("thady", topybins, topybins, "|y(t_{h})|");
	response_ps.AddMatrix("tleppt", topptbins, topptbins, "p_{T}(t_{l}) [GeV]");
	response_ps.AddMatrix("tlepy", topybins, topybins, "|y(t_{l})|");
	response_ps.AddMatrix("ttm", ttmbins, ttmbins, "m(t#bar{t}) [GeV]");
	response_ps.AddMatrix("ttpt", ttptbins, ttptbins, "p_{T}(t#bar{t}) [GeV]");
	response_ps.AddMatrix("tty", ttybins, ttybins, "|y(t#bar{t})|");

	if(PDFTEST)
	{
		pdfunc->Add1dHist("pdfunc_thardpt", topptbins, "hard p_{T}(t) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_tsoftpt", topptbins, "soft p_{T}(t) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_thadpt", topptbins, "p_{T}(t_{h}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_thady", topybins, "|y(t_{h})|", "Events");
		pdfunc->Add1dHist("pdfunc_tleppt", topptbins, "p_{T}(t_{l}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_tlepy", topybins, "|y(t_{l})|", "Events");
		pdfunc->Add1dHist("pdfunc_ttm", ttmbins, "M(t#bar{t}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_tty", ttybins, "|y(t#bar{t})|", "Events");
		pdfunc->Add1dHist("pdfunc_ttpt", ttptbins, "p_{T}(t#bar{t}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_njet", jetbins, "n-jets", "Events");
		pdfunc->Add1dHist("pdfunc_nobin", nobins, "total", "Events");
		pdfunc->Add1dHist("pdfunc_dymp", dybins, "y(t)-y(#bar{t})", "Events");
		pdfunc->Add1dHist("pdfunc_dy", dybins, "|y(t)|-|y(#bar{t})|", "Events");
		pdfunc->Add1dHist("pdfunc_ht", htbins, "ht", "Events");
		pdfunc->Add1dHist("pdfunc_evtmass", evtmassbins, "evtmass", "Events");
		pdfunc->Add1dHist("pdfunc_njet+thadpt", response2dvar.GetNBins("njet+thadpt"), 0., response2dvar.GetNBins("njet+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_njet+ttpt", response2dvar.GetNBins("njet+ttpt"), 0., response2dvar.GetNBins("njet+ttpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_njet+ttm", response2dvar.GetNBins("njet+ttm"), 0., response2dvar.GetNBins("njet+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_thady+thadpt", response2dvar.GetNBins("thady+thadpt"), 0., response2dvar.GetNBins("thady+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_ttpt+thadpt", response2dvar.GetNBins("ttpt+thadpt"), 0., response2dvar.GetNBins("ttpt+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_thadpt+ttm", response2dvar.GetNBins("thadpt+ttm"), 0., response2dvar.GetNBins("thadpt+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_ttpt+ttm", response2dvar.GetNBins("ttpt+ttm"), 0., response2dvar.GetNBins("ttpt+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_ttm+tty", response2dvar.GetNBins("ttm+tty"), 0., response2dvar.GetNBins("ttm+tty") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_ttm+dy", response2dvar.GetNBins("ttm+dy"), 0., response2dvar.GetNBins("ttm+dy") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_jet+jetpt", response2dvar.GetNBins("jet+jetpt"), 0., response2dvar.GetNBins("jet+jetpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_jet+jetptOF", response2dvar.GetNBins("jet+jetptOF"), 0., response2dvar.GetNBins("jet+jetptOF") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_jet+jeteta", response2dvar.GetNBins("jet+jeteta"), 0., response2dvar.GetNBins("jet+jeteta") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_jet+jetdr", response2dvar.GetNBins("jet+jetdr"), 0., response2dvar.GetNBins("jet+jetdr") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_jet+jetdrtop", response2dvar.GetNBins("jet+jetdrtop"), 0., response2dvar.GetNBins("jet+jetdrtop") , "bin", "Events");

		pdfunc->Add1dHist("pdfunc_reco_thardpt", topptbins, "hard p_{T}(t) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_tsoftpt", topptbins, "soft p_{T}(t) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_thadpt", topptbins, "p_{T}(t_{h}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_thady", topybins, "|y(t_{h})|", "Events");
		pdfunc->Add1dHist("pdfunc_reco_tleppt", topptbins, "p_{T}(t_{l}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_tlepy", topybins, "|y(t_{l})|", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttm", ttmbins, "M(t#bar{t}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_tty", ttybins, "|y(t#bar{t})|", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttpt", ttptbins, "p_{T}(t#bar{t}) [GeV]", "Events");
		pdfunc->Add1dHist("pdfunc_reco_njet", jetbins, "n-jets", "Events");
		pdfunc->Add1dHist("pdfunc_reco_nobin", nobins, "total", "Events");
		pdfunc->Add1dHist("pdfunc_reco_dymp", dybins, "y(t)-y(#bar{t})", "Events");
		pdfunc->Add1dHist("pdfunc_reco_dy", dybins, "|y(t)|-|y(#bar{t})|", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ht", htbins, "ht", "Events");
		pdfunc->Add1dHist("pdfunc_reco_evtmass", evtmassbins, "evtmass", "Events");
		pdfunc->Add1dHist("pdfunc_reco_njet+thadpt", response2dvar.GetNBins("njet+thadpt"), 0., response2dvar.GetNBins("njet+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_njet+ttpt", response2dvar.GetNBins("njet+ttpt"), 0., response2dvar.GetNBins("njet+ttpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_njet+ttm", response2dvar.GetNBins("njet+ttm"), 0., response2dvar.GetNBins("njet+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_thady+thadpt", response2dvar.GetNBins("thady+thadpt"), 0., response2dvar.GetNBins("thady+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttpt+thadpt", response2dvar.GetNBins("ttpt+thadpt"), 0., response2dvar.GetNBins("ttpt+thadpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_thadpt+ttm", response2dvar.GetNBins("thadpt+ttm"), 0., response2dvar.GetNBins("thadpt+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttpt+ttm", response2dvar.GetNBins("ttpt+ttm"), 0., response2dvar.GetNBins("ttpt+ttm") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttm+tty", response2dvar.GetNBins("ttm+tty"), 0., response2dvar.GetNBins("ttm+tty") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_ttm+dy", response2dvar.GetNBins("ttm+dy"), 0., response2dvar.GetNBins("ttm+dy") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_jet+jetpt", response2dvar.GetNBins("jet+jetpt"), 0., response2dvar.GetNBins("jet+jetpt") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_jet+jetptOF", response2dvar.GetNBins("jet+jetptOF"), 0., response2dvar.GetNBins("jet+jetptOF") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_jet+jeteta", response2dvar.GetNBins("jet+jeteta"), 0., response2dvar.GetNBins("jet+jeteta") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_jet+jetdr", response2dvar.GetNBins("jet+jetdr"), 0., response2dvar.GetNBins("jet+jetdr") , "bin", "Events");
		pdfunc->Add1dHist("pdfunc_reco_jet+jetdrtop", response2dvar.GetNBins("jet+jetdrtop"), 0., response2dvar.GetNBins("jet+jetdrtop") , "bin", "Events");
	}

	ttp_truth.Init(this);
	ttp_right.Init(this);
	ttp_wrong.Init(this);
	ttp_semi.Init(this);
	ttp_other.Init(this);

//	ttp_hadjets_right.Init(this);
//	ttp_hadjets_wrong.Init(this);
//	ttp_jets_right.Init(this);
//	ttp_jets_wrong.Init(this);
//	ttp_blep_right.Init(this);
//	ttp_blep_wrong.Init(this);
//	ttp_whad_right.Init(this);
//	ttp_whad_wrong.Init(this);
//	ttp_jetspos_right.Init(this);
//	ttp_jetspos_wrong.Init(this);

	ttp_tlepthad_right.Init(this);
	ttp_tlep_right.Init(this);
	ttp_thad_right.Init(this);
	ttp_nn_right.Init(this);
	ttp_nsemi_right.Init(this);
	ttp_alljets_right.Init(this);
	ttp_alljets_wrong.Init(this);
	for(int j = 0 ; j < 8 ; ++j)
	{
		stringstream ss;
		ss << "jet_" << j << "_";
		truth1d.AddHist(ss.str()+"pt_fb_right", jetptbins[j], "p_{T}("+jetnames[j]+") [GeV]", "Events");
		truth1d.AddHist(ss.str()+"eta_fb_right", jetetabins[j], "|#eta("+jetnames[j]+")|", "Events");
		truth1d.AddHist(ss.str()+"dr_fb_right", jetdrbins[j], "#DeltaR_{j_{t}}("+jetnames[j]+")", "Events");
		truth1d.AddHist(ss.str()+"drtop_fb_right", jetdrtopbins[j], "#DeltaR_{t}("+jetnames[j]+")", "Events");
		truth1d.AddHist(ss.str()+"pt_fb_wrong", jetptbins[j], "p_{T}("+jetnames[j]+") [GeV]", "Events");
		truth1d.AddHist(ss.str()+"eta_fb_wrong", jetetabins[j], "|#eta("+jetnames[j]+")|", "Events");
		truth1d.AddHist(ss.str()+"dr_fb_wrong", jetdrbins[j], "#DeltaR_{j_{t}}("+jetnames[j]+")", "Events");
		truth1d.AddHist(ss.str()+"drtop_fb_wrong", jetdrtopbins[j], "#DeltaR_{t}("+jetnames[j]+")", "Events");
	}

	double bkgcutmin = 0.0;
	double bkgcutmax = 0.6;
	if(cnbtag == "bkgl") {bkgcutmin = 0.0; bkgcutmax = 0.3;}
	if(cnbtag == "bkgh") {bkgcutmin = 0.4; bkgcutmax = 0.7;}
	btageff.Init(B_MEDIUM, bkgcutmin, bkgcutmax);

	TDirectory* dir_reco = outFile_.mkdir("RECO");
	dir_reco->cd();
	reco1d.AddHist("counter", 20, 0., 20., "counter", "Events");
	reco1d.AddHist("counter_test", 20, 0., 20., "counter", "Events");
	reco2d.AddHist("Wmasshad_tmasshad", 500, 0., 500., 500, 0., 500, "M(W) [GeV]", "M(t) [GeV]");
	reco1d.AddHist("btag_high", 100, 0., 1., "btag", "Events");
	reco1d.AddHist("btag_low", 100, 0., 1., "btag", "Events");
	reco1d.AddHist("bjetmulti", 10, 0., 10., "b-jets", "Events");
	reco1d.AddHist("jetmulti", 20, 0., 20., "jets", "Events");
	reco1d.AddHist("bjetmultiW", 10, 0., 10., "b-jets", "Events");
	reco1d.AddHist("btagtest", 1000, -100, 100., "-Log(p) btag-test", "Events");
	reco1d.AddHist("masstest", 1000, -100, 100., "-Log(p) mass-test", "Events");
	reco1d.AddHist("nstest", 200, 0, 20., "neutrino-test", "Events");
	reco1d.AddHist("MuIsolation", 200, 0, 2, "rel. iso #mu", "Events");
	reco1d.AddHist("ElIsolation", 200, 0, 2, "rel. iso el", "Events");
	reco1d.AddHist("NumVertices", 200, 0, 200, "num vertices", "Events");
	reco1d.AddHist("NumVerticesWeighted", 200, 0, 200, "num vertices", "Events");
	reco1d.AddHist("Rho", 100, 0, 100, "#rho", "Events");
	reco1d.AddHist("RhoWeighted", 100, 0, 100, "#rho", "Events");
	reco1d.AddHist("drmu", 600, 0, 6, "#DeltaR(#mu#mu)", "Events");
	reco1d.AddHist("drel", 600, 0, 6, "#DeltaR(#mu#mu)", "Events");
	reco1d.AddHist("pupt_tight", 400, -2, 2, "pu fraction", "Events");
	reco1d.AddHist("pupt_loose", 400, -2, 2, "pu fraction", "Events");
	reco1d.AddHist("allpt_tight", 400, -2, 2, "pu fraction", "Events");
	reco1d.AddHist("allpt_loose", 400, -2, 2, "pu fraction", "Events");
	//reco1d.AddHist("Mt_W", 200, 0, 200, "M_{t}(W) [GeV]", "Events");
	ttp_all.Init(this);

	jetscaler.Init(cJetEnergyUncertainty, cjecuncertainty);
	jetscaler.Init(cJetEnergyUncertainty, "FlavorPureBottom");
	jetscaler.Init(cJetEnergyUncertainty, "FlavorQCD");
	//jetscaler.Init("Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt", cjecuncertainty);
	jetscaler.InitResolution(cJetResolution, cJetResolutionSF);
	jetscaler.InitMCrescale(this, "jetrescale.root");
	string probfilename("LH_parton.root");
	if(PSEUDOTOP){probfilename = "LH_pseudo.root";}
	if(BTAGMODE)
	{
		ttsolver.Init(PSEUDOTOP, probfilename, false, true, true);//for btag
	}
	if(JETSCALEMODE)
	{
		ttsolver.Init(PSEUDOTOP, probfilename, true, true, false);//don't use mass info
	}
	else
	{
		ttsolver.Init(PSEUDOTOP, probfilename, false, true, true);
		//ttsolver.Init(probfilename, false, true);
	}
	btagweight.Init(this, cBTaggingSF, cBTaggingEff, cbtagunc, cltagunc);
	string puhistname("pu_central");
	if(cpileup == -1) puhistname = "pu_minus";
	if(cpileup == 1) puhistname = "pu_plus";

	TFile* f = TFile::Open("PUweight.root");
	puhist = dynamic_cast<TH1D*>(f->Get(puhistname.c_str()));
	TFile* fl = TFile::Open(cLeptonScaleFactor.c_str());
	musfhist = dynamic_cast<TH2D*>(fl->Get("MuSF"));
	elsfhist = dynamic_cast<TH2D*>(fl->Get("ElSF"));

	bdecayweights.Init(cbdecay);
	bfragweights.Init("bfragweights.root", cbfrag);

}

ttbar::~ttbar()
{
	if(STUDENT)
	{
		stud_tf->Write();
		stud_tf->Close();
	}
}

void ttbar::SelectGenParticles(URStreamer& event)
{
	SEMILEP = false;
	SEMILEPACC = false;
	FULLLEP = false;
	FULLHAD = false;

	const vector<Ttgen>& ttgens = event.TTGens();
	if(ttgens.size() != 8) { return; }
	for(size_t n = 0 ; n < 8 ; ++n) { gps[n] = ttgens[n]; }
	gen1d["TYP"]->Fill(0.5, weight);
	gentq = gps[0];
	gentqbar = gps[4];

	weight *= 1.+cttptweight*((gentq + gentqbar).Pt()-200.)/2000.;
	weight *= 1.+ctopptweight*(gentq.Pt()-200.)/1000.;
	weight *= 1.+ctoprapweight*(0.2-0.2*Abs(gentq.Rapidity()));
	bool tlep = false;
	bool tbarlep = false;
	bool tem = false;
	bool tbarem = false;
	if(gps[2].pdgId() == -11 || gps[2].pdgId() == -13){tem = true;}
	if(gps[7].pdgId() == 11 || gps[7].pdgId() == 13){tbarem = true;}
	if(gps[2].pdgId() == -11 || gps[2].pdgId() == -13 || gps[2].pdgId() == -15){tlep = true;}
	if(gps[7].pdgId() == 11 || gps[7].pdgId() == 13 || gps[7].pdgId() == 15){tbarlep = true;}

	if(STUDENT)
	{
		if(tlep || tbarlep) num_gen = 2;
		gen_px[0] = gentq.Px();
		gen_py[0] = gentq.Py();
		gen_pz[0] = gentq.Pz();
		gen_e[0] = gentq.E();
		gen_type[0] = gps[2].pdgId();
		gen_px[1] = gentqbar.Px();
		gen_py[1] = gentqbar.Py();
		gen_pz[1] = gentqbar.Pz();
		gen_e[1] = gentqbar.E();
		gen_type[1] = gps[7].pdgId();
	}

	if(tem && !tbarlep)
	{
		SEMILEPACC = true;
		gen1d["TYP"]->Fill(4.5, weight);
		//genallper.Init(&gps[6], &gps[7], &gps[5], &gps[1], &gps[2], gps[2].pdgId(), gps[3]);
		gentqhad = gentqbar;
		gentqlep = gentq;

	}
	else if(!tlep && tbarem)
	{
		SEMILEPACC = true;
		gen1d["TYP"]->Fill(4.5, weight);
		//genallper.Init(&gps[2], &gps[3], &gps[1], &gps[5], &gps[7], gps[7].pdgId(), gps[6]);
		gentqhad = gentq;
		gentqlep = gentqbar;
	}

	if(tlep && tbarlep)
	{
		FULLLEP = true;
		gen1d["TYP"]->Fill(3.5, weight);
	}
	else if(!tlep && !tbarlep)
	{
		FULLHAD = true;
		gen1d["TYP"]->Fill(1.5, weight);
	}
	else 
	{
		SEMILEP = true;
		gen1d["TYP"]->Fill(2.5, weight);
	}

}

void ttbar::SelectRivetPS(URStreamer& event)
{
	const vector<Rivetpl>& rivgens = event.RivetPLs();
	if(rivgens.size() == 0) { return; }
	for(size_t n = 0 ; n < 8 ; ++n){rivetobjs[n] = rivgens[n];}

	if(rivetobjs[2].pdgId() == -11 || rivetobjs[2].pdgId() == -13)
	{
		psper.Init(&rivetobjs[6], &rivetobjs[7], &rivetobjs[5], &rivetobjs[1], &rivetobjs[2], rivetobjs[2].pdgId(), rivetobjs[3]);
	}
	else if(rivetobjs[6].pdgId() == 11 || rivetobjs[6].pdgId() == 13)
	{
		psper.Init(&rivetobjs[2], &rivetobjs[3], &rivetobjs[1], &rivetobjs[5], &rivetobjs[6], rivetobjs[6].pdgId(), rivetobjs[7]);
	}
	vector<TLorentzVector*> genaddjets;
	for(size_t j = 8 ; j < rivgens.size() ;  ++j)
	{
		const Rivetpl& rp = rivgens[j]; 	
		if(abs(rp.pdgId()) < 6)
		{
			sgenparticles.push_back(rp);
			genaddjets.push_back(&(sgenparticles.back()));
		}
	}
	psper.SetAdditionalJets(genaddjets, [](TLorentzVector* j){return true;});
}

void ttbar::SelectPseudoTopLHC(URStreamer& event)
{
	const vector<Pl>& pls = event.PLs();
	GenObject* lepton = nullptr;
	int lc = 0;
	TLorentzVector nu;
	vector<TLorentzVector*> pstbjets;
	vector<TLorentzVector*> pstljets;
	vector<TLorentzVector*> pstphotons;
	for(const Pl& pl : pls)
	{
		if(abs(pl.pdgId()) < 6) continue;

		if(abs(pl.pdgId()) == 12 || abs(pl.pdgId()) == 14 || abs(pl.pdgId()) == 16)
		{
			nu += pl;
			continue;
		}

		if((abs(pl.pdgId()) == 11 || abs(pl.pdgId()) == 13) && pl.isoR3() < 0.2)
		{
			if(Abs(pl.Eta()) < 2.4 && pl.Pt() > 15)
			{
				lc++;
				if(lc == 2) {return;}
			}
			if(Abs(pl.Eta()) < cpletamax && pl.Pt() > cplptmin)
			{
				sgenparticles.push_back(pl);
				lepton = &(sgenparticles.back());
			}
		}
	}

	if(lepton == nullptr) {return;}

	for(const Pl& pl : pls)
	{
		if(abs(pl.pdgId()) > 6) continue;
		if(pl.Pt() < min({cpbjetptsoft, cpwjetptsoft}) || Abs(pl.Eta()) > cpjetetamax) continue;
		if(lepton->DeltaR(pl) < 0.4) {continue;}
		sgenparticles.push_back(pl);
		if(abs(pl.pdgId()) == 5)
		{
			pstbjets.push_back(&(sgenparticles.back()));	
		}
		else
		{
			pstljets.push_back(&(sgenparticles.back()));
		}
	}

	if(pstbjets.size() < 2 || pstljets.size() < 2){return;}

	double mt = sqrt(2*lepton->Pt()*nu.Pt()*(1.-cos(lepton->DeltaPhi(nu))));
	
	gen2d["met_mt"]->Fill(nu.Pt(), mt);
	//if(mt < 30. || nu.Pt() < 30.) {return;}

	sort(pstbjets.begin(), pstbjets.end(), sortpt<TLorentzVector>);
	sort(pstljets.begin(), pstljets.end(), sortpt<TLorentzVector>);
	double Mw = 80.4;
	double A = Mw*Mw*0.5 + nu.Px()*lepton->Px() + nu.Py()*lepton->Py();
	double N = lepton->Pz()*lepton->Pz() - lepton->E()*lepton->E();
	double p = -1.*A*lepton->Pz()/N;
	double q = (A*A - lepton->E()*lepton->E()*nu.Pt()*nu.Pt())/N;

	double pz = p;
	if(p*p - q > 0.)
	{
		if(p > 0) {pz += sqrt(p*p - q);}
		else {pz -= sqrt(p*p - q);}
	}
	nu.SetXYZM(nu.Px(), nu.Py(), pz, 0.);

	int bh = 0;
	int bl = 1;
	if(lepton->DeltaR(*pstbjets[0]) < lepton->DeltaR(*pstbjets[1])) {bl = 0; bh = 1;}
	psper.Init(pstljets[0], pstljets[1], pstbjets[bh], pstbjets[bl], lepton, lepton->pdgId(), nu);
	
	vector<TLorentzVector*> genaddjets;
	for(TLorentzVector* gj : pstbjets)
	{
		if(psper.IsJetIn(gj) == -1)
		{
			genaddjets.push_back(gj);
		}
	}
	for(TLorentzVector* gj : pstljets)
	{
		if(psper.IsJetIn(gj) == -1)
		{
			genaddjets.push_back(gj);
		}
	}
	psper.SetAdditionalJets(genaddjets, [](TLorentzVector* j){return j->Pt() > 30.;});
}

void ttbar::SelectPseudoTop(URStreamer& event)
{
	// const vector<Pl>& pls = event.PLs();

	// GenObject* lepton = nullptr;
	// int lc = 0;
	// TLorentzVector nu;
	// vector<GenObject*> pstbjets;
	// vector<GenObject*> pstljets;
	// vector<GenObject*> pstphotons;
	// for(const Pl& pl : pls)
	// {
	// 	if(abs(pl.pdgId()) < 6) continue;

	// 	if(abs(pl.pdgId()) == 12 || abs(pl.pdgId()) == 14 || abs(pl.pdgId()) == 16)
	// 	{
	// 		nu += pl;
	// 		continue;
	// 	}

	// 	if(pl.pdgId() == 22)
	// 	{
	// 		if(Abs(pl.Eta()) > 2.4 || pl.Pt() < 15.) continue;
	// 		if(pl.isoR3() < 0.2)
	// 		{
	// 			sgenparticles.push_back(pl);
	// 			pstphotons.push_back(&(sgenparticles.back()));
	// 		}
	// 		continue;
	// 	}

	// 	if(abs(pl.pdgId()) == 13 || abs(pl.pdgId()) == 11)
	// 	{
	// 		if(Abs(pl.Eta()) > 2.4 || pl.Pt() < 15) continue;

	// 		if(pl.isoR3() < 0.35)
	// 		{
	// 			lc++;
	// 			if(Abs(pl.Eta()) < cpletamax && pl.Pt() > cplptmin)
	// 			{
	// 				sgenparticles.push_back(pl);
	// 				lepton = &(sgenparticles.back());
	// 			}
	// 		}
	// 		continue;
	// 	}
	// }

	// if(lc != 1 || lepton == nullptr) {return;}

	// for(const Pl& pl : pls)
	// {
	// 	if(abs(pl.pdgId()) > 6) continue;
	// 	if(pl.Pt() < min({cpbjetptsoft, cpwjetptsoft}) || Abs(pl.Eta()) > cpjetetamax) continue;
	// 	if(lepton->DeltaR(pl) < 0.4) {continue;}
	// 	bool isph = false;
	// 	for(TLorentzVector* ph : pstphotons) {if(ph->DeltaR(pl) < 0.4) {isph = true; break;}}
	// 	if(isph) continue;
	// 	sgenparticles.push_back(pl);
	// 	if(abs(pl.pdgId()) == 5)
	// 	{
	// 		pstbjets.push_back(&(sgenparticles.back()));	
	// 	}	
	// 	else
	// 	{
	// 		pstljets.push_back(&(sgenparticles.back()));	
	// 	}
	// }

	// if(pstbjets.size() < 2 || pstljets.size() < 2){return;}
	
	// TLorentzVector wl(*lepton + nu);
	// double chi2min = 1.E100;
	// double MW = 80.4;
	// double Mt = 172.5;
	// for(int wa = 0 ; wa < pstljets.size() ; wa++)
	// {
	// 	for(int wb = 0 ; wb < wa ; wb++)
	// 	{
	// 		if((pstljets[wa]->Pt() < cpwjetpthard && pstljets[wa]->Pt() < cpwjetpthard)){continue;}
	// 		TLorentzVector wh(*pstljets[wa] + *pstljets[wb]);
	// 		for(int bl = 0 ; bl < pstbjets.size() ; bl++)
	// 		{
	// 			if(pstbjets[bl] == pstljets[wa] || pstbjets[bl] == pstljets[wb]) {continue;}
	// 			for(int bh = 0 ; bh < pstbjets.size() ; bh++)
	// 			{
	// 				if(pstbjets[bh] == pstljets[wa] || pstbjets[bh] == pstljets[wb] || pstbjets[bh] == pstbjets[bl]){continue;}
	// 				if(pstbjets[bh]->Pt() < cpbjetpthard && pstbjets[bl]->Pt() < cpbjetpthard) {continue;}
	// 				TLorentzVector th(wh + *pstbjets[bh]);
	// 				TLorentzVector tl(wl + *pstbjets[bl]);
	// 				double chi2 = Power(th.M() - Mt, 2) + Power(tl.M() - Mt, 2) + Power(wh.M() - MW, 2);
	// 				if(chi2 < chi2min)
	// 				{
	// 					chi2min = chi2;
	// 					psper.Init(pstljets[wa], pstljets[wb], pstbjets[bh], pstbjets[bl], lepton, lepton->pdgId(), nu);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// vector<TLorentzVector*> genaddjets;
	// for(GenObject* gj : pstbjets)
	// {
	// 	if(psper.IsJetIn(gj) == -1)
	// 	{
	// 		genaddjets.push_back(gj);
	// 	}
	// }
	// for(GenObject* gj : pstljets)
	// {
	// 	if(psper.IsJetIn(gj) == -1)
	// 	{
	// 		genaddjets.push_back(gj);
	// 	}
	// }
	// psper.SetAdditionalJets(genaddjets, [](TLorentzVector* j){return j->Pt() > 30.;});
}

void ttbar::AddGenJetSelection(URStreamer& event)
{
	const vector<Pl>& pls = event.PLs();
	int addbs = 0;	
	int addcs = 0;	
	double bweights = 1.;
	for(const Pl& gj : pls)
	{
		if(abs(gj.pdgId()) == 13 || abs(gj.pdgId()) == 11)
		{
			if(Abs(gj.Eta()) > 2.4 || gj.Pt() < 30) continue;
			sgenparticles.push_back(gj);
			plleptons.push_back(&(sgenparticles.back()));
		}
		else if(gj.pdgId() == 22)
		{
			if(Abs(gj.Eta()) > 2.4 || gj.Pt() < 15) continue;
			sgenparticles.push_back(gj);
			plphotons.push_back(&(sgenparticles.back()));
		}

		if(abs(gj.pdgId()) > 5) {continue;}
		sgenparticles.push_back(gj);
		genalljets.push_back(&(sgenparticles.back()));
		if(gj.pdgId() == 5)
		{
			genbjets.push_back(&(sgenparticles.back()));
			bweights *= bfragweights.Weight(gj);
			bweights *= bdecayweights.Weight(gj);
			if(genallper.IsComplete() && genallper.DRminTTjets(&gj) > 0.4) {addbs++;}
		}
		else if(gj.pdgId() == 4)
		{
			gencjets.push_back(&(sgenparticles.back()));
			if(genallper.IsComplete() && genallper.DRminTTjets(&gj) > 0.4) {addcs++;}
		}
		else
		{
			genljets.push_back(&(sgenparticles.back()));
		}
	}
	weight *= bweights;
	for(const Pl& gj : pls)
	{
		if(abs(gj.pdgId()) == 4)
		{
			truth1d["xc"]->Fill(gj.isoR4(), weight);
		}
		else if(abs(gj.pdgId()) == 5)
		{
			truth1d["xb"]->Fill(gj.isoR4(), weight);
			float pdgid = gj.isoR3();
			if(fmod(pdgid,1) != 0.)
			{
				truth1d["bids_lep"]->Fill(abs(pdgid-0.5), weight);
			}
			else
			{
				truth1d["bids_had"]->Fill(abs(pdgid), weight);
			}
		}
	}
	weight *= pow(cbsplitting, addbs);
	weight *= pow(ccsplitting, addbs);
}

void ttbar::SelectRecoParticles(URStreamer& event)
{
//Muons
	const vector<Muon>& muons = event.muons();
	for(vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon)
	{
		IDMuon mu(*muon);
		if(genper->IsComplete() && mu.DeltaR(*genper->L()) < 0.3)
		{
			double eta = Abs(mu.Eta());
			if(eta < 0.8){truth2d["Murho_iso_1"]->Fill(event.rho().value(), mu.pfNeutralIso04() + mu.pfPhotonIso04(), weight);}
			else if(eta < 1.5){truth2d["Murho_iso_2"]->Fill(event.rho().value(), mu.pfNeutralIso04() + mu.pfPhotonIso04(), weight);}
			else if(eta < 2.1){truth2d["Murho_iso_3"]->Fill(event.rho().value(), mu.pfNeutralIso04() + mu.pfPhotonIso04(), weight);}
			else if(eta < 2.4){truth2d["Murho_iso_4"]->Fill(event.rho().value(), mu.pfNeutralIso04() + mu.pfPhotonIso04(), weight);}
		}
		if(mu.ID(IDMuon::LOOSE_16) && mu.Pt() > 15.)
		{
			smuons.push_back(mu);
			loosemuons.push_back(&(smuons.back()));
			if(MUONS && mu.ID(IDMuon::TIGHT_16) && mu.Pt() > clptmin && Abs(mu.Eta()) < cletamax)
			{
				tightmuons.push_back(&(smuons.back()));
			}
		}
	}
	if(loosemuons.size() == 2) { reco1d["drmu"]->Fill(loosemuons[0]->DeltaR(*loosemuons[1]));}

	for(GenObject* pllep : plleptons)
	{
		if(abs(pllep->pdgId()) != 13) continue;
		gen2d["plmuiso"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		if(pllep->fromhadron() == false) gen2d["plmuisoprompt"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		if(TTMC && ((abs(gps[2].pdgId()) == 13 && pllep->DeltaR(gps[2]) < 0.1) || (abs(gps[7].pdgId()) == 13 && pllep->DeltaR(gps[7]) < 0.1))) 
		{
			gen2d["plmuisott"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		}
		for(IDMuon* mu : tightmuons)
		{
			if(mu->DeltaR(*pllep) < 0.3)
			{
				gen2d["plmuisoreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				if(pllep->fromhadron() == false) gen2d["plmuisopromptreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				if(TTMC && ((abs(gps[2].pdgId()) == 13 && pllep->DeltaR(gps[2]) < 0.1) || (abs(gps[7].pdgId()) == 13 && pllep->DeltaR(gps[7]) < 0.1))) gen2d["plmuisottreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				break;
			}
		}
	}

//Electrons
	const vector<Electron>& electrons = event.electrons();
	for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron)
	{
		IDElectron el(*electron);
		if(genper->IsComplete() && el.DeltaR(*genper->L()) < 0.3)
		{
			double eta = abs(el.SCeta());
			if(eta < 0.8){truth2d["Elrho_iso_1"]->Fill(event.rho().value(), el.pfNeutralIso() + el.pfPhotonIso(), weight);}
			else if(eta < 1.4442){truth2d["Elrho_iso_2"]->Fill(event.rho().value(), el.pfNeutralIso() + el.pfPhotonIso(), weight);}
			else if(eta < 1.566){}
			else if(eta < 2.1){truth2d["Elrho_iso_3"]->Fill(event.rho().value(), el.pfNeutralIso() + el.pfPhotonIso(), weight);}
			else if(eta < 2.5){truth2d["Elrho_iso_4"]->Fill(event.rho().value(), el.pfNeutralIso() + el.pfPhotonIso(), weight);}
		}
		if(el.ID(IDElectron::LOOSE_16) && el.Pt() > 15.)
		{
			selectrons.push_back(el);
			looseelectrons.push_back(&(selectrons.back()));
			if(ELECTRONS && el.ID(IDElectron::TIGHT_16) && el.Pt() > clptmin && Abs(el.Eta()) < cletamax)
			{
				mediumelectrons.push_back(&(selectrons.back()));
			}
		}
	}
	if(looseelectrons.size() == 2) {reco1d["drel"]->Fill(looseelectrons[0]->DeltaR(*looseelectrons[1]));}

	for(GenObject* pllep : plleptons)
	{
		if(abs(pllep->pdgId()) != 11) continue;
		gen2d["pleliso"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		if(pllep->fromhadron() == false) gen2d["plelisoprompt"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		if(TTMC && ((abs(gps[2].pdgId()) == 11 && pllep->DeltaR(gps[2]) < 0.1) || (abs(gps[7].pdgId()) == 11 && pllep->DeltaR(gps[7]) < 0.1))) gen2d["plelisott"]->Fill(pllep->Pt(), pllep->isolation(), weight);
		for(IDElectron* el : mediumelectrons)
		{
			if(el->DeltaR(*pllep) < 0.3)
			{
				gen2d["plelisoreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				if(pllep->fromhadron() == false) gen2d["plelisopromptreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				if(TTMC && ((abs(gps[2].pdgId()) == 11 && pllep->DeltaR(gps[2]) < 0.1) || (abs(gps[7].pdgId()) == 11 && pllep->DeltaR(gps[7]) < 0.1))) gen2d["plelisottreco"]->Fill(pllep->Pt(), pllep->isolation(), weight);
				break;
			}
		}
	}

	TVector2 metcorr(0.,0.);
	const vector<Jet>& jets = event.jets();
	for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
	{
		IDJet jet(*jetit);
		if(!jet.ID() || !jet.Clean(loosemuons, looseelectrons)) {continue;}
		double sfres = jetscaler.GetRes(jet, event.rho().value(), cjetres);
		double sfgen = sfres * jetscaler.GetScale(jet, csigmajet);
		double sfwj = 0.;
		if(csigmajetwj != 1000.) {sfwj = sfres * jetscaler.GetScale(jet, csigmajetwj);}
		else{ sfwj = sfres * (1. + cscalejetwj);}
		jet.SetSF(1, sfres, metcorr);
		jet.SetSF(2, sfgen, metcorr);
		jet.SetSF(3, sfwj, metcorr);
		if(jet.Pt() < jetptmin || Abs(jet.Eta()) > cjetetamax) {continue;}

		sjets.push_back(jet);
		cleanedjets.push_back(&(sjets.back()));
	}

	for(IDJet* j:cleanedjets)
	{
		j->ApplySF(2, metcorr);
	}
	//const vector<Nohfmet>& mets = event.NoHFMETs();
	const vector<Met>& mets = event.METs();
	if(mets.size() == 1)
	{
		met = mets[0];
		met.SetPx(met.Px() + csigmamet*met.pxunc());
		met.SetPy(met.Py() + csigmamet*met.pyunc());
		met.Update(metcorr);
	}

	for(GenObject* plph : plphotons)
	{
		gen2d["plphiso"]->Fill(plph->Pt(), plph->isolation(), weight);
		if(!plph->fromhadron()){gen2d["plphisoprompt"]->Fill(plph->Pt(), plph->isolation(), weight);}
		for(IDJet* j:cleanedjets)
		{
			if(j->DeltaR(*plph) < 0.3)
			{
				gen2d["plphisoreco"]->Fill(plph->Pt(), plph->isolation(), weight);
				break;
			}
		}
	}

	for(IDJet* j:cleanedjets)
	{
		if(genalljets.size() == 0) continue;
		GenObject* genj = *min_element(genalljets.begin(), genalljets.end(), [&](TLorentzVector* A, TLorentzVector* B){return j->DeltaR(*A) < j->DeltaR(*B);});
		if(j->DeltaR(*genj) > 0.3) continue;
		double dpt = (j->Pt() - genj->Pt())/genj->Pt();
		if(genj->pdgId() == 5)
		{
			if(Abs(j->Eta()) < 1.5)
			{
				truth2d["jetscale_bB"]->Fill(j->Pt(), dpt, weight);
			}
			else
			{
				truth2d["jetscale_bE"]->Fill(j->Pt(), dpt, weight);
			}					
		}
		else
		{
			if(Abs(j->Eta()) < 1.5)
			{
				truth2d["jetscale_lB"]->Fill(j->Pt(), dpt, weight);
			}
			else
			{
				truth2d["jetscale_lE"]->Fill(j->Pt(), dpt, weight);
			}					
		}
	}

	if(SEMILEPACC)
	{
		rightper.MET(&met);
		for(IDElectron* el : mediumelectrons)
		{
			if(el->DeltaR(*genper->L()) < 0.2)
			{
				rightper.L(el, el->charge()*-13);
				truth1d["found"]->Fill(4.5, weight);
			}
		}

		for(IDMuon* mu : tightmuons)
		{
			if(mu->DeltaR(*genper->L()) < 0.2)
			{
				rightper.L(mu, mu->charge()*-11);
				truth1d["found"]->Fill(5.5, weight);
			}
		}

		metcorr = TVector2(0.,0.);
		bool wa = false;
		bool wb = false;
		bool ba = false;
		bool bb = false;
		double ptbhadmax = 0.;
		double ptblepmax = 0.;
		double wjaptmax = 0.;
		double wjbptmax = 0.;
		for(size_t j = 0 ; j < cleanedjets.size() ; ++j)
		{
			IDJet* jet = cleanedjets[j];
			if(jet->DeltaR(*genper->BHad()) < 0.2 && jet->Pt() > ptbhadmax)
			{
				ptbhadmax = jet->Pt();
				rightper.BHad(jet);
				if(!ba){ba = true; truth1d["found"]->Fill(0.5, weight);}
			}
			if(jet->DeltaR(*genper->BLep()) < 0.2 && jet->Pt() > ptblepmax)
			{
				ptblepmax = jet->Pt();
				rightper.BLep(jet);
				if(!bb){bb = true; truth1d["found"]->Fill(1.5, weight);}
			}
			if(jet->DeltaR(*genper->WJa()) < 0.2 && jet->Pt() > wjaptmax)
			{
				wjaptmax = jet->Pt();
				rightper.WJa(jet);
				if(!wa){wa = true; truth1d["found"]->Fill(2.5, weight);}
			}
			if(jet->DeltaR(*genper->WJb()) < 0.2 && jet->Pt() > wjbptmax)
			{
				wjbptmax = jet->Pt();
				rightper.WJb(jet);
				if(!wb){wb = true; truth1d["found"]->Fill(3.5, weight);}
			}
		}

		if(genallper.IsComplete() && rightper.IsComplete())
		{
			TVector2 metcorrwj(0.,0.);
			if(rightper.WJa()->DeltaR(*genallper.WJa()) < 0.4 || rightper.WJa()->DeltaR(*genallper.WJb()) < 0.4 ) {dynamic_cast<IDJet*>(rightper.WJa())->ApplySF(3, metcorrwj);}
			if(rightper.WJb()->DeltaR(*genallper.WJa()) < 0.4 || rightper.WJb()->DeltaR(*genallper.WJb()) < 0.4 ){dynamic_cast<IDJet*>(rightper.WJb())->ApplySF(3, metcorrwj);}
			//dynamic_cast<IDJet*>(rightper.WJa())->ApplySF(3, metcorrwj);
			//dynamic_cast<IDJet*>(rightper.WJb())->ApplySF(3, metcorrwj);
			met.Update(metcorrwj);
		} 
	}
}

void ttbar::ttanalysis(URStreamer& event)
{
// 	reco1d["counter"]->Fill(0.5, weight);

// 	if(SEMILEPACC && rightper.IsComplete()){ttp_truth.Fill(rightper, weight);}

// 	if(SEMILEPACC)
// 	{
// 		if(Abs(genper->LPDGId()) == 11) {truth2d["Ne_Nmu"]->Fill(mediumelectrons.size()+0.5, tightmuons.size()+0.5, weight);}
// 		if(Abs(genper->LPDGId()) == 13) {truth2d["Nmu_Ne"]->Fill(tightmuons.size()+0.5, mediumelectrons.size()+0.5, weight);}
// 	}

// 	//check for lepton:
// 	if(tightmuons.size() == 1 && loosemuons.size() == 1 && looseelectrons.size() == 0)
// 	{
// 		lep = dynamic_cast<TLorentzVector*>(tightmuons[0]);
// 		leppdgid = tightmuons[0]->charge()*-13;
// 	}
// 	if(loosemuons.size() == 0 && mediumelectrons.size() == 1 && looseelectrons.size() == 1)
// 	{
// 		lep = dynamic_cast<TLorentzVector*>(mediumelectrons[0]);
// 		leppdgid = mediumelectrons[0]->charge()*-11;
// 	}

// 	if(lep == nullptr){return;}
// 	weight *= lepeffweight(lep, event);

// 	if(cleanedjets.size() < 3){return;}

// 	reco1d["counter"]->Fill(1.5, weight);


// 	sort(cleanedjets.begin(), cleanedjets.end(), [](IDJet* A, IDJet* B){return(A->csvIncl() > B->csvIncl());});
// 	int nbjets = count_if(cleanedjets.begin(), cleanedjets.end(), [&](IDJet* A){return(A->csvIncl() > B_MEDIUM);});
// 	reco1d["bjetmulti"]->Fill(nbjets, weight);
// 	if(isDA == 0 && !BTAGMODE)
// 	{
// 		double btw = btagweight.SF(cleanedjets);
// //cout << weight << " " << btw << endl;
// 		weight *= btw;
// 	}
// 	reco1d["bjetmultiW"]->Fill(nbjets, weight);

// 	if(SEMILEPACC)
// 	{
// 		truth1d["counter"]->Fill(3.5, weight);
// 	}
// 	//jet number plots
// 	if(SEMILEPACC)
// 	{
// 		//if(lep != rightper.L()) {cout << "Wrong Lep" << endl;}
// 		truth2d["Jetstt_JetsAll"]->Fill(rightper.NumTTBarJets()+0.5, cleanedjets.size()+0.5, weight);
// 	}
// //cout << "NC: " << cleanedjets.size() << endl;
// 	//plot b-tag distribution
// 	if(rightper.IsComplete())
// 	{
// 		truth2d["btag2d_true"]->Fill(Min(dynamic_cast<IDJet*>(rightper.BHad())->csvIncl(), dynamic_cast<IDJet*>(rightper.BLep())->csvIncl()), Max(dynamic_cast<IDJet*>(rightper.BHad())->csvIncl(), dynamic_cast<IDJet*>(rightper.BLep())->csvIncl()), weight);
// 		truth1d["btag_true"]->Fill(dynamic_cast<IDJet*>(rightper.BHad())->csvIncl(), weight);
// 		truth1d["btag_true"]->Fill(dynamic_cast<IDJet*>(rightper.BLep())->csvIncl(), weight);
// 		truth1d["btag_wrong"]->Fill(dynamic_cast<IDJet*>(rightper.WJa())->csvIncl(), weight);
// 		truth1d["btag_wrong"]->Fill(dynamic_cast<IDJet*>(rightper.WJb())->csvIncl(), weight);
// 	}

// 	reco1d["jetmulti"]->Fill(cleanedjets.size(), weight);

// 	reco1d["counter"]->Fill(2.5, weight);
// 	if(cleanedjets.size() < 4){return;}
// 	if(BTAGMODE && cleanedjets.size() > 4){return;}
// 	reco1d["counter"]->Fill(3.5, weight);
// 	if(SEMILEPACC) truth1d["counter"]->Fill(4.5, weight);
// 	if(tightmuons.size() == 1)
// 	{
// 		reco1d["MuIsolation"]->Fill(tightmuons[0]->PFIsoDb()/tightmuons[0]->Pt() , weight);
// 	}
// 	else
// 	{
// 		reco1d["ElIsolation"]->Fill(mediumelectrons[0]->CorPFIsolation(), weight);
// 	}
// 	//double Mt_W = Sqrt(2.*met.Pt()*lep->Pt()-2.*(met.Px()*lep->Px() + met.Py()*lep->Py()));
// 	//reco1d["Mt_W"]->Fill(Mt_W, weight);
// 	//calculating btag eff.

// 	if(rightper.IsComplete())
// 	{
// 		if(dynamic_cast<IDJet*>(rightper.BLep())->csvIncl() > B_MEDIUM) truth1d["EffL_BpassingM"]->Fill(rightper.BLep()->Pt(), weight);
// 		if(dynamic_cast<IDJet*>(rightper.BLep())->csvIncl() > B_LOOSE) truth1d["EffL_BpassingL"]->Fill(rightper.BLep()->Pt(), weight);
// 		truth1d["EffL_Ball"]->Fill(rightper.BLep()->Pt(), weight);
// 		if(dynamic_cast<IDJet*>(rightper.BHad())->csvIncl() > B_MEDIUM) truth1d["EffH_BpassingM"]->Fill(rightper.BHad()->Pt(), weight);
// 		if(dynamic_cast<IDJet*>(rightper.BHad())->csvIncl() > B_LOOSE) truth1d["EffH_BpassingL"]->Fill(rightper.BHad()->Pt(), weight);
// 		truth1d["EffH_Ball"]->Fill(rightper.BHad()->Pt(), weight);
// 	}
// 	//check for b-jets
// 	reco1d["btag_high"]->Fill(cleanedjets[0]->csvIncl(), weight);
// 	reco1d["btag_low"]->Fill(cleanedjets[1]->csvIncl(), weight);
// 	if(cleanedjets[0]->csvIncl() > B_MEDIUM){ reco1d["counter"]->Fill(4.5, weight);}

// 	if(!BTAGMODE)
// 	{
// 		double bkgweight = 1.;
// 		if(cnbtag == "sigml") {if(cleanedjets[0]->csvIncl() < B_MEDIUM || cleanedjets[1]->csvIncl() < B_LOOSE){return;}}
// 		else if(cnbtag == "bkg") {bkgweight = 2.; if(cleanedjets[0]->csvIncl() < 0.0 || cleanedjets[0]->csvIncl() > 0.6){return;}}
// 		else if(cnbtag == "bkgh") {bkgweight = 2.; if(cleanedjets[0]->csvIncl() < 0.4 || cleanedjets[0]->csvIncl() > 0.7){return;}}
// 		else if(cnbtag == "bkgl") {bkgweight = 2.; if(cleanedjets[0]->csvIncl() < 0.0 || cleanedjets[0]->csvIncl() > 0.4){return;}}
// 		else if(cnbtag == "sigmm") {if(cleanedjets[0]->csvIncl() < B_MEDIUM || cleanedjets[1]->csvIncl() < B_MEDIUM){return;}}
// 		else {cout << "nbtag != sigml, sigmm, bkg, bkgl, bkgh" << endl;}
		
// 		if(bkgweight == 2. && isDA && cleanedjets.size() > 5) weight *= 2.5;
// 	}
// 	reco1d["counter"]->Fill(5.5, weight);
// 	if(SEMILEPACC) truth1d["counter"]->Fill(5.5, weight);
// 	double nvtx = event.vertexs().size();
// 	reco1d["NumVertices"]->Fill(nvtx , weight/puweight);
// 	reco1d["NumVerticesWeighted"]->Fill(nvtx , weight);
// 	reco1d["Rho"]->Fill(event.rho().value() , weight/puweight);
// 	reco1d["RhoWeighted"]->Fill(event.rho().value() , weight);

// 	//check what we have reconstructed
// 	if(SEMILEPACC)
// 	{
// 		truth2d["tt_jets"]->Fill(rightper.NumBJets()+0.5, rightper.NumWJets()+0.5, weight);
// 	}


// 	if(SEMILEPACC && rightper.IsComplete()) truth1d["counter"]->Fill(6.5, weight);

// 	if(rightper.IsComplete())
// 	{
// 		rightper.Solve(ttsolver);
// 		truth1d["nschi_right"]->Fill(ttsolver.NSChi2()/Sqrt(Abs(ttsolver.NSChi2())), weight);
// 		truth2d["Wmasshad_tmasshad_right"]->Fill(rightper.THad().M(), rightper.WHad().M(), weight);
// 		truth2d["Wmtlep_tmtlep_right"]->Fill(rightper.MttLep(), rightper.MtWLep(), weight);
// 	}

// 	//reconstruction
// 	if(LHCPS)
// 	{
// 		bestper.Reset();
// 		//double mt = sqrt(2*lep->Pt()*met.Pt()*(1.-cos(lep->DeltaPhi(met))));
// 		//if(mt < 30. || met.Pt() < 30.) {return;}
// 		vector<IDJet*> bjs;
// 		vector<IDJet*> ljs;
// 		for(IDJet* j : cleanedjets)
// 		{
// 			if(j->csvIncl() > B_MEDIUM)
// 			{
// 				bjs.push_back(j);
// 			}
// 			else
// 			{
// 				ljs.push_back(j);
// 			}

// 		}

// 		if(bjs.size() < 2 || ljs.size() < 2){return;}
// 		sort(bjs.begin(), bjs.end(), sortpt<IDJet>);
// 		sort(ljs.begin(), ljs.end(), sortpt<IDJet>);

// 		double MW = 80.4;
// 		double A = 0.5*MW*MW + lep->Px()*met.Px() + lep->Py()*met.Py();
// 		double N = lep->Pz()*lep->Pz() - lep->E()*lep->E();
// 		double p = -1.*A*lep->Pz()/N;
// 		double q = (A*A - met.Pt()*met.Pt()*lep->E()*lep->E())/N;

// 		double pz = p;
// 		if(p*p - q > 0.)
// 		{
// 			if(p > 0 ){pz += sqrt(p*p - q);}
// 			else{pz -= sqrt(p*p - q);}
// 		}
// 		TLorentzVector nu;
// 		nu.SetXYZM(met.Px(), met.Py(), pz, 0.);

// 		IDJet* bl = (lep->DeltaR(*bjs[0]) < lep->DeltaR(*bjs[1]) ? bjs[0] : bjs[1]);
// 		IDJet* bh = (bl == bjs[0] ? bjs[1] : bjs[0]);

// 		bestper.Init(ljs[0], ljs[1], bh, bl, lep, leppdgid, nu);
// 		//bestper.Solve(ttsolver, false);
// 	}
// 	else
// 	{
// 		reconstruction();
// 		if(bestper.Prob() > 1E9){return;}
// 	}
// 	if(bestper.IsComplete() == false){return;}
// 	bestper.SetAdditionalJets(cleanedjets, [&](IDJet* jet){return jet->Pt() >= jetptmin;});

// 	if(rightper.IsComplete() && cleanedjets.size() == 4)
// 	{
// 		truth1d["TTRECO"]->Fill(0.5);
// 		if(rightper.IsWHadCorrect(bestper)) {truth1d["TTRECO"]->Fill(1.5);}
// 		if(rightper.IsTHadCorrect(bestper)) {truth1d["TTRECO"]->Fill(2.5);}
// 		if(rightper.IsCorrect(bestper)) {truth1d["TTRECO"]->Fill(3.5);}
// 		if(!rightper.IsWHadCorrect(bestper) && rightper.IsBHadCorrect(bestper)) {truth1d["TTRECO"]->Fill(4.5);}
// 		if(!rightper.IsWHadCorrect(bestper) && rightper.IsBLepCorrect(bestper)) {truth1d["TTRECO"]->Fill(5.5);}
		
// 	}
// 	//if(!BTAGMODE && bestper.MassDiscr() > clikelihoodcut){return;}
// 	//if(!BTAGMODE && bestper.Prob() > clikelihoodcut){return;}
// 	if(STUDENT)
// 	{
// 		num_det = 2;
// 		det_px[0] = bestper.T().Px();
// 		det_py[0] = bestper.T().Py();
// 		det_pz[0] = bestper.T().Pz();
// 		det_e[0] = bestper.T().E();
// 		det_type[0] = (bestper.LPDGId() < 0 ? bestper.LPDGId() : -1);
// 		det_px[1] = bestper.Tb().Px();
// 		det_py[1] = bestper.Tb().Py();
// 		det_pz[1] = bestper.Tb().Pz();
// 		det_e[1] = bestper.Tb().E();
// 		det_type[1] = (bestper.LPDGId() > 0 ? bestper.LPDGId() : 1);
// 	}

// 	//bestper.Solve(ttsolver, true);
// 	reco1d["counter"]->Fill(6.5, weight);
// 	if(SEMILEPACC && rightper.IsComplete()) truth1d["counter"]->Fill(7.5, weight);
// 	if(BTAGMODE){btageff.Fill(bestper, nvtx, bestper.IsTHadCorrect(rightper), bestper.IsBLepCorrect(rightper), weight);}
// 	if(BTAGMODE && (dynamic_cast<IDJet*>(bestper.BLep())->csvIncl() < B_MEDIUM && dynamic_cast<IDJet*>(bestper.BHad())->csvIncl() < B_MEDIUM)) return;
// 	//Fill reconstructed hists
// 	ttp_all.Fill(bestper, weight);
// 	response.FillAll("thardpt", bestper.THard().Pt(), weight);
// 	response.FillAll("tsoftpt", bestper.TSoft().Pt(), weight);
// 	response.FillAll("nobin", bestper.THad().Pt(), weight);
// 	response.FillAll("thadpt", bestper.THad().Pt(), weight);
// 	response.FillAll("thady", Abs(bestper.THad().Rapidity()), weight);
// 	response.FillAll("tleppt", bestper.TLep().Pt(), weight);
// 	response.FillAll("tlepy", Abs(bestper.TLep().Rapidity()), weight);
// 	response.FillAll("ttm", bestper.TT().M(), weight);
// 	response.FillAll("ttpt", bestper.TT().Pt(), weight);
// 	response.FillAll("tty", Abs(bestper.TT().Rapidity()), weight);
// 	response.FillAll("njet", bestper.NAddJets(), weight);
// 	response.FillAll("dymp", bestper.T().Rapidity() - bestper.Tb().Rapidity(), weight);
// 	response.FillAll("dy", Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), weight);
// 	response.FillAll("ht", bestper.Ht(), weight);
// 	response.FillAll("evtmass", bestper.EvtMass(), weight);
// 	response2d.FillAll("njets_thadpt", bestper.THad().Pt(), bestper.NAddJets(), weight);
// 	response2d.FillAll("njets_ttpt", bestper.TT().Pt(), bestper.NAddJets(), weight);
// 	response2d.FillAll("njets_ttm", bestper.TT().M(), bestper.NAddJets(), weight);
// 	response2d.FillAll("thady_thadpt", bestper.THad().Pt(), Abs(bestper.THad().Rapidity()), weight);
// 	response2d.FillAll("ttpt_thadpt", bestper.THad().Pt(), bestper.TT().Pt(), weight);
// 	response2d.FillAll("thadpt_ttm", bestper.TT().M(), bestper.THad().Pt(), weight);
// 	response2d.FillAll("ttpt_ttm", bestper.TT().M(), bestper.TT().Pt(), weight);
// 	response2d.FillAll("ttm_tty", Abs(bestper.TT().Rapidity()), bestper.TT().M(), weight);
// 	response2d.FillAll("ttm_dy", Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), bestper.TT().M(), weight);

// 	response2dvar.FillAll("njet+thadpt", bestper.NAddJets(), bestper.THad().Pt(), weight);
// 	response2dvar.FillAll("njet+ttpt", bestper.NAddJets(), bestper.TT().Pt(), weight);
// 	response2dvar.FillAll("njet+ttm", bestper.NAddJets(), bestper.TT().M(), weight);
// 	response2dvar.FillAll("thady+thadpt", Abs(bestper.THad().Rapidity()), bestper.THad().Pt(), weight);
// 	response2dvar.FillAll("ttpt+thadpt", bestper.TT().Pt(), bestper.THad().Pt(), weight);
// 	response2dvar.FillAll("thadpt+ttm", bestper.THad().Pt(), bestper.TT().M(), weight);
// 	response2dvar.FillAll("ttpt+ttm", bestper.TT().Pt(), bestper.TT().M(), weight);
// 	response2dvar.FillAll("ttm+tty", bestper.TT().M(), Abs(bestper.TT().Rapidity()), weight);
// 	response2dvar.FillAll("ttm+dy", bestper.TT().M(), Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), weight);
// 	for(size_t n = 0 ; n < bestper.NJets() ; ++n)
// 	{
// 		response2dvar.FillAll("jet+jetpt", n, bestper.GetJet(n)->Pt(), weight);
// 		response2dvar.FillAll("jet+jetptOF", n, bestper.GetJet(n)->Pt(), weight);
// 		response2dvar.FillAll("jet+jeteta", n, abs(bestper.GetJet(n)->Eta()), weight);
// 		response2dvar.FillAll("jet+jetdr", n, bestper.DRminTTjets(bestper.GetJet(n)), weight);
// 		response2dvar.FillAll("jet+jetdrtop", n, bestper.DRminTop(bestper.GetJet(n)), weight);
// 	}


// 	if(SEMILEPACC)
// 	{
// 		if(PDFTEST)
// 		{
// 			pdfunc->Fill1d("pdfunc_reco_thardpt", genper->THard().Pt(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_tsoftpt", genper->TSoft().Pt() , weight);
// 			pdfunc->Fill1d("pdfunc_reco_nobin", genthad.Pt(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_thadpt", genthad.Pt(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_tleppt", gentlep.Pt(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_thady", Abs(genthad.Rapidity()), weight);
// 			pdfunc->Fill1d("pdfunc_reco_tlepy", Abs(gentlep.Rapidity()), weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttm", (gent+gentbar).M(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_tty", Abs((gent+gentbar).Rapidity()), weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttpt", (gent+gentbar).Pt(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_njet", genper->NAddJets(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_dymp", genper->T().Rapidity() - genper->Tb().Rapidity(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_ht", genper->Ht(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_evtmass", genper->EvtMass(), weight);
// 			pdfunc->Fill1d("pdfunc_reco_njet+thadpt", response2dvar.GetBin("njet+thadpt", genper->NAddJets(), genthad.Pt())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_njet+ttpt", response2dvar.GetBin("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_njet+ttm", response2dvar.GetBin("njet+ttm", genper->NAddJets(), (gent+gentbar).M())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_thady+thadpt", response2dvar.GetBin("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttpt+thadpt", response2dvar.GetBin("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_thadpt+ttm", response2dvar.GetBin("thadpt+ttm", genthad.Pt(), (gent+gentbar).M())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttpt+ttm", response2dvar.GetBin("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M())-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttm+tty", response2dvar.GetBin("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()))-0.5, weight);
// 			pdfunc->Fill1d("pdfunc_reco_ttm+dy", response2dvar.GetBin("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()))-0.5, weight);
// 			for(size_t n = 0 ; n < genper->NJets() ; ++n)
// 			{
// 				pdfunc->Fill1d("pdfunc_reco_jet+jetpt", response2dvar.GetBin("jet+jetpt", n, genper->GetJet(n)->Pt())-0.5, weight);
// 				pdfunc->Fill1d("pdfunc_reco_jet+jetptOF", response2dvar.GetBin("jet+jetptOF", n, genper->GetJet(n)->Pt())-0.5, weight);
// 				pdfunc->Fill1d("pdfunc_reco_jet+jeteta", response2dvar.GetBin("jet+jeteta", n, abs(genper->GetJet(n)->Eta()))-0.5, weight);
// 				pdfunc->Fill1d("pdfunc_reco_jet+jetdr", response2dvar.GetBin("jet+jetdr", n, genper->DRminTTjets(genper->GetJet(n)))-0.5, weight);
// 				pdfunc->Fill1d("pdfunc_reco_jet+jetdrtop", response2dvar.GetBin("jet+jetdrtop", n, genper->DRminTop(genper->GetJet(n)))-0.5, weight);
// 			}
// 		}
// 		response.FillTruthReco("thardpt", Max(gent.Pt(), gentbar.Pt()), bestper.THard().Pt(), weight);
// 		response.FillTruthReco("tsoftpt", Min(gent.Pt(), gentbar.Pt()), bestper.TSoft().Pt(), weight);
// 		response.FillTruthReco("nobin", genthad.Pt(), bestper.THad().Pt(), weight);
// 		response.FillTruthReco("thadpt", genthad.Pt(), bestper.THad().Pt(), weight);
// 		response.FillTruthReco("thady", Abs(genthad.Rapidity()), Abs(bestper.THad().Rapidity()), weight);
// 		response.FillTruthReco("tleppt", gentlep.Pt(), bestper.TLep().Pt(), weight);
// 		response.FillTruthReco("tlepy", Abs(gentlep.Rapidity()), Abs(bestper.TLep().Rapidity()), weight);
// 		response.FillTruthReco("ttm", (gent+gentbar).M(), bestper.TT().M(), weight);
// 		response.FillTruthReco("ttpt", (gent+gentbar).Pt(), bestper.TT().Pt(), weight);
// 		response.FillTruthReco("tty", Abs((gent+gentbar).Rapidity()), Abs(bestper.TT().Rapidity()), weight);
// 		response.FillTruthReco("njet", genper->NAddJets(), bestper.NAddJets(), weight);
// 		response.FillTruthReco("dymp", genper->T().Rapidity() - genper->Tb().Rapidity(), bestper.T().Rapidity() - bestper.Tb().Rapidity() , weight);
// 		response.FillTruthReco("dy", Abs(genper->T().Rapidity()) - Abs(genper->Tb().Rapidity()), Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()) , weight);
// 		response.FillTruthReco("ht", genper->Ht(), bestper.Ht() , weight);
// 		response.FillTruthReco("evtmass", genper->EvtMass(), bestper.EvtMass() , weight);
// 		response2d.FillTruthReco("njets_thadpt", genthad.Pt(), genper->NAddJets(), bestper.THad().Pt(), bestper.NAddJets(), weight);
// 		response2d.FillTruthReco("njets_ttpt", (gent+gentbar).Pt(), genper->NAddJets(), bestper.TT().Pt(), bestper.NAddJets(), weight);
// 		response2d.FillTruthReco("njets_ttm", (gent+gentbar).M(), genper->NAddJets(), bestper.TT().M(), bestper.NAddJets(), weight);
// 		response2d.FillTruthReco("thady_thadpt", genthad.Pt(), Abs(genthad.Rapidity()), bestper.THad().Pt(), Abs(bestper.THad().Rapidity()), weight);
// 		response2d.FillTruthReco("ttpt_thadpt", genthad.Pt(), (gent+gentbar).Pt(), bestper.THad().Pt(), bestper.TT().Pt(), weight);
// 		response2d.FillTruthReco("thadpt_ttm", (gent+gentbar).M(), genthad.Pt(), bestper.TT().M(), bestper.THad().Pt(), weight);
// 		response2d.FillTruthReco("ttpt_ttm", (gent+gentbar).M(), (gent+gentbar).Pt(), bestper.TT().M(), bestper.TT().Pt(), weight);
// 		response2d.FillTruthReco("ttm_tty", Abs((gent+gentbar).Rapidity()), (gent+gentbar).M(), Abs(bestper.TT().Rapidity()), bestper.TT().M(), weight);
// 		response2d.FillTruthReco("ttm_dy", Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), (gent+gentbar).M(), Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), bestper.TT().M(), weight);

// 		response2dvar.FillTruthReco("njet+thadpt", genper->NAddJets(), genthad.Pt(), bestper.NAddJets(), bestper.THad().Pt(), weight);
// 		response2dvar.FillTruthReco("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt(), bestper.NAddJets(), bestper.TT().Pt(), weight);
// 		response2dvar.FillTruthReco("njet+ttm", genper->NAddJets(), (gent+gentbar).M(), bestper.NAddJets(), bestper.TT().M(), weight);
// 		response2dvar.FillTruthReco("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt(), Abs(bestper.THad().Rapidity()), bestper.THad().Pt(), weight);
// 		response2dvar.FillTruthReco("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt(), bestper.TT().Pt(), bestper.THad().Pt(), weight);
// 		response2dvar.FillTruthReco("thadpt+ttm", genthad.Pt(), (gent+gentbar).M(), bestper.THad().Pt(), bestper.TT().M(), weight);
// 		response2dvar.FillTruthReco("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M(), bestper.TT().Pt(), bestper.TT().M(), weight);
// 		response2dvar.FillTruthReco("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()), bestper.TT().M(), Abs(bestper.TT().Rapidity()), weight);
// 		response2dvar.FillTruthReco("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), bestper.TT().M(), Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), weight);
// 		for(size_t n = 0 ; n < cleanedjets.size() ; ++n)
// 		{
// 			int jtrec = bestper.IsJetIn(cleanedjets[n]);
// 			int jttruth = genper->IsJetIn(cleanedjets[n]);
// 			if(jtrec != -1 && jttruth != -1)
// 			{
// 				const TLorentzVector* rjet = bestper.GetJet(jtrec);
// 				const TLorentzVector* gjet = genper->GetJet(jttruth);
// 				response2dvar.FillTruthReco("jet+jetpt", jttruth, gjet->Pt(), jtrec, rjet->Pt(), weight);
// 				response2dvar.FillTruthReco("jet+jetptOF", jttruth, gjet->Pt(), jtrec, rjet->Pt(), weight);
// 				response2dvar.FillTruthReco("jet+jeteta", jttruth, abs(gjet->Eta()), jtrec, abs(rjet->Eta()), weight);
// 				response2dvar.FillTruthReco("jet+jetdr", jttruth, genper->DRminTTjets(gjet), jtrec, bestper.DRminTTjets(rjet), weight);
// 				response2dvar.FillTruthReco("jet+jetdrtop", jttruth, genper->DRminTop(gjet), jtrec, bestper.DRminTop(rjet), weight);
// 			}
// 		}
// 	}


// 	//Fill reconstructed hists with matching information
// 	if(rightper.IsCorrect(bestper))
// 	{
// 		ttp_right.Fill(bestper, weight);
// 		truth1d["counter"]->Fill(8.5, weight);
// 		response.FillRes("thardpt", Max(gent.Pt(), gentbar.Pt()), bestper.THard().Pt(), weight);
// 		response.FillRes("tsoftpt", Min(gent.Pt(), gentbar.Pt()), bestper.TSoft().Pt(), weight);
// 		response.FillRes("nobin", genthad.Pt(), bestper.THad().Pt(), weight);
// 		response.FillRes("thadpt", genthad.Pt(), bestper.THad().Pt(), weight);
// 		response.FillRes("thady", Abs(genthad.Rapidity()), Abs(bestper.THad().Rapidity()), weight);
// 		response.FillRes("tleppt", gentlep.Pt(), bestper.TLep().Pt(), weight);
// 		response.FillRes("tlepy", Abs(gentlep.Rapidity()), Abs(bestper.TLep().Rapidity()), weight);
// 		response.FillRes("ttm", (gent+gentbar).M(), bestper.TT().M(), weight);
// 		response.FillRes("ttpt", (gent+gentbar).Pt(), bestper.TT().Pt(), weight);
// 		response.FillRes("tty", Abs((gent+gentbar).Rapidity()), Abs(bestper.TT().Rapidity()), weight);
// 		response.FillRes("njet", genper->NAddJets(), bestper.NAddJets(), weight);
// 		response.FillRes("dymp", gent.Rapidity() - gentbar.Rapidity(), bestper.T().Rapidity() - bestper.Tb().Rapidity() , weight);
// 		response.FillRes("ht", genper->Ht(), bestper.Ht() , weight);
// 		response.FillRes("evtmass", genper->EvtMass(), bestper.EvtMass() , weight);
// 		response2dvar.FillRes("njet+thadpt", genper->NAddJets(), genthad.Pt(), bestper.NAddJets(), bestper.THad().Pt(), weight);
// 		response2dvar.FillRes("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt(), bestper.NAddJets(), bestper.TT().Pt(), weight);
// 		response2dvar.FillRes("njet+ttm", genper->NAddJets(), (gent+gentbar).M(), bestper.NAddJets(), bestper.TT().M(), weight);
// 		response2dvar.FillRes("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt(), Abs(bestper.THad().Rapidity()), bestper.THad().Pt(), weight);
// 		response2dvar.FillRes("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt(), bestper.TT().Pt(), bestper.THad().Pt(), weight);
// 		response2dvar.FillRes("thadpt+ttm", genthad.Pt(), (gent+gentbar).M(), bestper.THad().Pt(), bestper.TT().M(), weight);
// 		response2dvar.FillRes("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M(), bestper.TT().Pt(), bestper.TT().M(), weight);
// 		response2dvar.FillRes("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()), bestper.TT().M(), Abs(bestper.TT().Rapidity()), weight);
// 		response2dvar.FillRes("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), bestper.TT().M(), Abs(bestper.T().Rapidity()) - Abs(bestper.Tb().Rapidity()), weight);
// 		truth1d["dRNu_right"]->Fill(bestper.Nu().DeltaR(genper->Nu()), weight);
// 		truth1d["dPtNu_right"]->Fill((bestper.Nu().Pt() - genper->Nu().Pt())/genper->Nu().Pt(), weight);
// 		truth1d["dPzNu_right"]->Fill((bestper.Nu().Pz() - genper->Nu().Pz())/genper->Nu().Pz(), weight);
// 		truth2d["dPzNu_dPhi_right"]->Fill((bestper.Nu().Pz() - genper->Nu().Pz())/genper->Nu().Pz(), Abs(bestper.Nu().DeltaPhi(met)), weight);
// 		truth2d["dPzNu_chi2_right"]->Fill((bestper.Nu().Pz() - genper->Nu().Pz())/genper->Nu().Pz(), Sqrt(bestper.NuChisq()), weight);
// 		truth2d["RES_Mtt_right"]->Fill(((bestper.THad() + bestper.TLep()).M() - (genthad + gentlep).M())/(genthad + gentlep).M(), (genthad + gentlep).M(), weight);
// 		truth2d["RES_dbeta_right"]->Fill(((bestper.THad().BoostVector() - bestper.TLep().BoostVector()).Mag() - (genthad.BoostVector() - gentlep.BoostVector()).Mag())/(genthad.BoostVector() - gentlep.BoostVector()).Mag(), (genthad.BoostVector() - gentlep.BoostVector()).Mag(), weight);
// 		truth1d["dRNuMet_right"]->Fill(met.DeltaR(genper->Nu()), weight);
// 		truth1d["dPtNuMet_right"]->Fill((met.Pt() - genper->Nu().Pt())/genper->Nu().Pt(), weight);
// 		truth2d["RES_Mtt_all"]->Fill((bestper.TT().M() - (gent+gentbar).M())/(gent+gentbar).M(), (gent+gentbar).M(), weight);
// 		truth2d["RES_dbeta_all"]->Fill(((bestper.THad().BoostVector() - bestper.TLep().BoostVector()).Mag() - (genthad.BoostVector() - gentlep.BoostVector()).Mag())/(genthad.BoostVector() - gentlep.BoostVector()).Mag(), (genthad.BoostVector() - gentlep.BoostVector()).Mag(), weight);
// 	}
// 	else if(rightper.IsComplete())
// 	{
// 		ttp_wrong.Fill(bestper, weight);
// 		truth1d["counter"]->Fill(9.5, weight);
// 	}
// 	else if(SEMILEPACC)
// 	{
// 		ttp_semi.Fill(bestper, weight);
// 	}
// 	else
// 	{
// //		cout << "RECOL:" << bestper.LPDGId() << " " << bestper.L()->Eta() << " " << bestper.L()->Pt() << endl;
// //		for(int j = 0 ; j < min(bestper.NJets(),8) ; ++j)
// //		{
// //			const TLorentzVector* jet = bestper.GetJet(j); 
// //			cout << "RECO: " << j << " " << jet->Eta() << " " << jet->Pt() << " " << jet->DeltaR(*bestper.L()) << endl;
// //		}

// 		ttp_other.Fill(bestper, weight);
// 	}

// 	if(bestper.AreHadJetsCorrect(rightper) && bestper.IsBLepCorrect(rightper))
// 	{
// 		ttp_tlepthad_right.Fill(bestper, weight);
// 	}
// 	else if(bestper.AreHadJetsCorrect(rightper))
// 	{
// 		ttp_thad_right.Fill(bestper, weight);
// 	}
// 	else if(bestper.IsBLepCorrect(rightper))
// 	{
// 		ttp_tlep_right.Fill(bestper, weight);
// 	}
// 	else if(SEMILEPACC)
// 	{
// 		ttp_nn_right.Fill(bestper, weight);
// 	}
// 	else
// 	{
// 		ttp_nsemi_right.Fill(bestper, weight);
// 	}

// 	if(SEMILEPACC)
// 	{
// 		for(size_t j = 0 ; j < min(bestper.NJets(), (size_t)8) ; ++j)
// 		{
// 			stringstream ss;
// 			ss << "jet_" << j << "_";
// 			const TLorentzVector* jet = bestper.GetJet(j); 
// 			if(genper->IsJetIn(jet) == j)
// 			{
// 				truth1d[ss.str()+"pt_fb_right"]->Fill(jet->Pt(),weight);
// 				truth1d[ss.str()+"eta_fb_right"]->Fill(Abs(jet->Eta()), weight);
// 				truth1d[ss.str()+"dr_fb_right"]->Fill(bestper.DRminTTjets(jet), weight);
// 				truth1d[ss.str()+"drtop_fb_right"]->Fill(bestper.DRminTop(jet), weight);
// 				const TLorentzVector* gjet = genper->GetJet(j); 
// 				response2dvar.FillRes("jet+jetpt", j, gjet->Pt(), j, jet->Pt(), weight);
// 				response2dvar.FillRes("jet+jetptOF", j, gjet->Pt(), j, jet->Pt(), weight);
// 				response2dvar.FillRes("jet+jeteta", j, abs(gjet->Eta()), j, abs(jet->Eta()), weight);
// 				response2dvar.FillRes("jet+jetdr", j, genper->DRminTTjets(gjet), j, bestper.DRminTTjets(jet), weight);
// 				response2dvar.FillRes("jet+jetdrtop", j, genper->DRminTop(gjet), j, bestper.DRminTop(jet), weight);
// 			}
// 			else
// 			{
// 				truth1d[ss.str()+"pt_fb_wrong"]->Fill(jet->Pt(), weight);
// 				truth1d[ss.str()+"eta_fb_wrong"]->Fill(Abs(jet->Eta()), weight);
// 				truth1d[ss.str()+"dr_fb_wrong"]->Fill(bestper.DRminTTjets(jet), weight);
// 				truth1d[ss.str()+"drtop_fb_wrong"]->Fill(bestper.DRminTTjets(jet), weight);
// 			}
// 		}

// 		if(bestper.AreJetsCorrect(*genper))
// 		{
// 			ttp_alljets_right.Fill(bestper, weight);
// 		}
// 		else
// 		{
// 			ttp_alljets_wrong.Fill(bestper, weight);
// 		}
// 	}

}

void ttbar::reconstruction()
{
	bestper.Reset();
	int nbtaglocal = 2;
	if(BTAGMODE)
	{
		nbtaglocal = 0;
	}
	int percount = 0;
	for(size_t i = nbtaglocal ; i < cleanedjets.size() ; ++i)
	{
		for(size_t j = nbtaglocal ; j < i ; ++j)
		{
			for(size_t k = 0 ; k < (nbtaglocal == 2 ? 2 : cleanedjets.size()) ; ++k)
			{
				if(i == k || j == k) continue;
				for(size_t l = 0 ; l < (nbtaglocal == 2 ? 2 : cleanedjets.size()) ; ++l)
				{
					if(l == i || l == j || l == k) continue;
					if(nbtaglocal == 1 && k != 0 && l != 0) continue;
					testper.Init(cleanedjets[i], cleanedjets[j], cleanedjets[k], cleanedjets[l], lep, leppdgid, &met);
					if(testper.WJa()->Pt() < cwjetpthard && testper.WJb()->Pt() < cwjetpthard) continue;
					if(testper.WJa()->Pt() < cwjetptsoft || testper.WJb()->Pt() < cwjetptsoft) continue;
					if(testper.BHad()->Pt() < cbjetpthard && testper.BLep()->Pt() < cbjetpthard) continue;
					if(testper.BHad()->Pt() < cbjetptsoft || testper.BLep()->Pt() < cbjetptsoft) continue;
					testper.Solve(ttsolver, false);

					//if(rightper.IsComplete()) cout << (percount++) << " " << rightper.IsCorrect(testper) << " "  << testper.WHad().M() << " " << testper.THad().M() << " " <<  testper.Prob() << " " << testper.MassDiscr() << " " << testper.NuChisq() << endl;

					reco2d["Wmasshad_tmasshad"]->Fill(testper.WHad().M(), testper.THad().M());
					if(rightper.IsComplete())
					{
						if(rightper.IsBLepCorrect(testper))
						{
							truth1d["nstest_right"]->Fill(ttsolver.NSRes(), weight);
						}
						else
						{
							truth1d["nstest_wrong"]->Fill(ttsolver.NSRes(), weight);
							truth1d["nschi_wrong"]->Fill(ttsolver.NSChi2()/Sqrt(Abs(ttsolver.NSChi2())), weight);
							truth2d["Wmtlep_tmtlep_wrong"]->Fill(testper.MttLep(), testper.MtWLep(), weight);
						}

						if(rightper.AreBsCorrect(testper))
						{
							truth1d["btagtest_right"]->Fill(ttsolver.BTagRes(), weight);
						}
						else
						{
							truth1d["btagtest_wrong"]->Fill(ttsolver.BTagRes(), weight);
						}

						if(rightper.IsCorrect(testper))
						{
							truth1d["masstest_right"]->Fill(ttsolver.MassRes(), weight);
							truth1d["comtest_right"]->Fill(ttsolver.Res(), weight);
						}
						else
						{
							if(rightper.IsWHadCorrect(testper))
							{
								truth2d["Wmasshad_tmasshad_rightw"]->Fill(testper.WHad().M(), testper.THad().M(), weight);
							}
							else
							{
								truth2d["Wmasshad_tmasshad_wrongw"]->Fill(testper.WHad().M(), testper.THad().M(), weight);
							}
							truth2d["Wmasshad_tmasshad_wrong"]->Fill(testper.THad().M(), testper.WHad().M(), weight);
							truth1d["masstest_wrong"]->Fill(ttsolver.MassRes(), weight);
							truth1d["comtest_wrong"]->Fill(ttsolver.Res(), weight);
						}

					}
					reco1d["btagtest"]->Fill(ttsolver.BTagRes(), weight);
					reco1d["masstest"]->Fill(ttsolver.MassRes(), weight);
					reco1d["nstest"]->Fill(ttsolver.NSRes(), weight);

					if(testper < bestper)
					{
						bestper = testper;
					}
				}
			}
		}
	}
}

double ttbar::lepeffweight(TLorentzVector* lep, URStreamer& event)
{
	//only do this on actual data!
	if(isDA != 0) return 1.;
	double lepw = 1.;
	double leperror = 0.02;
	if(tightmuons.size() == 1)
	{	
		int bx = musfhist->GetXaxis()->FindFixBin(lep->Eta());
		int by = musfhist->GetYaxis()->FindFixBin(Min(lep->Pt(), 170.));
		lepw = musfhist->GetBinContent(bx, by) + csigmalep * sqrt(pow(musfhist->GetBinError(bx, by),2) + pow(leperror, 2));
	}
	else if(mediumelectrons.size() == 1)
	{
		int bx = elsfhist->GetXaxis()->FindFixBin(dynamic_cast<IDElectron*>(lep)->SCeta());
		int by = elsfhist->GetYaxis()->FindFixBin(Min(lep->Pt(), 170.));
		lepw = elsfhist->GetBinContent(bx, by) + csigmalep * sqrt(pow(elsfhist->GetBinError(bx, by),2) + pow(leperror, 2));
		int l1ptmax = event.trigger().El27ptmax();
		if(l1ptmax != -1 && l1ptmax < 34)
		{
			lepw*=cel27eff[(l1ptmax-24)/2];
		}
	}
	return lepw;
}


//This method is called once every file, contains the event loop
//run your proper analysis here
void ttbar::analyze()
{
	int nevent = 0;
	URStreamer event(tree_);
	IDElectron::streamer = &event;
	IDMuon::streamer = &event;
	PDFuncertainty::streamer = &event;
	while(event.next())
	{
		nevent++;
		if(nevent % 10000 == 0)cout << "Event:" << nevent << " " << event.run << endl;
		num_det = 0;
		num_gen = 0;
		sgenparticles.clear();
		genalljets.clear();
		genbjets.clear();
		gencjets.clear();
		genljets.clear();
		plleptons.clear();
		plphotons.clear();

		genallper.Reset();
		psper.Reset();
		rightper.Reset();
		genper = 0;

		sjets.clear();
		cleanedjets.clear();
		smuons.clear();
		tightmuons.clear();
		loosemuons.clear();
		selectrons.clear();
		mediumelectrons.clear();
		looseelectrons.clear();

		lep = nullptr;
		leppdgid = 0;
		truth1d["counter"]->Fill(0.5);
		weight = 1.;	
		mcweight = 1.;	
		puweight = 1.;	
		if(isDA == 0)
		{
			const Geninfo& info = event.genInfo();
			mcweight = (info.weight() < 1. ? -1. : 1.);
			truth1d["counter"]->Fill(19.5, weight);
			if(TTMC)
			{
				const vector<Mcweight>& ws =  event.MCWeights();
				if(cfacscale == -1 && crenscale == -1) mcweight = ws[8].weights()/Abs(ws[0].weights());
				else if(cfacscale == 1 && crenscale == 1) mcweight = ws[4].weights()/Abs(ws[0].weights());
				else if(cfacscale == -1 && crenscale == 1) mcweight = ws[5].weights()/Abs(ws[0].weights());
				else if(cfacscale == 1 && crenscale == -1) mcweight = ws[7].weights()/Abs(ws[0].weights());
				else if(cfacscale == -1 && crenscale == 0) mcweight = ws[2].weights()/Abs(ws[0].weights());
				else if(cfacscale == 1 && crenscale == 0) mcweight = ws[1].weights()/Abs(ws[0].weights());
				else if(cfacscale == 0 && crenscale == -1) mcweight = ws[6].weights()/Abs(ws[0].weights());
				else if(cfacscale == 0 && crenscale == 1) mcweight = ws[3].weights()/Abs(ws[0].weights());
				if(chdamp == -1) mcweight = ws[227].weights()/Abs(ws[0].weights());
				else if(chdamp == 1) mcweight = ws[245].weights()/Abs(ws[0].weights());
			}
			truth1d["counter"]->Fill(18.5, mcweight);
			weight *= mcweight;
			double npu = event.PUInfos()[0].nInteractions();
			//cout << event.PUInfos()[0].nInteractions() << " " << event.PUInfos()[1].nInteractions() << endl;
			truth1d["Mu"]->Fill(npu, weight);
			puweight = puhist->GetBinContent(puhist->FindFixBin(npu));
			weight *= puweight;
			//cout << weight << " " << npu << endl;
			truth1d["MuWeighted"]->Fill(npu, weight);
			truth1d["counter"]->Fill(17.5, weight);
		}
		else
		{
			runinfo[event.run].insert(event.lumi);
		}

		if(TTMC)
		{
			SelectGenParticles(event);
			gen1d["tpt"]->Fill(gentq.Pt(), weight);
			gen1d["ty"]->Fill(Abs(gentq.Rapidity()), weight);
			gen1d["ttpt"]->Fill((gentq+gentqbar).Pt(), weight);
			gen1d["tty"]->Fill(Abs((gentq+gentqbar).Rapidity()), weight);
			if(LHCPS)
			{
				SelectPseudoTopLHC(event);
			}
			else
			{
				SelectPseudoTop(event);
				//SelectRivetPS(event);

			}
			if(psper.IsComplete() && psper.NAddJets() > 0)
			{
				double addjetpt = psper.GetJet(4)->Pt();
				weight *= 1 + cjetptweight*(addjetpt-250)/250 * 0.2;
			}
			
		}

		if(PSEUDOTOP)
		{
			if(genallper.IsComplete())
			{
				cout<<"?";
				response_ps.FillTruth("thardpt", Max(gentq.Pt(), gentqbar.Pt()), weight);
				response_ps.FillTruth("tsoftpt", Min(gentq.Pt(), gentqbar.Pt()), weight);
				response_ps.FillTruth("thadpt", gentqhad.Pt(), weight);
				response_ps.FillTruth("thady", Abs(gentqhad.Rapidity()), weight);
				response_ps.FillTruth("tleppt", gentqlep.Pt(), weight);
				response_ps.FillTruth("tlepy", Abs(gentqlep.Rapidity()), weight);
				response_ps.FillTruth("ttm", (gentq+gentqbar).M(), weight);
				response_ps.FillTruth("ttpt", (gentq+gentqbar).Pt(), weight);
				response_ps.FillTruth("tty", Abs((gentq+gentqbar).Rapidity()), weight);
			}
			if(psper.IsComplete())
			{
				response_ps.FillAll("thardpt", psper.THard().Pt(), weight);
				response_ps.FillAll("tsoftpt", psper.TSoft().Pt(), weight);
				response_ps.FillAll("thadpt", psper.THad().Pt(), weight);
				response_ps.FillAll("thady", Abs(psper.THad().Rapidity()), weight);
				response_ps.FillAll("tleppt", psper.TLep().Pt(), weight);
				response_ps.FillAll("tlepy", Abs(psper.TLep().Rapidity()), weight);
				response_ps.FillAll("ttm", psper.TT().M(), weight);
				response_ps.FillAll("ttpt", psper.TT().Pt(), weight);
				response_ps.FillAll("tty", Abs(psper.TT().Rapidity()), weight);
			}
			if(genallper.IsComplete() && psper.IsComplete())
			{
				response_ps.FillTruthReco("thardpt", Max(gentq.Pt(), gentqbar.Pt()), psper.THard().Pt(), weight);
				response_ps.FillTruthReco("tsoftpt", Min(gentq.Pt(), gentqbar.Pt()), psper.TSoft().Pt(), weight);
				response_ps.FillTruthReco("thadpt", gentqhad.Pt(), psper.THad().Pt(), weight);
				response_ps.FillTruthReco("thady", Abs(gentqhad.Rapidity()), Abs(psper.THad().Rapidity()), weight);
				response_ps.FillTruthReco("tleppt", gentqlep.Pt(), psper.TLep().Pt(), weight);
				response_ps.FillTruthReco("tlepy", Abs(gentqlep.Rapidity()), Abs(psper.TLep().Rapidity()), weight);
				response_ps.FillTruthReco("ttm", (gentq+gentqbar).M(), psper.TT().M(), weight);
				response_ps.FillTruthReco("ttpt", (gentq+gentqbar).Pt(), psper.TT().Pt(), weight);
				response_ps.FillTruthReco("tty", Abs((gentq+gentqbar).Rapidity()), Abs(psper.TT().Rapidity()), weight);
			}
			genper = &psper;
			SEMILEP = psper.IsComplete();
			SEMILEPACC = SEMILEP;
			if(psper.IsComplete())
			{
				gent = psper.T();
				gentbar = psper.Tb();
				gentlep = psper.TLep();
				genthad = psper.THad();
			}
		}
		else
		{
			genper = &genallper;
			gent = gentq;
			gentbar = gentqbar;
			gentlep = gentqlep;
			genthad = gentqhad;
		}

		if(isDA == 0) {AddGenJetSelection(event);}

		if(SEMILEPACC) 
		{
			truth1d["counter"]->Fill(1.5, weight);
			ttp_genall.Fill(*genper, weight);

			response.FillTruth("thardpt", Max(gent.Pt(), gentbar.Pt()), weight);
			response.FillTruth("tsoftpt", Min(gent.Pt(), gentbar.Pt()), weight);
			response.FillTruth("nobin", genthad.Pt(), weight);
			response.FillTruth("thadpt", genthad.Pt(), weight);
			response.FillTruth("thady", Abs(genthad.Rapidity()), weight);
			response.FillTruth("tleppt", gentlep.Pt(), weight);
			response.FillTruth("tlepy", Abs(gentlep.Rapidity()), weight);
			response.FillTruth("ttm", (gent+gentbar).M(), weight);
			response.FillTruth("ttpt", (gent+gentbar).Pt(), weight);
			response.FillTruth("tty", Abs((gent+gentbar).Rapidity()), weight);
			response.FillTruth("njet", genper->NAddJets(), weight);
			response.FillTruth("dymp", gent.Rapidity() - gentbar.Rapidity(), weight);
			response.FillTruth("dy", Abs(gent.Rapidity())- Abs(gentbar.Rapidity()), weight);
			response.FillTruth("ht", genper->Ht(), weight);
			response.FillTruth("evtmass", genper->EvtMass(), weight);
			response2d.FillTruth("njets_thadpt", genthad.Pt(), genper->NAddJets(), weight);
			response2d.FillTruth("njets_ttpt", (gent+gentbar).Pt(), genper->NAddJets(), weight);
			response2d.FillTruth("njets_ttm", (gent+gentbar).M(), genper->NAddJets(), weight);
			response2d.FillTruth("thady_thadpt", genthad.Pt(), Abs(genthad.Rapidity()), weight);
			response2d.FillTruth("ttpt_thadpt", genthad.Pt(), (gent+gentbar).Pt(), weight);
			response2d.FillTruth("thadpt_ttm", (gent+gentbar).M(), genthad.Pt(), weight);
			response2d.FillTruth("ttpt_ttm", (gent+gentbar).M(), (gent+gentbar).Pt(), weight);
			response2d.FillTruth("ttm_tty", Abs((gent+gentbar).Rapidity()), (gent+gentbar).M(), weight);
			response2d.FillTruth("ttm_dy", Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), (gent+gentbar).M(), weight);

			response2dvar.FillTruth("njet+thadpt", genper->NAddJets(), genthad.Pt(), weight);
			response2dvar.FillTruth("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt(), weight);
			response2dvar.FillTruth("njet+ttm", genper->NAddJets(), (gent+gentbar).M(), weight);
			response2dvar.FillTruth("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt(), weight);
			response2dvar.FillTruth("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt(), weight);
			response2dvar.FillTruth("thadpt+ttm", genthad.Pt(), (gent+gentbar).M(), weight);
			response2dvar.FillTruth("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M(), weight);
			response2dvar.FillTruth("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()), weight);
			response2dvar.FillTruth("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), weight);
			for(size_t n = 0 ; n < genper->NJets() ; ++n)
			{
				response2dvar.FillTruth("jet+jetpt", n, genper->GetJet(n)->Pt(), weight);
				response2dvar.FillTruth("jet+jetptOF", n, genper->GetJet(n)->Pt(), weight);
				response2dvar.FillTruth("jet+jeteta", n, abs(genper->GetJet(n)->Eta()), weight);
				response2dvar.FillTruth("jet+jetdr", n, genper->DRminTTjets(genper->GetJet(n)), weight);
				response2dvar.FillTruth("jet+jetdrtop", n, genper->DRminTop(genper->GetJet(n)), weight);
			}
			if(PDFTEST)
			{
				pdfunc->Fill1d("pdfunc_thardpt", Max(gent.Pt(), gentbar.Pt()), weight);
				pdfunc->Fill1d("pdfunc_tsoftpt", Min(gent.Pt(), gentbar.Pt()), weight);
				pdfunc->Fill1d("pdfunc_nobin", genthad.Pt(), weight);
				pdfunc->Fill1d("pdfunc_thadpt", genthad.Pt(), weight);
				pdfunc->Fill1d("pdfunc_tleppt", gentlep.Pt(), weight);
				pdfunc->Fill1d("pdfunc_thady", Abs(genthad.Rapidity()), weight);
				pdfunc->Fill1d("pdfunc_tlepy", Abs(gentlep.Rapidity()), weight);
				pdfunc->Fill1d("pdfunc_ttm", (gent+gentbar).M(), weight);
				pdfunc->Fill1d("pdfunc_tty", Abs((gent+gentbar).Rapidity()), weight);
				pdfunc->Fill1d("pdfunc_ttpt", (gent+gentbar).Pt(), weight);
				pdfunc->Fill1d("pdfunc_njet", genper->NAddJets(), weight);
				pdfunc->Fill1d("pdfunc_dymp", gent.Rapidity() - gentbar.Rapidity(), weight);
				pdfunc->Fill1d("pdfunc_dy", Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), weight);
				pdfunc->Fill1d("pdfunc_ht", genper->Ht(), weight);
				pdfunc->Fill1d("pdfunc_evtmass", genper->EvtMass(), weight);

				pdfunc->Fill1d("pdfunc_njet+thadpt", response2dvar.GetBin("njet+thadpt", genper->NAddJets(), genthad.Pt())-0.5, weight);
				pdfunc->Fill1d("pdfunc_njet+ttpt", response2dvar.GetBin("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt())-0.5, weight);
				pdfunc->Fill1d("pdfunc_njet+ttm", response2dvar.GetBin("njet+ttm", genper->NAddJets(), (gent+gentbar).M())-0.5, weight);
				pdfunc->Fill1d("pdfunc_thady+thadpt", response2dvar.GetBin("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt())-0.5, weight);
				pdfunc->Fill1d("pdfunc_ttpt+thadpt", response2dvar.GetBin("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt())-0.5, weight);
				pdfunc->Fill1d("pdfunc_thadpt+ttm", response2dvar.GetBin("thadpt+ttm", genthad.Pt(), (gent+gentbar).M())-0.5, weight);
				pdfunc->Fill1d("pdfunc_ttpt+ttm", response2dvar.GetBin("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M())-0.5, weight);
				pdfunc->Fill1d("pdfunc_ttm+tty", response2dvar.GetBin("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()))-0.5, weight);
				pdfunc->Fill1d("pdfunc_ttm+dy", response2dvar.GetBin("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()))-0.5, weight);
				for(size_t n = 0 ; n < genper->NJets() ; ++n)
				{
					pdfunc->Fill1d("pdfunc_jet+jetpt", response2dvar.GetBin("jet+jetpt", n, genper->GetJet(n)->Pt())-0.5, weight);
					pdfunc->Fill1d("pdfunc_jet+jetptOF", response2dvar.GetBin("jet+jetptOF", n, genper->GetJet(n)->Pt())-0.5, weight);
					pdfunc->Fill1d("pdfunc_jet+jeteta", response2dvar.GetBin("jet+jeteta", n, abs(genper->GetJet(n)->Eta()))-0.5, weight);
					pdfunc->Fill1d("pdfunc_jet+jetdr", response2dvar.GetBin("jet+jetdr", n, genper->DRminTTjets(genper->GetJet(n)))-0.5, weight);
					pdfunc->Fill1d("pdfunc_jet+jetdrtop", response2dvar.GetBin("jet+jetdrtop", n, genper->DRminTop(genper->GetJet(n)))-0.5, weight);
				}
			}
		}
		if(SEMILEPACC)
		{
			ttp_genacc.Fill(*genper, weight);
			truth1d["counter"]->Fill(2.5, weight);
		}

		SelectRecoParticles(event);

		if(isDA && Abs(event.trigger().HLT_IsoMu24()) != 1) {cout << "TRIGGER UNDEFINED IsoMu24:" << event.trigger().HLT_IsoMu24() << endl; }
		if(isDA && Abs(event.trigger().HLT_IsoTkMu24()) != 1) {cout << "TRIGGER UNDEFINED: TKMu24" << event.trigger().HLT_IsoTkMu24() << endl; }
		if(isDA && Abs(event.trigger().HLT_Ele27_WPTight_Gsf()) != 1) {cout << "TRIGGER UNDEFINED EL:" <<  event.trigger().HLT_Ele27_WPTight_Gsf() << endl; }
		if(
				(
				 isDA == 0
				 && (
				  event.trigger().HLT_IsoMu24() == 1 || event.trigger().HLT_IsoTkMu24() == 1
				  || (event.trigger().El27ptmax() != -1 && event.trigger().HLT_Ele27_WPTight_Gsf() == 1)
				 )
				) ||
				(
				 isDA == 13 &&
					(
						event.trigger().HLT_IsoMu24() == 1 || event.trigger().HLT_IsoTkMu24() == 1 //2016
					)
				 ) ||
				(
				 isDA == 11 &&
					(
						event.trigger().HLT_IsoMu24() == -1 && event.trigger().HLT_IsoTkMu24() == -1 && event.trigger().El27ptmax() != -1 && event.trigger().HLT_Ele27_WPTight_Gsf() == 1 //2016
				)
			)
		)
		{
			//cout << "sel " << event.run << " " << event.lumi << " " << event.evt << " " << -1 << " " << cleanedjets.size() << endl;
			//SelectRecoParticles(event);
			ttanalysis(event);
		}

if(STUDENT){stud_tr->Fill();}
	}

	outFile_.cd();
	TTree* lumitree = new TTree("runls", "runls");
	UInt_t run = 0;
	UInt_t ls = 0;
	lumitree->Branch("run",&run,"run/i");
	lumitree->Branch("lumisec",&ls,"lumisec/i");
	for(map<int, set<int> >::const_iterator ita = runinfo.begin() ; ita != runinfo.end() ; ++ita)
	{
		run = ita->first;
		for(set<int>::const_iterator itb = ita->second.begin() ; itb != ita->second.end() ; ++itb)
		{
			ls = *itb;
			lumitree->Fill();
		}
	}

}


//make it executable
int main(int argc, char *argv[])
{
	URParser &parser = URParser::instance(argc, argv);
	URDriver<ttbar> test;
	return test.run();
}
