#include "ttbarxsec.h"
#include "TRandom3.h"
#include <time.h>

#include "Permutation.h"
#include "Dilepton.h"
#include "PDFuncertainty.h"
#include "NeutrinoSolver.h"
#include "ConfigParser.h"
#include <TCanvas.h>
#include <TMarker.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/GSLMinimizer.h"
#include "DiNeutrinoSolver.h"

using namespace std;

ttbar::ttbar(const std::string output_filename):
	AnalyzerBase("ttbar", output_filename),
	FULLHAD(false),
	SEMILEP(false),
	DILEP(false),
	DILEPACC(false),
	SEMILEPACC(false),
	FULLLEP(false),
	gen1d("gen"),
	gen2d("gen"),
	nugraphs("nu"),
	right_bjets("rightb"),
	right_bjets_ave("rightb_ave"),
	right_bjets_lop("rightb_lop"),
//	right_bjets_smear1("rightb_s1"),
	right_bjets_smear2("rightb_s2"),
//	right_bjets_smear3("rightb_s3"),
	swapped_bjets("swappedb"),
	gen_test("gentest"),
	gen_MET("genMET"),
	nu_MET("nuMET"),
	gen_bl("genblep"),
	ttp_genall("genall"),
	ttp_genacc("genacc"),
	ttp_genjet("genjet"),
	reco1d("reco"),
	reco2d("reco"),
	truth1d("truth"),
	truth2d("truth"),
	ditruth1d("ditruth"),
	ditruth2d("ditruth"),
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
	cletamax(2.3),//max |eta| of leptons (max allowed value is 2.4)
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
	TUNEDOWN(false),
	MW(80.4),
	Mt(172.5),
	dilepton_event_counter(0),
	DILEP_count(0),
	DILEPACC_count(0),
	neutrino_graphs(0),
	neutrino_bins(200),
	smeargraphs(100),
	neutrino_step(twopi/neutrino_bins)
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
	ttmbins = {300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 1500.0, 2500.0};
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
	DiNeutrinoSolver::mini = ROOT::Math::Factory::CreateMinimizer("Minuit2");

	testcounter=0;
	outFile_.cd();
	TDirectory* dir_gen = outFile_.mkdir("GEN");
	dir_gen->cd();


	gen1d.AddHist("dileppttop", 500, 0., 1000, "p_{t}(t)+p_t(#bar{t}) [GeV]", "Events");
	gen1d.AddHist("tmass", 800, 0., 400, "M(t) [GeV]", "Events");
	gen1d.AddHist("dilep_tmass", 800, 0., 400, "M(t) [GeV]", "Events");
	gen1d.AddHist("dilepptb", 500, 0., 1000, "p_{t}(b)+p_t(#bar{b})  [GeV]", "Events");
	gen1d.AddHist("dilepptW1", 500, 0., 1000, "p_{t}(W1) [GeV]", "Events");
	gen1d.AddHist("dilepptW2", 500, 0., 1000, "p_{t}(W2) [GeV]", "Events");
	gen1d.AddHist("dilepptW", 500, 0., 1000, "p_{t}(W1)+p_t(W2})  [GeV]", "Events");
	gen1d.AddHist("neutrino_pt",500,0.,1000,"Neutrino Pt","Events");
	gen1d.AddHist("neutrino_momentum_vs_lep",500,0.,500,"Neutrino momentumDistance from Lep","Events");
	gen1d.AddHist("neutrino_momentum_vs_bj",500,0.,500,"Neutrino momentumDistance from bjet","Events");
	gen1d.AddHist("neutrino_dr_vs_lep",500,0.,10,"Neutrino dr from lep","Events");
	gen1d.AddHist("neutrino_dr_vs_bj",500,0.,10,"Neutrino dr from bjet","Events");
	gen1d.AddHist("lep_dr_vs_bj",500,0.,10,"lep dr from bjet","Events");
	gen2d.AddHist("neutrino_error_vs_IM",250,0.,250.,2,-0.5,1.5,"IM(b+lep)in [GeV]","0=no error, 1=error");
	gen1d.AddHist("neutrino_minnum",7,-0.5,6.5,"number of mins","events");
	gen1d.AddHist("neutrino_log_minchi2",7,-1.5,5.5,"log","events");

	gen1d.AddHist("TYP", 5, 0., 5., "Decay TYP", "Events");
	gen1d.AddHist("tpt", 500, 0, 1000, "p_{T}(t) [GeV]", "Events");
	gen1d.AddHist("ty", 500, 0, 10, "y(t)", "Events");
	gen1d.AddHist("ttpt", 500, 0, 1000, "p_{T}(t#bar{t}) [GeV]", "Events");
	gen1d.AddHist("tty", 500, 0, 10, "y(t#bar{t})", "Events");

	gen1d.AddHist("nu_Pt",50,0,300.,"Nu Pt","Events");
	gen1d.AddHist("nu_Pz",50,0,300.,"Nu Pt","Events");
	gen1d.AddHist("nu_P",50,0,500.,"Nu Pt","Events");


	gen2d.AddHist("MDY_bl", 200, 0., 1000,3,0.,3., "M bl", "DeltaY");
	gen2d.AddHist("MDY_tt", 200, 0., 1000,3,0.,3., "M ytt", "DeltaY");
	gen1d.AddHist("M_tt", 180, 100., 1000, "M(tt)", "Dilepton Events");
	gen1d.AddHist("DY_tt", 140,-3.5,3.5, "DeltaY(tt)", "Dilepton Events");
	gen1d.AddHist("M_bl", 180, 100., 1000, "M(bs and leps)", "Dilepton Events");
	gen1d.AddHist("DY_bl", 140,-3.5,3.5, "DeltaY(bs and leps)", "Dilepton Events");



	for (size_t i = 0; i < 6; i++) {
		gen2d.AddHist("MDY_bl_HATHOR_y"+to_string(i), 200, 0., 1000,3,0.,3., "Mbl [GeV] after reweighting", "DY bin");
		gen2d.AddHist("MDY_tt_HATHOR_y"+to_string(i), 200, 0., 1000,3,0.,3., "Mtt [GeV] after reweighting", "DY bin");

		gen1d.AddHist("M_tt_HATHOR_y"+to_string(i), 180,100., 1000, "M(tt) [GeV] after reweighting", "Dilepton Events");
		gen1d.AddHist("DY_tt_HATHOR_y"+to_string(i), 140, -3.5,3.5, "DeltaY(tt) after reweighting", "Dilepton Events");
		gen1d.AddHist("M_bl_HATHOR_y"+to_string(i), 180, 100., 1000, "M(bs and leps) [GeV] after reweighting", "Dilepton Events");
		gen1d.AddHist("DY_bl_HATHOR_y"+to_string(i), 140, -3.5,3.5, "DeltaY(bs and leps) after reweighting", "Dilepton Events");
	}

	gen2d.AddHist("ttpt_vs_IM",200,0,2000,200,0,600,"M_tt","ttpt");

    // ttp_genall.Init(this);


   TDirectory* dir_nu = outFile_.mkdir("NUGRAPHS");
	dir_nu->cd();
	for (int i = 0; i <= neutrino_graphs; ++i)
	{
		nugraphs.AddNeutrinoHist("neutrino_solver"+to_string(i)+"_gen",neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,"t_{1}","t_{2}");
		nugraphs.AddNeutrinoHist("neutrino_solver"+to_string(i)+"_right",neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,"t_{1}","t_{2}");
		nugraphs.AddNeutrinoHist("neutrino_solver"+to_string(i)+"_swapped",neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,neutrino_bins+1,0-neutrino_step/2,twopi+neutrino_step/2,"t_{1}","t_{2}");
		nugraphs["neutrino_solver"+to_string(i)+"_gen"]->SetTitle("(gen) Neutrino Solver #Delta P_{xy}^{2} (MET vs #nu#bar{#nu})");
	}

	TDirectory* dir_ditruth = outFile_.mkdir("DITRUTH");
	dir_ditruth->cd();
	//Evan's hists
	right_bjets.Init(this);
	right_bjets_ave.Init(this);
	right_bjets_lop.Init(this);
//	right_bjets_smear1.Init(this);
	right_bjets_smear2.Init(this);
//	right_bjets_smear3.Init(this);
	swapped_bjets.Init(this);
	gen_test.Init(this);
	gen_MET.Init(this);
	nu_MET.Init(this);
	gen_bl.Init(this);

	truth1d.AddHist("jet_number_dilepACC",7,1.5,8.5,"number of jets","Events");

	truth1d.AddHist("min_deltaPxPy_success_vs_fail",2,-0.5,1.5, "0 = wrong b pairing, 1 = right b pairing", "Events");
	truth1d.AddHist("num_mins_guess",3,-1.5,1.5, "-1: swapped has more mins 0 = same, 1 = right has more mins", "Events");
	truth1d.AddHist("gen_vs_reco_MET_deltaPxy",25,0,150,"delta P","Events");
	truth1d.AddHist("gen_MET_vs_Nus_deltaPxy",25,0,150,"delta P","Events");
	truth2d.AddHist("cleanedjets_vs_bestmin_deltaPxy",100,0,100,25,1.5,6.5,"delta Pxy", "total jets");
	truth2d.AddHist("cleanedjets_vs_bestmin_deltaP",100,0,400,25,1.5,6.5,"delta P", "total jets");

	truth1d.AddHist("bjet_match",4,-0.5,3.5,"number of jets matching gen bjets","Events");
	truth1d.AddHist("bjet_match_tag",4,-0.5,3.5,"number of b-tagged jets matching gen bjets","Events");
	truth1d.AddHist("bjet_match_tag_fraction",4,-0.5,3.5,"number of b-tagged jets matching gen bjets given that a matching pair exists","Events");
	truth1d.AddHist("bjet_ID",3,-0.5,2.5,"how many of the 2 reco bjets match gen data?","Events");
	truth1d.AddHist("bjet_match_PERCENT",3,-0.5,2.5,"best guess match ID","Events");
	truth2d.AddHist("jetmatching",500,0.,500.,30,-0.5,2.5,"p_{T} in [GeV]","How many jets match gen data after guess");
	truth2d.AddHist("jetIDvsPt",200,0.,200.,30,-0.5,2.5,"p_{T} in [GeV]","0=not found, 1=detected, 2=detected+tagged");
	truth1d.AddHist("Neutrino_errors_right_jetl",2,-0.5,1.5,"0 = no error, 1 = error","Events");
	truth1d.AddHist("Neutrino_errors_wrong_jetl",2,-0.5,1.5,"0 = no error, 1 = error","Events");
	truth1d.AddHist("best_mindist",100,0.,1500.,"DeltaP1 + deltaP2", "events");

	truth1d.AddHist("old_vs_new_nusolve",250,0.,50.,"delta_p(old,new)","solutions");
	truth2d.AddHist("old_vs_new_nusolve2d",250,0.,50.,5,-0.5,4.5,"delta_p(old,new)","nsolns");

	truth1d.AddHist("NUSMEAR1_deltaP",50,0,400,"smearing deltaP","Events");
	truth1d.AddHist("NUSMEAR2_deltaP",50,0,400,"smearing deltaP","Events");
	truth1d.AddHist("NUSMEAR3_deltaP",50,0,400,"smearing deltaP","Events");


	truth2d.AddHist("deltaP_vs_Pt",100,0.,250.,100,0.,500.,"Pt","deltaP from gen data");


	truth1d.AddHist("jetIDerr_genbPt",500, 0, 1000, "gen b p_{T} [GeV]", "Events");
	truth1d.AddHist("jetIDerr_recobPt",500, 0, 1000, "reco b p_{T} [GeV]", "Events");

	truth2d.AddHist("bjet_matching_mindeltaR",500,0.,10.,3,-0.5,2.5,"min delta R","# of correct assigned jets");
	truth2d.AddHist("bjet_matching_minchi2",80,-1.5,7.5,6,-0.5,2.5,"min chi 2","# of correct assigned jets");

	truth1d.AddHist("compare_sum_Pt",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_Pz",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_E",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_deltaR_lep",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_deltaR_lepandb",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_IMW",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_sum_blep_deltaR",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");
	truth1d.AddHist("compare_min_pxpy",2,-0.5,1.5,"0 = bigger for swapped, 1 = bigger for correct pairing","events");


	truth1d.AddHist("tt_IM_resolution",25,0,100,"IM difference fron GEN data","events");


	truth1d.AddHist("bj_p_mag",60,0,300,"just the momentum of the b jets","events");
	truth1d.AddHist("bj_resolution_theta",30,0,0.2,"DeltaTheta between GEN and RECO bjets","events");
	truth1d.AddHist("bj_resolution_px",50,-.5,.5,"DeltaPx between GEN and RECO bjets","events");
	truth1d.AddHist("bj_resolution_py",50,-.5,.5,"DeltaPy between GEN and RECO bjets","events");
	truth1d.AddHist("bj_resolution_pz",50,-.5,.5,"DeltaPz between GEN and RECO bjets","events");
	truth1d.AddHist("bj_resolution_E",50,-.6,.6,"DeltaPz between GEN and RECO bjets","events");
	truth1d.AddHist("bj_resolution_vs_pt",50,-1,1,"DeltaPx between GEN and RECO bjets","events");


	truth1d.AddHist("DeltaY_error_case_1",60,-.6,.6,"DeltaY error case 1 : naive reco","events");
	truth1d.AddHist("DeltaY_error_case_2",60,-.6,.6,"DeltaY error case 2 : reco nusolver","events");
	truth1d.AddHist("DeltaY_error_case_3",60,-.6,.6,"DeltaY error case 3 : nusolver runs on gen info","events");


	truth1d.AddHist("DeltaY_error_case_2_1s",60,-.6,.6,"DeltaY error case 2 : reco nusolver","events");
	truth1d.AddHist("DeltaY_error_case_2_2s",60,-.6,.6,"DeltaY error case 2 : reco nusolver","events");
	truth1d.AddHist("DeltaY_error_case_2_4s",60,-.6,.6,"DeltaY error case 2 : reco nusolver","events");

	truth2d.AddHist("Mtt_error_vs_pt",100,0,200,50,0,.6,"Pt (GeV)","Mtt Error");
	truth2d.AddHist("Mtt_error_vs_pz",100,0,200,50,0,.6,"Pz (GeV)","Mtt Error");
	truth2d.AddHist("dy_error_vs_pt",100,0,200,50,0,.6,"Pt (GeV)","delta Y Error");
	truth2d.AddHist("dy_error_vs_pz",100,0,200,50,0,.6,"Pz (GeV)","delta Y Error");

	truth2d.AddHist("reco_MDY_bl", 200, 0., 1000,3,0.,3., "M bl", "DY");
	truth2d.AddHist("reco_MDY_tt", 200, 0., 1000,3,0.,3., "Mtt", "DY");
	truth2d.AddHist("ideal_MDY_tt", 200, 0., 1000,3,0.,3., "Mtt", "DY");


	vector<double> MDYbins = {0,300,360,400,440,480,540,2000,  2300,2400,2440,2480,2520,2580,4000,  4300,4500,4560,4620,4700,6000};
	vector<double> MDYbins_bl = {0,100,200,250,300,400,520,2000,  2200,2260,2320,2400,2480,2600,4000,  4250,4350,4440,4520,4620,6000};

	vector<double> mbins = {0,300,360,380,420,460,520,2000};
	vector<double> mbins_bl = {0,100,200,250,300,400,520,2000};

	truth1d.AddHist("idealnu_M_tt", mbins, "M(tt)", "Dilepton Events");
	truth1d.AddHist("idealnu_DY_tt", 40, -4,4, "DeltaY(tt)", "Dilepton Events");
	truth1d.AddHist("reco_M_tt", mbins, "M(tt)", "Dilepton Events");
	truth1d.AddHist("reco_DY_tt", 40, -4,4, "DeltaY(tt)", "Dilepton Events");
	truth1d.AddHist("reco_M_blmet", mbins_bl, "M(tt) [GeV] after reweighting", "Dilepton Events");
	truth1d.AddHist("reco_M_bl", mbins_bl, "M(bs and leps)", "Dilepton Events");
	truth1d.AddHist("reco_DY_bl", 40, -4,4, "DeltaY(bs and leps)", "Dilepton Events");


	for (size_t i = 0; i < 6; i++) {


		truth2d.AddHist("reco_MDY_bl_HATHOR_y"+to_string(i), 200,0,2000,3,0,3, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth2d.AddHist("reco_MDY_tt_HATHOR_y"+to_string(i), 200,0,2000,3,0,3, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth2d.AddHist("smearnu_MDY_tt_HATHOR_y"+to_string(i), 200,0,2000,3,0,3, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth2d.AddHist("avenu_MDY_tt_HATHOR_y"+to_string(i), 200,0,2000,3,0,3, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth2d.AddHist("idealnu_MDY_tt_HATHOR_y"+to_string(i), 200,0,2000,3,0,3, "Mtt [GeV] after reweighting", "Dilepton Events");

		truth1d.AddHist("idealnu_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("idealnu_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu1_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu1_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu2_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu2_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu3_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("smearnu3_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("lonu_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("lonu_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("avenu_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "Mtt [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("avenu_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");

		truth1d.AddHist("reco_M_tt_HATHOR_y"+to_string(i), 200,0,2000, "M(tt) [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("reco_DY_tt_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(tt) after reweighting", "Dilepton Events");
		truth1d.AddHist("reco_M_blmet_HATHOR_y"+to_string(i),200,0,2000, "M(tt) [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("reco_M_bl_HATHOR_y"+to_string(i), 200,0,2000, "M(bs and leps) [GeV] after reweighting", "Dilepton Events");
		truth1d.AddHist("reco_DY_bl_HATHOR_y"+to_string(i), 40, -4,4, "DeltaY(bs and leps) after reweighting", "Dilepton Events");
	}

	truth1d.AddHist("Nu_solver_overshoot_1", 100,-2.5,2.5, "nu solver P() vs gen nu P()", "Dilepton Events");
	truth1d.AddHist("Nu_solver_overshoot_2",100,-2.5,2.5, "nu solver P() vs gen nu P()", "Dilepton Events");
	truth1d.AddHist("Nu_solver_overshoot_3",100,-2.5,2.5,"nu solver P() vs gen nu P()", "Dilepton Events");
	truth1d.AddHist("LowPNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPtNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPzNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPsNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPs2NuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPmaxNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowMttNuIndex",4,0,4,"low index after sorting","events");
	truth1d.AddHist("LowPNu_percentChange",50,0,2.5,"fractional difference: lowPnu vs bestNu","events");

	truth2d.AddHist("bestnu_IMerr_vs_IMgen",200,0,2000,100,0,3,"Mtt at gen level","relative error solver ");
	truth2d.AddHist("bestnu_IMerr_vs_IMreco",200,0,2000,100,0,3,"Mtt at gen level","relative error solver ");
	truth2d.AddHist("bestnu_IMdiff_vs_IMgen",200,0,2000,100,0,1000,"Mtt at gen level","relative error solver ");
	truth2d.AddHist("bestnu_IMdiff_vs_IMreco",200,0,2000,100,0,1000,"Mtt at gen level","relative error solver ");

	truth2d.AddHist("b_smear_Mtt_err_vs_stdev",200,0,300,100,0,3,"Standard Dev from smearing", "relative error");
	truth2d.AddHist("b_smear_Mtt_diff_vs_stdev",200,0,300,100,0,1200,"Standard Dev from smearing", "difference from truth");
	truth2d.AddHist("met_smear_Mtt_err_vs_stdev",200,0,300,100,0,3,"Standard Dev from met smearing", "relative error");
	truth2d.AddHist("met_smear_Mtt_diff_vs_stdev",200,0,300,100,0,1000,"Standard Dev from met smearing", "difference from truth");
	truth2d.AddHist("total_smear_Mtt_err_vs_stdev",200,0,300,100,0,3,"Standard Dev from met smearing", "relative error");
	truth2d.AddHist("total_smear_Mtt_diff_vs_stdev",200,0,300,100,0,1000,"Standard Dev from met smearing", "difference from truth");

	truth1d.AddHist("smear2_Mtt",200,0.,2000,"Mtt smear error","events");
	truth1d.AddHist("b_smear_Mtt",200,0.,2000,"Mtt smear error","events");
	truth1d.AddHist("met_smear_Mtt",200,0.,2000,"Mtt smear error","events");

	// for (size_t i = 0; i < smeargraphs+1; i++)
	// {
	// 	truth1d.AddHist("smear_Mtt"+to_string(i),200,0.,2000,"Mtt smear error","events");
	// 	truth1d.AddHist("smear_sumP"+to_string(i),200,0,1000,"Mtt smear error","events");
	// 	truth1d.AddHist("smear_sumDP"+to_string(i),200,0,500,"Mtt smear error","events");
	// }
	truth1d.AddHist("Mtterr_best_gen",200,-1000.,1000,"Mtt error best","events");
	truth1d.AddHist("NuPerr_best_gen",200,-600.,600,"sum Nu P","events");

	truth1d.AddHist("Mtterr_ave_gen",200,-1000.,1000,"Mtt error ave","events");
	truth1d.AddHist("Mtterr_low_gen",200,-1000.,1000,"Mtt error low stdev","events");



	truth2d.AddHist("Mtterr_vs_jeterror",200,0,250,100,0,1200,"jet error (deltaP)", "IM difference from truth");
	truth2d.AddHist("Mtterr_vs_meterror",200,0,200,100,0,1200,"met error (deltaP)", "IM difference from truth");
	truth2d.AddHist("Nuerr_vs_ttpt",200,0,300,100,0,300,"ttpt error (deltaP)", "IM difference from truth");
	truth1d.AddHist("Meterr",100,0,200,"met error (deltaP)","events");
	truth1d.AddHist("m_top_gen_postselection",100,150,200,"title","events");
	truth1d.AddHist("m_w_gen_postselection",100,65,95,"title","events");
	truth2d.AddHist("ptop_err_vs_topoffshell",100,0,75,100,0,3,"how far off shell","top delta P gen vs ideal rel err");
	truth2d.AddHist("ptop_err_vs_TWoffshell",100,0,100,100,0,1.0,"how far off shell T+W","top delta P gen vs ideal rel err");
	truth2d.AddHist("ptop_err_vs_TWoffshell_alt",100,0,100,100,0,300,"how far off shell","top delta P gen vs ideal rel err");
	truth2d.AddHist("pnu_err_vs_TWoffshell",100,0,100,100,0,3,"how far off shell T+W","top delta P gen vs ideal rel err");
	truth2d.AddHist("pnu_err_vs_TWoffshell_alt",100,0,100,100,0,300,"how far off shell T+W","top delta P gen vs ideal rel err");

	truth2d.AddHist("pnu_err_vs_noncon_E",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation E","relative error: nus P");
	truth2d.AddHist("pnu_err_vs_noncon_p",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation p","relative error: nus P");
	truth2d.AddHist("top_err_vs_noncon_E",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation E","relative error: tops P");
	truth2d.AddHist("top_err_vs_noncon_p",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation p","relative error: tops P");
	truth2d.AddHist("pnu_err_vs_noncon_EP",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation E+p","relative error: nus P");
	truth2d.AddHist("top_err_vs_noncon_EP",100,0,50,100,0,2.0,"ttbar 4-momentum non-conservation E+p","relative error: tops P");
	truth2d.AddHist("pnu_err_vs_gen_allsource",100,0,50,100,0,2.0,"all gen err sources","relative error: nus P");
	truth2d.AddHist("top_err_vs_gen_allsource",100,0,50,100,0,2.0,"all gen err sources","relative error: tops P");




	truth1d.AddHist("nu_test",100,0,300,"title","events");
	truth1d.AddHist("nu_testb",100,0,300,"title","events");
	truth1d.AddHist("TTPT_test",100,0,150,"title","events");

	truth1d.AddHist("MET_res_x",100,-2,2,"x res (fraction)","events");
	truth1d.AddHist("MET_res_y",100,-2,2,"y res (fraction)","events");
	truth1d.AddHist("MET_res2_x",100,-2,2,"x res (fraction)","events");
	truth1d.AddHist("MET_res2_y",100,-2,2,"y res (fraction)","events");
	truth1d.AddHist("MET_diff_x",200,-250,250,"x res (fraction)","events");
	truth1d.AddHist("MET_diff_y",200,-250,250,"y res (fraction)","events");
	truth1d.AddHist("MET_res_E",100,-2,2,"E res (fraction)","events");
	truth1d.AddHist("MET_res_theta",100,0,4,"theta res ","events");
	truth2d.AddHist("MET_res_y_vs_x",100,0,2,100,0,2,"x err ","y err");








	TDirectory* dir_truth = outFile_.mkdir("TRUTH");
	dir_truth->cd();


	truth1d.AddHist("counter", 20, 0., 20., "counter", "Events");
	truth1d.AddHist("Mu", 100, 0, 100, "#mu", "Events");
	truth1d.AddHist("MuWeighted", 100, 0, 100, "#mu", "Events");


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


	double bkgcutmin = 0.0;
	double bkgcutmax = 0.6;
	if(cnbtag == "bkgl") {bkgcutmin = 0.0; bkgcutmax = 0.3;}
	if(cnbtag == "bkgh") {bkgcutmin = 0.4; bkgcutmax = 0.7;}
	btageff.Init(B_MEDIUM, bkgcutmin, bkgcutmax);

	TDirectory* dir_reco = outFile_.mkdir("RECO");
	dir_reco->cd();
	reco1d.AddHist("counter", 20, 0., 20., "counter", "Events");

//	ttp_all.Init(this);

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
	TFile* fHAT_0 = TFile::Open("yukawa_reweighting0.0y.root");
	HAThists.push_back((TH2D*)fHAT_0->Get("EWtoLO"));
	TFile* fHAT_1 = TFile::Open("yukawa_reweighting1.0y.root");
	HAThists.push_back((TH2D*)fHAT_1->Get("EWtoLO"));
	TFile* fHAT_2 = TFile::Open("yukawa_reweighting2.0y.root");
	HAThists.push_back((TH2D*)fHAT_2->Get("EWtoLO"));
	TFile* fHAT_3 = TFile::Open("yukawa_reweighting3.0y.root");
	HAThists.push_back((TH2D*)fHAT_3->Get("EWtoLO"));
	TFile* fHAT_4 = TFile::Open("yukawa_reweighting4.0y.root");
	HAThists.push_back((TH2D*)fHAT_4->Get("EWtoLO"));
	TFile* fHAT_5 = TFile::Open("yukawa_reweighting5.0y.root");
	HAThists.push_back((TH2D*)fHAT_5->Get("EWtoLO"));

	//HAThist = dynamic_cast<TH2D*>(fHAT->Get("EWtoLO"));

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
	DILEP = false;
	FULLHAD = false;

	const vector<Ttgen>& ttgens = event.TTGens();
	if(ttgens.size() != 8) { return; }
	for(size_t n = 0 ; n < 8 ; ++n) { gps[n] = ttgens[n]; }
	gentq = gps[0];
	gentqbar = gps[4];

	weight *= 1.+cttptweight*((gentq + gentqbar).Pt()-200.)/2000.;
	weight *= 1.+ctopptweight*(gentq.Pt()-200.)/1000.;
	weight *= 1.+ctoprapweight*(0.2-0.2*Abs(gentq.Rapidity()));
	bool tlep = false;
	bool tbarlep = false;
	bool tem = false;
	bool tbarem = false;
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
		//SEMILEPACC = true;
		//genallper.Init(&gps[6], &gps[7], &gps[5], &gps[1], &gps[2], gps[2].pdgId(), gps[3]);
		// gentqhad = gentqbar;
		// gentqlep = gentq;

	}
	else if(!tlep && tbarem)
	{
		//SEMILEPACC = true;
		//genallper.Init(&gps[2], &gps[3], &gps[1], &gps[5], &gps[7], gps[7].pdgId(), gps[6]);
		// gentqhad = gentq;
		// gentqlep = gentqbar;
	}

	if(tlep && tbarlep)
	{
		if(abs(gps[1].Eta())>2.4 || abs(gps[5].Eta())>2.4 || gps[1].Pt() <35. || gps[5].Pt() < 35.){return;}
		if(abs(gps[2].Eta())>2.4 || abs(gps[7].Eta())>2.4 || gps[2].Pt() <25. || gps[7].Pt() < 25.){return;}


		DILEP=true;
		DILEP_count++;

		gendilep.InitGen(&gentq,&gentqbar,&gps[2],&gps[7],&gps[1],&gps[5],gps[2].pdgId(),gps[7].pdgId(),&gps[3],&gps[6]);

//		cout<<endl<<"{"<<gendilep.MET()->Pt()<<", "<<gendilep.NuMET()->Pt()<<"} ";


		idealdilep.Init(&gps[2],&gps[7],&gps[1],&gps[5],gps[2].pdgId(),gps[7].pdgId(),gendilep.NuMET());
		gen_bl_dilep.Init(&gps[2],&gps[7],&gps[1],&gps[5],gps[2].pdgId(),gps[7].pdgId(),&gendilep.NuNu());

		TLorentzVector W1 = gps[2] + gps[3];
		TLorentzVector W2 = gps[7] + gps[6];

		gen1d["dilep_tmass"]->Fill(gentq.M(), weight);
		gen1d["dilep_tmass"]->Fill(gentqbar.M(), weight);
		gen1d["dileppttop"]->Fill(gentq.Pt()+gentqbar.Pt(), weight);
		gen1d["dilepptW1"]->Fill(W1.Pt(), weight);
		gen1d["dilepptW2"]->Fill(W2.Pt(), weight);
		gen1d["dilepptW"]->Fill(W1.Pt()+W2.Pt(), weight);

		gen1d["nu_Pt"]->Fill(gendilep.Nut()->Pt(),weight);
		gen1d["nu_Pt"]->Fill(gendilep.Nutbar()->Pt(),weight);
		gen1d["nu_Pz"]->Fill(gendilep.Nut()->Pz(),weight);
		gen1d["nu_Pz"]->Fill(gendilep.Nutbar()->Pz(),weight);
		gen1d["nu_P"]->Fill(gendilep.Nut()->P(),weight);
		gen1d["nu_P"]->Fill(gendilep.Nutbar()->P(),weight);

		gen1d["neutrino_momentum_vs_lep"]->Fill(deltaP(*gendilep.Nut(),*gendilep.Lt()),weight);
		gen1d["neutrino_momentum_vs_bj"]->Fill(deltaP(*gendilep.Nut(),*gendilep.Bt()),weight);
		gen1d["neutrino_dr_vs_lep"]->Fill(gendilep.Nut()->DeltaR(*gendilep.Lt()),weight);
		gen1d["neutrino_dr_vs_bj"]->Fill(gendilep.Nut()->DeltaR(*gendilep.Bt()),weight);
		gen1d["lep_dr_vs_bj"]->Fill(gendilep.Lt()->DeltaR(*gendilep.Bt()),weight);

		NeutrinoSolver genntguess;
		genntguess.Build(gendilep.Lt(),gendilep.Bt(),MW,Mt);
		NeutrinoSolver genntbarguess;
		genntbarguess.Build(gendilep.Ltbar(),gendilep.Btbar(),MW,Mt);

		int nterr = 0;
		int ntbarerr = 0;
		if (genntguess.GetSolution(0)[0]==0)
		{
			nterr++;
		}
		if (genntbarguess.GetSolution(0)[0]==0)
		{
			ntbarerr++;
		}

		gen2d["neutrino_error_vs_IM"]->Fill((*gendilep.Lt()+*gendilep.Bt()).M(),nterr,weight);
		gen2d["neutrino_error_vs_IM"]->Fill((*gendilep.Ltbar()+*gendilep.Btbar()).M(),ntbarerr,weight);

	}
	else if(!tlep && !tbarlep)
	{
		FULLHAD = true;
	}
	else
	{
		SEMILEP = true;
	}

}

void ttbar::SelectRivetPS(URStreamer& event){}

void ttbar::SelectPseudoTopLHC(URStreamer& event){}

void ttbar::SelectPseudoTop(URStreamer& event){}



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


//Electrons
	const vector<Electron>& electrons = event.electrons();
	for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron)
	{
		IDElectron el(*electron);
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

//Jets
	TVector2 metcorr(0.,0.);
	const vector<Jet>& jets = event.jets();
	for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
	{
		IDJet jet(*jetit);
		if(!jet.ID() || !jet.Clean(loosemuons, looseelectrons)) {continue;}
		double sfres = jetscaler.GetRes(jet, event.rho().value(), cjetres);
		double sfgen = sfres * jetscaler.GetScale(jet, csigmajet);
		double sf_minus = sfres * jetscaler.GetScale(jet, -1);
		double sf_plus = sfres * jetscaler.GetScale(jet, 1);

		//cout<<"sf_minus="<<sf_minus<<", sf_gen="<<sfgen<<", sf_plus="<<sf_plus<<", metcorr="<<endl;

		double sfwj = 0.;
		if(csigmajetwj != 1000.) {sfwj = sfres * jetscaler.GetScale(jet, csigmajetwj);}
		else{ sfwj = sfres * (1. + cscalejetwj);}
		// jet.SetSF(1, sfres, metcorr);
		// jet.SetSF(2, sfgen, metcorr);
		// jet.SetSF(3, sfwj, metcorr);
		jet.SetSF(1, sf_minus, metcorr);
		jet.SetSF(2, sfgen, metcorr);
		jet.SetSF(3, sf_plus, metcorr);
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


}

void ttbar::ttanalysis(URStreamer& event)
{
	reco1d["counter"]->Fill(0.5, weight);

//##################
//check for leptons:
//##################

	if(!DILEP){return;}

	if(tightmuons.size() == 2 && loosemuons.size() == 2 && looseelectrons.size() == 0)
	{
		mumu = true;
		int mu0_charge = tightmuons[0]->charge();
		lt = (mu0_charge > 0 ? dynamic_cast<TLorentzVector*>(tightmuons[0]) : dynamic_cast<TLorentzVector*>(tightmuons[1]));
		ltpdgid = mu0_charge*-13;
		ltbar = (mu0_charge < 0 ? dynamic_cast<TLorentzVector*>(tightmuons[0]) : dynamic_cast<TLorentzVector*>(tightmuons[1]));
		ltbarpdgid = mu0_charge*13;
	}
	else if(tightmuons.size() == 1 && loosemuons.size() == 1 && mediumelectrons.size() == 1 && looseelectrons.size() == 1)
	{
		emu = true;
		int mu_charge = tightmuons[0]->charge();
		lt = (mu_charge > 0 ? dynamic_cast<TLorentzVector*>(tightmuons[0]) : dynamic_cast<TLorentzVector*>(mediumelectrons[0]));
		ltpdgid = mu_charge*-11;
		ltbar = (mu_charge < 0 ? dynamic_cast<TLorentzVector*>(tightmuons[0]) : dynamic_cast<TLorentzVector*>(mediumelectrons[0]));
		ltbarpdgid = mu_charge*11;
	}
	else if(loosemuons.size() == 0 && mediumelectrons.size() == 2 && looseelectrons.size() == 2)
	{
		ee = true;
		int e0_charge = mediumelectrons[0]->charge();
		lt = (e0_charge > 0 ? dynamic_cast<TLorentzVector*>(mediumelectrons[0]) : dynamic_cast<TLorentzVector*>(mediumelectrons[1]));
		ltpdgid = e0_charge*-11;
		ltbar = (e0_charge < 0 ? dynamic_cast<TLorentzVector*>(mediumelectrons[0]) : dynamic_cast<TLorentzVector*>(mediumelectrons[1]));
		ltbarpdgid = e0_charge*11;
	}

	RECOISDILEP = (ee || mumu || emu);
	if(!RECOISDILEP){return;}
	recoleps={lt,ltbar};

//##################
//check for b jets:
//##################

	int genbjsmatch=0;
	int btagsmatch=0;
	vector<int> jetsthatmatch;
	bool btfound = false;
	bool bttag = false;
	bool btbarfound = false;
	bool btbartag = false;

	for (size_t j = 0; j < cleanedjets.size(); ++j)
	{
		if (cleanedjets[j]->DeltaR(*gendilep.Bt())< .3)
		{
			jetsthatmatch.push_back(j);
			btfound=true;
			if (cleanedjets[j]->csvIncl() > B_MEDIUM) {btagsmatch++; bttag=true;}
		}
		if (cleanedjets[j]->DeltaR(*gendilep.Btbar()) < .3)
		{
			jetsthatmatch.push_back(j);
			btbarfound=true;
			if (cleanedjets[j]->csvIncl() > B_MEDIUM) {btagsmatch++; btbartag=true;}
		}
	}

	truth2d["jetIDvsPt"]->Fill(gendilep.Bt()->Pt(),btfound+bttag);
	truth2d["jetIDvsPt"]->Fill(gendilep.Btbar()->Pt(),btbarfound+btbartag);

	truth1d["bjet_match"]->Fill(jetsthatmatch.size(),weight);
   truth1d["bjet_match_tag"]->Fill(btagsmatch,weight);

    if (jetsthatmatch.size()>1)
    {
    	truth1d["bjet_match_tag_fraction"]->Fill(btagsmatch,weight);
    }

    for (IDJet* j: cleanedjets)
    {
    	if (j->csvIncl() > B_MEDIUM)
    	{
    		recobjets.push_back(j);
    	}
    	else
    	{
    		addjets.push_back(j);
    	}
    }

//2 b-jets only, please
	if(recobjets.size() !=2 ){return;}


	int recobsuccess =0;

	for (size_t j = 0; j < 2; ++j)
	{
		if (recobjets[j]->DeltaR(*gendilep.Bt())< .3)
		{
			recobsuccess++;
		}
		if (recobjets[j]->DeltaR(*gendilep.Btbar()) < .3)
		{
			recobsuccess++;
		}
	}
  	truth1d["bjet_ID"]->Fill(recobsuccess,weight);

//#################
//# b-jet pairing #
//#################

//use deltaR method to set initial pairing guess. The NeutrinoSolver will allow us to improve this guess later

	double min_delta_R =10000000;
	for (bool j : {0,1})
	{
		for (bool l : {0,1})
		{
			double current_delta_R =(recobjets[j])->DeltaR(*recoleps[l])+(recobjets[1-j])->DeltaR(*recoleps[1-l]);
			if(current_delta_R < min_delta_R) {min_delta_R = current_delta_R; bjetmatch={j,l};}
		}
	}

//initialize reco Dilepton and its evil twin, which has the bjets swapped
	TLorentzVector* metTLV = dynamic_cast<TLorentzVector*>(&met);


	ltjet = (bjetmatch[0]==0 ? recobjets[bjetmatch[1]] : recobjets[!bjetmatch[1]]);
	ltbarjet = (bjetmatch[0]==1 ? recobjets[bjetmatch[1]] : recobjets[!bjetmatch[1]]);
	recodilep.Init(lt,ltbar,ltjet,ltbarjet,ltpdgid,ltbarpdgid,metTLV);
	recodilep.SetAddJets(addjets);
	recodilep_bswap.Init(lt,ltbar,ltbarjet,ltjet,ltpdgid,ltbarpdgid,metTLV);
	recodilep_bswap.SetAddJets(addjets);







	recodilep.SetAddJets(addjets);
	recodilep_bswap.SetAddJets(addjets);


	int jetIDguesscorrect=recodilep.BJetsCorrect(gendilep);
	int jetIDguessswapped=recodilep.BJetsSwapped(gendilep);

// Fill some plots and exit if the jets were not found correctly

	if (jetIDguesscorrect+jetIDguessswapped != 2)
	{
		if (!gendilep.BtCorrect(recodilep) && !gendilep.BtCorrect(recodilep_bswap) )
		{
			truth1d["jetIDerr_genbPt"]->Fill(gendilep.Bt()->Pt(),weight);
		}
		if (!gendilep.BtbarCorrect(recodilep) && !gendilep.BtCorrect(recodilep_bswap) )
		{
			truth1d["jetIDerr_genbPt"]->Fill(gendilep.Btbar()->Pt(),weight);
		}
		return;
	}

// set correct and "swapped" pairing dilepton objects for truth-level plotting

	truth2d["bjet_matching_mindeltaR"]->Fill(min_delta_R,jetIDguesscorrect,weight);
	if (jetIDguesscorrect==2)
	{
		rightbdilep = &recodilep;
		swappedbdilep = &recodilep_bswap;

	}
	else if (jetIDguessswapped==2)
	{
		rightbdilep = &recodilep_bswap;
		swappedbdilep = &recodilep;
	}
	else
	{
		//the code thinks the jets were ID'd right, but they weren't? check your catch radius
		cout<<"something unexpected has happened, check your code"; return;
	}

//fill b-jet related plots before moving on

	truth1d["jet_number_dilepACC"]->Fill(cleanedjets.size(),weight);
	truth2d["jetmatching"]->Fill(gendilep.Bt()->Pt()+gendilep.Btbar()->Pt(),jetIDguesscorrect);

	TVector3 p_bt_gen = gendilep.Bt()->Vect();
	TVector3 p_btbar_gen = gendilep.Btbar()->Vect();
	TVector3 p_bt_reco = rightbdilep->Bt()->Vect();
	TVector3 p_btbar_reco = rightbdilep->Btbar()->Vect();

	truth1d["bj_p_mag"]->Fill(p_bt_gen.Mag());
	truth1d["bj_p_mag"]->Fill(p_btbar_gen.Mag());

	truth1d["bj_resolution_theta"]->Fill(p_bt_gen.Angle(p_bt_reco)*100/p_bt_gen.Mag(),weight);
	truth1d["bj_resolution_theta"]->Fill(p_btbar_gen.Angle(p_btbar_reco)*100/p_btbar_gen.Mag(),weight);

	truth1d["bj_resolution_px"]->Fill((p_bt_gen.X()-p_bt_reco.X())/abs(p_bt_gen.X()),weight);
	truth1d["bj_resolution_px"]->Fill((p_btbar_gen.X()-p_btbar_reco.X())/abs(p_btbar_gen.X()),weight);

	truth1d["bj_resolution_py"]->Fill((p_bt_gen.Y()-p_bt_reco.Y())/abs(p_bt_gen.Y()),weight);
	truth1d["bj_resolution_py"]->Fill((p_btbar_gen.Y()-p_btbar_reco.Y())/abs(p_btbar_gen.Y()),weight);

	truth1d["bj_resolution_pz"]->Fill((p_bt_gen.Z()-p_bt_reco.Z())/abs(p_bt_gen.Z()),weight);
	truth1d["bj_resolution_pz"]->Fill((p_btbar_gen.Z()-p_btbar_reco.Z())/abs(p_btbar_gen.Z()),weight);

	truth1d["bj_resolution_E"]->Fill((gendilep.Bt()->E()-rightbdilep->Bt()->E())/abs(gendilep.Bt()->E()),weight);
	truth1d["bj_resolution_E"]->Fill((gendilep.Btbar()->E()-rightbdilep->Btbar()->E())/abs(gendilep.Btbar()->E()),weight);

	truth1d["bj_resolution_vs_pt"]->Fill((p_bt_gen-p_bt_reco).Z()/gendilep.Bt()->Pt(),weight);
	truth1d["bj_resolution_vs_pt"]->Fill((p_btbar_gen-p_btbar_reco).Z()/gendilep.Btbar()->Pt(),weight);


//###################
//# Neutrino Solver #
//###################


// some preliminary MET graphs

	truth1d["gen_vs_reco_MET_deltaPxy"]->Fill(deltaPxy(gendilep.NuNu(),*rightbdilep->MET()),weight);
	truth1d["gen_MET_vs_Nus_deltaPxy"]->Fill(deltaPxy(gendilep.NuNu(),*gendilep.NuMET()),weight);

// run the Neutrino solver

	idealdilep.Solve();
	recodilep.Solve();

	recodilep_bswap.Solve();

	rightbdilep_ave = *rightbdilep;
	rightbdilep_ave.AverageNeutrinos();


	Dilepton rightbdilep_lop = *rightbdilep;
	rightbdilep_lop.Solve();


	double Mtt_recoav = 0.;
	double Mtt_recolo = 0.;

	if (rightbdilep->MinNum()>0) {
		int lowpi = rightbdilep_lop.LowPsIndex();
		rightbdilep_lop.SetNeutrinos(lowpi);
		Mtt_recoav = rightbdilep_ave.Mtt();
		Mtt_recolo = rightbdilep_lop.Mtt();
	}

//#########################################
//testerino
//#########################################


	if(idealdilep.MinNum()>1)
	{
		idealdilep.SortSolutions(gendilep);
		idealdilep.SetNeutrinos(0);
		TVector3 p_top= gendilep.Top().Vect();
		TVector3 p_l= idealdilep.Lt()->Vect();
		TVector3 p_b= idealdilep.Bt()->Vect();
		TVector3 p_nu= idealdilep.Nut()->Vect();
		TVector3 p_trunu = gendilep.Nut()->Vect();

		TVector3 p_top_nu = p_l + p_b + p_nu;

		TVector3 p_top_bar= gendilep.Tbar().Vect();
		TVector3 p_l_bar= idealdilep.Ltbar()->Vect();
		TVector3 p_b_bar= idealdilep.Btbar()->Vect();
		TVector3 p_nu_bar= idealdilep.Nutbar()->Vect();
		TVector3 p_top_nu_bar = p_l_bar + p_b_bar +p_nu_bar;
		TVector3 p_trunu_bar = gendilep.Nutbar()->Vect();

		double p_top_err = (p_top-p_top_nu).Mag();
		double m_top_os = gendilep.Tbar().M() - 172.5;
		double m_w_os = (*gendilep.Lt()+*gendilep.Nut()).M() - 80.4;
		double p_top_err_bar = (p_top_bar-p_top_nu_bar).Mag();
		double m_top_os_bar = gendilep.Top().M() - 172.5;
		double m_w_os_bar = (*gendilep.Ltbar()+ *gendilep.Nutbar()).M() - 80.4;
		double p_nu_err = (p_trunu-p_nu).Mag();
		double p_nu_err_bar = (p_trunu_bar-p_nu_bar).Mag();

		TLorentzVector top_by_add = *gendilep.Lt() + *gendilep.Bt() + *gendilep.Nut();
		TLorentzVector tbar_by_add = *gendilep.Ltbar() + *gendilep.Btbar() + *gendilep.Nutbar();

		double top_noncon_p = (p_top-top_by_add.Vect()).Mag();
		double tbar_noncon_p = (p_top_bar-tbar_by_add.Vect()).Mag();
		double top_noncon_E = abs((gendilep.Top()-top_by_add).E());
		double tbar_noncon_E = abs((gendilep.Tbar()-tbar_by_add).E());



		truth2d["ptop_err_vs_topoffshell"]->Fill(abs(m_top_os)+abs(m_top_os_bar),abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);
		truth2d["ptop_err_vs_TWoffshell"]->Fill(abs(m_top_os)+abs(m_top_os_bar)+abs(m_w_os)+abs(m_w_os_bar),abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);;

		truth2d["ptop_err_vs_TWoffshell_alt"]->Fill(abs(m_top_os)+abs(m_top_os_bar)+abs(m_w_os)+abs(m_w_os_bar),abs(p_top_err)+abs(p_top_err_bar),weight);;
		truth2d["pnu_err_vs_TWoffshell"]->Fill(abs(m_top_os)+abs(m_top_os_bar)+abs(m_w_os)+abs(m_w_os_bar),abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag()),weight);;
		truth2d["pnu_err_vs_TWoffshell_alt"]->Fill(abs(m_top_os)+abs(m_top_os_bar)+abs(m_w_os)+abs(m_w_os_bar),abs(p_nu_err)+abs(p_nu_err_bar),weight);;

		truth1d["m_top_gen_postselection"]->Fill(m_top_os+172.5,weight);
		truth1d["m_w_gen_postselection"]->Fill(m_w_os+80.4,weight);

		truth1d["nu_test"]->Fill(p_nu_err,weight);
		truth1d["nu_testb"]->Fill(p_nu_err_bar,weight);
		truth1d["TTPT_test"]->Fill(gendilep.TT().Pt(),weight);

		truth2d["pnu_err_vs_noncon_E"]->Fill(top_noncon_E+tbar_noncon_E,abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag()),weight);
		truth2d["pnu_err_vs_noncon_p"]->Fill(top_noncon_p+tbar_noncon_p,abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag()),weight);
		truth2d["top_err_vs_noncon_E"]->Fill(top_noncon_E+tbar_noncon_E,abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);
		truth2d["top_err_vs_noncon_p"]->Fill(top_noncon_p+tbar_noncon_p,abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);
		truth2d["pnu_err_vs_noncon_EP"]->Fill(top_noncon_E+tbar_noncon_E+top_noncon_p+tbar_noncon_p,abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag()),weight);
		truth2d["top_err_vs_noncon_EP"]->Fill(top_noncon_E+tbar_noncon_E+top_noncon_p+tbar_noncon_p,abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);

		double allgensources = top_noncon_E+tbar_noncon_E+top_noncon_p+tbar_noncon_p + abs(m_top_os)+abs(m_top_os_bar)+abs(m_w_os)+abs(m_w_os_bar);

		truth2d["pnu_err_vs_gen_allsource"]->Fill(allgensources,abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag()),weight);
		truth2d["top_err_vs_gen_allsource"]->Fill(allgensources,abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()),weight);




		if ( allgensources < 5.0 && abs(p_top_err/p_top.Mag())+abs(p_top_err_bar/p_top_bar.Mag()) > 1.) {
			cout<<endl<<"( "<<abs(m_top_os)<<", "<<abs(m_top_os_bar)<<","<<abs(m_w_os)<<", "<< abs(m_w_os_bar);
			cout<<", "<<abs(p_nu_err/p_trunu.Mag())+abs(p_nu_err_bar/p_trunu_bar.Mag())<<" )";
			NuDetailPrint(idealdilep,"i");
			NuDetailPrint(gendilep,"g");

			cout<<endl<<"~~~~~"<<endl;
			idealdilep.Solve();
			NuDetailPrint(idealdilep,"test");
		}

	}


	if (rightbdilep->AddJets().size() >0) {
		return;
	}

//###################
//# MET Res         #
//###################

	// Dilepton METdilep;
	// METdilep.Copy(rightbdilep);
	TLorentzVector MET_err = *recodilep.MET() - *gendilep.MET();
	truth1d["MET_res_y"]->Fill(MET_err.Y()/abs(gendilep.MET()->Y()),weight);
	truth1d["MET_res_x"]->Fill(MET_err.X()/abs(gendilep.MET()->X()),weight);
	truth1d["MET_res2_y"]->Fill(MET_err.Y()/abs(gendilep.MET()->P()),weight);
	truth1d["MET_res2_x"]->Fill(MET_err.X()/abs(gendilep.MET()->P()),weight);
	truth1d["MET_diff_x"]->Fill(MET_err.X(),weight);
	truth1d["MET_diff_y"]->Fill(MET_err.Y(),weight);
	truth1d["MET_res_E"]->Fill(MET_err.E()/abs(gendilep.MET()->E()),weight);
	truth1d["MET_res_theta"]->Fill(recodilep.MET()->Angle(gendilep.MET()->Vect()),weight);
	truth2d["MET_res_y_vs_x"]->Fill(abs(MET_err.X()/gendilep.MET()->X()),abs(MET_err.Y()/gendilep.MET()->Y()),weight);


//###################
//# MET Smearing    #
//###################

	double gen_Mtt = gendilep.Mtt();
	double gen_sumP = gendilep.Nut()->P()+gendilep.Nutbar()->P();


	Dilepton rightbdilep_smear;
	double met_smear_std = 0;
	int metAv = 200;
	int imet = 0;
	size_t count = 0;
	size_t maxcount = 10*metAv; //in case you want to run extra events when errors are occuring
	double met_sig = 0.7;
	truth1d["met_smear_Mtt"]->Reset();

	TLorentzVector smearnut0 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnutbar0 = TLorentzVector(0,0,0,0);

	while ( count < metAv )
	{
		rightbdilep_smear.Reset();
		rightbdilep_smear.Copy(rightbdilep);
		TLorentzVector smeared_met = METsmear(*rightbdilep->MET(),25.7);
		rightbdilep_smear.SetMET(&smeared_met);
		rightbdilep_smear.Solve();
		int n_solns = rightbdilep_smear.MinNum();

		if (n_solns>0) {
			imet++;
			int lowpi = rightbdilep_smear.LowPsIndex();
			rightbdilep_smear.SetNeutrinos(lowpi);
			smearnut0 = smearnut0+ rightbdilep_smear.NutSolns()[lowpi];
			smearnutbar0 = smearnutbar0+ rightbdilep_smear.NutbarSolns()[lowpi];
			double smear_Mtt = rightbdilep_smear.Mtt();
			truth1d["met_smear_Mtt"]->Fill( smear_Mtt, weight);
		}

		count++;
	}

	if (imet>0 && rightbdilep->MinNum()>0) {
		int lowpi = rightbdilep->LowPsIndex();
		rightbdilep->SetNeutrinos(lowpi);
		met_smear_std = truth1d["met_smear_Mtt"]->GetStdDev();
		truth2d["met_smear_Mtt_err_vs_stdev"]->Fill(met_smear_std, abs(Mtt_recolo-gen_Mtt )/gen_Mtt,weight);
		truth2d["met_smear_Mtt_diff_vs_stdev"]->Fill(met_smear_std, abs(Mtt_recolo-gen_Mtt ),weight);
	}


//###################
//# B-jet Smearing  #
//###################




//parameter choices
	int NuAv = 200;
	int iNu1 = 0;
	int iNu2 = 0;
	count = 0;
	maxcount = 10*NuAv; //in case you want to run extra events when errors are occuring
	double bj_sig = 0.18;
	TLorentzVector smearnut = TLorentzVector(0,0,0,0);
	TLorentzVector smearnutbar = TLorentzVector(0,0,0,0);
	double b_smear_std = 0;

//initialize things
//smear1: each solution has weight 1
//smear2: each smeared event has weight 1, solutions averaged
//smear3: each smeared event has weight 1, low P solution

	rightbdilep_smear1 = *rightbdilep;
	rightbdilep_smear2 = *rightbdilep;
	rightbdilep_smear3 = *rightbdilep;

	TLorentzVector smearnut1 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnutbar1 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnut2 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnutbar2 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnut3 = TLorentzVector(0,0,0,0);
	TLorentzVector smearnutbar3 = TLorentzVector(0,0,0,0);
	truth1d["b_smear_Mtt"]->Reset();



	while ( count < NuAv ) {
		rightbdilep_smear.Reset();
		rightbdilep_smear.Copy(rightbdilep);
		TLorentzVector bt0 = *rightbdilep_smear.Bt();
		TLorentzVector btbar0 = *rightbdilep_smear.Btbar();
		TLorentzVector met0 = *rightbdilep_smear.MET();

//smear and find met correction
		TLorentzVector bt_smear = ESmear(*rightbdilep->Bt(),bj_sig);
		TLorentzVector btbar_smear = ESmear(*rightbdilep->Btbar(),bj_sig);
		TLorentzVector bt_smear_corr = bt_smear - bt0;
		TLorentzVector btbar_smear_corr = btbar_smear - btbar0;
		TLorentzVector met_smear = met0 - bt_smear_corr - btbar_smear_corr;
		met_smear = TLorentzVector(met_smear.X(),met_smear.Y(),0,0);

//set smeared jets, adjust MET, solve for neutrinos
		rightbdilep_smear.SetBs(&bt_smear,&btbar_smear);
		rightbdilep_smear.SetMET(&met_smear);
		rightbdilep_smear.Solve();

//implement the 3 methods
		int n_solns = rightbdilep_smear.MinNum();
		for (int i = 0; i < n_solns; i++) {
			smearnut1 = smearnut1 + rightbdilep_smear.NutSolns()[i];
			smearnutbar1 = smearnutbar1 + rightbdilep_smear.NutbarSolns()[i];
			iNu1++;
		}
		if (n_solns>0) {
			iNu2++;

			int lowpi = rightbdilep_smear.LowPsIndex();
			smearnut3 = smearnut3+ rightbdilep_smear.NutSolns()[lowpi];
			smearnutbar3 = smearnutbar3+ rightbdilep_smear.NutbarSolns()[lowpi];


			rightbdilep_smear.SetNeutrinos(lowpi);
			double smear_Mtt = rightbdilep_smear.Mtt();
			truth1d["b_smear_Mtt"]->Fill( smear_Mtt, weight);

			rightbdilep_smear.AverageNeutrinos();
			smearnut2 = smearnut2 + rightbdilep_smear.NutSolns()[0];
			smearnutbar2 = smearnutbar2 + rightbdilep_smear.NutbarSolns()[0];

//plots

			double smear_sumDP = deltaP(*rightbdilep_smear.Nut(), *gendilep.Nut()) + deltaP(*rightbdilep_smear.Nutbar(), *gendilep.Nutbar());
			double smear_sumP = rightbdilep_smear.Nut()->P() + rightbdilep_smear.Nutbar()->P();


			if (DILEPACC_count < smeargraphs) {
				// truth1d["smear_Mtt"+to_string(DILEPACC_count)]->Fill( smear_Mtt, 1.);
				// truth1d["smear_sumP"+to_string(DILEPACC_count)]->Fill( smear_sumDP, 1.);
				// truth1d["smear_sumDP"+to_string(DILEPACC_count)]->Fill(smear_sumP, 1.);
			}
		}
		count++;
	}



	if (rightbdilep_lop.MinNum()>0) {
		b_smear_std = truth1d["b_smear_Mtt"]->GetStdDev();
		truth2d["b_smear_Mtt_err_vs_stdev"]->Fill(b_smear_std, abs(Mtt_recolo-gen_Mtt )/gen_Mtt,weight);
		truth2d["b_smear_Mtt_diff_vs_stdev"]->Fill(b_smear_std, abs(Mtt_recolo-gen_Mtt ),weight);
		truth1d["Mtterr_ave_gen"]->Fill(Mtt_recolo-gen_Mtt,weight);
		truth1d["Mtterr_low_gen"]->Fill(Mtt_recolo-b_smear_std-gen_Mtt,weight);
	}

	if (DILEPACC_count < smeargraphs) {
		//cout<<endl<<"{"<<DILEPACC_count<<", "<<smear_std<<"} ";
}

	if (b_smear_std >0 && met_smear_std >0 && Mtt_recolo>0) {
		double total_std = sqrt(pow(b_smear_std,2)+pow(met_smear_std,2));
		truth2d["total_smear_Mtt_err_vs_stdev"]->Fill(total_std, abs(Mtt_recolo-gen_Mtt )/gen_Mtt,weight);
		truth2d["total_smear_Mtt_diff_vs_stdev"]->Fill(total_std, abs(Mtt_recolo-gen_Mtt ),weight);
	}


//if it worked, normalize

	if (iNu2>NuAv/4) {

		double nunorm1 = 1/(double)iNu1;
		double nunorm2 = 1/(double)iNu2;
		smearnut1 *= nunorm1;
		smearnutbar1 *= nunorm1;
		smearnut2 *= nunorm2;
		smearnutbar2 *= nunorm2;
		smearnut3 *= nunorm2;
		smearnutbar3 *= nunorm2;

		truth1d["NUSMEAR1_deltaP"]->Fill( deltaP( *gendilep.Nut(), smearnut1) ,weight);
		truth1d["NUSMEAR1_deltaP"]->Fill( deltaP( *gendilep.Nutbar(), smearnutbar1) ,weight);
		truth1d["NUSMEAR2_deltaP"]->Fill( deltaP( *gendilep.Nut(), smearnut2) ,weight);
		truth1d["NUSMEAR2_deltaP"]->Fill( deltaP( *gendilep.Nutbar(), smearnutbar2) ,weight);
		truth1d["NUSMEAR3_deltaP"]->Fill( deltaP( *gendilep.Nut(), smearnut3) ,weight);
		truth1d["NUSMEAR3_deltaP"]->Fill( deltaP( *gendilep.Nutbar(), smearnutbar3) ,weight);
		truth1d["Nu_solver_overshoot_1"]->Fill((smearnut2.P()-gendilep.Nut()->P())/abs(gendilep.Nut()->P()));

		rightbdilep_smear1.SetNeutrinos(smearnut1, smearnutbar1);
		rightbdilep_smear2.SetNeutrinos(smearnut2, smearnutbar2);
		rightbdilep_smear3.SetNeutrinos(smearnut3, smearnutbar3);
	}
	else{
		cout<<"bad event!!";
	}







//############################
//# Mtt and Delta Y Analysis #
//############################

	if (rightbdilep->MinNum()>0
		//&&  rightbdilep->AddJets().size() ==0
	) {
		DILEPACC=true;
		DILEPACC_count++;

		rightbdilep->SortSolutions(gendilep);

		if (rightbdilep->MinNum()>1) {
			int lowpi = rightbdilep->LowPIndex();
			int lowpti = rightbdilep->LowPtIndex();
			int lowpzi = rightbdilep->LowPzIndex();
			int lowpsi = rightbdilep->LowPsIndex();
			truth1d["LowPNuIndex"]->Fill(lowpi,weight);
			truth1d["LowPtNuIndex"]->Fill(lowpti,weight);
			truth1d["LowPzNuIndex"]->Fill(lowpzi,weight);
			truth1d["LowPsNuIndex"]->Fill(lowpsi,weight);
			truth1d["LowPs2NuIndex"]->Fill(rightbdilep->LowPsIndex2(),weight);
			truth1d["LowPmaxNuIndex"]->Fill(rightbdilep->LowPmaxIndex(),weight);
			truth1d["LowMttNuIndex"]->Fill(rightbdilep->LowMttIndex(),weight);
			truth1d["LowPNu_percentChange"]->Fill(deltaP(rightbdilep->NutSolns()[lowpi],rightbdilep->NutSolns()[0])/rightbdilep->NutSolns()[0].P(),weight);
		}

		int ilo = rightbdilep->LowPsIndex();
		rightbdilep->SetNeutrinos(ilo);
		double Mtt_lo = rightbdilep->Mtt();
		double DY_lo = rightbdilep->DY();
		//(*rightbdilep->Bt()+*rightbdilep->Btbar()+ *rightbdilep->Lt()+ *rightbdilep->Ltbar() + rightbdilep->NutSolns()[ilo] + rightbdilep->NutbarSolns()[ilo]).M();
		//(*rightbdilep->Bt() + *rightbdilep->Lt() + rightbdilep->NutSolns()[ilo]).Rapidity() - (*rightbdilep->Btbar() + *rightbdilep->Ltbar() + rightbdilep->NutbarSolns()[ilo]).Rapidity();

		rightbdilep->SetNeutrinos(0);
		double Mtt_best = rightbdilep->Mtt();
		double DY_best = rightbdilep->DY();
		truth1d["Nu_solver_overshoot_2"]->Fill((rightbdilep_ave.Nut()->P()-gendilep.Nut()->P())/abs(gendilep.Nut()->P()));
		truth1d["Nu_solver_overshoot_3"]->Fill((rightbdilep->Nut()->P()-gendilep.Nut()->P())/abs(gendilep.Nut()->P()));

		double Mtt_ideal = 0;
		double DY_ideal = 0;

		if (idealdilep.MinNum()>0) {
			Dilepton idealnudilep = *rightbdilep;
			idealdilep.SortSolutions(gendilep);
			idealdilep.SetNeutrinos(0);
			idealnudilep.SetNeutrinos(*idealdilep.Nut(),*idealdilep.Nutbar());
			Mtt_ideal = idealnudilep.Mtt();
			DY_ideal = idealnudilep.DY();

			truth1d["idealnu_M_tt"]->Fill( Mtt_ideal , weight);
			truth1d["idealnu_DY_tt"]->Fill( DY_ideal , weight);

			for (size_t i = 0; i < 6; i++)
			{
				truth2d["idealnu_MDY_tt_HATHOR_y"+to_string(i)]->Fill(Mtt_ideal, DYbin(abs(DY_ideal))  ,weight*HATweights[i]);
				truth1d["idealnu_M_tt_HATHOR_y"+to_string(i)]->Fill( Mtt_ideal , weight*HATweights[i]);
				truth1d["idealnu_DY_tt_HATHOR_y"+to_string(i)]->Fill( DY_ideal ,weight*HATweights[i]);
			}
		}

		double TTpt = gendilep.TT().Pt();

		double Mtt_av = rightbdilep_ave.Mtt();
		double DY_av = rightbdilep_ave.DY();
		double M_bl = (*rightbdilep->Bt()+*rightbdilep->Btbar()+ *rightbdilep->Lt()+ *rightbdilep->Ltbar()).M();
		double DY_bl = (*rightbdilep->Bt()+ *rightbdilep->Lt() ).Rapidity() - (*rightbdilep->Btbar()+ *rightbdilep->Ltbar() ).Rapidity();
		double Mtt_gen = gendilep.Mtt();
		double DY_gen = gendilep.DY();

		double BJ_err = deltaP(*gendilep.Bt(), *rightbdilep->Bt()) + deltaP(*gendilep.Btbar(),*rightbdilep->Btbar());
		double MET_err = deltaPxy(*gendilep.NuMET(), *rightbdilep->MET()) ;


		truth2d["bestnu_IMerr_vs_IMgen"]->Fill(Mtt_gen,abs(Mtt_best-Mtt_gen)/Mtt_gen,weight);
		truth2d["bestnu_IMerr_vs_IMreco"]->Fill(Mtt_best,abs(Mtt_best-Mtt_gen)/Mtt_gen,weight);
		truth2d["bestnu_IMdiff_vs_IMgen"]->Fill(Mtt_gen,abs(Mtt_best-Mtt_gen),weight);
		truth2d["bestnu_IMdiff_vs_IMreco"]->Fill(Mtt_best,abs(Mtt_best-Mtt_gen),weight);
		truth1d["Mtterr_best_gen"]->Fill(Mtt_best-Mtt_gen,weight);
		truth1d["NuPerr_best_gen"]->Fill(gendilep.Nut()->P()-rightbdilep->Nut()->P() +gendilep.Nutbar()->P()- rightbdilep->Nutbar()->P(), weight);

		truth2d["Mtterr_vs_jeterror"]->Fill(BJ_err, Mtt_best-Mtt_gen, weight);
		truth2d["Mtterr_vs_meterror"]->Fill(MET_err, Mtt_best-Mtt_gen, weight);
		truth2d["Nuerr_vs_ttpt"]->Fill(TTpt, deltaP(*rightbdilep->Nut(),*gendilep.Nut())+ deltaP(*rightbdilep->Nutbar(),*gendilep.Nutbar()), weight);
		truth1d["Meterr"]->Fill(MET_err,  weight);


		// if (rightbdilep_smear1.MinNum()>0) {
		// 	for (size_t i = 0; i < 6; i++) {
		// 		truth1d["smearnu1_M_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear1.Mtt()  ,weight*HATweights[i]);
		// 		truth1d["smearnu1_DY_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear1.DY() ,weight*HATweights[i]);
		// 		truth1d["smearnu2_M_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear2.Mtt()  ,weight*HATweights[i]);
		// 		truth1d["smearnu2_DY_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear2.DY() ,weight*HATweights[i]);
		// 		truth1d["smearnu3_M_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear3.Mtt()  ,weight*HATweights[i]);
		// 		truth1d["smearnu3_DY_tt_HATHOR_y"+to_string(i)]->Fill( rightbdilep_smear3.DY() ,weight*HATweights[i]);
		// 	}
		// }



		truth1d["reco_M_bl"]->Fill( (*rightbdilep->Bt()+*rightbdilep->Btbar()+ *rightbdilep->Lt()+ *rightbdilep->Ltbar()).M() ,weight);
		truth1d["reco_DY_bl"]->Fill( (*rightbdilep->Bt()+ *rightbdilep->Lt() ).Rapidity() - (*rightbdilep->Btbar()+ *rightbdilep->Ltbar() ).Rapidity() ,weight);
		truth1d["reco_M_tt"]->Fill( Mtt_best , weight);
		truth1d["reco_DY_tt"]->Fill( DY_best , weight);

		for (size_t i = 0; i < 6; i++)
		{
			truth2d["reco_MDY_bl_HATHOR_y"+to_string(i)]->Fill(M_bl,  DYbin(abs(DY_bl)), weight*HATweights[i]);
			truth2d["reco_MDY_tt_HATHOR_y"+to_string(i)]->Fill( Mtt_best, DYbin(abs(DY_best)) ,weight*HATweights[i]);
			truth1d["reco_M_bl_HATHOR_y"+to_string(i)]->Fill( M_bl , weight*HATweights[i]);
			truth1d["reco_DY_bl_HATHOR_y"+to_string(i)]->Fill( DY_bl , weight*HATweights[i]);
			truth1d["reco_M_tt_HATHOR_y"+to_string(i)]->Fill( Mtt_best , weight*HATweights[i]);
			truth1d["reco_DY_tt_HATHOR_y"+to_string(i)]->Fill( DY_best ,weight*HATweights[i]);
			truth1d["avenu_M_tt_HATHOR_y"+to_string(i)]->Fill( Mtt_av ,weight*HATweights[i]);
			truth1d["avenu_DY_tt_HATHOR_y"+to_string(i)]->Fill( DY_av ,weight*HATweights[i]);
			truth1d["lonu_M_tt_HATHOR_y"+to_string(i)]->Fill( Mtt_lo ,weight*HATweights[i]);
			truth1d["lonu_DY_tt_HATHOR_y"+to_string(i)]->Fill( DY_lo ,weight*HATweights[i]);
		}
	}




//------------------------------------------------------
// NuSolver specific plots
//--------------------------------------------------------
// Fill the graphs in DilepPlots.cc

	right_bjets.Fill(*rightbdilep,weight);
	right_bjets.FillNuInfo(*rightbdilep,gendilep,weight);
	gen_test.Fill(idealdilep,weight);
	gen_test.FillNuInfo(idealdilep,gendilep,weight);
	right_bjets_lop.Fill(rightbdilep_lop,weight);
	right_bjets_lop.FillNuInfo(rightbdilep_lop,gendilep,weight);

	right_bjets_ave.Fill(rightbdilep_ave,weight);
	right_bjets_ave.FillNuInfo(rightbdilep_ave,gendilep,weight);
	// right_bjets_smear1.Fill(rightbdilep_smear1,weight);
	// right_bjets_smear1.FillNuInfo(rightbdilep_smear1,gendilep,weight);
	right_bjets_smear2.Fill(rightbdilep_smear2,weight);
	right_bjets_smear2.FillNuInfo(rightbdilep_smear2,gendilep,weight);
	// right_bjets_smear3.Fill(rightbdilep_smear3,weight);
	// right_bjets_smear3.FillNuInfo(rightbdilep_smear3,gendilep,weight);
	swapped_bjets.FillNuInfo(*swappedbdilep,gendilep,weight);
	swapped_bjets.Fill(*swappedbdilep,weight);



//is the MET resolution or b/lep measurements the source of our woes??

	Dilepton genmet_dilep = *rightbdilep;
	Dilepton numet_dilep = *rightbdilep;
	genmet_dilep.SetMET(gendilep.MET());
	genmet_dilep.Solve();
	numet_dilep.SetMET(gendilep.NuMET());
	numet_dilep.Solve();
	gen_MET.Fill(genmet_dilep,weight);
	gen_MET.FillNuInfo(genmet_dilep,gendilep,weight);
	nu_MET.Fill(numet_dilep,weight);
	nu_MET.FillNuInfo(numet_dilep,gendilep,weight);

	gen_bl_dilep.SetMET(rightbdilep->MET());
	gen_bl_dilep.Solve();
	gen_bl.Fill(gen_bl_dilep,weight);
	gen_bl.FillNuInfo(gen_bl_dilep,gendilep,weight);

// you might want to plot some gen-level neutrino solver stuff here


//------------------------------------------------------
// compare the two b-jet pairings using NuSolver information
//--------------------------------------------------------
//*** this section is kinda messy, should be cleaned.

	vector<TLorentzVector> rightnuts = rightbdilep->NutSolns();
	vector<TLorentzVector> rightnutbars = rightbdilep->NutbarSolns();
	vector<TLorentzVector> wrongnuts = swappedbdilep->NutSolns();
	vector<TLorentzVector> wrongnutbars = swappedbdilep->NutbarSolns();

	for(TLorentzVector nu : rightnuts)
	{
		truth2d["deltaP_vs_Pt"]->Fill(nu.Pt(),deltaP(nu, *gendilep.Nut()));
	}
	for(TLorentzVector nu : rightnutbars)
	{
		truth2d["deltaP_vs_Pt"]->Fill(nu.Pt(),deltaP(nu, *gendilep.Nutbar()));
	}
	for(TLorentzVector nu : wrongnuts)
	{

	}
	for(TLorentzVector nu : wrongnutbars)
	{

	}


	if (rightbdilep->MinNum()>0 && swappedbdilep->MinNum()>0)
	{


		truth1d["compare_sum_Pt"]->Fill( (rightnuts[0].Pt() + rightnutbars[0].Pt()) > (wrongnuts[0].Pt() + wrongnutbars[0].Pt()),weight);
		truth1d["compare_sum_Pz"]->Fill( (rightnuts[0].Pz() + rightnutbars[0].Pz()) > (wrongnuts[0].Pz() + wrongnutbars[0].Pz()),weight);
		truth1d["compare_sum_E"]->Fill( (rightnuts[0].E() + rightnutbars[0].E()) > (wrongnuts[0].E() + wrongnutbars[0].E()),weight);

		double right_sum_deltaR_lep = rightnuts[0].DeltaR(*rightbdilep->Lt()) + rightnutbars[0].DeltaR(*rightbdilep->Ltbar());


		double wrong_sum_deltaR_lep = wrongnuts[0].DeltaR(*rightbdilep->Lt()) + wrongnutbars[0].DeltaR(*rightbdilep->Ltbar());


		truth1d["compare_sum_deltaR_lep"]->Fill( right_sum_deltaR_lep > wrong_sum_deltaR_lep , weight);
		double right_sum_IMW =  (rightnuts[0]+*recodilep.Lt()).M()  +  (rightnutbars[0]+*rightbdilep->Ltbar()).M();
		double wrong_sum_IMW =  (wrongnuts[0]+*recodilep.Lt()).M()  +  (wrongnutbars[0]+*rightbdilep->Ltbar()).M();

		truth1d["compare_sum_IMW"]->Fill( right_sum_IMW > wrong_sum_IMW, weight);
		double right_sum_deltaR_B = right_sum_deltaR_lep + rightnuts[0].DeltaR(*rightbdilep->Bt()) + rightnutbars[0].DeltaR(*rightbdilep->Btbar());
		double wrong_sum_deltaR_B = wrong_sum_deltaR_lep  + wrongnuts[0].DeltaR(*rightbdilep->Bt()) + wrongnutbars[0].DeltaR(*rightbdilep->Btbar());
		truth1d["compare_sum_deltaR_lepandb"]->Fill( right_sum_deltaR_B > wrong_sum_deltaR_B , weight);
			truth1d["compare_sum_blep_deltaR"]->Fill( rightbdilep->Bt()->DeltaR(*rightbdilep->Lt()) + rightbdilep->Btbar()->DeltaR(*rightbdilep->Ltbar())
											> swappedbdilep->Bt()->DeltaR(*swappedbdilep->Lt()) + swappedbdilep->Btbar()->DeltaR(*swappedbdilep->Ltbar()), weight);
		truth1d["compare_min_pxpy"]->Fill(rightbdilep->DNS().MinEllipseDist()  > swappedbdilep->DNS().MinEllipseDist(),weight);

		double minpxpy_correct =  rightbdilep->DNS().MinEllipseDist();
		double minpxpy_swapped =  swappedbdilep->DNS().MinEllipseDist();
		double diff = minpxpy_correct - minpxpy_swapped;

		int rightbn = rightbdilep->MinNum();
		int swappedbn = swappedbdilep->MinNum();
		if (rightbn > swappedbn)
		{
			truth1d["num_mins_guess"]->Fill(1.,weight);
		}
		if (rightbn < swappedbn)
		{
			truth1d["num_mins_guess"]->Fill(-1.,weight);
		}
		if (rightbn = swappedbn)
		{
			truth1d["num_mins_guess"]->Fill(0.,weight);
		}
		if (diff<0)
		{
				truth1d["min_deltaPxPy_success_vs_fail"]->Fill(1.,weight);
		}
		if (diff>0)
		{
				truth1d["min_deltaPxPy_success_vs_fail"]->Fill(0.,weight);
		}

	}


//-----------------------------------------------
//Make a few of those colorful  plots if you want
//-----------------------------------------------

		if( dilepton_event_counter < neutrino_graphs )
	{

		gendilep.NuPlot(nugraphs, "gen", dilepton_event_counter, 400 );
		rightbdilep->NuPlot(nugraphs, "right", dilepton_event_counter, 400 );
		swappedbdilep->NuPlot(nugraphs, "swapped", dilepton_event_counter, 400 );

	// Print a lot of details if needed

		// cout.precision(20);
		// cout<<endl<<"pb= {";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Bt())[i]<<", ";
		// }
		// cout<<"} pl=";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Lt())[i]<<", ";
		// }
		// cout<<"} pnu= {";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Nut())[i]<<", ";
		// }
		// cout<<"}"<<endl;
		// cout<<"} pnunu= {";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(gendilep.NuNu())[i]<<", ";
		// }
		// cout<<endl<<"pbb= {";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Btbar())[i]<<", ";
		// }
		// cout<<"} plb=";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Ltbar())[i]<<", ";
		// }
		// cout<<"}"<<endl;
		// cout<<"} pnu= {";
		// for (int i = 0; i < 4; ++i)
		// {
		// 	cout<<(*gendilep.Nutbar())[i]<<", ";
		// }
		// cout<<"}"<<endl;
      //
		// cout.precision(6);
		// cout<<"}"<<endl;
		// cout<<"Mmu= "<<gendilep.Lt()->M()<<" E="<<gendilep.Bt()->E()<<" theta="<<gendilep.Bt()->Angle(gendilep.Lt()->Vect())<<" beta="<<gendilep.Bt()->Beta()<<" gammamu="<<gendilep.Lt()->Gamma()<< "Mb=" <<gendilep.Bt()->M()<<endl;
      //


	}
		dilepton_event_counter++;

//############################
//# b jet error analysis     #
//############################




//############################
//# b jet pairing adjustment #
//############################

//obviously we swap the pairing if there is an error


	bool reco_err = recodilep.Error();
	bool reco_bswap_err = recodilep_bswap.Error();


	if (reco_err)
	{
		if (reco_bswap_err)
		{
			return;
		}
		else
		{
			recodilep = recodilep_bswap;
		}
	}
	else if (!recodilep_bswap.Error())
	{
		// if neither pairing has an error, we use some info about the neutrino solutions to pick the more likely one.


		// if (recodilep.DNS().MinEllipseDist()>recodilep_bswap.DNS().MinEllipseDist())
		// {
		// 	cout<<"!!!";
		// 	return;
		// }
	}

//test how well we did

	jetIDguesscorrect=recodilep.BJetsCorrect(gendilep);
	truth1d["bjet_match_PERCENT"]->Fill(jetIDguesscorrect,weight);

}

void ttbar::reconstruction(){}

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

double ttbar::HATHORweight(Dilepton& gendi, size_t yuk)
{
	double Hweight = 1.;
	double Mttbar=(gendi.Tbar() +gendi.Top()).M();
	double DelY=gendi.Top().Rapidity()-gendi.Tbar().Rapidity();

	if (yuk>HAThists.size()) {
		cout<<"err";
		return 1.;
	}

	TH2D* HAThist = HAThists[yuk];

	int binx = HAThist->GetXaxis()->FindFixBin(Mttbar);
	int biny = HAThist->GetYaxis()->FindFixBin(DelY);
	Hweight += HAThist->GetBinContent(binx,biny);
	return Hweight;
}



//This method is called once every file, contains the event loop
//run your proper analysis here
void ttbar::analyze()
{
	clock_t begint = clock();
	int nevent = 0;
	URStreamer event(tree_);
	IDElectron::streamer = &event;
	IDMuon::streamer = &event;
	PDFuncertainty::streamer = &event;
	while(event.next())
	{
		nevent++;
		testcounter++;
		if(nevent % 10000 == 0)cout << "(Event:" << nevent << " " << event.run << ") mem useage:" << mem_usage() << endl;
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
		recodilep.Reset();
		recodilep_bswap.Reset();
		gendilep.Reset();
		idealdilep.Reset();
		rightbdilep_ave.Reset();
		rightbdilep_smear1.Reset();
		rightbdilep_smear2.Reset();
		rightbdilep_smear3.Reset();

		genper = 0;

		HATweights.clear();

		sjets.clear();
		cleanedjets.clear();
		recobjets.clear();
		addjets.clear();
		smuons.clear();
		tightmuons.clear();
		loosemuons.clear();
		selectrons.clear();
		mediumelectrons.clear();
		looseelectrons.clear();

		DILEP=false;
		DILEPACC=false;
		ee=0;
		emu=0;
		mumu=0;
		lt=nullptr;
		ltbar=nullptr;
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
			if (DILEP) {

				gen1d["M_tt"]->Fill( gendilep.TT().M() ,weight);
				gen1d["DY_tt"]->Fill( gendilep.Top().Rapidity() - gendilep.Tbar().Rapidity()  ,weight);
				gen1d["M_bl"]->Fill( (*gendilep.Bt()+*gendilep.Btbar()+ *gendilep.Lt()+ *gendilep.Ltbar()).M() ,weight);
				gen1d["DY_bl"]->Fill( (*gendilep.Bt()+ *gendilep.Lt() ).Rapidity() - (*gendilep.Btbar()+ *gendilep.Ltbar() ).Rapidity() ,weight);

				double DY_tt = gendilep.Top().Rapidity() - gendilep.Tbar().Rapidity();
				double DY_bl = (*gendilep.Bt()+ *gendilep.Lt() ).Rapidity() - (*gendilep.Btbar()+ *gendilep.Ltbar() ).Rapidity();

				for (size_t i = 0; i < 6; i++) {
					HATweights.push_back(HATHORweight(gendilep,i));
					gen2d["MDY_bl_HATHOR_y"+to_string(i)]->Fill( (*gendilep.Bt()+*gendilep.Btbar()+ *gendilep.Lt()+ *gendilep.Ltbar()).M(), DYbin(abs(DY_bl)) ,weight*HATweights[i]);
					gen2d["MDY_tt_HATHOR_y"+to_string(i)]->Fill( gendilep.TT().M(), DYbin(abs(DY_tt)) ,weight*HATweights[i]);
					gen1d["M_tt_HATHOR_y"+to_string(i)]->Fill( gendilep.TT().M() ,weight*HATweights[i]);
					gen1d["DY_tt_HATHOR_y"+to_string(i)]->Fill( DY_tt ,weight*HATweights[i]);
					gen1d["M_bl_HATHOR_y"+to_string(i)]->Fill( (*gendilep.Bt()+*gendilep.Btbar()+ *gendilep.Lt()+ *gendilep.Ltbar()).M() ,weight*HATweights[i]);
					gen1d["DY_bl_HATHOR_y"+to_string(i)]->Fill( DY_bl ,weight*HATweights[i]);
				}


			}
			else continue;

			gen1d["ttpt"]->Fill((gentq+gentqbar).Pt(), weight);
			gen2d["ttpt_vs_IM"]->Fill(gendilep.TT().M(),(gentq+gentqbar).Pt(), weight);
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

		// if(PSEUDOTOP)
		// {
		// 	if(genallper.IsComplete())
		// 	{
		// 		response_ps.FillTruth("thardpt", Max(gentq.Pt(), gentqbar.Pt()), weight);
		// 		response_ps.FillTruth("tsoftpt", Min(gentq.Pt(), gentqbar.Pt()), weight);
		// 		response_ps.FillTruth("thadpt", gentqhad.Pt(), weight);
		// 		response_ps.FillTruth("thady", Abs(gentqhad.Rapidity()), weight);
		// 		response_ps.FillTruth("tleppt", gentqlep.Pt(), weight);
		// 		response_ps.FillTruth("tlepy", Abs(gentqlep.Rapidity()), weight);
		// 		response_ps.FillTruth("ttm", (gentq+gentqbar).M(), weight);
		// 		response_ps.FillTruth("ttpt", (gentq+gentqbar).Pt(), weight);
		// 		response_ps.FillTruth("tty", Abs((gentq+gentqbar).Rapidity()), weight);
		// 	}
		// 	if(psper.IsComplete())
		// 	{
		// 		response_ps.FillAll("thardpt", psper.THard().Pt(), weight);
		// 		response_ps.FillAll("tsoftpt", psper.TSoft().Pt(), weight);
		// 		response_ps.FillAll("thadpt", psper.THad().Pt(), weight);
		// 		response_ps.FillAll("thady", Abs(psper.THad().Rapidity()), weight);
		// 		response_ps.FillAll("tleppt", psper.TLep().Pt(), weight);
		// 		response_ps.FillAll("tlepy", Abs(psper.TLep().Rapidity()), weight);
		// 		response_ps.FillAll("ttm", psper.TT().M(), weight);
		// 		response_ps.FillAll("ttpt", psper.TT().Pt(), weight);
		// 		response_ps.FillAll("tty", Abs(psper.TT().Rapidity()), weight);
		// 	}
		// 	if(genallper.IsComplete() && psper.IsComplete())
		// 	{
		// 		response_ps.FillTruthReco("thardpt", Max(gentq.Pt(), gentqbar.Pt()), psper.THard().Pt(), weight);
		// 		response_ps.FillTruthReco("tsoftpt", Min(gentq.Pt(), gentqbar.Pt()), psper.TSoft().Pt(), weight);
		// 		response_ps.FillTruthReco("thadpt", gentqhad.Pt(), psper.THad().Pt(), weight);
		// 		response_ps.FillTruthReco("thady", Abs(gentqhad.Rapidity()), Abs(psper.THad().Rapidity()), weight);
		// 		response_ps.FillTruthReco("tleppt", gentqlep.Pt(), psper.TLep().Pt(), weight);
		// 		response_ps.FillTruthReco("tlepy", Abs(gentqlep.Rapidity()), Abs(psper.TLep().Rapidity()), weight);
		// 		response_ps.FillTruthReco("ttm", (gentq+gentqbar).M(), psper.TT().M(), weight);
		// 		response_ps.FillTruthReco("ttpt", (gentq+gentqbar).Pt(), psper.TT().Pt(), weight);
		// 		response_ps.FillTruthReco("tty", Abs((gentq+gentqbar).Rapidity()), Abs(psper.TT().Rapidity()), weight);
		// 	}
		// 	genper = &psper;
		// 	SEMILEP = psper.IsComplete();
		// 	SEMILEPACC = SEMILEP;
		// 	if(psper.IsComplete())
		// 	{
		// 		gent = psper.T();
		// 		gentbar = psper.Tb();
		// 		gentlep = psper.TLep();
		// 		genthad = psper.THad();
		// 	}
		// }
		// else
		// {
		// 	genper = &genallper;
		// 	gent = gentq;
		// 	gentbar = gentqbar;
		// 	gentlep = gentqlep;
		// 	genthad = gentqhad;
		// }

		if(isDA == 0) {AddGenJetSelection(event);}

		// if(SEMILEPACC)
		// {
		// 	truth1d["counter"]->Fill(1.5, weight);
		// 	ttp_genall.Fill(*genper, weight);

		// 	response.FillTruth("thardpt", Max(gent.Pt(), gentbar.Pt()), weight);
		// 	response.FillTruth("tsoftpt", Min(gent.Pt(), gentbar.Pt()), weight);
		// 	response.FillTruth("nobin", genthad.Pt(), weight);
		// 	response.FillTruth("thadpt", genthad.Pt(), weight);
		// 	response.FillTruth("thady", Abs(genthad.Rapidity()), weight);
		// 	response.FillTruth("tleppt", gentlep.Pt(), weight);
		// 	response.FillTruth("tlepy", Abs(gentlep.Rapidity()), weight);
		// 	response.FillTruth("ttm", (gent+gentbar).M(), weight);
		// 	response.FillTruth("ttpt", (gent+gentbar).Pt(), weight);
		// 	response.FillTruth("tty", Abs((gent+gentbar).Rapidity()), weight);
		// 	response.FillTruth("njet", genper->NAddJets(), weight);
		// 	response.FillTruth("dymp", gent.Rapidity() - gentbar.Rapidity(), weight);
		// 	response.FillTruth("dy", Abs(gent.Rapidity())- Abs(gentbar.Rapidity()), weight);
		// 	response.FillTruth("ht", genper->Ht(), weight);
		// 	response.FillTruth("evtmass", genper->EvtMass(), weight);
		// 	response2d.FillTruth("njets_thadpt", genthad.Pt(), genper->NAddJets(), weight);
		// 	response2d.FillTruth("njets_ttpt", (gent+gentbar).Pt(), genper->NAddJets(), weight);
		// 	response2d.FillTruth("njets_ttm", (gent+gentbar).M(), genper->NAddJets(), weight);
		// 	response2d.FillTruth("thady_thadpt", genthad.Pt(), Abs(genthad.Rapidity()), weight);
		// 	response2d.FillTruth("ttpt_thadpt", genthad.Pt(), (gent+gentbar).Pt(), weight);
		// 	response2d.FillTruth("thadpt_ttm", (gent+gentbar).M(), genthad.Pt(), weight);
		// 	response2d.FillTruth("ttpt_ttm", (gent+gentbar).M(), (gent+gentbar).Pt(), weight);
		// 	response2d.FillTruth("ttm_tty", Abs((gent+gentbar).Rapidity()), (gent+gentbar).M(), weight);
		// 	response2d.FillTruth("ttm_dy", Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), (gent+gentbar).M(), weight);

		// 	response2dvar.FillTruth("njet+thadpt", genper->NAddJets(), genthad.Pt(), weight);
		// 	response2dvar.FillTruth("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt(), weight);
		// 	response2dvar.FillTruth("njet+ttm", genper->NAddJets(), (gent+gentbar).M(), weight);
		// 	response2dvar.FillTruth("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt(), weight);
		// 	response2dvar.FillTruth("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt(), weight);
		// 	response2dvar.FillTruth("thadpt+ttm", genthad.Pt(), (gent+gentbar).M(), weight);
		// 	response2dvar.FillTruth("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M(), weight);
		// 	response2dvar.FillTruth("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()), weight);
		// 	response2dvar.FillTruth("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), weight);
		// 	for(size_t n = 0 ; n < genper->NJets() ; ++n)
		// 	{
		// 		response2dvar.FillTruth("jet+jetpt", n, genper->GetJet(n)->Pt(), weight);
		// 		response2dvar.FillTruth("jet+jetptOF", n, genper->GetJet(n)->Pt(), weight);
		// 		response2dvar.FillTruth("jet+jeteta", n, abs(genper->GetJet(n)->Eta()), weight);
		// 		response2dvar.FillTruth("jet+jetdr", n, genper->DRminTTjets(genper->GetJet(n)), weight);
		// 		response2dvar.FillTruth("jet+jetdrtop", n, genper->DRminTop(genper->GetJet(n)), weight);
		// 	}
		// 	if(PDFTEST)
		// 	{
		// 		pdfunc->Fill1d("pdfunc_thardpt", Max(gent.Pt(), gentbar.Pt()), weight);
		// 		pdfunc->Fill1d("pdfunc_tsoftpt", Min(gent.Pt(), gentbar.Pt()), weight);
		// 		pdfunc->Fill1d("pdfunc_nobin", genthad.Pt(), weight);
		// 		pdfunc->Fill1d("pdfunc_thadpt", genthad.Pt(), weight);
		// 		pdfunc->Fill1d("pdfunc_tleppt", gentlep.Pt(), weight);
		// 		pdfunc->Fill1d("pdfunc_thady", Abs(genthad.Rapidity()), weight);
		// 		pdfunc->Fill1d("pdfunc_tlepy", Abs(gentlep.Rapidity()), weight);
		// 		pdfunc->Fill1d("pdfunc_ttm", (gent+gentbar).M(), weight);
		// 		pdfunc->Fill1d("pdfunc_tty", Abs((gent+gentbar).Rapidity()), weight);
		// 		pdfunc->Fill1d("pdfunc_ttpt", (gent+gentbar).Pt(), weight);
		// 		pdfunc->Fill1d("pdfunc_njet", genper->NAddJets(), weight);
		// 		pdfunc->Fill1d("pdfunc_dymp", gent.Rapidity() - gentbar.Rapidity(), weight);
		// 		pdfunc->Fill1d("pdfunc_dy", Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()), weight);
		// 		pdfunc->Fill1d("pdfunc_ht", genper->Ht(), weight);
		// 		pdfunc->Fill1d("pdfunc_evtmass", genper->EvtMass(), weight);

		// 		pdfunc->Fill1d("pdfunc_njet+thadpt", response2dvar.GetBin("njet+thadpt", genper->NAddJets(), genthad.Pt())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_njet+ttpt", response2dvar.GetBin("njet+ttpt", genper->NAddJets(), (gent+gentbar).Pt())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_njet+ttm", response2dvar.GetBin("njet+ttm", genper->NAddJets(), (gent+gentbar).M())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_thady+thadpt", response2dvar.GetBin("thady+thadpt", Abs(genthad.Rapidity()), genthad.Pt())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_ttpt+thadpt", response2dvar.GetBin("ttpt+thadpt", (gent+gentbar).Pt(), genthad.Pt())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_thadpt+ttm", response2dvar.GetBin("thadpt+ttm", genthad.Pt(), (gent+gentbar).M())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_ttpt+ttm", response2dvar.GetBin("ttpt+ttm", (gent+gentbar).Pt(), (gent+gentbar).M())-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_ttm+tty", response2dvar.GetBin("ttm+tty", (gent+gentbar).M(), Abs((gent+gentbar).Rapidity()))-0.5, weight);
		// 		pdfunc->Fill1d("pdfunc_ttm+dy", response2dvar.GetBin("ttm+dy", (gent+gentbar).M(), Abs(gent.Rapidity()) - Abs(gentbar.Rapidity()))-0.5, weight);
		// 		for(size_t n = 0 ; n < genper->NJets() ; ++n)
		// 		{
		// 			pdfunc->Fill1d("pdfunc_jet+jetpt", response2dvar.GetBin("jet+jetpt", n, genper->GetJet(n)->Pt())-0.5, weight);
		// 			pdfunc->Fill1d("pdfunc_jet+jetptOF", response2dvar.GetBin("jet+jetptOF", n, genper->GetJet(n)->Pt())-0.5, weight);
		// 			pdfunc->Fill1d("pdfunc_jet+jeteta", response2dvar.GetBin("jet+jeteta", n, abs(genper->GetJet(n)->Eta()))-0.5, weight);
		// 			pdfunc->Fill1d("pdfunc_jet+jetdr", response2dvar.GetBin("jet+jetdr", n, genper->DRminTTjets(genper->GetJet(n)))-0.5, weight);
		// 			pdfunc->Fill1d("pdfunc_jet+jetdrtop", response2dvar.GetBin("jet+jetdrtop", n, genper->DRminTop(genper->GetJet(n)))-0.5, weight);
		// 		}
		// 	}
		// }
		// if(SEMILEPACC)
		// {
		// 	ttp_genacc.Fill(*genper, weight);
		// 	truth1d["counter"]->Fill(2.5, weight);
		// }

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
	clock_t endt = clock();
	double timespent = (endt-begint)/CLOCKS_PER_SEC;
	cout<<"analysis took about "<< timespent <<"seconds"<<endl;

	double bj2 = truth1d["bjet_match_PERCENT"]->GetBinContent(3);
	double bj1 = truth1d["bjet_match_PERCENT"]->GetBinContent(2);
	double bj0 = truth1d["bjet_match_PERCENT"]->GetBinContent(1);
	double bjetmatchfrac = bj2/(bj0+bj2);
	cout<<(bj0+bj2);

	cout<<endl<<"DeltaRlep:::"<<truth1d["compare_sum_deltaR_lep"]->GetBinContent(2)/(truth1d["compare_sum_deltaR_lep"]->GetBinContent(1)+truth1d["compare_sum_deltaR_lep"]->GetBinContent(2));
	cout<<endl<<"PT:::"<<truth1d["compare_sum_Pt"]->GetBinContent(2)/(truth1d["compare_sum_Pt"]->GetBinContent(1)+truth1d["compare_sum_Pt"]->GetBinContent(2));
	cout<<endl<<"E:::"<<truth1d["compare_sum_E"]->GetBinContent(2)/(truth1d["compare_sum_E"]->GetBinContent(1)+truth1d["compare_sum_E"]->GetBinContent(2));
	cout<<endl<<"pxpy:::"<<truth1d["compare_min_pxpy"]->GetBinContent(2)/(truth1d["compare_min_pxpy"]->GetBinContent(1)+truth1d["compare_min_pxpy"]->GetBinContent(2));
	cout<<endl<<"dR:::"<<truth1d["compare_sum_blep_deltaR"]->GetBinContent(2)/(truth1d["compare_sum_blep_deltaR"]->GetBinContent(1)+truth1d["compare_sum_blep_deltaR"]->GetBinContent(2));
	cout<<endl<<"bjet efficiency:  "<<bjetmatchfrac<<endl;
	cout<<endl<<"nevents = "<<testcounter;
	cout<<endl<<"DILEP_count = "<<DILEP_count<<endl;
	cout<<endl<<"DILEPACC_count = "<<DILEPACC_count<<endl;
	double eff = ((double)DILEPACC_count) / ((double)DILEP_count);
	cout<<endl<<"efficiency = "<<eff<<endl;




}

Long64_t ttbar::mem_usage()
{
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    Long64_t rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
        >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
        >> utime >> stime >> cutime >> cstime >> priority >> nice
        >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();
    return(rss*sysconf(_SC_PAGE_SIZE));
}

//make it executable
int main(int argc, char *argv[])
{
	URParser &parser = URParser::instance(argc, argv);
	URDriver<ttbar> test;
	return test.run();
}
