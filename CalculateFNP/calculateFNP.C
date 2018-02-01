#include <iostream>

using namespace std;

//-----------------------
// Attributes
//-----------------------

//Number of photonic primaries (pi, eta, direct gamma)
const int NPRIM = 3;

//Primaries
TH1F *h_prim_pT[NPRIM];

//Electrons accepted and rejected by the veto cut
TH1F *h_elec_pT_accepted[NPRIM];
TH1F *h_elec_pT_rejected[NPRIM];

//Combined total photonic electrons from all three species
TH1F *h_elec_pT_photonic_accepted;
TH1F *h_elec_pT_photonic_rejected;

//Photonic survival rate
TH1F *h_photonic_surv;

//Electrons from data
TH1F *h_elec_pT_data_accepted;
TH1F *h_elec_pT_data_rejected;
TH1F *h_elec_pT_data;

//Hadrons from data
TH1F *h_hadron_pT_data;
TH1F *h_hadron_pT_data_accepted;

//Hadronic survival rate
TH1F *h_hadronic_surv;

//FNP
TH1F *h_P;
TH1F *h_NP;
TH1F *h_FNP;
TGraphAsymmErrors *g_FNP;

//FNP with conversion veto taken into account
TH1F *h_P_conv;
TH1F *h_NP_conv;
TH1F *h_FNP_conv;
TGraphAsymmErrors *g_FNP_conv;

//Additional NP cocktail components
TH1F *h_elec_pT_accepted_ke3;
TH1F *h_elec_pT_accepted_hf;
TH1F *h_elec_pT_accepted_jpsi;

//Relative fractions
TH1F *h_total_photonic;
TH1F *h_total_nonphotonic;

TH1F *h_frac_jpsi;
TH1F *h_frac_ke3;

TH1F *h_frac_piz;
TH1F *h_frac_eta;
TH1F *h_frac_gam;

//Use updated electrons from Tim's Ntuple
bool useUpdatedElectrons = true;

//Window for the conversion veto cut
//Default 1x
//Twice 2x
string window = "2x";

//-----------------------
// Functions
//-----------------------

void readSimFiles()
{
	TFile *finPiz = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/Sims012318/pizeros.root");
	TFile *finEta = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/Sims012318/etas.root");
	TFile *finGam = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/Sims012318/photons.root");

	h_prim_pT[0]          = (TH1F*) finPiz->Get("h_prim_pT_rw");
	h_elec_pT_accepted[0] = (TH1F*) finPiz->Get("h_elec_pT_rw");
	h_elec_pT_rejected[0] = (TH1F*) finPiz->Get("h_elec_pT_rw_vetorejected");

	h_prim_pT[1]          = (TH1F*) finEta->Get("h_prim_pT_rw");
	h_elec_pT_accepted[1] = (TH1F*) finEta->Get("h_elec_pT_rw");
	h_elec_pT_rejected[1] = (TH1F*) finEta->Get("h_elec_pT_rw_vetorejected");

	h_prim_pT[2]          = (TH1F*) finGam->Get("h_prim_pT_rw");
	h_elec_pT_accepted[2] = (TH1F*) finGam->Get("h_elec_pT_rw");
	h_elec_pT_rejected[2] = (TH1F*) finGam->Get("h_elec_pT_rw_vetorejected");
}


void readDataFiles()
{
	TFile *fData;

	if (!useUpdatedElectrons)
	{
		fData = new TFile("11568/calcsurvivalrate.root");

		h_elec_pT_data_accepted = (TH1F*) fData->Get("h_elec_pT_notrejected_pTcut");
		h_elec_pT_data_rejected = (TH1F*) fData->Get("h_elec_pT_rejected_pTcut");
		h_elec_pT_data          = (TH1F*) fData->Get("h_elec_pT");

		h_hadron_pT_data_accepted = (TH1F*) fData->Get("h_hadron_pT_notrejected_pTcut");
		h_hadron_pT_data = (TH1F*) fData->Get("h_hadron_pT");
	}
	else
	{
		fData = new TFile("012918/updated_electron_hadron.root");

		h_elec_pT_data_accepted = (TH1F*) fData->Get(Form("h_elec_pT_notrejected_%s", window.c_str()));
		h_elec_pT_data_rejected = (TH1F*) fData->Get(Form("h_elec_pT_rejected_%s", window.c_str()));
		h_elec_pT_data          = (TH1F*) fData->Get("h_elec_pT");

		h_hadron_pT_data_accepted = (TH1F*) fData->Get("hadnotconversions");
		
		TH1F *h_hadron_pT_data_rejected = (TH1F*) fData->Get("hadconversions");
		
		h_hadron_pT_data = (TH1F*) h_hadron_pT_data_rejected->Clone("h_hadron_pT_data");	
		h_hadron_pT_data->Add(h_hadron_pT_data_accepted);
	}
}


void readAdditionalHFComponents()
{
	//My own Jpsis
	float hadron_density = 7.3E-4;
	float branching_ratio =  0.0602;
	float sigma_pp = 23.0;

	TFile *finJpsi = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/Sims012318/jpsis.root");
	TH1F *h_prim_pT_jpsi = (TH1F*) finJpsi->Get("h_prim_pT_rw");
	h_elec_pT_accepted_jpsi = (TH1F*) finJpsi->Get("h_elec_pT_rw");
	h_elec_pT_accepted_jpsi->Scale(hadron_density);
	h_elec_pT_accepted_jpsi->Scale(branching_ratio);
	h_elec_pT_accepted_jpsi->Scale(1.0 / sigma_pp);
	h_elec_pT_accepted_jpsi->Scale(1.0 / h_prim_pT_jpsi->Integral(1, h_prim_pT_jpsi->GetNbinsX()));

	const int NBINS = 9;
	double bins[NBINS + 1] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
	h_elec_pT_accepted_jpsi = (TH1F*) h_elec_pT_accepted_jpsi->Rebin(NBINS, "h_elec_pT_accepted_jpsi_rb", bins);

	//Ke3s from Tim, already with the right binning
	TFile *fKe3 = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/tim_cocktail.root");
	h_elec_pT_accepted_ke3 = (TH1F*) fKe3->Get("h_ke3_pt");
	h_elec_pT_accepted_ke3->Scale(1.0 / 23.0);

	//Combined HF yield from PPG077
	//Data from tables in published version of PPG077
	float pT[28] = {0.35,
	                0.45,
	                0.55,
	                0.65,
	                0.75,
	                0.85,
	                0.95,
	                1.10,
	                1.30,
	                1.50,
	                1.70,
	                1.90,
	                2.10,
	                2.30,
	                2.50,
	                2.70,
	                2.90,
	                3.10,
	                3.30,
	                3.50,
	                3.70,
	                3.90,
	                4.25,
	                4.75,
	                5.50,
	                6.50,
	                7.50,
	                8.50
	               };

	float yield[28] = {1.36E-2,
	                   6.05E-3,
	                   3.36E-3,
	                   1.81E-3,
	                   1.30E-3,
	                   7.50E-4,
	                   5.61E-4,
	                   3.18E-4,
	                   1.26E-4,
	                   6.58E-5,
	                   3.83E-5,
	                   1.89E-5,
	                   1.04E-5,
	                   6.08E-6,
	                   3.42E-6,
	                   2.06E-6,
	                   1.40E-6,
	                   8.64E-7,
	                   5.91E-7,
	                   4.11E-7,
	                   2.83E-7,
	                   1.91E-7,
	                   1.08E-7,
	                   4.84E-8,
	                   1.59E-8,
	                   5.30E-9,
	                   1.27E-9,
	                   8.19E-10
	                  };

	float errx[28] = {0.0};

	float erry[28] = {4.45E-3,
	                  1.70E-3,
	                  7.56E-4,
	                  3.98E-4,
	                  2.28E-4,
	                  1.35E-4,
	                  9.14E-5,
	                  3.77E-5,
	                  2.10E-5,
	                  1.25E-5,
	                  2.07E-6,
	                  1.18E-6,
	                  7.44E-7,
	                  4.96E-7,
	                  3.63E-7,
	                  6.49E-8,
	                  4.77E-8,
	                  3.53E-8,
	                  2.73E-8,
	                  2.15E-8,
	                  1.70E-8,
	                  1.36E-8,
	                  5.95E-9,
	                  3.73E-9,
	                  1.79E-9,
	                  9.26E-10,
	                  5.73E-10,
	                  5.16E-10
	                 };

	TGraphErrors *g_pT_hf_ppg077 = new TGraphErrors(28, pT, yield, errx, erry);

	//Fit the HF electron yield from PPG077 with a modified Hagedorn function
	TF1 *f_hf_fit = new TF1("f_hf_fit", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0.0, 10.0);
	f_hf_fit->SetParameter(0, 4.41211e-02);
	f_hf_fit->SetParameter(1, 0.96);
	f_hf_fit->SetParameter(2, -9.56027e-02);
	f_hf_fit->SetParameter(3, 7.10173e-01);
	f_hf_fit->SetParameter(4, 7.11577e+00);

	g_pT_hf_ppg077->Fit(f_hf_fit, "Q0R");

	//Use TGraph to fill histogram with the same binning as for J/psi and Ke3
	h_elec_pT_accepted_hf = (TH1F*) h_elec_pT_accepted_ke3->Clone("h_elec_pT_accepted_hf");
	h_elec_pT_accepted_hf->Reset();

	for (int i = 1; i <= h_elec_pT_accepted_hf->GetNbinsX(); i++)
	{
		float pTcenter = h_elec_pT_accepted_hf->GetBinCenter(i);
		float yield = f_hf_fit->Eval(pTcenter);

		h_elec_pT_accepted_hf->SetBinContent(i, yield);
	}
}


void rebinHistograms()
{
	const int NBINS = 9;
	double bins[NBINS + 1] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

	//Rebin simulated electrons
	for (int i = 0; i < NPRIM; i++)
	{
		h_elec_pT_accepted[i] = (TH1F*) h_elec_pT_accepted[i]->Rebin(NBINS, Form("h_elec_pT_accepted_rb_%i", i), bins);
		h_elec_pT_rejected[i] = (TH1F*) h_elec_pT_rejected[i]->Rebin(NBINS, Form("h_elec_pT_rejected_rb_%i", i), bins);
	}

	//Rebin data electrons and hadrons
	h_elec_pT_data = (TH1F*) h_elec_pT_data->Rebin(NBINS, "h_elec_pT_data_rb", bins);
	h_elec_pT_data_accepted = (TH1F*) h_elec_pT_data_accepted->Rebin(NBINS, "h_elec_pT_data_accepted_rb", bins);
	h_elec_pT_data_rejected = (TH1F*) h_elec_pT_data_rejected->Rebin(NBINS, "h_elec_pT_data_rejected_rb", bins);

	h_hadron_pT_data_accepted = (TH1F*) h_hadron_pT_data_accepted->Rebin(NBINS, "h_hadron_pT_data_accepted_rb", bins);
	h_hadron_pT_data = (TH1F*) h_hadron_pT_data->Rebin(NBINS, "h_hadron_pT_data_rb", bins);
}


void normalizeSpectra()
{
	float hadron_density[NPRIM]  = {43.5, 5.1, 0.10117};
	float branching_ratio[NPRIM] = {1.0, 0.4437, 1.0};
	float sigma_pp = 23.0;

	for (int i = 0; i < NPRIM; i++)
	{
		float numPrim = h_prim_pT[i]->Integral(0, h_prim_pT[i]->GetNbinsX());

		h_elec_pT_accepted[i]->Scale(1.0 / sigma_pp);
		h_elec_pT_accepted[i]->Scale(1.0 / numPrim);
		h_elec_pT_accepted[i]->Scale(hadron_density[i]);
		h_elec_pT_accepted[i]->Scale(branching_ratio[i]);

		h_elec_pT_rejected[i]->Scale(1.0 / sigma_pp);
		h_elec_pT_rejected[i]->Scale(1.0 / numPrim);
		h_elec_pT_rejected[i]->Scale(hadron_density[i]);
		h_elec_pT_rejected[i]->Scale(branching_ratio[i]);

		for (int j = 1; j <= h_elec_pT_rejected[i]->GetNbinsX(); j++)
		{
			h_elec_pT_rejected[i]->SetBinError(j, 0);
		}
	}
}


void calculatePhotonicSurvival()
{
	h_elec_pT_photonic_rejected = (TH1F*) h_elec_pT_rejected[0]->Clone("h_elec_pT_photonic_rejected");
	h_elec_pT_photonic_accepted = (TH1F*) h_elec_pT_accepted[0]->Clone("h_elec_pT_photonic_accepted");

	for (int i = 1; i < NPRIM; i++)
	{
		h_elec_pT_photonic_accepted->Add(h_elec_pT_accepted[i]);
		h_elec_pT_photonic_rejected->Add(h_elec_pT_rejected[i]);
	}

	TH1F *h_photonic_total = (TH1F*) h_elec_pT_photonic_accepted->Clone("h_photonic_total");
	h_photonic_total->Add(h_elec_pT_photonic_rejected);

	h_photonic_surv = (TH1F*) h_elec_pT_photonic_accepted->Clone("h_photonic_surv");
	h_photonic_surv->Divide(h_photonic_total);
}


void calculateHadronicSurvival()
{

	h_hadronic_surv = (TH1F*) h_hadron_pT_data_accepted->Clone("h_hadronic_surv");
	h_hadronic_surv->Divide(h_hadron_pT_data);
}


void computeFNP()
{
	//Histogram filled with all ones
	TH1F *h_one = (TH1F*) h_elec_pT_data->Clone("h_one");
	for (int i = 1; i <= h_one->GetNbinsX(); i++)
	{
		h_one->SetBinContent(i, 1.0);
		h_one->SetBinError(i, 0.0);
	}

	//Calculate P and NP without conv veto
	TH1F *h_aux1 = (TH1F*) h_elec_pT_data->Clone("h_aux1");
	h_aux1->Multiply(h_hadronic_surv);

	TH1F *h_aux2 = (TH1F*) h_elec_pT_data_accepted->Clone("h_aux2");
	h_aux2->Add(h_aux1, -1.0);

	TH1F *h_aux3 = (TH1F*) h_photonic_surv->Clone("h_aux3");
	h_aux3->Add(h_one, -1.0);

	TH1F *h_aux4 = (TH1F*) h_hadronic_surv->Clone("h_aux4");
	h_aux4->Multiply(h_aux3);

	h_P = (TH1F*) h_aux2->Clone("h_P");
	h_P->Divide(h_aux4);

	h_NP = (TH1F*) h_elec_pT_data->Clone("h_NP");
	h_NP->Add(h_P, -1.0);

	//Calculate FNP without conv veto
	TH1F *h_aux5 = (TH1F*) h_NP->Clone("h_aux5");
	h_aux5->Add(h_P);

	h_FNP = (TH1F*) h_NP->Clone("h_FNP");
	h_FNP->Divide(h_aux5);

	g_FNP = new TGraphAsymmErrors();
	g_FNP->Divide(h_NP, h_aux5, "cl=0.683 b(1,1) mode");

	//Calculate FNP with conversion veto taken into account
	//Like so: FNP = e_r*NP / (e_r*NP + e_r*e_i*P)
	h_P_conv = (TH1F*) h_P->Clone("h_P_conv");
	h_NP_conv = (TH1F*) h_NP->Clone("h_NP_conv");

	h_P_conv->Multiply(h_photonic_surv);
	h_P_conv->Multiply(h_hadronic_surv);
	h_NP_conv->Multiply(h_hadronic_surv);

	TH1F *h_aux6 = (TH1F*) h_P_conv->Clone("h_aux6");
	h_aux6->Add(h_NP_conv);

	h_FNP_conv = (TH1F*) h_NP_conv->Clone("h_FNP_conv");
	h_FNP_conv->Divide(h_aux6);

	g_FNP_conv = new TGraphAsymmErrors();
	g_FNP_conv->Divide(h_NP_conv, h_aux6, "cl=0.683 b(1,1) mode");

	g_FNP_conv->Draw();
}


void readTimsFNP()
{
	TFile *fFNP = new TFile("~/Documents/CUBoulder/JNLab/Spring18/ElectronCocktail/tim_cocktail.root");
	h_FNP_conv = (TH1F*) fFNP->Get("h_oldfnp_conveto");
}


void calculateFractions()
{
	//Add up total photonic contributions
	h_total_photonic = (TH1F*) h_elec_pT_photonic_accepted->Clone("h_total_photonic");

	//One minus FNP
	TH1F *h_FNP_conv_minus_one = (TH1F*) h_FNP_conv->Clone("h_FNP_conv_minus_one");
	for (int i = 1; i <= h_FNP_conv_minus_one->GetNbinsX(); i++)
	{
		float cont = h_FNP_conv_minus_one->GetBinContent(i);
		h_FNP_conv_minus_one->SetBinContent(i, 1.0 - cont);
	}

	//Photonic
	h_frac_piz = (TH1F*) h_FNP_conv_minus_one->Clone("h_frac_piz");
	h_frac_piz->Multiply(h_elec_pT_accepted[0]);
	h_frac_piz->Divide(h_total_photonic);

	h_frac_eta = (TH1F*) h_FNP_conv_minus_one->Clone("h_frac_eta");
	h_frac_eta->Multiply(h_elec_pT_accepted[1]);
	h_frac_eta->Divide(h_total_photonic);

	h_frac_gam = (TH1F*) h_FNP_conv_minus_one->Clone("h_frac_gam");
	h_frac_gam->Multiply(h_elec_pT_accepted[2]);
	h_frac_gam->Divide(h_total_photonic);

	//Non-photonic
	h_frac_jpsi = (TH1F*) h_frac_piz->Clone("h_frac_jpsi");
	h_frac_jpsi->Multiply(h_elec_pT_accepted_jpsi);
	h_frac_jpsi->Divide(h_elec_pT_accepted[0]);

	h_frac_ke3 = (TH1F*) h_frac_piz->Clone("h_frac_ke3");
	h_frac_ke3->Multiply(h_elec_pT_accepted_ke3);
	h_frac_ke3->Divide(h_elec_pT_accepted[0]);
}


void plotFractions()
{
	gStyle->SetOptStat(0);

	TCanvas *cFractions = new TCanvas("cFractions", "cFractions", 700, 700);
	TH1F *hTemplate = new TH1F("hTemplate", "hTemplate", 100, 1, 10);
	hTemplate->GetYaxis()->SetRangeUser(1E-6, 1);

	cFractions->SetLogy();
	cFractions->SetTicky();
	cFractions->SetTickx();

	hTemplate->SetTitle("");
	hTemplate->GetXaxis()->SetTitleFont(62);
	hTemplate->GetXaxis()->SetLabelFont(62);
	hTemplate->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	hTemplate->GetYaxis()->SetTitleFont(62);
	hTemplate->GetYaxis()->SetLabelFont(62);
	hTemplate->GetYaxis()->SetTitle("Normalization Factor");

	hTemplate->Draw();
	h_frac_piz->Draw("P,same");
	h_frac_eta->Draw("P,same");
	h_frac_gam->Draw("P,same");
	h_frac_jpsi->Draw("P,same");
	h_frac_ke3->Draw("P,same");

	for (int i = 1; i <= h_frac_piz->GetNbinsX(); i++)
	{
		h_frac_piz->SetBinError(i, 0.0);
		h_frac_eta->SetBinError(i, 0.0);
		h_frac_gam->SetBinError(i, 0.0);
		h_frac_jpsi->SetBinError(i, 0.0);
		h_frac_ke3->SetBinError(i, 0.0);
	}

	h_frac_piz->SetLineColor(kBlue);
	h_frac_piz->SetMarkerColor(kBlue);
	h_frac_piz->SetMarkerStyle(20);

	h_frac_eta->SetLineColor(kRed);
	h_frac_eta->SetMarkerColor(kRed);
	h_frac_eta->SetMarkerStyle(20);

	h_frac_gam->SetLineColor(kGreen + 3);
	h_frac_gam->SetMarkerColor(kGreen + 3);
	h_frac_gam->SetMarkerStyle(20);

	h_frac_jpsi->SetLineColor(kBlack);
	h_frac_jpsi->SetMarkerColor(kBlack);
	h_frac_jpsi->SetMarkerStyle(20);

	h_frac_ke3->SetLineColor(kOrange - 3);
	h_frac_ke3->SetMarkerColor(kOrange - 3);
	h_frac_ke3->SetMarkerStyle(20);

	TLegend *leg = new TLegend(0.15, 0.15, 0.35, 0.4);
	leg->SetLineColor(kWhite);
	leg->AddEntry(h_frac_piz, "#pi^{0} Electrons", "P");
	leg->AddEntry(h_frac_eta, "#eta Electrons", "P");
	leg->AddEntry(h_frac_gam, "Direct #gamma Electrons", "P");
	leg->AddEntry(h_frac_jpsi, "J/#psi Electrons", "P");
	leg->AddEntry(h_frac_ke3, "Ke3 Electrons", "P");
	leg->Draw("same");
}

void writeNormalizations()
{
	TFile *fout = new TFile(Form("jav_normalization_%s.root", window.c_str()), "RECREATE");
	h_frac_piz->Write();
	h_frac_eta->Write();
	h_frac_gam->Write();
	h_frac_jpsi->Write();
	h_frac_ke3->Write();
	fout->Close();
}


void calculateFNP()
{
	readSimFiles();
	readDataFiles();
	rebinHistograms();
	normalizeSpectra();
	calculatePhotonicSurvival();
	calculateHadronicSurvival();
	computeFNP();
	//readTimsFNP();
	readAdditionalHFComponents();
	calculateFractions();
	plotFractions();
	writeNormalizations();
}