//--------------------------------------------------------------------
// Take 2D DCA v. pT histograms and produce DCA distributions for
// each pT bin
//
// 1.0 - 1.5
// 1.5 - 2.0
// 2.0 - 2.5
// 2.5 - 3.0
// 3.0 - 4.0
// 4.0 - 5.0
// 5.0 - 6.0
// 6.0 - 8.0
// 8.0 - 10.0
//--------------------------------------------------------------------

#include <iostream>

using namespace std;

//-----------------------
// Variables
//-----------------------

//Number of pT bins
const int NPTBINS = 9;

//Number of background sources (pi, eta, gamma, jpsi, ke3, hadrons)
const int NBACKG = 6;

//Edges for pT bins
float binLow[NPTBINS] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0};
float binHigh[NPTBINS] = {1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

//DCA v pT for data electrons from input file
TH2D *h_elec_data_dca_pT;

//DCA v pT for simulated electrons from input file
TH2D *h_elec_back_dca_pT[NBACKG];

//DCA from data electrons in every pT bin
TH1D *h_dca_elec_data[NPTBINS];

//DCA in every pT bin for every background species
TH1D *h_dca_elec_back[NBACKG][NPTBINS];

//Normalized DCA in every pT bin for every background species
TH1D *h_dca_elec_back_norm[NBACKG][NPTBINS];

//Normalization for every background species
TH1D *h_normalization[NBACKG];

//Sum of normalized DCA templates for background
TH1D *h_dca_elec_back_norm_sum[NPTBINS];

//Colors for each background source
int colors[NBACKG] = {kBlue, kRed, kGreen + 3, kBlack, kOrange - 3, kViolet};

//-----------------------
// Functions
//-----------------------

void readDataFiles()
{
	TFile *fin = new TFile("electron_out.root");
	h_elec_data_dca_pT = (TH2D*) fin->Get("h_elec_dca_pT_notrejected_1x");
}


void readBackgroundInformation()
{
	//Read in normalizations for electron background
	TFile *fNormBackg = new TFile("jav_normalization_1x.root");
	h_normalization[0] = (TH1D*) fNormBackg->Get("h_frac_piz");
	h_normalization[1] = (TH1D*) fNormBackg->Get("h_frac_eta");
	h_normalization[2] = (TH1D*) fNormBackg->Get("h_frac_gam");
	h_normalization[3] = (TH1D*) fNormBackg->Get("h_frac_jpsi");
	h_normalization[4] = (TH1D*) fNormBackg->Get("h_frac_ke3");

	//Read in normalizations for hadron contamination
	TFile *fNormHadrons = new TFile("tim_normalization_1x.root");
	h_normalization[5] = (TH1D*) fNormHadrons->Get("h_hadcontam_n");

	//Read in the DCA v pT for each electron background source
	TFile *fDCABackg = new TFile("elecbkgdca.root");
	h_elec_back_dca_pT[0] = (TH2D*) fDCABackg->Get("h_pT_sdca_piz");
	h_elec_back_dca_pT[1] = (TH2D*) fDCABackg->Get("h_pT_sdca_eta");
	h_elec_back_dca_pT[2] = (TH2D*) fDCABackg->Get("h_pT_sdca_photon");
	h_elec_back_dca_pT[3] = (TH2D*) fDCABackg->Get("h_pT_sdca_jpsi");
	h_elec_back_dca_pT[4] = (TH2D*) fDCABackg->Get("h_pT_sdca_ke3");

	//Read in the DCA for every pT bin for hadron contamination
	TFile *fDCAHadron = new TFile("hadroncontdca.root");
	for (int i = 0; i < NPTBINS; i++)
	{
		h_dca_elec_back[5][i] = (TH1D*) fDCAHadron->Get(Form("h_pioncon_%i", i + 1));
	}
}


void projectIndividualDCA()
{
	//Project the 2D DCA v pT histograms to get DCA for data electrons
	for (int i = 0; i < NPTBINS; i++)
	{
		int pTbinLo = h_elec_data_dca_pT->GetYaxis()->FindBin(binLow[i]);
		int pTbinHi = h_elec_data_dca_pT->GetYaxis()->FindBin(binHigh[i]);

		h_dca_elec_data[i] = (TH1D*) h_elec_data_dca_pT->ProjectionX(Form("h_dca_elec_data_%i", i), pTbinLo, pTbinHi);
	}

	//Project the 2D DCA v pT histograms to get DCA for electron background sources
	for (int i = 0; i < NBACKG - 1; i++)
	{
		for (int j = 0; j < NPTBINS; j++)
		{
			int pTbinLo = h_elec_back_dca_pT[i]->GetXaxis()->FindBin(binLow[j]);
			int pTbinHi = h_elec_back_dca_pT[i]->GetXaxis()->FindBin(binHigh[j]);

			h_dca_elec_back[i][j] = (TH1D*) h_elec_back_dca_pT[i]->ProjectionY(Form("h_dca_elec_back_%i_%i", i, j), pTbinLo, pTbinHi);
		}
	}
}


void normalizeDCAs()
{
	for (int i = 0; i < NPTBINS; i++)
	{
		//Get the integral of the data electron DCA for pT bin at hand
		float dataArea = h_dca_elec_data[i]->Integral();

		cout << "pT BIN " << i << endl;

		for (int j = 0; j < NBACKG; j++)
		{
			float desiredArea = dataArea * h_normalization[j]->GetBinContent(i + 1);
			h_dca_elec_back_norm[j][i] = (TH1D*) h_dca_elec_back[j][i]->Clone(Form("h_dca_elec_back_norm_%i_%i", j, i));
			h_dca_elec_back_norm[j][i]->Scale(desiredArea / h_dca_elec_back_norm[j][i]->Integral());

			cout << "   ---> Source " << j << ": " << h_dca_elec_back_norm[j][i]->Integral() / dataArea << endl;
		}

		cout << endl;
	}
}


void combineDCAs()
{
	for (int i = 0; i < NPTBINS; i++)
	{
		h_dca_elec_back_norm_sum[i] = (TH1D*) h_dca_elec_back_norm[0][i]->Clone(Form("h_dca_elec_back_norm_sum_%i", i));

		for (int j = 1; j < NBACKG; j++)
		{
			h_dca_elec_back_norm_sum[i]->Add(h_dca_elec_back_norm[j][i]);
		}
	}
}


void plot()
{
	gStyle->SetOptStat(0);

	TCanvas *c[NPTBINS];

	for (int i = 0; i < NPTBINS; i++)
	{
		c[i] = new TCanvas(Form("c%i", i), Form("c%i", i), 600, 600);
		c[i]->SetLogy();

		TH1D *hTemplate = new TH1D(Form("hTemplate%i", i), ";DCA_{T} [cm]", 100, -0.2, 0.2);
		hTemplate->GetXaxis()->SetTitleFont(62);
		hTemplate->GetXaxis()->SetLabelFont(62);
		hTemplate->GetYaxis()->SetTitleFont(62);
		hTemplate->GetYaxis()->SetLabelFont(62);
		hTemplate->GetYaxis()->SetRangeUser(1E-2, 1E3);
		hTemplate->Draw();

		h_dca_elec_data[i]->Draw("same");
		h_dca_elec_data[i]->SetLineColor(kBlack);

		for (int j = 0; j < NBACKG; j++)
		{
			h_dca_elec_back_norm[j][i]->SetLineColor(colors[j]);
			h_dca_elec_back_norm[j][i]->Draw("hist,same");
		}
	}
}


void plotCombined()
{
	gStyle->SetOptStat(0);
	TCanvas *c[NPTBINS];

	for (int i = 0; i < NPTBINS; i++)
	{
		c[i] = new TCanvas(Form("cCombined%i", i), Form("cCombined%i", i), 600, 600);
		c[i]->SetLogy();

		TH1D *hTemplate = new TH1D(Form("hTemplate%i", i), ";DCA_{T} [cm]", 100, -0.2, 0.2);
		hTemplate->GetXaxis()->SetTitleFont(62);
		hTemplate->GetXaxis()->SetLabelFont(62);
		hTemplate->GetYaxis()->SetTitleFont(62);
		hTemplate->GetYaxis()->SetLabelFont(62);
		hTemplate->GetYaxis()->SetRangeUser(1E-2, 5E4);
		hTemplate->Draw();

		h_dca_elec_back_norm_sum[i]->Draw("hist,same");

		h_dca_elec_data[i]->Draw("same");
		h_dca_elec_data[i]->SetLineColor(kBlack);
	}
}


void writeResult()
{
	//The names of electron DCA should be run15ppdca<number>
	//The names of background DCA should be run15ppbkg<number>
	TFile *fout = new TFile("dca_for_unfold.root", "RECREATE");

	//01-31-18: Removing 1.0-1.5 bins from unfold, so that the 1.5-2.0 bin should have index 1
	for (int i = 1; i < NPTBINS-1; i++)
	{
		h_dca_elec_back_norm_sum[i]->SetName(Form("run15ppbkg%i", i));
		h_dca_elec_data[i]->SetName(Form("run15ppdca%i", i));

		h_dca_elec_back_norm_sum[i]->Write();
		h_dca_elec_data[i]->Write();
	}

	fout->Close();
}


void processDataElectrons()
{
	//Read in all input files
	readDataFiles();
	readBackgroundInformation();

	//Project 2D DCA v pT histos into individual DCA for each pT bin, and every background electron source
	projectIndividualDCA();

	//Normalize DCAs
	normalizeDCAs();

	//Add up DCA from all background sources, for every pT bin
	combineDCAs();

	//Plot, with a canvas for every pT bin showing the data + background sources
	//plot();
	//plotCombined();

	//Write histograms to provide to unfolding code
	writeResult();
}