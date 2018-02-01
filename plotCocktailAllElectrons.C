#include <iostream>

using namespace std;

//------------------------------------
// Variables
//------------------------------------

const int NPRIM = 4;

TH1F *h_elec_pT_nm[NPRIM];
TH1F *h_elec_pT_rw[NPRIM];
TH1F *h_prim_pT_rw[NPRIM];
TH1F *h_elec_pT_dalitz_rw[NPRIM];
TH1F *h_elec_pT_dalitz_nm[NPRIM];

TH1F *h_jpsi_dalitzeta_ratio;
TH1F *h_jpsi_dalitzeta_ratio_077;
TH1F *h_jpsi_dalitzpiz_ratio;
TH1F *h_jpsi_dalitzpiz_ratio_077;
TH1F *h_ratio_eta_dalitz;
TH1F *h_ratio_piz_dalitz;

TH1F *h_surv_piz;
TH1F *h_surv_dalitz_piz;
TH1F *h_surv_eta;
TH1F *h_surv_eta_piz;
TH1F *h_surv_jpsi;
TH1F *h_surv_gam;

TH1F *h_total;
TH1F *h_ratio_to_total[NPRIM];

float hadron_density[NPRIM]  = {43.5, 5.1, 0.10117, 7.3E-4};
float branching_ratio[NPRIM] = {1.0, 0.4437, 1.0, 0.0602};
float sigma_pp = 23.0;

float colors[NPRIM] = {kBlue, kRed, kGreen + 3, kBlack};

//Include all electrons, regardless of convveto status?
bool ingnoreVeto = true;

//------------------------------------
// Functions
//------------------------------------

void getDalitzJpsiRatio077()
{
	TFile *fin = new TFile("cocktail_ppg077.root");
	TH1F *h_jpsi = (TH1F*) fin->Get("pteJPsi");
	TH1F *h_eta  = (TH1F*) fin->Get("pteEta");
	TH1F *h_piz  = (TH1F*) fin->Get("ptePion");

	double bins[21] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
	h_jpsi = (TH1F*) h_jpsi->Rebin(20, "h_jpsi_rb", bins);
	h_eta = (TH1F*) h_eta->Rebin(20, "h_eta_rb", bins);
	h_piz = (TH1F*) h_piz->Rebin(20, "h_piz_rb", bins);

	h_jpsi_dalitzeta_ratio_077 = (TH1F*) h_jpsi->Clone("h_jpsi_dalitzeta_ratio_077");
	h_jpsi_dalitzeta_ratio_077->Divide(h_eta);

	h_jpsi_dalitzpiz_ratio_077 = (TH1F*) h_jpsi->Clone("h_jpsi_dalitzpiz_ratio_077");
	h_jpsi_dalitzpiz_ratio_077->Divide(h_piz);

	h_jpsi_dalitzeta_ratio_077->SetMarkerColor(kBlue);
	h_jpsi_dalitzeta_ratio_077->SetMarkerSize(0.8);
	h_jpsi_dalitzeta_ratio_077->SetMarkerStyle(20);
	h_jpsi_dalitzeta_ratio_077->SetLineColor(kBlue);

	h_jpsi_dalitzpiz_ratio_077->SetMarkerColor(kBlue);
	h_jpsi_dalitzpiz_ratio_077->SetMarkerSize(0.8);
	h_jpsi_dalitzpiz_ratio_077->SetMarkerStyle(20);
	h_jpsi_dalitzpiz_ratio_077->SetLineColor(kBlue);
}

void plotCocktailAllElectrons()
{
	//Read files
	TFile *fin0 = new TFile("Sims012518/pizeros.root");
	TFile *fin1 = new TFile("Sims012518/etas.root");
	TFile *fin2 = new TFile("Sims012518/photons.root");
	TFile *fin3 = new TFile("Sims012518/jpsis.root");

	h_elec_pT_rw[0] = (TH1F*) fin0->Get("h_elec_pT_rw");
	h_prim_pT_rw[0] = (TH1F*) fin0->Get("h_prim_pT_rw");
	h_elec_pT_dalitz_rw[0] = (TH1F*) fin0->Get("h_elec_pT_dalitz_rw");

	h_elec_pT_rw[1] = (TH1F*) fin1->Get("h_elec_pT_rw");
	h_prim_pT_rw[1] = (TH1F*) fin1->Get("h_prim_pT_rw");
	h_elec_pT_dalitz_rw[1] = (TH1F*) fin1->Get("h_elec_pT_dalitz_rw");

	h_elec_pT_rw[2] = (TH1F*) fin2->Get("h_elec_pT_rw");
	h_prim_pT_rw[2] = (TH1F*) fin2->Get("h_prim_pT_rw");
	h_elec_pT_dalitz_rw[2] = (TH1F*) fin2->Get("h_elec_pT_dalitz_rw");

	h_elec_pT_rw[3] = (TH1F*) fin3->Get("h_elec_pT_rw");
	h_prim_pT_rw[3] = (TH1F*) fin3->Get("h_prim_pT_rw");
	h_elec_pT_dalitz_rw[3] = (TH1F*) fin3->Get("h_elec_pT_dalitz_rw");

	//Add electrons rejected by the conversion veto as well
	TH1F *h_elec_pT_rw_rej_piz  = (TH1F*) fin0->Get("h_elec_pT_rw_vetorejected");
	TH1F *h_elec_pT_rw_rej_eta  = (TH1F*) fin1->Get("h_elec_pT_rw_vetorejected");
	TH1F *h_elec_pT_rw_rej_gam  = (TH1F*) fin2->Get("h_elec_pT_rw_vetorejected");
	TH1F *h_elec_pT_rw_rej_jpsi = (TH1F*) fin3->Get("h_elec_pT_rw_vetorejected");

	TH1F *h_elec_pT_dalitz_rw_rej_piz  = (TH1F*) fin0->Get("h_elec_pT_dalitz_rw_vetorejected");
	TH1F *h_elec_pT_dalitz_rw_rej_eta  = (TH1F*) fin1->Get("h_elec_pT_dalitz_rw_vetorejected");
	TH1F *h_elec_pT_dalitz_rw_rej_gam  = (TH1F*) fin2->Get("h_elec_pT_dalitz_rw_vetorejected");
	TH1F *h_elec_pT_dalitz_rw_rej_jpsi = (TH1F*) fin3->Get("h_elec_pT_dalitz_rw_vetorejected");

	for (int i = 1; i <= h_elec_pT_rw_rej_piz->GetNbinsX(); i++)
		{
			h_elec_pT_rw_rej_piz->SetBinError(i, 0);
			h_elec_pT_rw_rej_eta->SetBinError(i, 0);
			h_elec_pT_rw_rej_gam->SetBinError(i, 0);
			h_elec_pT_rw_rej_jpsi->SetBinError(i, 0);

			h_elec_pT_dalitz_rw_rej_piz->SetBinError(i, 0);
			h_elec_pT_dalitz_rw_rej_eta->SetBinError(i, 0);
			h_elec_pT_dalitz_rw_rej_gam->SetBinError(i, 0);
			h_elec_pT_dalitz_rw_rej_jpsi->SetBinError(i, 0);
		}

	if (!ingnoreVeto)
	{
		h_elec_pT_rw[0]->Add(h_elec_pT_rw_rej_piz);
		h_elec_pT_rw[1]->Add(h_elec_pT_rw_rej_eta);
		h_elec_pT_rw[2]->Add(h_elec_pT_rw_rej_gam);
		h_elec_pT_rw[3]->Add(h_elec_pT_rw_rej_jpsi);

		h_elec_pT_dalitz_rw[0]->Add(h_elec_pT_dalitz_rw_rej_piz);
		h_elec_pT_dalitz_rw[1]->Add(h_elec_pT_dalitz_rw_rej_eta);
		h_elec_pT_dalitz_rw[2]->Add(h_elec_pT_dalitz_rw_rej_gam);
		h_elec_pT_dalitz_rw[3]->Add(h_elec_pT_dalitz_rw_rej_jpsi);
	}

	//Calculate veto survival rate
	h_surv_piz = (TH1F*) h_elec_pT_rw[0]->Clone("h_surv_piz");
	TH1F *h_aux_piz = (TH1F*) h_elec_pT_rw[0]->Clone("h_aux_piz");
	h_aux_piz->Add(h_elec_pT_rw_rej_piz);
	h_surv_piz->Divide(h_aux_piz);

	h_surv_dalitz_piz = (TH1F*) h_elec_pT_dalitz_rw[0]->Clone("h_surv_dalitz_piz");
	TH1F *h_aux_dalitz_piz = (TH1F*) h_elec_pT_dalitz_rw[0]->Clone("h_aux_dalitz_piz");
	h_aux_dalitz_piz->Add(h_elec_pT_dalitz_rw_rej_piz);
	h_surv_dalitz_piz->Divide(h_aux_dalitz_piz);

	h_surv_eta = (TH1F*) h_elec_pT_rw[1]->Clone("h_surv_eta");
	TH1F *h_aux_eta = (TH1F*) h_elec_pT_rw[1]->Clone("h_aux_eta");
	h_aux_eta->Add(h_elec_pT_rw_rej_eta);
	h_surv_eta->Divide(h_aux_eta);

	h_surv_dalitz_eta = (TH1F*) h_elec_pT_dalitz_rw[1]->Clone("h_surv_dalitz_eta");
	TH1F *h_aux_dalitz_eta = (TH1F*) h_elec_pT_dalitz_rw[1]->Clone("h_aux_dalitz_eta");
	h_aux_dalitz_eta->Add(h_elec_pT_dalitz_rw_rej_eta);
	h_surv_dalitz_eta->Divide(h_aux_dalitz_eta);

	h_surv_gam = (TH1F*) h_elec_pT_rw[2]->Clone("h_surv_gam");
	TH1F *h_aux_gam = (TH1F*) h_elec_pT_rw[2]->Clone("h_aux_gam");
	h_aux_gam->Add(h_elec_pT_rw_rej_gam);
	h_surv_gam->Divide(h_aux_gam);

	h_surv_jpsi = (TH1F*) h_elec_pT_rw[3]->Clone("h_surv_jpsi");
	TH1F *h_aux_jpsi = (TH1F*) h_elec_pT_rw[3]->Clone("h_aux_jpsi");
	h_aux_jpsi->Add(h_elec_pT_rw_rej_jpsi);
	h_surv_jpsi->Divide(h_aux_jpsi);

	for (int ispecies = 0; ispecies < NPRIM; ispecies++)
	{
		//Determine number of primary particles
		float numPrim = h_prim_pT_rw[ispecies]->Integral(1, h_prim_pT_rw[ispecies]->GetNbinsX());

		//Normalize reweighted electrons
		h_elec_pT_nm[ispecies] = (TH1F*) h_elec_pT_rw[ispecies]->Clone(Form("h_elec_pT_nm_%i", ispecies));
		h_elec_pT_dalitz_nm[ispecies] = (TH1F*) h_elec_pT_dalitz_rw[ispecies]->Clone(Form("h_elec_pT_dalitz_nm_%i", ispecies));

		//Phase space factor
		for (int i = 1; i <= h_elec_pT_nm[ispecies]->GetNbinsX(); i++)
		{
			float pT = h_elec_pT_nm[ispecies]->GetBinCenter(i);
			float cont = h_elec_pT_nm[ispecies]->GetBinContent(i) / (2 * TMath::Pi() * pT);
			//h_elec_pT_nm[ispecies]->SetBinContent(i, cont);

			cont = h_elec_pT_dalitz_nm[ispecies]->GetBinContent(i) / (2 * TMath::Pi() * pT);
			//h_elec_pT_dalitz_nm[ispecies]->SetBinContent(i, cont);
		}

		h_elec_pT_nm[ispecies]->Scale(1.0 / sigma_pp);
		h_elec_pT_nm[ispecies]->Scale(1.0 / numPrim);
		h_elec_pT_nm[ispecies]->Scale(branching_ratio[ispecies]);
		h_elec_pT_nm[ispecies]->Scale(hadron_density[ispecies]);

		h_elec_pT_nm[ispecies]->SetMarkerColor(colors[ispecies]);
		h_elec_pT_nm[ispecies]->SetMarkerSize(0.8);
		h_elec_pT_nm[ispecies]->SetMarkerStyle(20);
		h_elec_pT_nm[ispecies]->SetLineColor(colors[ispecies]);

		h_elec_pT_dalitz_nm[ispecies]->Scale(1.0 / sigma_pp);
		h_elec_pT_dalitz_nm[ispecies]->Scale(1.0 / numPrim);
		h_elec_pT_dalitz_nm[ispecies]->Scale(branching_ratio[ispecies]);
		h_elec_pT_dalitz_nm[ispecies]->Scale(hadron_density[ispecies]);

		h_elec_pT_dalitz_nm[ispecies]->SetMarkerColor(colors[ispecies]);
		h_elec_pT_dalitz_nm[ispecies]->SetMarkerSize(0.8);
		h_elec_pT_dalitz_nm[ispecies]->SetMarkerStyle(20);
		h_elec_pT_dalitz_nm[ispecies]->SetLineColor(colors[ispecies]);
	}

	//Add up the full cocktail
	h_total = (TH1F*) h_elec_pT_nm[0]->Clone("h_total");
	for (int i = 1; i < NPRIM; i++)
	{
		h_total->Add(h_elec_pT_nm[i]);
	}

	for (int i = 0; i < NPRIM; i++)
	{
		h_ratio_to_total[i] = (TH1F*) h_elec_pT_nm[i]->Clone(Form("h_ratio_to_total_%i", i));
		h_ratio_to_total[i]->Divide(h_total);
	}

	//Create J/psi to Dalitz eta ratio
	getDalitzJpsiRatio077();

	h_jpsi_dalitzeta_ratio = (TH1F*) h_elec_pT_nm[3]->Clone("h_jpsi_dalitzeta_ratio");
	h_jpsi_dalitzeta_ratio->Divide(h_elec_pT_dalitz_nm[1]);
	h_jpsi_dalitzeta_ratio->SetMarkerColor(kBlack);
	h_jpsi_dalitzeta_ratio->SetMarkerSize(0.8);
	h_jpsi_dalitzeta_ratio->SetMarkerStyle(20);
	h_jpsi_dalitzeta_ratio->SetLineColor(kBlack);

	h_ratio_eta_dalitz = (TH1F*) h_jpsi_dalitzeta_ratio->Clone("h_ratio_eta_dalitz");
	h_ratio_eta_dalitz->Divide(h_jpsi_dalitzeta_ratio_077);

	h_jpsi_dalitzpiz_ratio = (TH1F*) h_elec_pT_nm[3]->Clone("h_jpsi_dalitzpiz_ratio");
	h_jpsi_dalitzpiz_ratio->Divide(h_elec_pT_dalitz_nm[0]);
	h_jpsi_dalitzpiz_ratio->SetMarkerColor(kBlack);
	h_jpsi_dalitzpiz_ratio->SetMarkerSize(0.8);
	h_jpsi_dalitzpiz_ratio->SetMarkerStyle(20);
	h_jpsi_dalitzpiz_ratio->SetLineColor(kBlack);

	h_ratio_piz_dalitz = (TH1F*) h_jpsi_dalitzpiz_ratio->Clone("h_ratio_piz_dalitz");
	h_ratio_piz_dalitz->Divide(h_jpsi_dalitzpiz_ratio_077);

	h_ratio_piz_dalitz->GetXaxis()->SetRangeUser(1, 10);
	h_ratio_piz_dalitz->GetYaxis()->SetRangeUser(1, 10);

	h_jpsi_dalitzeta_ratio->GetXaxis()->SetRangeUser(1, 10);
	h_jpsi_dalitzeta_ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_jpsi_dalitzeta_ratio->GetYaxis()->SetTitle("(J/#psi) / #eta Dalitz Electron Ratio");

	h_jpsi_dalitzpiz_ratio->GetXaxis()->SetRangeUser(1, 10);
	h_jpsi_dalitzpiz_ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_jpsi_dalitzpiz_ratio->GetYaxis()->SetTitle("(J/#psi) / #pi^{0} Dalitz Electron Ratio");

	///////////////////////////////////////////////////////
	gStyle->SetErrorX(0);
	TCanvas *cJpsiEtaRatio = new TCanvas("cJpsiEtaRatio", "Jpsi to Dalitz Eta", 700, 900);
	TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1);
	pad11->SetLogy();
	pad11->SetTickx();
	pad11->SetTicky();
	pad11->SetBottomMargin(0);
	pad11->Draw();
	pad11->cd();

	TH1F *hTemplate5 = new TH1F("hTemplate5", "hTemplate5", 100, 0, 20);
	hTemplate5->SetTitle("");
	hTemplate5->GetXaxis()->SetTitleFont(62);
	hTemplate5->GetXaxis()->SetLabelFont(62);
	hTemplate5->GetXaxis()->SetRangeUser(1, 10);
	hTemplate5->GetYaxis()->SetTitleFont(62);
	hTemplate5->GetYaxis()->SetLabelFont(62);
	hTemplate5->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate5->GetYaxis()->SetTitle("(J/#psi) / Dalitz #eta Electrons");
	hTemplate5->GetYaxis()->SetTitleOffset(1.3);
	hTemplate5->GetYaxis()->SetRangeUser(0.5E-2, 40);
	hTemplate5->Draw();

	h_jpsi_dalitzeta_ratio->Draw("same");
	h_jpsi_dalitzeta_ratio_077->Draw("same");

	TLatex latex2;
	latex2.SetNDC();
	latex2.SetTextSize(0.04);
	latex2.DrawLatex(.15, .83, "J/#psi Electron to #eta Dalitz Electron Ratio");
	latex2.SetTextSize(0.03);
	if (ingnoreVeto)
	{
		latex2.DrawLatex(.15, .78, "- Conversion Veto Applied");
	}
	else
	{
		latex2.DrawLatex(.15, .78, "- No Conversion Veto Applied");
	}


	TLegend *legendJpsiEtaRatio = new TLegend(0.45, 0.3, 0.88, 0.5);
	legendJpsiEtaRatio->AddEntry(h_jpsi_dalitzeta_ratio, "Run15 p+p Ana", "PL");
	legendJpsiEtaRatio->AddEntry(h_jpsi_dalitzeta_ratio_077, "PPG077", "PL");
	legendJpsiEtaRatio->SetFillStyle(0.0);
	legendJpsiEtaRatio->SetLineColor(kWhite);
	legendJpsiEtaRatio->Draw("same");

	cJpsiEtaRatio->cd();
	TPad *pad4 = new TPad("pad4", "pad4", 0, 0.0, 1, 0.3);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->Draw();
	pad4->cd();
	pad4->SetTickx();
	pad4->SetTicky();

	TH1F *hTemplate4 = new TH1F("hTemplate4", "hTemplate4", 100, 0, 20);
	hTemplate4->SetTitle("");
	hTemplate4->GetYaxis()->CenterTitle();
	hTemplate4->GetXaxis()->SetTitleFont(62);
	hTemplate4->GetXaxis()->SetLabelFont(62);
	hTemplate4->GetXaxis()->SetRangeUser(1, 10);
	hTemplate4->GetYaxis()->SetTitleFont(62);
	hTemplate4->GetYaxis()->SetLabelFont(62);
	hTemplate4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate4->GetYaxis()->SetTitle("Run15 / PPG077");
	hTemplate4->GetYaxis()->SetTitleOffset(0.6);
	hTemplate4->GetYaxis()->SetRangeUser(0.0, 5.0);
	hTemplate4->GetYaxis()->SetLabelSize(0.08);
	hTemplate4->GetXaxis()->SetTitleSize(0.08);
	hTemplate4->GetYaxis()->SetTitleSize(0.08);
	hTemplate4->GetXaxis()->SetTitleOffset(1.2);
	hTemplate4->GetXaxis()->SetLabelSize(0.08);
	hTemplate4->Draw();

	h_ratio_eta_dalitz->GetXaxis()->SetRangeUser(1, 10);
	h_ratio_eta_dalitz->Draw("same");

	TLine *tl = new TLine(1.0, 1.0, 10.0, 1.0);
	tl->SetLineStyle(7);
	tl->Draw("same");

	cJpsiEtaRatio->cd();

	////////////////////////////////////////////

	gStyle->SetErrorX(0);
	TCanvas *cJpsiPizRatio = new TCanvas("cJpsiPizRatio", "Jpsi to Dalitz Eta", 700, 900);
	TPad *pad12 = new TPad("pad12", "pad12", 0, 0.3, 1, 1);
	pad12->SetLogy();
	pad12->SetTickx();
	pad12->SetTicky();
	pad12->SetBottomMargin(0);
	pad12->Draw();
	pad12->cd();

	TH1F *hTemplate6 = new TH1F("hTemplate6", "hTemplate6", 100, 0, 20);
	hTemplate6->SetTitle("");
	hTemplate6->GetXaxis()->SetTitleFont(62);
	hTemplate6->GetXaxis()->SetLabelFont(62);
	hTemplate6->GetXaxis()->SetRangeUser(1, 10);
	hTemplate6->GetYaxis()->SetTitleFont(62);
	hTemplate6->GetYaxis()->SetLabelFont(62);
	hTemplate6->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate6->GetYaxis()->SetTitle("(J/#psi) / Dalitz #pi^{0} Electrons");
	hTemplate6->GetYaxis()->SetTitleOffset(1.3);
	hTemplate6->GetYaxis()->SetRangeUser(0.5E-2, 40);
	hTemplate6->Draw();

	TLatex latex1;
	latex1.SetNDC();
	latex1.SetTextSize(0.04);
	latex1.DrawLatex(.15, .83, "J/#psi Electron to #pi^{0} Dalitz Electron Ratio");
	latex1.SetTextSize(0.03);
	if (ingnoreVeto)
	{
		latex1.DrawLatex(.15, .78, "- Conversion Veto Applied");
	}
	else
	{
		latex1.DrawLatex(.15, .78, "- No Conversion Veto Applied");
	}

	h_jpsi_dalitzpiz_ratio->Draw("same");
	h_jpsi_dalitzpiz_ratio_077->Draw("same");

	TLegend *legendJpsiPizeroRatio = new TLegend(0.45, 0.3, 0.88, 0.5);
	legendJpsiPizeroRatio->AddEntry(h_jpsi_dalitzpiz_ratio, "Run15 p+p Ana", "PL");
	legendJpsiPizeroRatio->AddEntry(h_jpsi_dalitzpiz_ratio_077, "PPG077", "PL");
	legendJpsiPizeroRatio->SetFillStyle(0.0);
	legendJpsiPizeroRatio->SetLineColor(kWhite);
	legendJpsiPizeroRatio->Draw("same");

	cJpsiPizRatio->cd();
	TPad *pad6 = new TPad("pad6", "pad6", 0, 0.0, 1, 0.3);
	pad6->SetTopMargin(0);
	pad6->SetBottomMargin(0.3);
	pad6->Draw();
	pad6->cd();
	pad6->SetTickx();
	pad6->SetTicky();

	TH1F *hTemplate7 = new TH1F("hTemplate7", "hTemplate7", 100, 0, 20);
	hTemplate7->SetTitle("");
	hTemplate7->GetYaxis()->CenterTitle();
	hTemplate7->GetXaxis()->SetTitleFont(62);
	hTemplate7->GetXaxis()->SetLabelFont(62);
	hTemplate7->GetXaxis()->SetRangeUser(1, 10);
	hTemplate7->GetYaxis()->SetTitleFont(62);
	hTemplate7->GetYaxis()->SetLabelFont(62);
	hTemplate7->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate7->GetYaxis()->SetTitle("Run15 / PPG077");
	hTemplate7->GetYaxis()->SetTitleOffset(0.6);
	hTemplate7->GetYaxis()->SetRangeUser(0.0, 5.0);
	hTemplate7->GetYaxis()->SetLabelSize(0.08);
	hTemplate7->GetXaxis()->SetTitleSize(0.08);
	hTemplate7->GetYaxis()->SetTitleSize(0.08);
	hTemplate7->GetXaxis()->SetTitleOffset(1.2);
	hTemplate7->GetXaxis()->SetLabelSize(0.08);
	hTemplate7->Draw();

	h_ratio_piz_dalitz->GetXaxis()->SetRangeUser(1, 10);
	h_ratio_piz_dalitz->Draw("same");
	tl->Draw("same");

	cJpsiPizRatio->cd();

	///////////////////////////////////////////////////////

	//Plot cocktail
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0);
	TCanvas *cCocktail = new TCanvas("cCocktail", "PPG077 Cocktail", 700, 900);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetLogy();
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();

	TH1F *hTemplate2 = new TH1F("hTemplate2", "hTemplate2", 100, 0, 20);
	hTemplate2->SetTitle("");
	hTemplate2->GetXaxis()->SetTitleFont(62);
	hTemplate2->GetXaxis()->SetLabelFont(62);
	hTemplate2->GetXaxis()->SetRangeUser(0, 10);
	hTemplate2->GetYaxis()->SetTitleFont(62);
	hTemplate2->GetYaxis()->SetLabelFont(62);
	hTemplate2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate2->GetYaxis()->SetTitle("d^{2}N / dp^{T}dy");
	hTemplate2->GetYaxis()->SetTitleOffset(1.3);
	hTemplate2->GetYaxis()->SetRangeUser(0.5e-13, 1E-1);
	hTemplate2->Draw();

	for (int i = 0; i < NPRIM; i++)
	{
		h_elec_pT_nm[i]->Draw("same");
	}

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.DrawLatex(.15, .83, "Reconstructed Electron Yield from Simulations");
	latex.SetTextSize(0.03);
	latex.DrawLatex(.15, .78, "- Kinematic and Quality Cuts Applied");
	if (ingnoreVeto)
	{
		latex.DrawLatex(.15, .73, "- Conversion Veto Applied");
	}
	else
	{
		latex.DrawLatex(.15, .73, "- No Conversion Veto Applied");
	}

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(h_elec_pT_nm[0], "#pi^{0} Electrons", "PL");
	legend->AddEntry(h_elec_pT_nm[1], "#eta Electrons", "PL");
	legend->AddEntry(h_elec_pT_nm[2], "Direct #gamma Electrons", "PL");
	legend->AddEntry(h_elec_pT_nm[3], "J/#psi Electrons", "PL");
	legend->SetFillStyle(0.0);
	legend->SetLineColor(kWhite);
	legend->Draw("same");

	cCocktail->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();

	//gPad->SetLogy();

	TH1F *hTemplate3 = new TH1F("hTemplate3", "hTemplate3", 100, 0, 20);
	hTemplate3->SetTitle("");
	hTemplate3->GetYaxis()->CenterTitle();
	hTemplate3->GetXaxis()->SetTitleFont(62);
	hTemplate3->GetXaxis()->SetLabelFont(62);
	hTemplate3->GetXaxis()->SetRangeUser(0, 10);
	hTemplate3->GetYaxis()->SetTitleFont(62);
	hTemplate3->GetYaxis()->SetLabelFont(62);
	hTemplate3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate3->GetYaxis()->SetTitle("Ratio to Total");
	hTemplate3->GetYaxis()->SetTitleOffset(0.6);
	hTemplate3->GetYaxis()->SetRangeUser(1E-3, 1.5);
	hTemplate3->GetYaxis()->SetLabelSize(0.08);
	hTemplate3->GetXaxis()->SetTitleSize(0.08);
	hTemplate3->GetYaxis()->SetTitleSize(0.08);
	hTemplate3->GetXaxis()->SetTitleOffset(1.2);
	hTemplate3->GetXaxis()->SetLabelSize(0.08);
	hTemplate3->Draw();

	for (int i = 0; i < NPRIM; i++)
	{
		h_ratio_to_total[i]->Draw("same");
	}

	cCocktail->cd();

	//Plot survival rates
	if(!ingnoreVeto) return;
	TCanvas *cSurvPiz = new TCanvas("cSurvPiz", "cSurvPiz", 600, 600);
	h_surv_piz->GetYaxis()->SetTitle("Veto Survival Rate");
	h_surv_piz->GetYaxis()->SetTitleOffset(1.2);
	h_surv_piz->GetYaxis()->SetTitleFont(62);
	h_surv_piz->GetYaxis()->SetLabelFont(62);
	h_surv_piz->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
	h_surv_piz->GetXaxis()->SetTitleOffset(1.2);
	h_surv_piz->GetXaxis()->SetTitleFont(62);
	h_surv_piz->GetXaxis()->SetLabelFont(62);
	h_surv_piz->GetYaxis()->SetRangeUser(0, 1);
	h_surv_piz->SetMarkerColor(kBlack);
	h_surv_piz->SetLineColor(kBlack);
	h_surv_piz->SetMarkerStyle(20);
	h_surv_piz->SetMarkerSize(0.8);

	h_surv_dalitz_piz->SetMarkerColor(kBlue);
	h_surv_dalitz_piz->SetLineColor(kBlue);
	h_surv_dalitz_piz->SetMarkerStyle(20);
	h_surv_dalitz_piz->SetMarkerSize(0.8);
	h_surv_piz->Draw();
	h_surv_dalitz_piz->Draw("same");

	TLatex latexPiz;
	latexPiz.SetNDC();
	latexPiz.SetTextSize(0.04);
	latexPiz.DrawLatex(.15, .83, "#pi^{0} Electrons");

	TCanvas *cSurvEta = new TCanvas("cSurvEta", "cSurvEta", 600, 600);
	h_surv_eta->GetYaxis()->SetTitle("Veto Survival Rate");
	h_surv_eta->GetYaxis()->SetTitleOffset(1.2);
	h_surv_eta->GetYaxis()->SetTitleFont(62);
	h_surv_eta->GetYaxis()->SetLabelFont(62);
	h_surv_eta->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
	h_surv_eta->GetXaxis()->SetTitleOffset(1.2);
	h_surv_eta->GetXaxis()->SetTitleFont(62);
	h_surv_eta->GetXaxis()->SetLabelFont(62);
	h_surv_eta->GetYaxis()->SetRangeUser(0, 1);
	h_surv_eta->SetMarkerColor(kBlack);
	h_surv_eta->SetLineColor(kBlack);
	h_surv_eta->SetMarkerStyle(20);
	h_surv_eta->SetMarkerSize(0.8);
	h_surv_dalitz_eta->SetMarkerColor(kBlue);
	h_surv_dalitz_eta->SetLineColor(kBlue);
	h_surv_dalitz_eta->SetMarkerStyle(20);
	h_surv_dalitz_eta->SetMarkerSize(0.8);
	h_surv_eta->Draw();
	h_surv_dalitz_eta->Draw("same");

	TLatex latexEta;
	latexEta.SetNDC();
	latexEta.SetTextSize(0.04);
	latexEta.DrawLatex(.15, .83, "#eta Electrons");

	TCanvas *cSurvJpsi = new TCanvas("cSurvJpsi", "cSurvJpsi", 600, 600);
	h_surv_jpsi->GetYaxis()->SetTitle("Veto Survival Rate");
	h_surv_jpsi->GetYaxis()->SetTitleOffset(1.2);
	h_surv_jpsi->GetYaxis()->SetTitleFont(62);
	h_surv_jpsi->GetYaxis()->SetLabelFont(62);
	h_surv_jpsi->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
	h_surv_jpsi->GetXaxis()->SetTitleOffset(1.2);
	h_surv_jpsi->GetXaxis()->SetTitleFont(62);
	h_surv_jpsi->GetXaxis()->SetLabelFont(62);
	h_surv_jpsi->GetYaxis()->SetRangeUser(0, 1);
	h_surv_jpsi->SetMarkerColor(kBlack);
	h_surv_jpsi->SetLineColor(kBlack);
	h_surv_jpsi->SetMarkerStyle(20);
	h_surv_jpsi->SetMarkerSize(0.8);
	h_surv_jpsi->Draw();

	TLatex latexJpsi;
	latexJpsi.SetNDC();
	latexJpsi.SetTextSize(0.04);
	latexJpsi.DrawLatex(.15, .83, "J/#psi Electrons");

	TCanvas *cSurvGamma = new TCanvas("cSurvGamma", "cSurvGamma", 600, 600);
	h_surv_gam->GetYaxis()->SetTitle("Veto Survival Rate");
	h_surv_gam->GetYaxis()->SetTitleOffset(1.2);
	h_surv_gam->GetYaxis()->SetTitleFont(62);
	h_surv_gam->GetYaxis()->SetLabelFont(62);
	h_surv_gam->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
	h_surv_gam->GetXaxis()->SetTitleOffset(1.2);
	h_surv_gam->GetXaxis()->SetTitleFont(62);
	h_surv_gam->GetXaxis()->SetLabelFont(62);
	h_surv_gam->GetYaxis()->SetRangeUser(0, 1);
	h_surv_gam->SetMarkerColor(kBlack);
	h_surv_gam->SetLineColor(kBlack);
	h_surv_gam->SetMarkerStyle(20);
	h_surv_gam->SetMarkerSize(0.8);
	h_surv_gam->Draw();

	TLatex latexGamma;
	latexGamma.SetNDC();
	latexGamma.SetTextSize(0.04);
	latexGamma.DrawLatex(.15, .83, "Direct #gamma Electrons");
}