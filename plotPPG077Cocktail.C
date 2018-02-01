//-------------------------------------------------------------
// Plot the p+p electron cocktail from PPG077
//
// JDOK
// 01-24-18
//-------------------------------------------------------------

//--------------------------------------
// Variables
//--------------------------------------

//Electrons from each source
TH1F *h_jpsi;
TH1F *h_pizero;
TH1F *h_eta;
TH1F *h_ke3;
TH1F *h_upsilon;
TH1F *h_total;

TH1F *h_jpsi_ratio;
TH1F *h_pizero_ratio;
TH1F *h_eta_ratio;
TH1F *h_ke3_ratio;
TH1F *h_upsilon_ratio;

//Colors
int colors[5] = {kBlue, kRed, kBlack, kOrange - 3, kViolet - 7};

//--------------------------------------
// Functions
//--------------------------------------

void plot()
{
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
	hTemplate2->GetYaxis()->SetTitle("E d^{3}N / dp^{3} [GeV^{-2}/c^{3}]");
	hTemplate2->GetYaxis()->SetTitleOffset(1.3);
	hTemplate2->GetYaxis()->SetRangeUser(0.5e-12, 1E2);
	hTemplate2->Draw();

	h_jpsi->Draw("same");
	h_eta->Draw("same");
	h_pizero->Draw("same");
	h_ke3->Draw("same");
	h_upsilon->Draw("same");

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.DrawLatex(.15, .83, "Electron Background Cocktail");
	latex.DrawLatex(.15, .78, "PPG077");

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(h_pizero, "#pi^{0} Electrons", "PL");
	legend->AddEntry(h_eta, "#eta Electrons", "PL");
	legend->AddEntry(h_ke3, "Ke3 Electrons", "PL");
	legend->AddEntry(h_jpsi, "J/#psi Electrons", "PL");
	legend->AddEntry(h_upsilon, "#Upsilon Electrons", "PL");
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

	h_pizero_ratio->Draw("same");
	h_eta_ratio->Draw("same");
	h_ke3_ratio->Draw("same");
	h_jpsi_ratio->Draw("same");
	h_upsilon_ratio->Draw("same");

	cCocktail->cd();
}

void plotPPG077Cocktail()
{
	//Read histograms in from file
	TFile *fin = new TFile("cocktail_ppg077.root");

	h_jpsi    = (TH1F*) fin->Get("pteJPsi");
	h_pizero  = (TH1F*) fin->Get("ptePion");
	h_eta     = (TH1F*) fin->Get("pteEta");
	h_ke3     = (TH1F*) fin->Get("pteKe3");
	h_upsilon = (TH1F*) fin->Get("pteUpsilon");

	//Add up total
	h_total = (TH1F*) h_pizero->Clone("h_total");
	h_total->Add(h_eta);
	h_total->Add(h_jpsi);
	h_total->Add(h_ke3);
	h_total->Add(h_upsilon);

	//Get ratios to total
	h_jpsi_ratio = (TH1F*) h_jpsi->Clone("h_jpsi_ratio");
	h_pizero_ratio = (TH1F*) h_pizero->Clone("h_pizero_ratio");
	h_eta_ratio = (TH1F*) h_eta->Clone("h_eta_ratio");
	h_ke3_ratio = (TH1F*) h_ke3->Clone("h_ke3_ratio");
	h_upsilon_ratio = (TH1F*) h_upsilon->Clone("h_upsilon_ratio");

	h_jpsi_ratio->Divide(h_total);
	h_pizero_ratio->Divide(h_total);
	h_eta_ratio->Divide(h_total);
	h_ke3_ratio->Divide(h_total);
	h_upsilon_ratio->Divide(h_total);

	//Format each species
	h_total->SetLineColor(kBlack);
	h_total->SetMarkerColor(kBlack);
	h_total->SetMarkerColor(kBlack);
	h_total->SetMarkerStyle(20);
	h_total->SetMarkerSize(0.4);

	h_pizero->SetLineColor(colors[0]);
	h_eta->SetLineColor(colors[1]);
	h_jpsi->SetLineColor(colors[2]);
	h_ke3->SetLineColor(colors[3]);
	h_upsilon->SetLineColor(colors[4]);

	h_pizero->SetMarkerColor(colors[0]);
	h_eta->SetMarkerColor(colors[1]);
	h_jpsi->SetMarkerColor(colors[2]);
	h_ke3->SetMarkerColor(colors[3]);
	h_upsilon->SetMarkerColor(colors[4]);

	h_pizero->SetMarkerStyle(21);
	h_eta->SetMarkerStyle(21);
	h_jpsi->SetMarkerStyle(21);
	h_ke3->SetMarkerStyle(21);
	h_upsilon->SetMarkerStyle(21);

	h_pizero->SetMarkerSize(0.4);
	h_eta->SetMarkerSize(0.4);
	h_jpsi->SetMarkerSize(0.4);
	h_ke3->SetMarkerSize(0.4);
	h_upsilon->SetMarkerSize(0.4);

	h_pizero_ratio->SetLineColor(colors[0]);
	h_eta_ratio->SetLineColor(colors[1]);
	h_jpsi_ratio->SetLineColor(colors[2]);
	h_ke3_ratio->SetLineColor(colors[3]);
	h_upsilon_ratio->SetLineColor(colors[4]);

	h_pizero_ratio->SetMarkerColor(colors[0]);
	h_eta_ratio->SetMarkerColor(colors[1]);
	h_jpsi_ratio->SetMarkerColor(colors[2]);
	h_ke3_ratio->SetMarkerColor(colors[3]);
	h_upsilon_ratio->SetMarkerColor(colors[4]);

	h_pizero_ratio->SetMarkerStyle(21);
	h_eta_ratio->SetMarkerStyle(21);
	h_jpsi_ratio->SetMarkerStyle(21);
	h_ke3_ratio->SetMarkerStyle(21);
	h_upsilon_ratio->SetMarkerStyle(21);

	h_pizero_ratio->SetMarkerSize(0.4);
	h_eta_ratio->SetMarkerSize(0.4);
	h_jpsi_ratio->SetMarkerSize(0.4);
	h_ke3_ratio->SetMarkerSize(0.4);
	h_upsilon_ratio->SetMarkerSize(0.4);

	//Plot things
	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c", "c", 700, 700);

	c->SetLeftMargin(1.0);
	c->SetLogy();

	TH1F *hTemplate = new TH1F("hTemplate", "hTemplate", 100, 0, 10);
	hTemplate->SetTitle("");
	hTemplate->GetXaxis()->SetTitleOffset(1.3);
	hTemplate->GetXaxis()->SetRangeUser(0, 10);
	hTemplate->GetXaxis()->SetTitleFont(62);
	hTemplate->GetXaxis()->SetLabelFont(62);
	hTemplate->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	hTemplate->GetYaxis()->SetTitleOffset(1.4);
	hTemplate->GetYaxis()->SetRangeUser(1E-12, 100);
	hTemplate->GetYaxis()->SetTitleFont(62);
	hTemplate->GetYaxis()->SetLabelFont(62);
	hTemplate->GetYaxis()->SetTitle("E d^{3}N/dp^{3} [GeV^{-2}/c^{3}]");

	hTemplate->Draw();
	h_pizero->Draw("same");
	h_eta->Draw("same");
	h_jpsi->Draw("same");
	h_ke3->Draw("same");
	h_upsilon->Draw("same");
	h_total->Draw("same");

	plot();
}