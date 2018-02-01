//---------------------------------------------
// Fit J/psi yield in p+p from PPG104
//
// JDOK
// 11-14-17
//---------------------------------------------

#include <iostream>

using namespace std;

//------------------------
// Variables
//------------------------

const int NPOINTS = 24;

//Published spectrum and its fit
TGraphErrors *g_spectrum;
TF1 *f_published_fit;
TBox b_syst[NPOINTS];

//Fit by mT scaling
TF1 *f_mtscaled_fit;

//Ratio
TGraphErrors *g_ratio_fit_data;
TGraphErrors *g_ratio_mtscaled_data;

//------------------------
// Functions
//------------------------

void constructPublishedSpectrum()
{
	float x_vals[NPOINTS] = {0.125,
	                         0.375,
	                         0.625,
	                         0.875,
	                         1.125,
	                         1.375,
	                         1.625,
	                         1.875,
	                         2.125,
	                         2.375,
	                         2.625,
	                         2.875,
	                         3.125,
	                         3.375,
	                         3.625,
	                         3.875,
	                         4.125,
	                         4.375,
	                         4.625,
	                         4.875,
	                         5.5,
	                         6.5,
	                         7.5,
	                         8.5
	                        };

	float y_vals[NPOINTS] = {4.9E-6,
	                         3.9E-6,
	                         3.9E-6,
	                         3.5E-6,
	                         3.19E-6,
	                         2.21E-6,
	                         1.69E-6,
	                         1.42E-6,
	                         0.95E-6,
	                         0.66E-6,
	                         0.56E-6,
	                         0.37E-6,
	                         0.29E-6,
	                         0.199E-6,
	                         0.136E-6,
	                         0.106E-6,
	                         0.0689E-6,
	                         0.044E-6,
	                         0.047E-6,
	                         0.031E-6,
	                         0.0135E-6,
	                         0.0041E-6,
	                         0.0016E-6,
	                         0.00037E-6
	                        };

	float err_y[NPOINTS] = {0.6E-6,
	                        0.5E-6,
	                        0.5E-6,
	                        0.4E-6,
	                        0.38E-6,
	                        0.27E-6,
	                        0.2E-6,
	                        0.17E-6,
	                        0.12E-6,
	                        0.08E-6,
	                        0.07E-6,
	                        0.05E-6,
	                        0.04E-6,
	                        0.024E-6,
	                        0.017E-6,
	                        0.013E-6,
	                        0.008E-6,
	                        0.005E-6,
	                        0.006E-6,
	                        0.004E-6,
	                        0.00175E-6,
	                        0.00055E-6,
	                        0.0002E-6,
	                        0.00055E-6
	                       };

	float err_x[NPOINTS] = {0.0};

	float err_syst[NPOINTS] = {0.5E-6,
	                           0.3E-6,
	                           0.2E-6,
	                           0.2E-6,
	                           0.17E-6,
	                           0.14E-6,
	                           0.12E-6,
	                           0.1E-6,
	                           0.08E-6,
	                           0.07E-6,
	                           0.06E-6,
	                           0.05E-6,
	                           0.04E-6,
	                           0.0355E-6,
	                           0.028E-6,
	                           0.025E-6,
	                           0.0225E-6,
	                           0.0165E-6,
	                           0.0145E-6,
	                           0.012E-6,
	                           0.00315E-6,
	                           0.0013E-6,
	                           0.0007E-6,
	                           0.000295E-6
	                          };

	g_spectrum = new TGraphErrors(NPOINTS, x_vals, y_vals, err_x, err_y);

	//Construct systematics
	for (int i = 0; i < NPOINTS; i++)
	{
		b_syst[i] = TBox(x_vals[i] - 0.1, y_vals[i] - err_syst[i] / 2.0, x_vals[i] + 0.1, y_vals[i] + err_syst[i] / 2.0);
		b_syst[i].SetLineColor(kBlack);
		b_syst[i].SetFillStyle(0);
	}

	//Fit to published data
	f_published_fit = new TF1("f_published_fit", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_fit->SetParameters(4.51297e-06, 3.74161e-01, 2.62854e-03, 2.76577e+00, 7.83726e+00);
	g_spectrum->Fit(f_published_fit, "R");

	//Fit by mT scaling
	//Here we put in the branching ratio B = 5.94% pm 0.06 from PPG104
	f_mtscaled_fit = new TF1("f_mtscaled_fit", "0.0594*0.054*(([0]/TMath::Power((TMath::Exp(-[1]*TMath::Sqrt(3.096*3.096-0.135*0.135+x*x)-[2]*TMath::Sqrt(3.096*3.096-0.135*0.135+x*x)*TMath::Sqrt(3.096*3.096-0.135*0.135+x*x))+(TMath::Sqrt(3.096*3.096-0.135*0.135+x*x)/[3])),[4])))", 0, 18);
	f_mtscaled_fit->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_mtscaled_fit->SetLineColor(kBlue);

	g_spectrum->SetMarkerStyle(20);
	g_spectrum->SetMarkerSize(0.7);
	g_spectrum->SetMarkerColor(kBlack);
	g_spectrum->SetLineColor(kBlack);
}


void getRatio()
{
	g_ratio_fit_data = new TGraphErrors();

	for (int i = 0; i < NPOINTS; i++)
	{
		double x, y, ex, ey;
		g_spectrum->GetPoint(i, x, y);
		ex = 0;
		ey = g_spectrum->GetErrorY(i);

		g_ratio_fit_data->SetPoint(i, x, y / f_published_fit->Eval(x));
		g_ratio_fit_data->SetPointError(i, 0.0, ey / f_published_fit->Eval(x));
	}

	g_ratio_mtscaled_data = new TGraphErrors();

	for (int i = 0; i < NPOINTS; i++)
	{
		double x, y, ex, ey;
		g_spectrum->GetPoint(i, x, y);
		ex = 0;
		ey = g_spectrum->GetErrorY(i);

		g_ratio_mtscaled_data->SetPoint(i, x, y / f_mtscaled_fit->Eval(x));
		g_ratio_mtscaled_data->SetPointError(i, 0.0, ey / f_mtscaled_fit->Eval(x));
	}
}


void plot()
{
	TCanvas *cNP = new TCanvas("cNP", "Stacked Representation of Fit", 700, 900);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1);
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
	hTemplate2->GetXaxis()->SetRangeUser(0, 18);
	hTemplate2->GetYaxis()->SetTitleFont(62);
	hTemplate2->GetYaxis()->SetLabelFont(62);
	hTemplate2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate2->GetYaxis()->SetTitle("B_{ee} #times E d^{3}#sigma / dp^{3}");
	hTemplate2->GetYaxis()->SetTitleOffset(1.3);
	hTemplate2->GetYaxis()->SetRangeUser(3E-12, 1E-4);
	hTemplate2->Draw();

	g_spectrum->SetTitle("");
	g_spectrum->GetXaxis()->SetTitleFont(62);
	g_spectrum->GetXaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetRangeUser(0, 18);
	g_spectrum->GetYaxis()->SetTitleFont(62);
	g_spectrum->GetYaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum->GetYaxis()->SetTitle("B_{ee} #times E d^{3}#sigma / dp^{3}");
	g_spectrum->GetYaxis()->SetRangeUser(1e-8, 150);
	g_spectrum->SetMarkerStyle(20);
	g_spectrum->SetMarkerSize(0.8);
	g_spectrum->SetMarkerColor(kBlack);
	g_spectrum->Draw("P,same");
	f_published_fit->Draw("same");

	f_mtscaled_fit->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		b_syst[i].Draw("same");
	}

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.DrawLatex(.15, .83, "PHENIX J/#psi Spectrum");
	latex.DrawLatex(.15, .78, "PPG104");

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(f_published_fit, "Modified Hagedorn Fit", "l");
	legend->AddEntry(f_mtscaled_fit, "m_{T} Scaling of #pi^{0} Fit", "l");
	legend->SetFillStyle(0.0);
	legend->SetLineColor(kWhite);
	legend->Draw("same");

	cNP->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.23, 1, 0.4);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();

	TH1F *hTemplate3 = new TH1F("hTemplate3", "hTemplate3", 100, 0, 20);
	hTemplate3->SetTitle("");
	hTemplate3->GetYaxis()->CenterTitle();
	hTemplate3->GetXaxis()->SetTitleFont(62);
	hTemplate3->GetXaxis()->SetLabelFont(62);
	hTemplate3->GetXaxis()->SetRangeUser(0, 18);
	hTemplate3->GetYaxis()->SetTitleFont(62);
	hTemplate3->GetYaxis()->SetLabelFont(62);
	hTemplate3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate3->GetYaxis()->SetTitle("Data / Fit");
	hTemplate3->GetYaxis()->SetTitleOffset(0.4);
	hTemplate3->GetYaxis()->SetRangeUser(0.7, 1.28);
	hTemplate3->GetYaxis()->SetLabelSize(0.11);
	hTemplate3->GetXaxis()->SetTitleSize(0.11);
	hTemplate3->GetYaxis()->SetTitleSize(0.11);
	hTemplate3->GetXaxis()->SetTitleOffset(1.2);
	hTemplate3->GetXaxis()->SetLabelSize(0.11);
	hTemplate3->Draw();

	g_ratio_fit_data->SetMarkerStyle(20);
	g_ratio_fit_data->SetMarkerSize(0.7);
	g_ratio_fit_data->GetXaxis()->SetTitleFont(62);
	g_ratio_fit_data->GetXaxis()->SetLabelFont(62);
	g_ratio_fit_data->GetXaxis()->SetRangeUser(0, 18);
	g_ratio_fit_data->GetYaxis()->SetTitleFont(62);
	g_ratio_fit_data->GetYaxis()->SetLabelFont(62);
	g_ratio_fit_data->SetTitle("");
	g_ratio_fit_data->GetYaxis()->CenterTitle();
	g_ratio_fit_data->GetYaxis()->SetRangeUser(0.7, 1.28);
	g_ratio_fit_data->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_ratio_fit_data->GetYaxis()->SetTitle("Data / Fit");
	g_ratio_fit_data->GetYaxis()->SetTitleSize(0.14);
	g_ratio_fit_data->GetYaxis()->SetTitleOffset(0.4);
	g_ratio_fit_data->GetYaxis()->SetLabelSize(0.12);
	g_ratio_fit_data->GetXaxis()->SetTitleSize(0.12);
	g_ratio_fit_data->GetXaxis()->SetTitleOffset(1.2);
	g_ratio_fit_data->GetXaxis()->SetLabelSize(0.14);
	g_ratio_fit_data->SetMarkerColor(kRed);
	g_ratio_fit_data->SetLineColor(kRed);
	g_ratio_fit_data->Draw("P,same");

	TLine *lref1 = new TLine(0.0, 1.0, 18.0, 1.0);
	lref1->SetLineStyle(7);
	lref1->Draw("same");
	cNP->cd();

	cNP->cd();
	TPad *pad3 = new TPad("pad3", "pad3", 0, 0, 1, 0.23);
	pad3->SetTopMargin(0);
	pad3->SetBottomMargin(0.3);
	pad3->Draw();
	pad3->cd();
	pad3->SetTickx();
	pad3->SetTicky();

	TH1F *hTemplate4 = new TH1F("hTemplate4", "hTemplate4", 100, 0, 20);
	hTemplate4->SetTitle("");
	hTemplate4->GetYaxis()->CenterTitle();
	hTemplate4->GetXaxis()->SetTitleFont(62);
	hTemplate4->GetXaxis()->SetLabelFont(62);
	hTemplate4->GetXaxis()->SetRangeUser(0, 18);
	hTemplate4->GetYaxis()->SetTitleFont(62);
	hTemplate4->GetYaxis()->SetLabelFont(62);
	hTemplate4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate4->GetYaxis()->SetTitle("Data / Fit");
	hTemplate4->GetYaxis()->SetTitleOffset(0.4);
	hTemplate4->GetYaxis()->SetRangeUser(0.7, 1.28);
	hTemplate4->GetYaxis()->SetLabelSize(0.09);
	hTemplate4->GetXaxis()->SetTitleSize(0.11);
	hTemplate4->GetYaxis()->SetTitleSize(0.09);
	hTemplate4->GetXaxis()->SetTitleOffset(1.2);
	hTemplate4->GetXaxis()->SetLabelSize(0.09);
	hTemplate4->Draw();

	g_ratio_mtscaled_data->SetMarkerStyle(20);
	g_ratio_mtscaled_data->SetMarkerSize(0.7);
	g_ratio_mtscaled_data->GetXaxis()->SetTitleFont(62);
	g_ratio_mtscaled_data->GetXaxis()->SetLabelFont(62);
	g_ratio_mtscaled_data->GetXaxis()->SetRangeUser(0, 18);
	g_ratio_mtscaled_data->GetYaxis()->SetTitleFont(62);
	g_ratio_mtscaled_data->GetYaxis()->SetLabelFont(62);
	g_ratio_mtscaled_data->SetTitle("");
	g_ratio_mtscaled_data->GetYaxis()->CenterTitle();
	g_ratio_mtscaled_data->GetYaxis()->SetRangeUser(0.7, 1.28);
	g_ratio_mtscaled_data->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_ratio_mtscaled_data->GetYaxis()->SetTitle("Data / Fit");
	g_ratio_mtscaled_data->GetYaxis()->SetTitleSize(0.14);
	g_ratio_mtscaled_data->GetYaxis()->SetTitleOffset(0.4);
	g_ratio_mtscaled_data->GetYaxis()->SetLabelSize(0.09);
	g_ratio_mtscaled_data->GetXaxis()->SetTitleSize(0.1);
	g_ratio_mtscaled_data->GetXaxis()->SetTitleOffset(1.2);
	g_ratio_mtscaled_data->GetXaxis()->SetLabelSize(0.14);
	g_ratio_mtscaled_data->SetMarkerColor(kBlue);
	g_ratio_mtscaled_data->SetLineColor(kBlue);
	g_ratio_mtscaled_data->Draw("P,same");

	lref1->Draw("same");
	cNP->cd();
}


void fitJpsi()
{
	gStyle->SetOptStat(0);

	constructPublishedSpectrum();
	getRatio();
	plot();
}