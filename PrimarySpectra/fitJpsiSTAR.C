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

const int NPOINTS_STAR = 11;
const int NPOINTS_PHENIX = 24;
const int NPOINTS_COMBINED = 26;

//Published spectrum and its fit
TGraphErrors *g_spectrum_star;
TBox b_syst_star[NPOINTS_STAR];
TF1 *f_published_fit_star;

//PHENIX published data
TGraphErrors *g_spectrum_phenix;
TBox b_syst_phenix[NPOINTS_PHENIX];
TF1 *f_published_fit_phenix;

//PHENIX and STAR combined published data
TGraphErrors *g_spectrum_combined;
TBox b_syst_combined[NPOINTS_COMBINED];
TF1 *f_published_fit_combined;

//Ratio of combined data to its fit
TGraphErrors *g_ratio_fit_data_combined;

//Ratio of STAR data to its fit
TGraphErrors *g_ratio_fit_data_star;

//Ratio of PHENIX data to its fit
TGraphErrors *g_ratio_fit_data_phenix;

//Fit used in the current version of the cocktail analysis
TF1 *f_jpsi_current;

//------------------------
// Functions
//------------------------

void constructCurrentFit()
{
	f_jpsi_current = new TF1("f_jpsi_current", "0.0594*0.0527311*(([0]/TMath::Power((TMath::Exp(-[1]*TMath::Sqrt(x**2+3.09691**2-0.134976**2)-[2]*TMath::Sqrt(x**2+3.09691**2-0.134976**2)*TMath::Sqrt(x**2+3.09691**2-0.134976**2))+(TMath::Sqrt(x**2+3.09691**2-0.134976**2)/[3])),[4])))", 0.0, 20.0);
	f_jpsi_current->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_jpsi_current->SetLineColor(kAzure - 5);
}

void constructPublishedSpectrumCombined()
{
	//Combine  using PHENIX values below 6 GeV, and STAR points above that

	/////////////////////// PHX ///////////////////////
	float x_vals_phx[NPOINTS_PHENIX] = {0.125,
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

	float y_vals_phx[NPOINTS_PHENIX] = {4.9E-6,
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

	float err_y_phx[NPOINTS_PHENIX] = {0.6E-6,
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

	float err_syst_phx[NPOINTS_PHENIX] = {0.5E-6,
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
	///////////////////////////////////////////////////

	//////////////////////// STR //////////////////////
	float x_vals_str[NPOINTS_STAR] = {2.25,
	                                  2.75,
	                                  3.25,
	                                  3.75,
	                                  4.5,
	                                  5.5,
	                                  6.5,
	                                  7.5,
	                                  9,
	                                  11,
	                                  13
	                                 };

	float y_vals_str[NPOINTS_STAR] = {0.68,
	                                  0.318,
	                                  0.187,
	                                  0.1032,
	                                  0.0334,
	                                  0.00905,
	                                  0.00154,
	                                  0.00084,
	                                  2.01E-04,
	                                  4.55E-05,
	                                  9.71E-06
	                                 };

	float err_y_str[NPOINTS_STAR] = {0.17,
	                                 0.071,
	                                 0.032,
	                                 0.0189,
	                                 0.0058,
	                                 0.00183,
	                                 0.00052,
	                                 0.00029,
	                                 0.0000747,
	                                 0.0000186,
	                                 0.00000487
	                                };

	float err_syst_str[NPOINTS_STAR] = {0.03,
	                                    0.021,
	                                    0.009,
	                                    0.0071,
	                                    0.0024,
	                                    0.00066,
	                                    0.0001,
	                                    0.00012,
	                                    0.0000505,
	                                    0.0000114,
	                                    0.00000246
	                                   };
	///////////////////////////////////////////////////

	float x_vals[NPOINTS_COMBINED] = {0.0};
	float y_vals[NPOINTS_COMBINED] = {0.0};
	float err_y[NPOINTS_COMBINED] = {0.0};
	float err_syst[NPOINTS_COMBINED] = {0.0};
	float err_x[NPOINTS_COMBINED] = {0.0};

	for (int i = 0; i <= 21; i++)
	{
		x_vals[i] = x_vals_phx[i];
		y_vals[i] = y_vals_phx[i];
		err_y[i] = err_y_phx[i];
		err_syst[i] = err_syst_phx[i];
	}

	for (int i = 22; i <= 25; i++)
	{
		x_vals[i] = x_vals_str[7 + (i - 22)];
		y_vals[i] = y_vals_str[7 + (i - 22)];
		err_y[i] = err_y_str[7 + (i - 22)];
		err_syst[i] = err_syst_str[7 + (i - 22)];
	}

	//Set the units to mb for the STAR points
	for (int i = 22; i <= 25; i++)
	{
		y_vals[i] = (float) y_vals[i] / 1E6;
		err_y[i] = (float) err_y[i] / 1E6;
		err_syst[i] = (float) err_syst[i] / 1E6;
	}

	g_spectrum_combined = new TGraphErrors(NPOINTS_COMBINED, x_vals, y_vals, err_x, err_y);

	//Construct systematics
	for (int i = 0; i < NPOINTS_COMBINED; i++)
	{
		b_syst_combined[i] = TBox(x_vals[i] - 0.1, y_vals[i] - err_syst[i] / 2.0, x_vals[i] + 0.1, y_vals[i] + err_syst[i] / 2.0);
		b_syst_combined[i].SetLineColor(kBlack);
		b_syst_combined[i].SetFillStyle(0);
	}

	//Fit to published data
	f_published_fit_combined = new TF1("f_published_fit_combined", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 20.0);
	f_published_fit_combined->SetParameters(4.51297e-06, 3.74161e-01, 2.62854e-03, 2.76577e+00, 7.83726e+00);
	g_spectrum_combined->Fit(f_published_fit_combined, "Q0R");

	g_spectrum_combined->SetMarkerStyle(20);
	g_spectrum_combined->SetMarkerSize(0.7);
	g_spectrum_combined->SetMarkerColor(kBlack);
	g_spectrum_combined->SetLineColor(kBlack);
}


void constructPublishedSpectrumPHENIX()
{
	float x_vals[NPOINTS_PHENIX] = {0.125,
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

	float y_vals[NPOINTS_PHENIX] = {4.9E-6,
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

	float err_y[NPOINTS_PHENIX] = {0.6E-6,
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

	float err_x[NPOINTS_PHENIX] = {0.0};

	float err_syst[NPOINTS_PHENIX] = {0.5E-6,
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

	g_spectrum_phenix = new TGraphErrors(NPOINTS_PHENIX, x_vals, y_vals, err_x, err_y);

	//Construct systematics
	for (int i = 0; i < NPOINTS_PHENIX; i++)
	{
		b_syst_phenix[i] = TBox(x_vals[i] - 0.1, y_vals[i] - err_syst[i] / 2.0, x_vals[i] + 0.1, y_vals[i] + err_syst[i] / 2.0);
		b_syst_phenix[i].SetLineColor(kBlack);
		b_syst_phenix[i].SetFillStyle(0);
	}

	//Fit to published data
	f_published_fit_phenix = new TF1("f_published_fit_phenix", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_fit_phenix->SetParameters(-6.00000e+02, -1.63212e-04, -3.26764e+00, 1.62532e+00, 8.19851e+00);
	g_spectrum_phenix->Fit(f_published_fit_phenix, "Q0R");

	g_spectrum_phenix->SetMarkerStyle(20);
	g_spectrum_phenix->SetMarkerSize(0.7);
	g_spectrum_phenix->SetMarkerColor(kBlack);
	g_spectrum_phenix->SetLineColor(kBlack);
}


void constructPublishedSpectrumSTAR()
{
	float x_vals[NPOINTS_STAR] = {2.25,
	                              2.75,
	                              3.25,
	                              3.75,
	                              4.5,
	                              5.5,
	                              6.5,
	                              7.5,
	                              9,
	                              11,
	                              13
	                             };

	float y_vals[NPOINTS_STAR] = {0.68,
	                              0.318,
	                              0.187,
	                              0.1032,
	                              0.0334,
	                              0.00905,
	                              0.00154,
	                              0.00084,
	                              2.01E-04,
	                              4.55E-05,
	                              9.71E-06
	                             };

	float err_y[NPOINTS_STAR] = {0.17,
	                             0.071,
	                             0.032,
	                             0.0189,
	                             0.0058,
	                             0.00183,
	                             0.00052,
	                             0.00029,
	                             0.0000747,
	                             0.0000186,
	                             0.00000487
	                            };

	float err_x[NPOINTS_STAR] = {0.0};

	float err_syst[NPOINTS_STAR] = {0.03,
	                                0.021,
	                                0.009,
	                                0.0071,
	                                0.0024,
	                                0.00066,
	                                0.0001,
	                                0.00012,
	                                0.0000505,
	                                0.0000114,
	                                0.00000246
	                               };

	//Set the units to mb
	for (int i = 0; i < NPOINTS_STAR; i++)
	{
		y_vals[i] = (float) y_vals[i] / 1E6;
		err_y[i] = (float) err_y[i] / 1E6;
		err_syst[i] = (float) err_syst[i] / 1E6;
	}

	g_spectrum_star = new TGraphErrors(NPOINTS_STAR, x_vals, y_vals, err_x, err_y);

	//Construct systematics
	for (int i = 0; i < NPOINTS_STAR; i++)
	{
		b_syst_star[i] = TBox(x_vals[i] - 0.1, y_vals[i] - err_syst[i] / 2.0, x_vals[i] + 0.1, y_vals[i] + err_syst[i] / 2.0);
		b_syst_star[i].SetLineColor(kBlack);
		b_syst_star[i].SetFillStyle(0);
	}

	//Fit to published data
	f_published_fit_star = new TF1("f_published_fit_star", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_fit_star->SetParameters(4.91496e-06, 3.17564e-01, 1.90922e-02, 2.92862e+00, 8.99895e+00);
	g_spectrum_star->Fit(f_published_fit_star, "Q0R");
	f_published_fit_star->SetLineColor(kSpring - 6);
	g_spectrum_star->SetMarkerStyle(20);
	g_spectrum_star->SetMarkerSize(0.7);
	g_spectrum_star->SetMarkerColor(kBlack);
	g_spectrum_star->SetLineColor(kBlack);
}


void getRatio()
{
	//Combined
	g_ratio_fit_data_combined = new TGraphErrors();

	for (int i = 0; i < NPOINTS_COMBINED; i++)
	{
		double x, y, ex, ey;
		g_spectrum_combined->GetPoint(i, x, y);
		ex = 0;
		ey = g_spectrum_combined->GetErrorY(i);

		g_ratio_fit_data_combined->SetPoint(i, x, y / f_published_fit_combined->Eval(x));
		g_ratio_fit_data_combined->SetPointError(i, 0.0, ey / f_published_fit_combined->Eval(x));
	}

	//STAR
	g_ratio_fit_data_star = new TGraphErrors();

	for (int i = 0; i < NPOINTS_STAR; i++)
	{
		double x, y, ex, ey;
		g_spectrum_star->GetPoint(i, x, y);
		ex = 0;
		ey = g_spectrum_star->GetErrorY(i);

		g_ratio_fit_data_star->SetPoint(i, x, y / f_published_fit_star->Eval(x));
		g_ratio_fit_data_star->SetPointError(i, 0.0, ey / f_published_fit_star->Eval(x));
	}

	//PHENIX
	g_ratio_fit_data_phenix = new TGraphErrors();

	for (int i = 0; i < NPOINTS_PHENIX; i++)
	{
		double x, y, ex, ey;
		g_spectrum_phenix->GetPoint(i, x, y);
		ex = 0;
		ey = g_spectrum_phenix->GetErrorY(i);

		g_ratio_fit_data_phenix->SetPoint(i, x, y / f_jpsi_current->Eval(x));
		g_ratio_fit_data_phenix->SetPointError(i, 0.0, ey / f_jpsi_current->Eval(x));
	}
}


void plotCombined()
{
	TCanvas *cNP = new TCanvas("cNP", "JPsi Spectrum", 700, 900);
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

	g_spectrum_phenix->SetTitle("");
	g_spectrum_phenix->GetXaxis()->SetTitleFont(62);
	g_spectrum_phenix->GetXaxis()->SetLabelFont(62);
	g_spectrum_phenix->GetXaxis()->SetRangeUser(0, 18);
	g_spectrum_phenix->GetYaxis()->SetTitleFont(62);
	g_spectrum_phenix->GetYaxis()->SetLabelFont(62);
	g_spectrum_phenix->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum_phenix->GetYaxis()->SetTitle("B_{ee} #times E d^{3}#sigma / dp^{3}");
	g_spectrum_phenix->GetYaxis()->SetRangeUser(1e-8, 150);
	g_spectrum_phenix->SetMarkerStyle(20);
	g_spectrum_phenix->SetMarkerSize(0.8);
	g_spectrum_phenix->SetMarkerColor(kBlack);
	g_spectrum_phenix->Draw("P,same");

	g_spectrum_star->SetMarkerStyle(24);
	g_spectrum_star->SetMarkerSize(0.8);
	g_spectrum_star->SetMarkerColor(kBlack);
	g_spectrum_star->Draw("P,same");

	f_published_fit_combined->Draw("same");
	f_jpsi_current->Draw("same");
	//f_published_fit_star->Draw("same");

	for (int i = 0; i < NPOINTS_STAR; i++)
	{
		b_syst_combined[i].Draw("same");
	}

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.DrawLatex(.15, .83, "PHENIX and STAR J/#psi Spectra #times B_{ee}");
	latex.DrawLatex(.15, .78, "PPG104 + PLB 722 55-62");

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(g_spectrum_phenix, "PHENIX PPG104 Data", "P");
	legend->AddEntry(g_spectrum_star, "STAR PLB 722 55-62 Data", "P");
	legend->AddEntry(f_published_fit_combined, "Modified Hagedorn Fit to Combined Data*", "l");
	//legend->AddEntry(f_published_fit_star, "Modified Hagedorn Fit to STAR Data", "l");
	legend->AddEntry(f_jpsi_current, "m_{T} Scaled Fit to PHENIX Data", "l");
	legend->SetFillStyle(0.0);
	legend->SetLineColor(kWhite);
	legend->Draw("same");

	cNP->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.24, 1, 0.4);
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
	hTemplate3->GetYaxis()->SetTitle("Combined Data / Fit");
	hTemplate3->GetYaxis()->SetTitleOffset(0.4);
	hTemplate3->GetYaxis()->SetRangeUser(0.75, 1.55);
	hTemplate3->GetYaxis()->SetLabelSize(0.11);
	hTemplate3->GetXaxis()->SetTitleSize(0.11);
	hTemplate3->GetYaxis()->SetTitleSize(0.11);
	hTemplate3->GetXaxis()->SetTitleOffset(1.2);
	hTemplate3->GetXaxis()->SetLabelSize(0.11);
	hTemplate3->Draw();

	g_ratio_fit_data_combined->SetMarkerStyle(20);
	g_ratio_fit_data_combined->SetMarkerSize(0.7);
	g_ratio_fit_data_combined->GetXaxis()->SetTitleFont(62);
	g_ratio_fit_data_combined->GetXaxis()->SetLabelFont(62);
	g_ratio_fit_data_combined->GetXaxis()->SetRangeUser(0, 18);
	g_ratio_fit_data_combined->GetYaxis()->SetTitleFont(62);
	g_ratio_fit_data_combined->GetYaxis()->SetLabelFont(62);
	g_ratio_fit_data_combined->SetTitle("");
	g_ratio_fit_data_combined->GetYaxis()->CenterTitle();
	g_ratio_fit_data_combined->GetYaxis()->SetRangeUser(0.75, 1.55);
	g_ratio_fit_data_combined->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_ratio_fit_data_combined->GetYaxis()->SetTitle("Data / Fit");
	g_ratio_fit_data_combined->GetYaxis()->SetTitleSize(0.14);
	g_ratio_fit_data_combined->GetYaxis()->SetTitleOffset(0.4);
	g_ratio_fit_data_combined->GetYaxis()->SetLabelSize(0.12);
	g_ratio_fit_data_combined->GetXaxis()->SetTitleSize(0.12);
	g_ratio_fit_data_combined->GetXaxis()->SetTitleOffset(1.2);
	g_ratio_fit_data_combined->GetXaxis()->SetLabelSize(0.14);
	g_ratio_fit_data_combined->SetMarkerColor(kRed);
	g_ratio_fit_data_combined->SetLineColor(kRed);
	g_ratio_fit_data_combined->Draw("P,same");

	TLine *lref1 = new TLine(0.0, 1.0, 18.0, 1.0);
	lref1->SetLineStyle(7);
	lref1->Draw("same");
	cNP->cd();

	cNP->cd();
	TPad *pad4 = new TPad("pad4", "pad4", 0, 0, 1, 0.24);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->Draw();
	pad4->cd();
	pad4->SetTickx();
	pad4->SetTicky();

	TH1F *hTemplate5 = new TH1F("hTemplate5", "hTemplate5", 100, 0, 20);
	hTemplate5->SetTitle("");
	hTemplate5->GetYaxis()->CenterTitle();
	hTemplate5->GetXaxis()->SetTitleFont(62);
	hTemplate5->GetXaxis()->SetLabelFont(62);
	hTemplate5->GetXaxis()->SetRangeUser(0, 18);
	hTemplate5->GetYaxis()->SetTitleFont(62);
	hTemplate5->GetYaxis()->SetLabelFont(62);
	hTemplate5->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate5->GetYaxis()->SetTitle("PHENIX Data / Fit");
	hTemplate5->GetYaxis()->SetTitleOffset(0.5);
	hTemplate5->GetYaxis()->SetRangeUser(0.75, 1.55);
	hTemplate5->GetYaxis()->SetLabelSize(0.09);
	hTemplate5->GetXaxis()->SetTitleSize(0.11);
	hTemplate5->GetYaxis()->SetTitleSize(0.09);
	hTemplate5->GetXaxis()->SetTitleOffset(1.2);
	hTemplate5->GetXaxis()->SetLabelSize(0.09);
	hTemplate5->Draw();

	g_ratio_fit_data_phenix->SetMarkerColor(kAzure - 5);
	g_ratio_fit_data_phenix->SetLineColor(kAzure - 5);
	g_ratio_fit_data_phenix->SetMarkerStyle(20);
	g_ratio_fit_data_phenix->SetMarkerSize(0.7);
	g_ratio_fit_data_phenix->Draw("P,same");
	lref1->Draw("same");

	cNP->cd();
}


void printFunctions()
{
	//Combined Fit
	cout << "*********** Combined Fit **************" << endl;
	cout << "f_hadron = new TF1(\"f_hadron\",\"([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))\",0,20);" << endl;
	cout << Form("f_hadron->SetParameters(%f,%f,%f,%f,%f)", f_published_fit_combined->GetParameter(0), f_published_fit_combined->GetParameter(1), f_published_fit_combined->GetParameter(2), f_published_fit_combined->GetParameter(3), f_published_fit_combined->GetParameter(4)) << endl << endl;
}


void fitJpsiSTAR()
{
	gStyle->SetOptStat(0);

	constructCurrentFit();
	constructPublishedSpectrumSTAR();
	constructPublishedSpectrumPHENIX();
	constructPublishedSpectrumCombined();
	getRatio();
	printFunctions();
	//plot();
	plotCombined();
}