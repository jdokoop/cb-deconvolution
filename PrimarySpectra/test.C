void test()
{
	TF1 *f_jpsi_current = new TF1("f_jpsi_current", "(2*TMath::Pi()*x)*0.0527311*(([0]/TMath::Power((TMath::Exp(-[1]*TMath::Sqrt(x**2+3.09691**2-0.134976**2)-[2]*TMath::Sqrt(x**2+3.09691**2-0.134976**2)*TMath::Sqrt(x**2+3.09691**2-0.134976**2))+(TMath::Sqrt(x**2+3.09691**2-0.134976**2)/[3])),[4])))", 0.0, 20.0);
	f_jpsi_current->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_jpsi_current->SetLineColor(kAzure - 5);

	TF1 *f_hadron = new TF1("f_hadron","(2*TMath::Pi()*x)*(0.888985/0.0594)*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))",0,20);
	f_hadron->SetParameters(0.000005,0.266872,-0.012166,3.740257,10.077482);

	cout << "mt scaled " << 7.3E-4/f_jpsi_current->Integral(0,20) << endl;
	cout << "hagedorn  " << 7.3E-4/f_hadron->Integral(0,20) << endl;

	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c", "c", 600, 600);
	c->SetLogy();
	f_jpsi_current->Draw();
	f_hadron->Draw("same");
}