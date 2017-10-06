#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::endl;

Double_t par[9];

Double_t fun_tripleGauss(Double_t *x, Double_t *par)
{
	Double_t cte = 2.5066297;	

        Double_t G;
	Double_t N1, Function1, G1;
	Double_t N2, Function2, G2;
	Double_t N3, Function3, G3;

	N1 = 1.0/(cte*par[2]);
	Function1 = N1*(TMath::Exp(- (par[1]-x[0]) * (par[1] - x[0]) / (2.0 * par[2] * par[2]) ) );
	G1 = par[0]*Function1;

	N2 = 1.0/(cte*par[5]);
	Function2 = N2*(TMath::Exp(- (par[4]-x[0]) * (par[4] - x[0]) / (2.0 * par[5] * par[5]) ) );
	G2 = par[3]*Function2;

	N3 = 1.0/(cte*par[8]);
	Function3 = N3*(TMath::Exp(- (par[7]-x[0]) * (par[7] - x[0]) / (2.0 * par[8] * par[8]) ) );
	G3 = par[6]*Function3;

        G = G1+G2+G3; 
 
	return G;

}

Double_t fun_Gauss(Double_t *x, Double_t *par)
{
	Double_t cte = 2.5066297;	

        Double_t G;
	Double_t N, Function;

	N = 1.0/(cte*par[2]);
	Function = N*(TMath::Exp(- (par[1]-x[0]) * (par[1] - x[0]) / (2.0 * par[2] * par[2]) ) );
	G = par[0]*Function;
 
	return G;

}

void Measurement_3pts(TString const& filename1,TString const& filename2,TString const& filename3,double resbin=300, double xmin=0, double xmax=15000){
	
        //double res       = 300.; // ps
        //double hist_xmin =  0.;
        //double hist_xmax =  15000.;
        int n_bins = (xmax - xmin )/resbin;

        Float_t position[3]  = {0.597, 1.343, 1.994}; // meters
        Float_t eposition[3] = {0.001, 0.001, 0.001};

	TBranch *b_event;
   	TBranch *b_tdc;
        int event, tdc;

        int n_sigmas = 5;
        double mean, sigma, chi2, x_min, x_max;
        int ndf;
        double par0, err_par0, par1, err_par1, par2, err_par2;
        double par3, err_par3, par4, err_par4, par5, err_par5;
        double par6, err_par6, par7, err_par7, par8, err_par8;

        Float_t time[3];
        Float_t etime[3];

        TFile *file;
        TTree *tree;

        Long64_t nentries;

        //////// first point //////////////////

	file = new TFile(filename1,"read");
   	tree = (TTree*)file->Get("t1");
	tree->SetBranchAddress("event", &event, &b_event);
  	tree->SetBranchAddress("tdc", &tdc, &b_tdc);

	TH1F* hist1 = new TH1F("hist1", "Events X TDC", n_bins, xmin, xmax);
        hist1->Sumw2(true);

	nentries = tree->GetEntriesFast();
	for (Long64_t i=0; i < nentries; ++i) {
    	   tree->GetEntry(i);
    	   hist1->Fill(tdc);
        }

        mean  = hist1->GetMean();
        sigma = hist1->GetStdDev();
        x_min = mean - n_sigmas*sigma; 
        x_max = mean + n_sigmas*sigma;

        cout << " debug = " << x_min << " " << x_max << endl;
 
        //TF1* fit_function1 = new TF1("fit_function1","gaus",x_min,x_max);
        TF1* fit_function1 = new TF1("fit_function1",fun_tripleGauss,x_min,x_max,9);
        TF1* fun1_g0 = new TF1("fun1_g0",fun_Gauss,x_min,x_max,3);
        TF1* fun1_g1 = new TF1("fun1_g1",fun_Gauss,x_min,x_max,3);
        TF1* fun1_g2 = new TF1("fun1_g2",fun_Gauss,x_min,x_max,3);

	fit_function1->SetParameter(0,10000000);
	fit_function1->SetParameter(1,3000);
	fit_function1->SetParameter(2,1000);
	fit_function1->SetParameter(3,1000000);
	fit_function1->SetParameter(4,5000);
	fit_function1->SetParameter(5,500);
	fit_function1->SetParameter(6,10000);
	fit_function1->SetParameter(7,7000);
	fit_function1->SetParameter(8,10000);

	fit_function1->SetParLimits(0,0,1000000000);
	fit_function1->SetParLimits(1,0,10000);
	fit_function1->SetParLimits(2,0,1000000000);
	fit_function1->SetParLimits(3,0,1000000000);
	fit_function1->SetParLimits(4,0,10000);
	fit_function1->SetParLimits(5,0,1000000000);
	fit_function1->SetParLimits(6,0,1000000000);
	fit_function1->SetParLimits(7,0,10000000);
	fit_function1->SetParLimits(8,0,1000000000);

	fit_function1->FixParameter(6,0);

	hist1->Fit(fit_function1,"R","");
        
        chi2     = fit_function1->GetChisquare(); 
        ndf      = fit_function1->GetNDF();

        par0     = fit_function1->GetParameter(0);
        err_par0 = fit_function1->GetParError(0);
        par1     = fit_function1->GetParameter(1);
        err_par1 = fit_function1->GetParError(1);
        par2     = fit_function1->GetParameter(2);
        err_par2 = fit_function1->GetParError(2);

        par3     = fit_function1->GetParameter(3);
        err_par3 = fit_function1->GetParError(3);
        par4     = fit_function1->GetParameter(4);
        err_par4 = fit_function1->GetParError(4);
        par5     = fit_function1->GetParameter(5);
        err_par5 = fit_function1->GetParError(5);

        par6     = fit_function1->GetParameter(6);
        err_par6 = fit_function1->GetParError(6);
        par7     = fit_function1->GetParameter(7);
        err_par7 = fit_function1->GetParError(7);
        par8     = fit_function1->GetParameter(8);
        err_par8 = fit_function1->GetParError(8);

        fun1_g0->SetParameter(0,par0); 	
        fun1_g0->SetParameter(1,par1); 	
        fun1_g0->SetParameter(2,par2); 	
        fun1_g1->SetParameter(0,par3); 	
        fun1_g1->SetParameter(1,par4); 	
        fun1_g1->SetParameter(2,par5); 	
        fun1_g2->SetParameter(0,par6); 	
        fun1_g2->SetParameter(1,par7); 	
        fun1_g2->SetParameter(2,par8); 	

        time[0]  = par1;
        etime[0] = err_par1;

        //////// second point //////////////////

	file = new TFile(filename2,"read");
   	tree = (TTree*)file->Get("t1");
	tree->SetBranchAddress("event", &event, &b_event);
  	tree->SetBranchAddress("tdc", &tdc, &b_tdc);

	TH1F* hist2 = new TH1F("hist2", "Events X TDC", n_bins, xmin, xmax);
        hist2->Sumw2(true);

	nentries = tree->GetEntriesFast();
	for (Long64_t i=0; i < nentries; ++i) {
    	   tree->GetEntry(i);
    	   hist2->Fill(tdc);
        }

        mean  = hist2->GetMean();
        sigma = hist2->GetStdDev();
        x_min = mean - n_sigmas*sigma; 
        x_max = mean + n_sigmas*sigma;
 
        //TF1* fit_function2 = new TF1("fit_function2","gaus",x_min,x_max);

        TF1* fit_function2 = new TF1("fit_function2",fun_tripleGauss,x_min,x_max,9);
        TF1* fun2_g0 = new TF1("fun2_g0",fun_Gauss,x_min,x_max,3);
        TF1* fun2_g1 = new TF1("fun2_g1",fun_Gauss,x_min,x_max,3);
        TF1* fun2_g2 = new TF1("fun2_g2",fun_Gauss,x_min,x_max,3);

        /*
	fit_function2->SetParameter(0,10000000);
	fit_function2->SetParameter(1,3000);
	fit_function2->SetParameter(2,1000);
	fit_function2->SetParameter(3,1000000);
	fit_function2->SetParameter(4,5000);
	fit_function2->SetParameter(5,500);
	fit_function2->SetParameter(6,1000000);
	fit_function2->SetParameter(7,7000);
	fit_function2->SetParameter(8,10000);
        */ 

	fit_function2->SetParLimits(0,0,1000000000);
	fit_function2->SetParLimits(1,0,10000);
	fit_function2->SetParLimits(2,0,1000000000);
	fit_function2->SetParLimits(3,0,1000000000);
	fit_function2->SetParLimits(4,0,10000);
	fit_function2->SetParLimits(5,0,1000000000);
	fit_function2->SetParLimits(6,0,1000000000);
	fit_function2->SetParLimits(7,0,10000);
	fit_function2->SetParLimits(8,0,1000000000);

	fit_function2->FixParameter(6,0);

	hist2->Fit(fit_function2,"R","");
        
        chi2     = fit_function2->GetChisquare(); 
        ndf      = fit_function2->GetNDF();

        par0     = fit_function2->GetParameter(3);
        err_par0 = fit_function2->GetParError(3);
        par1     = fit_function2->GetParameter(4);
        err_par1 = fit_function2->GetParError(4);
        par2     = fit_function2->GetParameter(5);
        err_par2 = fit_function2->GetParError(5);

        par3     = fit_function2->GetParameter(0);
        err_par3 = fit_function2->GetParError(0);
        par4     = fit_function2->GetParameter(1);
        err_par4 = fit_function2->GetParError(1);
        par5     = fit_function2->GetParameter(2);
        err_par5 = fit_function2->GetParError(2);

        par6     = fit_function2->GetParameter(6);
        err_par6 = fit_function2->GetParError(6);
        par7     = fit_function2->GetParameter(7);
        err_par7 = fit_function2->GetParError(7);
        par8     = fit_function2->GetParameter(8);
        err_par8 = fit_function2->GetParError(8); 

        fun2_g0->SetParameter(0,par0); 	
        fun2_g0->SetParameter(1,par1); 	
        fun2_g0->SetParameter(2,par2); 	
        fun2_g1->SetParameter(0,par3); 	
        fun2_g1->SetParameter(1,par4); 	
        fun2_g1->SetParameter(2,par5); 	
        fun2_g2->SetParameter(0,par6); 	
        fun2_g2->SetParameter(1,par7); 	
        fun2_g2->SetParameter(2,par8); 		
 	
        time[1]  = par1;
        etime[1] = err_par1;

        //////// third point //////////////////

	file = new TFile(filename3,"read");
   	tree = (TTree*)file->Get("t1");
	tree->SetBranchAddress("event", &event, &b_event);
  	tree->SetBranchAddress("tdc", &tdc, &b_tdc);

	TH1F* hist3 = new TH1F("hist3", "Events X TDC", n_bins, xmin, xmax);
        hist3->Sumw2(true);

	nentries = tree->GetEntriesFast();
	for (Long64_t i=0; i < nentries; ++i) {
    	   tree->GetEntry(i);
    	   hist3->Fill(tdc);
        }

        mean  = hist3->GetMean();
        sigma = hist3->GetStdDev();
        x_min = mean - n_sigmas*sigma; 
        x_max = mean + n_sigmas*sigma;
 
        TF1* fit_function3 = new TF1("fit_function3",fun_tripleGauss,x_min,x_max,9);
        TF1* fun3_g0 = new TF1("fun3_g0",fun_Gauss,x_min,x_max,3);
        TF1* fun3_g1 = new TF1("fun3_g1",fun_Gauss,x_min,x_max,3);
        TF1* fun3_g2 = new TF1("fun3_g2",fun_Gauss,x_min,x_max,3);

        /*
	fit_function3->SetParameter(0,10000000);
	fit_function3->SetParameter(1,3000);
	fit_function3->SetParameter(2,1000);
	fit_function3->SetParameter(3,1000000);
	fit_function3->SetParameter(4,5000);
	fit_function3->SetParameter(5,500);
	fit_function3->SetParameter(6,10000000);
	fit_function3->SetParameter(7,7000);
	fit_function3->SetParameter(8,10000);
        */

	fit_function3->SetParLimits(0,0,1000000000);
	fit_function3->SetParLimits(1,0,10000);
	fit_function3->SetParLimits(2,0,1000000000);
	fit_function3->SetParLimits(3,0,1000000000);
	fit_function3->SetParLimits(4,0,10000);
	fit_function3->SetParLimits(5,0,1000000000);
	fit_function3->SetParLimits(6,0,1000000000);
	fit_function3->SetParLimits(7,0,10000000);
	fit_function3->SetParLimits(8,0,1000000000);

	fit_function3->FixParameter(6,0);

	hist3->Fit(fit_function3,"R","");
        
        chi2     = fit_function3->GetChisquare(); 
        ndf      = fit_function3->GetNDF();

        par0     = fit_function3->GetParameter(0);
        err_par0 = fit_function3->GetParError(0);
        par1     = fit_function3->GetParameter(1);
        err_par1 = fit_function3->GetParError(1);
        par2     = fit_function3->GetParameter(2);
        err_par2 = fit_function3->GetParError(2);

        par3     = fit_function3->GetParameter(3);
        err_par3 = fit_function3->GetParError(3);
        par4     = fit_function3->GetParameter(4);
        err_par4 = fit_function3->GetParError(4);
        par5     = fit_function3->GetParameter(5);
        err_par5 = fit_function3->GetParError(5);

        par6     = fit_function3->GetParameter(6);
        err_par6 = fit_function3->GetParError(6);
        par7     = fit_function3->GetParameter(7);
        err_par7 = fit_function3->GetParError(7);
        par8     = fit_function3->GetParameter(8);
        err_par8 = fit_function3->GetParError(8); 

        fun3_g0->SetParameter(0,par0); 	
        fun3_g0->SetParameter(1,par1); 	
        fun3_g0->SetParameter(2,par2); 	
        fun3_g1->SetParameter(0,par3); 	
        fun3_g1->SetParameter(1,par4); 	
        fun3_g1->SetParameter(2,par5); 	
        fun3_g2->SetParameter(0,par6); 	
        fun3_g2->SetParameter(1,par7); 	
        fun3_g2->SetParameter(2,par8); 		
 	
        time[2]  = par1;
        etime[2] = err_par1;

        /////////// Results ////////////////

	TCanvas *c0 = new TCanvas("c0","",100,9,1200,500);
   	//c0->SetFillColor(42);
   	//c0->SetGrid();
   	//c0->GetFrame()->SetFillColor(21);
   	//c0->GetFrame()->SetBorderSize(12);
        c0->Divide(3,1);
    
        c0->cd(1);
 	hist1->GetYaxis()->SetTitle("Events");
   	hist1->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist1->SetTitle("Events X TDC");	
   	hist1->Draw(); 
        fun1_g0->SetLineColor(4);
        fun1_g0->Draw("SAME");
        fun1_g1->SetLineColor(6);
        fun1_g1->Draw("SAME");
        fun1_g2->SetLineColor(3);
        fun1_g2->Draw("SAME");

        c0->cd(2);
 	hist2->GetYaxis()->SetTitle("Events");
   	hist2->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist2->SetTitle("Events X TDC");	
   	hist2->Draw(); 
        fun2_g0->SetLineColor(4);
        fun2_g0->Draw("SAME");
        fun2_g1->SetLineColor(6);
        fun2_g1->Draw("SAME");
        fun2_g2->SetLineColor(3);
        fun2_g2->Draw("SAME");

        c0->cd(3);
 	hist3->GetYaxis()->SetTitle("Events");
   	hist3->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist3->SetTitle("Events X TDC");	
   	hist3->Draw(); 
        fun3_g0->SetLineColor(4);
        fun3_g0->Draw("SAME");
        fun3_g1->SetLineColor(6);
        fun3_g1->Draw("SAME");
        fun3_g2->SetLineColor(3);
        fun3_g2->Draw("SAME");

	TCanvas *c1 = new TCanvas("c1","",100,9,700,700);
   	//c1->SetFillColor(42);
   	//c1->SetGrid();
   	//c1->GetFrame()->SetFillColor(21);
   	//c1->GetFrame()->SetBorderSize(12);
 
        // /*
        Int_t num = 3;
        TGraphErrors *gr = new TGraphErrors(num, position, time, eposition, etime);
        gr->SetTitle("muon speed");
 	gr->GetYaxis()->SetTitle("#Delta time (ps)");
   	gr->GetXaxis()->SetTitle("#Delta position (m)");
        gr->SetMarkerColor(4);
        gr->SetMarkerStyle(21);
        gr->Draw("AP");

        TF1* fit = new TF1("fit","pol1");
        gr->Fit(fit);
        c1->Update();

        chi2     = fit->GetChisquare(); 
        ndf      = fit->GetNDF();
        par0     = fit->GetParameter(0);
        err_par0 = fit->GetParError(0);
        par1     = fit->GetParameter(1);
        err_par1 = fit->GetParError(1);

        cout << " --------------- "<< endl;
        cout << " muon velocity (10^8 m/s) = " << 10000.0/par1 << " +- " << ErrDiv(10000,0,par1,err_par1) << endl;   

	//TCanvas *c2 = new TCanvas("c2","",100,9,700,700);
        //fit_functiont ->Draw("");

        // */      
}

Float_t ErrDiv(Float_t a, Float_t da, Float_t b, Float_t db)
{

  /// SIGX =SQRT( { D_A**2/A**2 + D_B**2/B**2)} * X**2)
  //Float_t Final = a/b;
  //Float_t dRes = sqrt(da*da+db*db);
  
  return ErrDiv = (1/b)*(1/b)*TMath::Sqrt((b*b*da*da) + (a*a*db*db));

}


