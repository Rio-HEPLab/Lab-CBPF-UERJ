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

        int n_sigmas = 3;
        double mean, sigma, chi2, x_min, x_max;
        int ndf;
        double par0, err_par0, par1, err_par1, par2, err_par2;

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
 
        TF1* fit_function1 = new TF1("fit_function1","gaus",x_min,x_max);
	hist1->Fit(fit_function1,"R","");
        
        chi2     = fit_function1->GetChisquare(); 
        ndf      = fit_function1->GetNDF();
        par0     = fit_function1->GetParameter(0);
        err_par0 = fit_function1->GetParError(0);
        par1     = fit_function1->GetParameter(1);
        err_par1 = fit_function1->GetParError(1);
        par2     = fit_function1->GetParameter(2);
        err_par2 = fit_function1->GetParError(2);
 	
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
 
        TF1* fit_function2 = new TF1("fit_function2","gaus",x_min,x_max);
	hist2->Fit(fit_function2,"R","");
        
        chi2     = fit_function2->GetChisquare(); 
        ndf      = fit_function2->GetNDF();
        par0     = fit_function2->GetParameter(0);
        err_par0 = fit_function2->GetParError(0);
        par1     = fit_function2->GetParameter(1);
        err_par1 = fit_function2->GetParError(1);
        par2     = fit_function2->GetParameter(2);
        err_par2 = fit_function2->GetParError(2);
 	
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
 
        TF1* fit_function3 = new TF1("fit_function3","gaus",x_min,x_max);
	hist3->Fit(fit_function3,"R","");
        
        chi2     = fit_function3->GetChisquare(); 
        ndf      = fit_function3->GetNDF();
        par0     = fit_function3->GetParameter(0);
        err_par0 = fit_function3->GetParError(0);
        par1     = fit_function3->GetParameter(1);
        err_par1 = fit_function3->GetParError(1);
        par2     = fit_function3->GetParameter(2);
        err_par2 = fit_function3->GetParError(2);
 	
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

        c0->cd(2);
 	hist2->GetYaxis()->SetTitle("Events");
   	hist2->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist2->SetTitle("Events X TDC");	
   	hist2->Draw(); 

        c0->cd(3);
 	hist3->GetYaxis()->SetTitle("Events");
   	hist3->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist3->SetTitle("Events X TDC");	
   	hist3->Draw(); 

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
        cout << " muon velocity (10^8 m/s) = " << 10000.0/par1 << " +- " << 0.0 << endl;   

        // */

         
}

