#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TSpectrum.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::endl;

TH1F* analysis(TString const& filename, TString const& label){
	
	TFile *file = new TFile(filename,"read");
   	TTree* t1 = (TTree*)file->Get("t1");

	TBranch *b_event;
   	//TBranch *b_counts;
   	TBranch *b_charge;

        int event, counts;
        double charge;
	t1->SetBranchAddress("event", &event, &b_event);
  	//t1->SetBranchAddress("counts", &counts, &b_counts);
  	t1->SetBranchAddress("charge", &charge, &b_charge);

        double res = 0.500; // pC
        double hist_xmin =  0.;
        double hist_xmax =  400.; // pC
        int n_bins = ( hist_xmax - hist_xmin )/res;
        TString hist_name = "histQDC"; hist_name += "_"; hist_name += label;
	TH1F* hist = new TH1F(hist_name, "", n_bins, hist_xmin, hist_xmax);
        hist->Sumw2(true);

	Long64_t nentries = t1->GetEntriesFast();

        int n_events = 0;
	for (Long64_t i=0; i < nentries; ++i) {
    	   t1->GetEntry(i);
    	   hist->Fill(charge);
           ++n_events;
        }

        hist->Scale( 1./n_events );

        /*int n_sigmas = 3;
        double mean = hist->GetMean();
        double sigma = hist->GetStdDev();
        double x_min = mean - n_sigmas*sigma; 
        double x_max = mean + n_sigmas*sigma;
 

        TF1* fit_function0 = new TF1("fit_function","gaus",x_min,x_max);
	hist->Fit(fit_function0,"R","");
        
        double chi2 = fit_function0->GetChisquare(); 
        int ndf = fit_function0->GetNDF();
        double par0 = fit_function0->GetParameter(0);
        double err_par0 = fit_function0->GetParError(0);
        double par1 = fit_function0->GetParameter(1);
        double err_par1 = fit_function0->GetParError(1);
        double par2 = fit_function0->GetParameter(2);
        double err_par2 = fit_function0->GetParError(2);*/

        TString canvas_name = "canvas"; canvas_name += "_"; canvas_name += label;
        TCanvas* canvas = new TCanvas(canvas_name);
 	hist->GetYaxis()->SetTitle("Normalized");
   	hist->GetXaxis()->SetTitle("charge (pC)");
   	hist->Draw(); 

        TSpectrum spectrum;
        spectrum.Search(hist,10.,"",0.20);

        Int_t n_peaks = spectrum.GetNPeaks();  
        Float_t* peaks = spectrum.GetPositionX();
        for(int idx = 0; idx < n_peaks; ++idx) cout << "Peak " << idx << ": " << peaks[idx] << endl;       

	//TF1* fit_function0 = new TF1("fit_function","gaus",62,72);
	//hist->Fit(fit_function0,"R","");

	return hist;
}

void analysis(){

   TString filename_pedestal = "analysis_QDC_pedestal_241016-01.root";
   TString filename_data_1 = "analysis_QDC_241016-01.root";
   TString filename_data_2 = "analysis_QDC_241016-02.root";

   TH1F* hist_pedestal = analysis(filename_pedestal,"pedestal");
   TH1F* hist_data_1 = analysis(filename_data_1,"data_1");
   TH1F* hist_data_2 = analysis(filename_data_2,"data_2");

   TCanvas* canvas = new TCanvas("canvas_all");
   canvas->cd();

   hist_pedestal->SetLineWidth(3);
   hist_pedestal->SetLineStyle(2);
   hist_pedestal->SetLineColor(kRed);
   hist_data_1->SetLineWidth(3);
   hist_data_1->SetLineStyle(1);
   hist_data_1->SetLineColor(kBlack);
   hist_data_2->SetLineWidth(3);
   hist_data_2->SetLineStyle(3);
   hist_data_2->SetLineColor(56);
  
   hist_data_1->Draw("HISTE1");
   hist_data_2->Draw("HISTE1SAME");
   hist_pedestal->Draw("HISTE1SAME");
} 

