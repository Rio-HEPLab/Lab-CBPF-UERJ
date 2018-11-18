#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TSpectrum.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void fit_Gauss(TString const& filename, TString const& label){

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

   double res = 1.000; // pC
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
   hist->GetYaxis()->SetTitle("Normalized");
   hist->GetXaxis()->SetTitle("charge (pC)");
   hist->SetMarkerStyle( 20 );
   hist->SetMarkerSize( 0.8 );

   // Modify fit ranges
   TF1* fit_function0 = new TF1("fit_function0","gaus",0,50);	
   TF1* fit_function1 = new TF1("fit_function1","gaus",50,115);
   TF1* fit_function2 = new TF1("fit_function2","gaus",115,175);       	
   TF1* fit_function3 = new TF1("fit_function3","gaus",175,240);
   TF1* fit_function4 = new TF1("fit_function4","gaus",240,300);
   fit_function0->SetLineColor(2);
   fit_function1->SetLineColor(2);
   fit_function2->SetLineColor(2);
   fit_function3->SetLineColor(2);
   fit_function4->SetLineColor(2);
   fit_function0->SetLineStyle(2);
   fit_function1->SetLineStyle(2);
   fit_function2->SetLineStyle(2);
   fit_function3->SetLineStyle(2);
   fit_function4->SetLineStyle(2);

   TF1* total = new TF1("Multi Gauss","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+pol2(15)",0,300);    
   total->SetLineColor(2);
   total->SetLineWidth(2);

   TString canvas_name = "canvas"; canvas_name += "_"; canvas_name += label;
   TCanvas* canvas = new TCanvas(canvas_name);

   hist->Draw("E1"); 

   hist->Fit(fit_function0,"R","");
   hist->Fit(fit_function1,"R+","");
   hist->Fit(fit_function2,"R+","");
   hist->Fit(fit_function3,"R+","");
   hist->Fit(fit_function4,"R+","");

   Double_t par[18];
   for(size_t i = 0; i < 18; ++i) par[i] = 0.; 
   fit_function0->GetParameters( &par[0]  );
   fit_function1->GetParameters( &par[3]  );
   fit_function2->GetParameters( &par[6]  );
   fit_function3->GetParameters( &par[9]  );
   fit_function4->GetParameters( &par[12] );
   //par[15] = 0.;
   //par[16] = 0.00001;
   //par[17] = 0.;

   total->SetParameters(par);
   hist->Fit(total,"R+"); 

   total->GetParameters( par );
   Double_t const* par_errors = total->GetParErrors();

   TF1* fit_function_pol = new TF1("fit_function_pol","pol2",0,300);
   fit_function_pol->SetParameters( &par[15] );
   fit_function_pol->SetLineColor(56);
   fit_function_pol->SetLineStyle(2);
   fit_function_pol->Draw("SAME");

   Double_t mean[5];
   mean[0] = par[1]; 
   mean[1] = par[4]; 
   mean[2] = par[7]; 
   mean[3] = par[10]; 
   mean[4] = par[13];
   Double_t error_mean[5];
   error_mean[0] = par_errors[1];
   error_mean[1] = par_errors[4];
   error_mean[2] = par_errors[7];
   error_mean[3] = par_errors[10];
   error_mean[4] = par_errors[13];
   for(size_t i = 0; i < 5; ++i){
      std::cout << "mean " << i << ": " << mean[i] << " +/- " << error_mean[i] << std::endl;
      if( i > 0 ){
         Double_t gain_from_previous = mean[i] - mean[i-1];
         Double_t err_gain_from_previous = TMath::Sqrt( error_mean[i]*error_mean[i] + error_mean[i-1]*error_mean[i-1] );
         std::cout << "  gain wrt previous: " << gain_from_previous << " +/- " << err_gain_from_previous << std::endl;
      }
   }
    
}
