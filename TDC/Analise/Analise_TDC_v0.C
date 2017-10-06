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

void MakeTree(TString const& input, TString const& output, TString const& output_root = "TreeTDC.root"){
	
	FILE *data;
	
        data = fopen(input,"rt");
	if(data == NULL){
		printf("Arquivo input não encontrado");
		exit(0);
	}
	
	char urlNull[]="Null.txt";

	FILE *arqNull;
	arqNull = fopen(urlNull,"w");
	if(arqNull == NULL){
		printf("Não foi possivel abrir o arquivo NULL\n");
		exit(0);
	}	

	FILE *arq;
        arq = fopen(output,"w"); 
	if(arq == NULL){
		printf("Não foi possivel abrir o arquivo output\n");
		exit(0);
	}
	
	TFile *file = new TFile(output_root,"recreate");
   	TTree *t1 = new TTree("t1","Tree com dados TDC");
	
	bool use_all_channels = true;
		
	int event, prev_event, eventN, channel0, channel1, channel2, tdc0, tdc1, tdc2, tdc;
        eventN = 0; 
        const int NSEARCHMAX = 100;
        event = -1;
        prev_event = -1;
        bool valid = true;
        fpos_t pos;

	t1->Branch("event",&event,"event/I");
   	t1->Branch("tdc",&tdc,"tdc/I");
	
	while( !feof(data) ){

                int ret = -1;
                
                if ( !valid ) {
                   cout << "Find correct start line." << endl;
                   for(int i=0; i < NSEARCHMAX; ++i) {
                      cout << "  Iteration " << (i+1) << endl;
                      cout << "  Prev. event " << prev_event << endl;
                      ret = fscanf(data,"%d", &event);
                      if ( event == (prev_event + 1) ) { valid = true; break; }
                   }
                   if ( !valid ) continue;
                } else{
                   prev_event = event;
                   fgetpos( data, &pos ); 
   		   ret = fscanf(data,"%d", &event);
                   if( feof(data) ) break;

                   if( event == 0 ) prev_event = 0;
                   else if ( event != (prev_event + 1) ) { valid = false; fsetpos( data, &pos ); continue; }
                }

                fgetpos( data, &pos );
		ret = fscanf(data,"%d %d", &channel0, &tdc0);
                if ( channel0 != 1 ) { valid = false; fsetpos( data, &pos ); prev_event = event; continue; }

                fgetpos( data, &pos ); 
		ret = fscanf(data,"%d %d", &channel1, &tdc1);
                if ( channel1 != 2 && channel1 != 1 ) { valid = false; fsetpos( data, &pos ); prev_event = event; continue; }

                if( use_all_channels ){
                   fgetpos( data, &pos ); 
		   ret = fscanf(data,"%d %d", &channel2, &tdc2);
                   if ( channel2 != 0 && channel2 != 2 ) { valid = false; fsetpos( data, &pos ); prev_event = event; continue; }
                }

		//tdc = ( use_all_channels ) ? abs( (tdc2-tdc1)*25 ) : (tdc1-tdc0)*25;
		//tdc = ( use_all_channels ) ? abs( (tdc2-tdc0)*25 ) : (tdc1-tdc0)*25; //conf channels		
		tdc = ( use_all_channels ) ? abs( (tdc1-tdc0)*25 ) : (tdc1-tdc0)*25; //conf channels		

		if( tdc0 == 0 || tdc1 == 0 || tdc2 == 0 || tdc < 0 ){
			(eventN++);
			fprintf(arqNull, "Eventos Nulo: %d Event: %d TDC0: %d TDC1: %d\n",eventN,event,tdc0,tdc1);
			 }
		 else{
			fprintf(arq,"Evento: %d TDC: %d \n", event, tdc);
			t1->Fill();
	        }

        }
	t1->Write();
   	t1->Print();
   	file->Close();

        fclose(data);
	fclose(arq);
	fclose(arqNull);
	return ;

}

void ReadTree(TString const& filename){
	

	TFile *file = new TFile(filename,"read");
   	TTree* t1 = (TTree*)file->Get("t1");

	TBranch *b_event;
   	TBranch *b_tdc;

        int event, tdc;
	t1->SetBranchAddress("event", &event, &b_event);
  	t1->SetBranchAddress("tdc", &tdc, &b_tdc);

        double res = 300.; // ps
        double hist_xmin =  0.;
        double hist_xmax =  50000.;
        int n_bins = ( hist_xmax - hist_xmin )/res;
	TH1F* hist = new TH1F("Histogram", "Events X TDC", n_bins, hist_xmin, hist_xmax);
        hist->Sumw2(true);

	Long64_t nentries = t1->GetEntriesFast();

	for (Long64_t i=0; i < nentries; ++i) {
    	   t1->GetEntry(i);
    	   hist->Fill(tdc);
        }

        int n_sigmas = 3;
        double mean = hist->GetMean();
        double sigma = hist->GetStdDev();
        double x_min = mean - n_sigmas*sigma; 
        double x_max = mean + n_sigmas*sigma;
 

        TF1* fit_function0 = new TF1("fit_function","gaus",x_min,x_max);
	hist->Fit(fit_function0,"R","");
        
	
	//TF1* fit_function1 = new TF1("fit_function","landau",118000,126000);  
        //hist->Fit(fit_function1,"R","");

        double chi2 = fit_function0->GetChisquare(); 
        int ndf = fit_function0->GetNDF();
        double par0 = fit_function0->GetParameter(0);
        double err_par0 = fit_function0->GetParError(0);
        double par1 = fit_function0->GetParameter(1);
        double err_par1 = fit_function0->GetParError(1);
        double par2 = fit_function0->GetParameter(2);
        double err_par2 = fit_function0->GetParError(2);
 	
	/*double chi2 = fit_function1->GetChisquare(); 
        int ndf = fit_function1->GetNDF();
        double par0 = fit_function1->GetParameter(0);
        double err_par0 = fit_function1->GetParError(0);
        double par1 = fit_function1->GetParameter(1);
        double err_par1 = fit_function1->GetParError(1);
        double par2 = fit_function1->GetParameter(2);
        double err_par2 = fit_function1->GetParError(2);*/
 	

        //cout << "par 0 = " << std::setprecision(2) << par0 << " +/- " << err_par0 << endl;
        //cout << "par 1 = " << std::setprecision(2) << par1 << " +/- " << err_par1 << endl;
        //cout << "par 2 = " << std::setprecision(2) << par2 << " +/- " << err_par2 << endl;
        //cout << "chi2/ndf = " << std::setprecision(2) << chi2/ndf << endl;

 	hist->GetYaxis()->SetTitle("Events");
   	hist->GetXaxis()->SetTitle("#Delta t (ps)");
   	hist->SetTitle("Events X TDC");	
   	hist->Draw(); 
         
}

