///// Fast simulation:
///// muons cosmics passing through set of detectors positioned one on the top of each other

///// author: André Massafferri 
///// date: april 2016

/////////////// configured parameters /////////////////////////////
///// all dimensions in cm !!!!!!!!!!!!

// TRandom3

//-------------------------------------------
/// Detectors of the setup
const int NumDetectors = 2;

// vertical distance from origin
Float_t z0[NumDetectors] = {0,140};
//Float_t z0[NumDetectors] = {0,163.1};
//Float_t z0[NumDetectors] = {0,199.4};

// type refer to the geometry  
Int_t type[NumDetectors] = {0,0};

// number of noise hits
Int_t Noise[NumDetectors] = {0,0};

//--------------------------------------------
const int NumTypesOfDetectors = 1;

Int_t InitNum = 11110012;
TRandom3 *r1=new TRandom3();
r1->SetSeed(InitNum);

/// geometries of detectors with rectangles: 
/// vertices in x: x0 + icellsx * deltax, icellsx = 0 up to ncellsx
/// vertices in y: y0 + icellsy * deltay, icellsy = 0 up to ncellsy
Double_t x0[NumTypesOfDetectors]     = {0.0};   
Double_t y0[NumTypesOfDetectors]     = {0.0};   
Int_t ncellsx[NumTypesOfDetectors]   = {1};
Int_t ncellsy[NumTypesOfDetectors]   = {1};
Double_t deltax[NumTypesOfDetectors] = {40};
Double_t deltay[NumTypesOfDetectors] = {40};

//-------------------------------------------- 
/// Muon production parameters
/// geometry: defining rectangle in z=0 and area with vertices (0.0,0.0) (xprod,yprod) 
//  model: ModelProd 
Int_t ModelProd = 0;      // 0 = cos^2(theta), 1 = cos(2*theta)
Double_t xprod  = 40;     // (reference origin x = 0)
Double_t yprod  = 40;     // (reference origin y = 0)

//-------------------------------------------- 
/// Timing parameters
Double_t c = 0.0299792458; // v=deltaL(cm)/deltat(ps)
Float_t spreadDist_InsideDet = 0.0; // in cm
Float_t nrefr = 1.2;

//--------------------------------------------
/// global variables 
Double_t pi     = 3.141592;
Int_t linewidth = 3; // line width of histogram
Char_t name1[20],name2[20],name3[20];
Int_t chx[NumDetectors],chy[NumDetectors]; //output channels 
Int_t OutAcceptance[NumDetectors];
Int_t InAcceptance[NumDetectors];
Float_t Acceptance[NumDetectors];
Double_t dist[NumDetectors];
Int_t ntracks;
Double_t x_gen,y_gen,theta_gen,phi_gen;

//--------------------------------------------
// histograms 
TH2F *hch[NumDetectors];
TH2F *hxy[NumDetectors];
TH1F *hxprod = new TH1F("hxprod","h x prod",50,0,xprod);
TH1F *hyprod = new TH1F("hyprod","h y prod",50,0,yprod);
TH1F *hphiprod = new TH1F("hphiprod","h phi prod",50,-pi/2,pi/2);
TH1F *hthetaprod = new TH1F("hthetaprod","h theta prod",50,0,pi/2);
TH1F *htime[NumDetectors];

//--------------------------------------------
// main function 
Int_t Run(Int_t nEvents = 100) 
{

  Initialization();

  for (Int_t i = 0; i < nEvents; i++) 
  {
    Int_t nmuons = 1; 

    for  (Int_t j = 0; j < nmuons; j++)
    {
        /// production of one muon given by the Model -> function of theta, the angle wrt vertical axis
        /// in a rectangle in z=0 varying uniformly 0->xprod and 0->yprod
        GenerateMuon(ModelProd);// output : x_gen, y_gen, theta_gen, phi_gen  

        hxprod->Fill(x_gen);
        hyprod->Fill(y_gen);
        hthetaprod->Fill(theta_gen);
        hphiprod->Fill(phi_gen); 
                        
        ntracks++;                 

        for (Int_t idet = 0; idet<NumDetectors; idet++)
        {
          /// propagation of the muon track in each detector
          /// identify the cells chx and chy in which the muon hit (100% efficiently, without noise)
          /// output: chx, chy (if no hit is identified -> chx or chy = 1000)
          MuonInDetector(idet,type[idet],chx[idet],chy[idet],dist[idet]);

          if ((chx[idet]==1000) || (chy[idet]==1000))
          {
            OutAcceptance[idet]++;
          } else {
            hch[idet]->Fill(chx[idet],chy[idet]);
            InAcceptance[idet]++;

            Float_t time_arrival = dist[idet]/c; //assuming vmu=c  
            //time_arrival = z0[idet]/c; 
            Float_t speedLight_InsideDet = c/nrefr;         
            Float_t t_conv = r1->Gaus(time_arrival,spreadDist_InsideDet/speedLight_InsideDet);
            htime[idet]->Fill(t_conv);
          }

          for (Int_t inoise=0; inoise<(Noise[idet]); inoise++)
          {
             NoiseInDetector(idet,type[idet],chx[idet],chy[idet]);
             hch[idet]->Fill(chx[idet],chy[idet]);
          }
        } 
    }
  }

  Results();
  return 0;
}

void Initialization(){
	
   if (ModelProd==0) cout << " Production model: (costheta)^2 " << endl;
   if (ModelProd==1) cout << " Production model: cos(2*theta)  " << endl;
 
   for (Int_t idet=0; idet<NumDetectors; idet++)
   {
 
    	cout << " DETECTOR " << idet << endl;

      	if (type[idet]==0) cout << " type: single " << endl;
      	if (type[idet]==1) cout << " type: creat " << endl;

      	cout << "     location x (nchannels =" << ncellsx[type[idet]] << ") = " << x0[type[idet]] << " " << x0[type[idet]] + ncellsx[type[idet]]*deltax[type[idet]] << endl;
      	cout << "     location y (nchannels =" << ncellsy[type[idet]] << ") = " << y0[type[idet]] << " " << y0[type[idet]] + ncellsy[type[idet]]*deltay[type[idet]] << endl;
      	cout << "     location z = " << z0[idet] << endl;

      	if ((ncellsx[type[idet]]*deltax[type[idet]]) > xprod) cout << "WARNING: [x] detector larger than production area " << endl; 
      	if ((ncellsy[type[idet]]*deltay[type[idet]]) > yprod) cout << "WARNING: [y] detector larger than production area " << endl; 

        cout << " -------- " << endl;   

        OutAcceptance[idet] = 0;
        InAcceptance[idet]  = 0;
        ntracks             = 0;

        sprintf(name1,"hch_%d",idet);
        sprintf(name2,"htime_%d",idet);
        //sprintf(name3,"line_%d",idet);

        hch[idet] = new TH2F(name1,name1,ncellsx[type[idet]],0,Double_t(ncellsx[type[idet]]),ncellsy[type[idet]],0,Double_t(ncellsy[type[idet]]));
        htime[idet] = new TH1F(name2,name2,500,0.1*z0[idet]/c,2*z0[idet]/c);
    }
}

void Results(){

   TCanvas *c1 = new TCanvas("c1","Muon Production details",200,10,700,500);
   c1->SetFillColor(0);
   c1->GetFrame()->SetFillColor(6);
   c1->GetFrame()->SetBorderSize(10);
   c1->GetFrame()->SetBorderMode(3);
   c1->Divide(2,2);

   hxprod->SetLineWidth(linewidth);
   hyprod->SetLineWidth(linewidth);
   hthetaprod->SetLineWidth(linewidth);
   hphiprod->SetLineWidth(linewidth);

   c1->cd(1);
   hxprod->SetMinimum(0); 
   hxprod->Draw();
   c1->cd(2);
   hyprod->SetMinimum(0); 
   hyprod->Draw();
   c1->cd(3); 
   hthetaprod->SetMinimum(0);
   hthetaprod->Draw();
   c1->cd(4); 
   hphiprod->SetMinimum(0);
   hphiprod->Draw();

   TCanvas *c2 = new TCanvas("c2","Detector details CH",200,10,700,500);
   c2->SetFillColor(0);
   c2->GetFrame()->SetFillColor(6);
   c2->GetFrame()->SetBorderSize(10);
   c2->GetFrame()->SetBorderMode(3);
   Float_t Temp = TMath::Sqrt(NumDetectors);
   Int_t num1 = Temp;
   Int_t num2 = Temp+1;
   c2->Divide(num1,num2);

   for (Int_t idet=0; idet<NumDetectors; idet++)
   {
       c2->cd(idet+1); 
       hch[idet]->SetLineWidth(linewidth);
       hch[idet]->SetMinimum(0);
       hch[idet]->Draw("LEGO");

       Acceptance[idet] = Float_t(InAcceptance[idet])/Float_t(ntracks);
       cout << "Eff geom (det=" << idet << ") over production plane = " << Acceptance[idet] << endl; 
   }
 
   cout << " ntracks produced = " << ntracks << endl;
  
    TCanvas *c3 = new TCanvas("c3","timing plots",200,10,700,500);
    c3->SetFillColor(0);
    c3->GetFrame()->SetFillColor(6);
    c3->GetFrame()->SetBorderSize(10);
    c3->GetFrame()->SetBorderMode(3);
    c3->Divide(num1,num2);

    for (Int_t idet=0; idet<NumDetectors; idet++)
    {
        c3->cd(1+idet);
        htime[idet]->SetLineWidth(linewidth);
        htime[idet]->Draw();
        //Float_t Correction_time = z0[idet]/c - htime[idet]->GetMean();
        cout << "Timing values = " << z0[idet]/c << " - " << htime[idet]->GetMean() << endl;
        cout << "Correction (cm) = " << (TMath::Abs(z0[idet]/c-htime[idet]->GetMean()))*c << endl;
    }
}

void MuonInDetector(Int_t detector, Int_t Type, Int_t& channel_x, Int_t& channel_y, Double_t& distance){
        
      /// linear propagation from area of production to the detector       
      Double_t radii = tan(theta_gen)*z0[detector];
      Double_t x     = radii*cos(phi_gen) + x_gen;
      Double_t y     = radii*sin(phi_gen) + y_gen;
    
      distance = TMath::Sqrt(z0[detector]*z0[detector] + radii*radii); // in cm

      channel_x = SearchChannelX(Type,x);
      channel_y = SearchChannelY(Type,y);

}

void NoiseInDetector(Int_t detector, Int_t Type, Int_t& channel_x, Int_t& channel_y){ 
        
      gRandom->SetSeed();
 
      /// generating noise uniformly in the rectangle and z0
      Double_t x = (x0[Type]+ncellsx[Type]*deltax[Type])*(RooRandom::uniform());
      Double_t y = (y0[Type]+ncellsy[Type]*deltay[Type])*(RooRandom::uniform());

      channel_x = SearchChannelX(Type,x);
      channel_y = SearchChannelY(Type,y);

}

Int_t SearchChannelX(Int_t Type, Double_t x){

      // identifying the cell in x
      Int_t icell = 0;
      while (1)
      {
           if ( ( (x > (x0[Type] + icell*deltax[Type]) ) && (x <= (x0[Type] + (icell+1)*deltax[Type]) ) ) || (icell==ncellsx[Type]) )  break;
          icell++;
      }
      if (icell < ncellsx[Type])
      {
             return = icell;
      } else 
      {
             return = 1000; 
      }
}

Int_t SearchChannelY(Int_t Type, Double_t y){

      /// identifying the cell in y
      Int_t icell = 0;

      while (1)
      {
           if ( ( (y > (y0[Type] + icell*deltay[Type]) ) && (y <= (y0[Type] + (icell+1)*deltay[Type]) ) ) || (icell==ncellsy[Type]) )  break;
           icell++;
      }
      if (icell < ncellsy[Type])
      {
             return = icell;
      } else 
      {
             return = 1000; 
      }
}

void GenerateMuon(Int_t model){
                   
      gRandom->SetSeed(InitNum);
 
      /// generating muons uniformly in the rectangle
      /// with phi angle also uniform
        x_gen   = xprod*(RooRandom::uniform());
    	y_gen   = yprod*(RooRandom::uniform());
    	phi_gen = pi*(2.0*(RooRandom::uniform())-1.);        //// -pi to pi
 
    	Double_t r    = 1.;		
    	Double_t prod = 0.;

        /// generating the angle wrt the normal according the model
        /// 0->cos^2(theta), 1->cos(2*theta)
    	while (r>prod)
    	{
       		r = 1.1*RooRandom::uniform();
       		theta_gen = (pi/2)*(RooRandom::uniform());  //// 0 to pi/2
                if (model==0) prod = pow(cos(theta_gen),2.);  
      		if (model==1) prod = cos(2.0*theta_gen);                            
    	}

}

Double_t funModifiedGauss(Double_t *x, Double_t *param)
{

  Double_t G;
  Double_t func_mod, func;
  func_mod = TMath::Exp(param[3]*(x[0]-param[1]) );
  func = (TMath::Exp(- func_mod * func_mod * (param[1]-x[0]) * (param[1] - x[0]) / (2.0 * param[2] * param[2]) ) );
  G = param[0]*func;

  return G;

}
