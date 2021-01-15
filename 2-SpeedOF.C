#include "OF.h"

//--------------------------------------------------------------------------------------------- VIEW IN FULL SCREEN FOR BEST EXPERIENCE!-----------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------COMPLEX TO REAL FFT (INVERSE FOURIER TRANSFORM)-------------------------------------------
void ComplexToRealIFFT(const vector<TComplex>& inComp, vector<double>& outRe)
{
  int n = inComp.size();
  double *re_pvector = new double[inComp.size()];
  double *im_pvector = new double[inComp.size()];

  double re,im;
  int i;

  //copying std vector back into arrays to pass to IFFT
  for (i=0; i < n; i++){
    re_pvector[i] = inComp[i].Re();
    im_pvector[i] = inComp[i].Im();
  }

  TVirtualFFT *ifftc2r = TVirtualFFT::FFT(1,&n,"C2R ES"); //this should be the opposite of FFT
  ifftc2r->SetPointsComplex(re_pvector, im_pvector);
  ifftc2r->Transform();
  ifftc2r->GetPointComplex(0,re,im);
  outRe.push_back(re/sqrt((double)n)); //DC component
  
  for(i=1; i<n ;i++)
    {
      ifftc2r->GetPointComplex(i,re, im);
      outRe.push_back(re/sqrt((double)n));
    }
  
  delete[] re_pvector;
  delete[] im_pvector;
  delete ifftc2r;

  return;
}
//------------------------------------------------------------------------------------------------------------------------REAL TO COMPLEX FFT(FOURIER TRANSFORM)------------------------------------------
void RealToComplexFFT(const vector<double>& pulsevector, vector<TComplex>& outComp)
{
  int n = pulsevector.size();
  
  double *pvector = new double[pulsevector.size()];
  double re,im;
  int i;

  //copying std vector back into array to pass to FFT (to mirror what is done in PulseTools)
  for (i=0;i<(int)pulsevector.size();i++){
    pvector[i] = pulsevector[i];
  }
  
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1,&n,"R2C ES");
  fftr2c->SetPoints(pvector);
  fftr2c->Transform();
  fftr2c->GetPointComplex(0,re,im);
  TComplex tempCopy(re/sqrt((double)n), im/sqrt((double)n));
  outComp.push_back(tempCopy);
  for(i=1; i<n ;i++)
    {
      fftr2c->GetPointComplex(i,re,im);
      TComplex tempCopy(re/sqrt((double)n), im/sqrt((double)n));
      outComp.push_back(tempCopy);
    }
  
   //done! so delete new'd objects
  delete[] pvector;
  delete fftr2c;
  
  return;

}

//---------------------------------------------------------------------------------------------------------------------GENERATE TEMPLATE------------------------------------------------------------------
void generate_template(vector<double> &template_vec_32768, double dt){
  //now make the template
  double time_offset_us = 0.0;
  double A_f = 0.9;
  double T_f = 100;
  double R_f = 25;
  double F_f = 50;
  double A_s = 0.4;
  double T_s = 100;
  double R_s = 100.;
  double F_s = 3;
  TF1 *myPhononPulse = new TF1("myPhononPulse", "[0]*(1-TMath::Exp(-(x-[1])/[2]))*TMath::Exp(-(x-[1])/[3])+[4]*(1-TMath::Exp(-(x-[5])/[6]))*TMath::Exp(-(x-[5])/[7])", 0, 16552);
  myPhononPulse->SetParameter(0,A_f);
  myPhononPulse->SetParameter(1,T_f);
  myPhononPulse->SetParameter(2,R_f);
  myPhononPulse->SetParameter(3,F_f);
  myPhononPulse->SetParameter(4,A_s);
  myPhononPulse->SetParameter(5,T_s);
  myPhononPulse->SetParameter(6,R_s);
  myPhononPulse->SetParameter(7,F_s);

  for(int i = 0; i < 16384; i++){
    template_vec_32768.push_back(0);
  }

  for(int i = 0; i < 16384; i++){
    if (i*1.6+dt  < T_f) template_vec_32768.push_back(0);
    else template_vec_32768.push_back(myPhononPulse->Eval(i*1.6+dt));
  }
  delete myPhononPulse;
}
//-------------------------------------------------------------------------------------------------------------------SHIFT_TEMPLATE (SHIFT A TEMPLATE)---------------------------------------------------
void shift_template(vector<double> &template_vec_32768, vector<double> &shifted_template, double delay){
  if(delay < 0){
    for(int i=0; i<32768; i++){
      if (i<16384+delay){
	shifted_template.push_back(template_vec_32768[i]);
      }
      if (i >= 16384+delay && i<32768+delay){
	shifted_template.push_back(template_vec_32768[i-delay]);
      } 
      if(i >= 32768+delay){
	shifted_template.push_back(0);
      }
    }
  }
  else{
    for(int i=0; i<32768; i++){
      if (i<16384+delay){
	shifted_template.push_back(0);
      }
      if (i >= 16384+delay && i<32768+delay){
	shifted_template.push_back(template_vec_32768[i-delay]);
      } 
    }
  }
}

//-------------------------------------------------------------------------------------------------------------------GET NOISE TIME(GENERATES NOISE)------------------------------------------------------
void getNoise_time( double *noise32768, double *noisePSD){

  TRandom3 gRandom2(0);

  double freqres = 19.07348633;
  for(int i = 0; i < 32768; i++){
    noise32768[i] = 0.0;
  }

  vector<TComplex> noiseFFT;
  for(int i = 0; i < 16384;i++){
    if(i==0){
      //baseline
      double A = 0.;
      double B = 0.;
      TComplex complex(A,B);
      noiseFFT.push_back(complex);

    }
    else if(i<16384){
      double A = gRandom2.Gaus(0,sqrt(noisePSD[i]*625/2.));
      double B = gRandom2.Gaus(0,sqrt(noisePSD[i]*625/2.));
      TComplex complex(A,B);
      noiseFFT.push_back(complex);
    }
  }
  for(int i = 0; i < 16384-1; i++){
    TComplex complex(noiseFFT[16383-i].Re(),-noiseFFT[16383-i].Im());
    noiseFFT.push_back(complex);
  }
  vector<double> noise_time;
  ComplexToRealIFFT(noiseFFT, noise_time);
  for(int i=0; i<32768; i++){
    noise32768[i] = noise_time[i];
  }
  noiseFFT.clear();
  noise_time.clear();
}

//---------------------------------------------------------------------------------------------------------------COMPRESS TIME (HYBRIDIZES THE INPUT ARRAY)-----------------------------------------------
void compressTime(double *time32768, double *time3008){
  double sum;
  double average;
  int arg;
  for(int i = 0; i < 1024; i++){
    sum = 0.0;
    average = 0.0;
    for(int j = 0; j < 16; j++){
      arg = i*16+j;
      sum += time32768[arg];
    }
    average = (sum)/16.0;
    time3008[i] = average;
  }
  for(int i = 0; i < 1024; i++){
    sum = 0.0;
    //16*1024=16384
    arg = i+16384;
    sum = time32768[arg];
    time3008[i+1024] = sum;
  }
  for(int i = 0; i < 960; i++){
    sum = 0;
    average = 0;
    for(int j = 0; j < 16; j++){
      arg = (i*16+j)+(1024*16+1024);
      sum += time32768[arg];
    }
    average = (sum)/16.0;
    time3008[i+1024+1024] = average;
  }
}

//-----------------------------------------------------------------------------------------------------------FAST TO SLOW COMPRESS (DOWNSAMPLES FAST SAMPLED REGION)--------------------------------------
void FastToSlowCompress(double *signal2048, double *signal3008){
  int i=1024, b=1040, j=0;
  double sum;
  double downsampled[64];
  while (i < 2048){
    sum=0;
    while (i<b){
	sum += signal3008[i];
	i++;
      }
    i=b;
    b=b+16;
    downsampled[j]=sum/16;
    j++;
  }
  
  i=0;
  while (i<1024){
    signal2048[i]=signal3008[i];
    i++;
  }
  while (i>=1024 && i<1088){
    signal2048[i]=downsampled[i-1024];
    i++;
  }
  while (i<2048){
    signal2048[i]=signal3008[i+960];
    i++;
  }
}

//-------------------------------------------------------------------------------------------------------DOWNSAMPLE_TEMPLATE (DOWNSAMPLE SHIFTED TEMPLATE)------------------------------------------------
void downsample_template(double *template32768, double *template2048){
  int i=0, b=16, j=0;
  while( i<32768){
    double sum=0;
    while(i<b){
      sum += template32768[i];
      i++;
    }
    i=b;
    b+=16;
    template2048[j]=sum/16;
    j++;
  }
}

//-------------------------------------------------------------------------------------------------------ARRAY TO VECTOR (CONVERTS AN ARRAY TO A VECTOR)--------------------------------------------------
void ArrayToVector(double *array, vector<double> &vector, int length){

  for(int i=0; i<length; i++){
    vector.push_back(array[i]);
  }
  
}

//----------------------------------------------------------------------------------------------------TIME2PSD (CONVERTS TIME DATA TO PSD)---------------------------------------------------------------
void Time2PSD(const vector<double>& pulseVector, const double& fs, vector<double>& outPSD)
{
   int nBins = pulseVector.size();
   if(nBins == 0)
   {
      cerr <<"PulseTools::Time2PSD - ERROR! empty pulse passed into this function" << endl;
      exit(1);
   }
   int nBinsPSD=nBins/2+1;
   
   // Do FFT
   vector<TComplex> outCFFT;
   RealToComplexFFT(pulseVector,outCFFT);
   int i;


 
   // We  want the 2-sided distribution in pulseNorm/sqrt(Hz), exploiting symmtery for +/- frequency components
   // being carefull that 
   //    1. DC compoment should not be doubled for +/- frequency
   //    2. Max frequency (Nyquist) should not be doubled for +/- frequency (case trace length even)


   // First let's calculate the squared norm
   
   // DC component...
   outPSD.push_back(outCFFT[0].Re()*outCFFT[0].Re()); // DC component
   
   // .. and the rest
   for (i=1;i<nBinsPSD;i++)
     outPSD.push_back(outCFFT[i].Re()*outCFFT[i].Re() 
		      + outCFFT[nBins-i].Im()*outCFFT[nBins-i].Im());
   
   
   // Then normalize properly in pulseNorm/sqrt(Hz)
   for(uint i=0; i<outPSD.size(); i++)
     {
       // DC AND  Nyquist component for even length case-> no factor 2
       if (i==0 || (i==outPSD.size()-1) ) {
	 outPSD[i] = sqrt(outPSD[i]/fs);
       }else
	 outPSD[i] = sqrt(2.0*outPSD[i]/fs);
     }
   outCFFT.clear();
   return;
}

//----------------------------------------------------------------------------------------------------------PSD2FFT (CONVERTS PSD DATA TO FFT)------------------------------------------------------------
void PSD2FFT(const vector<double>& psdvector, vector<double>& outComp)
{
   int nBins = psdvector.size();
  
   if(nBins == 0)
   {
      cerr <<"PulseTools::PSD2FFT - ERROR! empty pulse passed into this function" << endl;
      exit(1);
   }

   
   int i;
   vector<double> outTemp;
   TComplex atemp;

       
   // Build a temporary vector  (similar than PSD2FFT with phase information)
   // Normalize with sqrt(2.0) but this will have to be taken out for DC and 
   // Nyquist component (case even trace length)
   for(i=0;i<nBins;i++){
     atemp = psdvector[i]/sqrt(2.0); // complex number not really required but using same code...
     outTemp.push_back(atemp);
   }
   
   // Now Build output vector
   
   //even number of elements in original pulse
       
   // 1. DC component  (removing sqrt(2.0))
   outComp.push_back(sqrt(2.0)*outTemp[0]);
   
   // 2. Rest of First half without Nyquist
   for(i=1;i<nBins-1;i++)
     outComp.push_back(outTemp[i]);
   
   // 3. Nyquist (removing sqrt(2.0))
   outComp.push_back(sqrt(2.0)*outTemp[nBins-1]); 
   
   // 4. second half 
   for(i=nBins-2;i>0;i--){
     outComp.push_back(outTemp[i]);
   }
   outTemp.clear();
   return;
   
}
//-----------------------------------------------------------------------------------------------------------DO OF (FUNCTION THAT PERFORMS THE OF)--------------------------------------------------------
double DoOptimalFilter(double *signal, double *templ,const vector<double>& noise, int length, double df, int fwindow1, int fwindow2, double dt, double *result, double f){
  vector<TComplex> pProd;
  vector<double> signalVector;
  vector<double> templateVector;
  vector<TComplex> signalFFT;
  vector<TComplex> templateFFT;
  //double sqrtdT=dt;
  double normalization=0;
  double amp=0;
  ArrayToVector(signal,signalVector, length);
  ArrayToVector(templ,templateVector,length);

  //---------------------FFT Calculations----------------------------------- 
  RealToComplexFFT(signalVector,signalFFT);
  RealToComplexFFT(templateVector,templateFFT);

  //----------------------------------OF---------------------------
  int nBins=length;
  double temp=0;
  for(int binItr=0; binItr < nBins; binItr++) {
    //signalFFT[binItr] *= sqrt(1.6*1e-6);   //This line was present in the cdmsbats code, but the number made no sense to Emanuele or me, so we avoided it
    pProd.push_back(signalFFT[binItr]*(TComplex::Conjugate(templateFFT[binItr]))/(noise[binItr]*noise[binItr]));
 
   if(binItr != 0){
     amp += pProd[binItr].Re();
     temp += pow(TComplex::Abs(templateFFT[binItr]),2)/pow(TComplex::Abs(noise[binItr]),2);
   }
   
   }
  //cout << sqrt((1/temp)*length*df)/0.003 <<endl;
   amp = amp/temp;
   double NormFFT=temp/sqrt(noise.size());
   TComplex comp_zero(0.,0.);
   vector<double> p_prod_ifftRe;
   //ignoring DC component;
   pProd[0] = comp_zero; 
    //===== 3. get max amplitude and delay =====

   //vector<double> p_ahat; //vector of all amplitudes
   double maxAmp = -1*numeric_limits<double>::infinity(); //maxAmp can be a negative number
   //double maxAmp = 0; //maxAmp can be a negative number
   double finalAmp = 0.0; 
   int idelay = 0;

   ComplexToRealIFFT(pProd, p_prod_ifftRe);  //C2R - no imaginary component!

   //pick out maximum amplitude and corresponding delay     
   for(int binItr=0; binItr < nBins; binItr++) {
     double amp0 =  p_prod_ifftRe[binItr];
     if(amp0 > maxAmp) { 
       maxAmp = amp0; 
       finalAmp = p_prod_ifftRe[binItr]/NormFFT;  
       idelay = binItr;
     }
     
     //search up to fwindow1, and after fwindow2
     if(binItr == fwindow1 - 1)  
       binItr = fwindow2 - 2; //-2 b/c for loop increments by 1 before next round, and c++ array convention adds -1
     
   }
   
   //compute phase factor with this delay
   int delay;
   delay = (idelay < nBins/2 ? idelay : (idelay - nBins));
   //ignoring DC component
   double chisq=0; 
       
   for(int binItr=1; binItr < nBins; binItr++)
     {
       double theta = 2.0*TMath::Pi()*((double)binItr/(double)nBins)*(double)delay;

       TComplex phase_factor(cos(theta), sin(theta));
       TComplex fit_fft( finalAmp*templateFFT[binItr]/phase_factor );
  
       double chisqBin =  pow(TComplex::Abs(signalFFT[binItr] - fit_fft), 2)/pow(TComplex::Abs(noise[binItr]),2);
       chisq += chisqBin;
     }
   signalVector.clear();
   templateVector.clear();
   signalFFT.clear();
   templateFFT.clear();
   pProd.clear();
   p_prod_ifftRe.clear();
   result[0]=finalAmp;
   //cout<<finalAmp<<endl;
   result[1]=delay*dt;
   result[2]=chisq/(f*f*(length-2));
   return 0;
}  


//--------------------------------------------------------------------------------------------------------------- MAIN BEGINS-----------------------------------------------------------------------------
int main(){
  
  //get the noise PSD
  TFile f_noise("TestPD2_Filter_23200104_0014.root");
  TH1D *h_noise = (TH1D*)f_noise.Get("zip1/PES1NoisePSD");
  double noisePSD[16384];
  for (int i=0; i<h_noise->GetNbinsX()-1; i++){
    noisePSD[i] = h_noise->GetBinContent(i+1);
  }

  /* This part is just to generate data for the python code, ignore
  ofstream filepsd;
  filepsd.open("PD2 PSD.txt");
  for(int l=0; l<16384 ;l++){
    filepsd << noisePSD[l]  << endl;
  }
  */
  //define the template in time domain
  //32768 is the number of time bins for the whole trace
  //3008 is the number of time bins for the hybrid trace 
  vector<double> template_time_32768;
  double template_time32768[32768];


  double signal_amp = 3e-3;

   //run a single time or multiple
  for (int k=0; k<1; k++){

    if (k%100==0) cout<<k<<endl;
    //generate the template -- look at the fuction
    generate_template(template_time_32768, 0);
    double max = *max_element( template_time_32768.begin(), template_time_32768.end());
    for (int i=0; i<32768; i++){
      template_time32768[i]=template_time_32768.at(i)/max;
    }
    //compressTime is simulating the hybrid readout

    //downsample_template(shifted_template32768, template_time2048);
    // make 'p' noise traces in a p by 32768 matrix. Then make fast and slow sampled PSDs from each trace (one trace = one row in this matrix). Average the PSDs (average over each column, because (column)_i =ith trace element.
    int p=1000;
    vector <vector<double>> noisetemp32768;
    for (int k=0; k<p;k++){
	double tempnoise[32768]={0.0};
	getNoise_time(tempnoise, noisePSD);
	vector <double> tempvector;
	ArrayToVector(tempnoise,tempvector,32768);
	noisetemp32768.push_back(tempvector);
      }
  
    //remove the baseline from the noise
      double prepulse_average1[1000] = {0.0};
      for(int k=0; k < p; k++){
	for(int i = 0; i < 16384; i++){
	  prepulse_average1[k] += noisetemp32768[k][i];
	}
	prepulse_average1[k] /= 16384.;
      }
      
      for(int k=0; k<p; k++){
	for(int i = 0; i < 32768; i++){
	  noisetemp32768[k][i] -= prepulse_average1[k];
	}
      }

      //making the fast and slow sampled noise traces
      vector <vector<double>> slowtemp2048;
      vector <vector<double>> fasttemp1024;
      for(int k=0; k<p; k++){
	double temp3008[3008];
	double slowtemp[2048];
	double fasttemp[1024];
	double temparray[32768];
	for(int j=0; j<32768; j++){
	  temparray[j] = noisetemp32768[k][j];
	}
	compressTime(temparray,temp3008);
	FastToSlowCompress(slowtemp,temp3008);
	vector <double> slowtempvector;
	ArrayToVector(slowtemp,slowtempvector,2048);
	slowtemp2048.push_back(slowtempvector);
	for(int i=0; i<1024;i++){
	  fasttemp[i] = temp3008[i+1024];
	}
	vector <double> fasttempvector;
	ArrayToVector(fasttemp,fasttempvector,1024);
	fasttemp1024.push_back(fasttempvector);
      }
      // Making data for python code, ignore
      vector <vector<double>> noise32768forfile;
      for(int i=0; i<p;i++){
	vector <double> noisetempPSD32768;
	vector <double> noisetempFFT32768;
	Time2PSD(noisetemp32768[i],625.000,noisetempPSD32768);
	PSD2FFT(noisetempPSD32768,noisetempFFT32768);
	noise32768forfile.push_back(noisetempFFT32768);
      }
      
      noisetemp32768.clear();
      // Making data for python code, ignore
      vector <double> noisePSDtwosided;
      for(int k=0;k<32768;k++){
	double sum=0;
	for(int j=0; j<p;j++){
	  sum +=pow(noise32768forfile[j][k],2);
	}
	noisePSDtwosided.push_back(sqrt(sum)/(double)p);
      }
      noise32768forfile.clear();
      
      
      //Find FFT of each trace
      vector <vector<double>> noise2048PSD;
      vector <vector<double>> noise1024PSD;
      for(int i=0; i<p; i++){
	vector<double> tempPSD2048;
	vector<double>tempFFT2048;
	vector<double> tempPSD1024;
	vector<double>tempFFT1024;
	Time2PSD(slowtemp2048[i],39.0625,tempPSD2048);
	PSD2FFT(tempPSD2048,tempFFT2048);
	noise2048PSD.push_back(tempFFT2048);

	Time2PSD(fasttemp1024[i],625.000,tempPSD1024);
	PSD2FFT(tempPSD1024,tempFFT1024);
	noise1024PSD.push_back(tempFFT1024);
      }
      //take average of FFT
      vector <double> noiseFFT2048;
      vector <double> noiseFFT1024;
      for(int k=0;k<2048;k++){
	double sum=0;
	for(int j=0; j<p;j++){
	  sum +=pow(noise2048PSD[j][k],2);
	}
	noiseFFT2048.push_back(sqrt(sum)/(double)p);
      }
      noise2048PSD.clear();

      for(int k=0;k<1024;k++){
	double sum=0;
	for(int j=0; j<p;j++){
	  sum +=pow(noise1024PSD[j][k],2);
	}
	noiseFFT1024.push_back(sqrt(sum)/(double)p);
      }
      noise1024PSD.clear();
      /*
      //open output files
      ofstream file1;
      ofstream file2;
      ofstream file3;
      ofstream file4;
      ofstream file5;
      ofstream file6;

      file1.open ("A1_delay.txt");
      file2.open ("A2_delay.txt");
      file3.open ("t1_delay.txt");
      file4.open ("A32768_delay.txt");
      file5.open ("t32768_delay.txt");
      file6.open ("t2_delay.txt");
      */
 //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     for (int trial=0; trial<1; trial++){
       vector<double> shifted_template_32768;
       double shifted_template32768[32768]={0.0};
       double template_time3008[3008]={0.0};
       double template_time2048[2048]={0.0};

       TRandom3 gRandom2(0);
       //double t = gRandom2.Gaus(0,12.8);
       double t = gRandom2.Uniform(-25.6,25.6);
       double t0=t;
       //if(t>0) t0=25.6;
       //else t0=-25.6;
       //double t0 = t;
       //int t0=100;
       //cout<<"delay="<<t0<<endl;

       //shift_template(template_time_32768,shifted_template_32768,t0);
       generate_template(shifted_template_32768, t0);
       double max_shifted = *max_element( shifted_template_32768.begin(), shifted_template_32768.end());
       for (int i=0; i<32768; i++){
	 shifted_template32768[i]=shifted_template_32768.at(i)/max_shifted;
       }
       shifted_template_32768.clear();
       // downsample_template(shifted_template32768, template_time2048);
       double signalnoise[32768]={0.0};
    //generate the noise spectrum in time domain, getNoise_time gives you the time domain noise given the PSD
     getNoise_time(signalnoise, noisePSD);
     //----------------------------------------------------remove baseline-----------------------------------------------------
     double prepulse_averagesig=0.0;
     for(int i = 0; i < 16384; i++){
	  prepulse_averagesig += signalnoise[i];
     }
     prepulse_averagesig /= 16384.;
      
     for(int i = 0; i < 32768; i++){
       signalnoise[i] -= prepulse_averagesig;
     }


     //----------------------------------------initialize variables---------------------------------------
     double signal3008[3008]={0.0};
     double signal32768[32768]={0.0};
     double signal2048[2048]={0.0};
     double FastSignal[1024]={0.0};
     double template_timeFast[1024]={0.0};
    //-----------------------------create your signal as template + noise---------------------------------
      for(int i=0;i<32768;i++){
	signal32768[i]  = signalnoise[i];
	signal32768[i] += signal_amp*template_time32768[i];
      }
      /*
      //-------------------------------------------------------------------------------------------Plots-------------------------------
      double  time32768[32768];
      for (int i=0; i<32768; i++){
	time32768[i]=i;
      }
     TCanvas* c=new TCanvas("c","Uniformly downsampled",800,600);
     TMultiGraph *mg = new TMultiGraph();
     TGraph *gr = new TGraph(32768, time32768, template_time32768);
     gr->SetMarkerStyle(4);
     gr->SetMarkerSize(0.35);
     gr->SetMarkerColor(kBlue);
     //     gr->GetXaxis()->SetRangeUser(16380,18000);
     mg->Add(gr);
     
     TGraph *grr = new TGraph(32768, time32768, shifted_template32768);
     grr->SetMarkerStyle(4);
     grr->SetMarkerSize(0.35);
     grr->SetMarkerColor(kRed);
     //grr->GetXaxis()->SetRangeUser(16380,18000);
     mg->Add(grr);
     mg->Draw("AP");

     
      double  time3008[3008];
      for (int i=0; i<3008; i++){
	time3008[i]=i;
      }
      TCanvas* ch=new TCanvas("ch","Hybrid",800,600);
      TMultiGraph *mgh = new TMultiGraph();
     TGraph *grrh = new TGraph(3008, time3008, template_time3008);
     grrh->SetMarkerStyle(4);
     grrh->SetMarkerSize(0.35);
     grrh->SetMarkerColor(kRed);
     //grr->GetXaxis()->SetRangeUser(16380,18000);
     mgh->Add(grrh);
     TGraph *grh = new TGraph(3008, time3008, template_time3008);
     grh->SetMarkerStyle(4);
     grh->SetMarkerSize(0.35);
     grh->SetMarkerColor(kBlue);
     //     gr->GetXaxis()->SetRangeUser(16380,18000);
     mgh->Add(grh);
     mgh->Draw("AP");
     */
      
//-----------------evaluate downsampled and fast sampled signals, template and noise----------------------------------------
      compressTime(signal32768, signal3008);
      FastToSlowCompress(signal2048,signal3008);
 
      for (int i=16384; i<17408; i++){
	FastSignal[i-16384]=signal32768[i];
	template_timeFast[i-16384]=shifted_template32768[i];
      }
      
          //---------------------------- DO the OF!-------------------------------
       
      double temp1[3]={0.0};
      double temp2[3]={0.0};
      double temp3[3]={0.0};
          
      // DoOptimalFilter(signal2048, template_time2048, noiseFFT2048, 2048, 19.07348633, 994, 1054,25.6, temp1);
      DoOptimalFilter(FastSignal, template_timeFast, noiseFFT1024, 1024, 610.3515625, 0, 1024, 1.6, temp2,625);
      DoOptimalFilter(signal32768, shifted_template32768, noisePSDtwosided, 32768, 19.07348633, 16354, 16484, 1.6, temp3,625);
      vector <double> true_template;
      vector <double> old_template;
      double true_temparr[32768];
      int delayBin= temp2[1]/1.6;
      ArrayToVector(shifted_template32768, old_template, 32768);
      shift_template(old_template, true_template, delayBin);
 
      for(int i=0; i<32768; i++){
	true_temparr[i]=true_template[i];
      }
      //old_template.clear();
      //true_template.clear();
      compressTime(true_temparr, template_time3008);
      FastToSlowCompress(template_time2048,template_time3008);
      DoOptimalFilter(signal2048, template_time2048, noiseFFT2048, 2048, 19.07348633, 994, 1054,25.6, temp1,39.0625);
      //--------------------------------------------------------------------------Finding best delay------------------------------------------
 
     
      //file1 << temp1[0] << endl;
      //file2 << temp2[0] << endl;
      //file3 << temp1[1] << endl;
      //file4 << temp3[0] << endl;
      //file5 << temp3[1] << endl;
      //file6 << temp2[1] << endl;
     }
      //-----------------------------------------------------------PLOTS-------------------------------------------------
     /*       
     double noise[32768];
     for (int k=0; k<32768; k++){
       noise[k]=noisePSDtwosided.at(k);
     }
     double  time32768[32768];
     for (int i=0; i<32768; i++){
       time32768[i]=i*19.07;
     }
     TCanvas* cn=new TCanvas("cn","Net",800,600);
     TMultiGraph *mg = new TMultiGraph();
     TGraph *PSDn = new TGraph(16384, time32768, noise);
     PSDn->SetMarkerStyle(7);
     PSDn->SetMarkerColor(kBlue);
     mg->Add(PSDn);
     TGraph *PSDog = new TGraph(16384, time32768, noisePSD);
     PSDog->SetMarkerColor(kRed);
     mg->Add(PSDog);
     //PSDn->SetTitle("Net");
     mg->Draw("AP");
     //PSDn->GetXaxis()->SetRangeUser(0,16384*610.35);
     cn->SetLogy(1);
     cn->SetLogx(1);
     
      double slowPSDarray[2048];
      for (int k=0; k<2048; k++){
      slowPSDarray[k]=noiseFFT2048.at(k);
      }
      double  time_2048[2048];
      for (int i=0; i<2048; i++){
	time_2048[i]=i*19.07;
      }
       TCanvas* c2=new TCanvas("c2","Uniformly downsampled",800,600);
       TGraph *PSD = new TGraph(1024, time_2048, slowPSDarray);
       PSD->SetMarkerStyle(7);
       PSD->SetTitle("Uniformly downsampled");
       PSD->Draw("AP");
       PSD->GetXaxis()->SetRangeUser(15,1024*19.07);
       c2->SetLogy(1);
       c2->SetLogx(1);

      double  time_Fast[1024];
      for (int i=0; i<1024; i++){
	time_Fast[i]=i*610.35;
      }
      
      double fastPSDarray[1024];
      for (int k=0; k<1024; k++){
      fastPSDarray[k]=noiseFFT1024.at(k);
      }
       TCanvas* c3=new TCanvas("c3","Fast sampled only",800,600);
       TGraph *PSD1 = new TGraph(512, time_Fast, fastPSDarray);
       PSD1->Draw("AP");
       PSD1->SetTitle("Fast sampled only");
       PSD1->SetMarkerStyle(7);
       PSD1->GetXaxis()->SetRangeUser(610.35,512*610.35);
       c3->SetLogy(1);
       c3->SetLogx(1);
     */
     return 0;
  }
}
