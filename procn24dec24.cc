#include <Riostream.h>
#include <stdio.h>
#include <stdint.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TApplication.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#define MaxDataLength 1024
#define HeaderLength 6

int Butterworth(double ind[1024], double dT, double C, double oud[1024]);

using namespace std;
struct data_t
{
        Int_t numev;
        Int_t number;
        float amp[1024];
        Int_t ti[1024];
        Float_t sum;
};

struct datar_t
{
        Int_t numev;
        Float_t ped;
        Int_t min;
        Int_t tmin;
        Int_t max;
        Int_t tmax;
        Int_t pedmax;
        Int_t pedmin;
        Int_t t1;
        Int_t t2;
};
//--------------------------------------------------------------------------
// This function returns the data filtered. Converted to C# 2 July 2014.
// Original source written in VBA for Microsoft Excel, 2000 by Sam Van
// Wassenbergh (University of Antwerp), 6 june 2007.
//--------------------------------------------------------------------------

//public static double[]
int Butterworth(double indata[1024], double deltaTimeinsec, double CutOff, double data[1024])
{
    if (sizeof(indata) == 0) {cout<<"Zero input data"<<endl; exit(0);};
    
    double Samplingrate = 1 / deltaTimeinsec;
    int dF2 = 1024; //sizeof(indata)/sizeof(indata[0]) - 1;        // The data range is set with dF2

    double* Dat2 = new double[dF2 + 4]; // Array with 4 extra points front and back
    //double data[dF2];
    for (int i=0; i<dF2; i++){data[i] = indata[i];} // Ptr., changes passed data

    if (CutOff == 0) {cout<<"Nothing to do with data"<<endl; exit(0);};
    
    // Copy indata to Dat2
    for (int r = 0; r < dF2; r++) {
        Dat2[2 + r] = indata[r];
    }
    Dat2[1] = Dat2[0] = indata[0];
    Dat2[dF2 + 3] = Dat2[dF2 + 2] = indata[dF2];
    
    const double pi = 3.14159265358979;
    double wc = tan(CutOff * pi / Samplingrate);
//    cout<<"wc="<<wc<<" =tan("<<(CutOff * pi / Samplingrate)<<")"<<" dF2="<<dF2<<endl;
    double k1 = 1.414213562 * wc; // Sqrt(2) * wc
    double k2 = wc * wc;
    double a = k2 / (1 + k1 + k2);
    double b = 2 * a;
    double c = a;
    double k3 = b / k2;
    double d = -2 * a + k3;
    double e = 1 - (2 * a) - k3;
    
    // RECURSIVE TRIGGERS - ENABLE filter is performed (first, last points constant)
    double* DatYt = new double[dF2 + 4];
    DatYt[1] = DatYt[0] = indata[0];
    for (int s = 2; s < dF2 + 2; s++) {
        DatYt[s] = a * Dat2[s] + b * Dat2[s - 1] + c * Dat2[s - 2]
        + d * DatYt[s - 1] + e * DatYt[s - 2];
    }
    DatYt[dF2 + 3] = DatYt[dF2 + 2] = DatYt[dF2 + 1];
    
    // FORWARD filter
    double* DatZt = new double[dF2 + 2];
    DatZt[dF2] = DatYt[dF2 + 2];
    DatZt[dF2 + 1] = DatYt[dF2 + 3];
    for (int t = -dF2 + 1; t <= 0; t++) {
        DatZt[-t] = a * DatYt[-t + 2] + b * DatYt[-t + 3] + c * DatYt[-t + 4]
        + d * DatZt[-t + 1] + e * DatZt[-t + 2];
    }
    
    // Calculated points copied for return
    for (int p = 0; p < dF2; p++) {
        data[p] = DatZt[p];
    }
    
    return sizeof(data)/sizeof(data[0]);
}

int main(int argc, char* argv[])
{
    int nch = 0;	
    int total = -1000, tot = 2;

    TString file[64];
    TString fname = "test.root", rfile;

    cout<<argc<<endl;

    if(argc > 1)
    {
    	nch = atoi(argv[1]);	
	for(int j=0; j<nch; j++ )	
    	{
           file[j]=argv[j+2];
	   cout<<nch<<"\t"<<file[j]<<endl;
    	}
	cout<<nch<<"\t"<<fname<<"\t"<<total<<"\t"<<tot<<endl;
	if( argc > (2 + nch) ) fname = argv[nch+2];
	if( argc > (2 + nch + 1) ) total = atoi(argv[nch+3]);
        if( argc > (2 + nch + 2) ) tot = atoi(argv[nch+4]);  

	cout<<nch<<"\t"<<fname<<"\t"<<total<<"\t"<<tot<<endl;
    }
    else{
    	cout<< "Please, set number of channels and list of waveforms: procn24dec24 7 wave_0.dat wave_2.dat wave_4.dat wave_6.dat wave_9.dat wave_11.dat wave_13.dat test.root -100 2 " <<endl;
	exit(0);
    }    
    cout<<"total2:"<< total<<" ToT2:"<<tot << endl;

    //--- CutOff friquency for Butterworth filter---//
    double fcutoff=500e6;

    uint32_t BinHeader[nch][HeaderLength];
    float Data[nch][MaxDataLength];
    int DataLength[nch];

    float chan=0;
    uint64_t tstamp=0;
    uint64_t k[nch]={0}; 
    uint64_t tstamp0=0;
    uint64_t tstprev=0;

    data_t ch[nch];
    datar_t chr[nch];
    TFile* fout = new TFile(fname,"RECREATE","exbeam_data");
    TTree* data=new TTree("exbeamdata","testbeam data");
 
    char branchname1[nch];
    char branchname2[nch];
    for( int j=0; j<nch; j++)
    {	    
        sprintf(branchname1,"ch%d",j+1); 
    	data->Branch(branchname1, &ch[j].numev, "nevent/i:number/i:amp[1024]/F:ti[1024]/i:sum/f");
	sprintf(branchname2,"chr%d",j+1);
	data->Branch(branchname2, &chr[j].numev, "nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
    }


    TH1F *htrg,*h1,*h3,*h5,*ttsimcp1,*ttsimcp5,*ttsmcpt;
    htrg=new TH1F("htrg","Amplitude TRG_MCP",4096,-99.5,3995.5);
    h1=new TH1F("htest1","Amplitude iMCP(ch1)",4096,-99.5,3995.5);
    h3=new TH1F("hmcp","Amplitude MCP PMT",4096,-99.5,3995.5);
    h5=new TH1F("htest2","Amplitude iMCP(ch5)",4096,-99.5,3995.5);
    ttsimcp1=new TH1F("ttsimcp1","TTS iMCP(T_iMCP-T_MCPPMT",1023,-99.5,99.5);
    ttsimcp5=new TH1F("ttsimcp5","TTS iMCP(T_iMCP-T_MCPPMT",1023,-99.5,99.5);
    ttsmcpt=new TH1F("ttsmcpt","TTS MCP PMT(T_MCPPMT-T_MCPPMTtrg",1023,-99.5,99.5);
    //-------- End TTree anouncement -----//

    char* strr1=new char [257]; char* strr2=new char [257]; char* strr3=new char [257]; char* strr4=new char [257];
    char ev[5],num[7];
    float ampl1, ampl2, ampl3, ampl5, amp1[1024], amp2[1024], amp3[1024], amp5[1024];
    int ti[1024];
    for (int i=0; i<1024; i++)
    {
        amp1[i]=0; amp2[i]=0; amp3[i]=0; amp5[i]=0;
    }
    cout << "Before ifstream" << endl;

    if (file[0]!=NULL)
    {
	FILE* fd[nch];
	for( int j=0; j<nch; j++ )
	{
	   fd[j] = fopen(file[j],"r");
	   cout<<"File is opened!!!: fd:"<<file[j]<<"\t"<<endl;
        }
	
   	int n=0;
    	float sum1[nch]={0};
    	float ped1[nch]={0};
    	float pedmin1[nch]={0};
    	float pedmax1[nch]={0};
    	float min1[nch]={4096};
    	float max1[nch]={0};
    	int tmin1[nch];
    	int tmax1[nch];

        int l=0;
        	
	while (!feof(fd[0]))
        {
		for( int j=0; j<nch; j++ )
		{	
	    		if( fread(BinHeader[j], sizeof(*BinHeader[j]), HeaderLength, fd[j]) != HeaderLength)
			{
                		fprintf(stderr,"Some fread() error during header reading!\n");
                		break;
            		}
			DataLength[j]=(BinHeader[j][0]-HeaderLength*sizeof(*BinHeader[j]))/sizeof(*Data[j]);

			if( DataLength[j] > MaxDataLength )
                	{
                    		fprintf(stderr,"Record is too long (>%d)! %d/\n",MaxDataLength,DataLength[j]);
                   	 	break;
                	}

			if( fread(Data[j], sizeof(*Data[j]), DataLength[j], fd[j]) != DataLength[j] )
            		{
                   		if( feof(fd[j]) )
                   		{
                     			fprintf(stderr,"EOF: Some fread() error during data reading!\n");
                     			break;
                   		}
                   		else
                   		{
                     			fprintf(stderr,"Some fread() error during data reading!\n");
                     			break;
                   		}
            		}

            		if (BinHeader[j][5]<tstprev)
            		{
                    		k[j]++;
            		}
	
	    		tstamp=(BinHeader[j][5]+k[j]*(pow(2,31)-1))-tstamp0;
            		tstprev=BinHeader[j][5];

            		Int_t number=BinHeader[j][4];

            		ped1[j]=0;
            		min1[j]=4096;
            		max1[j]=0; 
            		int t11=0; int t12=1024;
            		int flmin1=0;
           
			ch[j].numev=n; ch[j].number=number;    
			chr[j].numev=n;  

	    		//double din[DataLength[0]];
            		//double* dou = new double[DataLength[0]];
            
	    		//for(int i=0; i<DataLength[0]; i++) {din[i] = Data[1][i]; }
                
	    		//int res=Butterworth(din, 2e-10, fcutoff, dou);
                
	    		for( int i=0; i<DataLength[j]; i++ )
            		{
				ch[j].amp[i]=Data[j][i]; 
				ch[j].ti[i]=i;
                
                   		sum1[j]+=Data[j][i]; 

                		if( i < 100 && i>10 )
                		{
                   			ped1[j]+=Data[j][i]/89;
                   			pedmin1[j]=min1[j]; 
					pedmax1[j]=max1[j];
          			}
                    
				if( Data[j][i]<min1[j] ) { min1[j]=Data[j][i]; tmin1[j]=i; }
				if( Data[j][i]>max1[j] ) { max1[j]=Data[j][i]; tmax1[j]=i; } 

				// ------ For POSITIVE TRG and SIGNAL Pulses------//
				if (tot<0)
				{
                    			if(i>100 && (Data[j][i]-ped1[j])>tot*(-1) && flmin1==0 && t11==0)   {t11=i; flmin1=1;}
                    			if(i>100 && (Data[j][i]-ped1[j])<tot*(-1) && flmin1==1 && t12==1024){t12=i; flmin1=0;}
				}

				// ------ For NEGATIVE TRG and SIGNAL Pulses------//
				if( tot>0 )
				{
            	    			if(i>100 && (ped1[j]-Data[j][i])>tot && flmin1==0 && t11==0)   {t11=i; flmin1=1;}
            	    			if(i>100 && (ped1[j]-Data[j][i])<tot && flmin1==1 && t12==1024){t12=i; flmin1=0;}
				}
	    		}

			ch[j].sum=sum1[j]/DataLength[j];
            		chr[j].ped=ped1[j];
            		chr[j].pedmin=pedmin1[j];
            		chr[j].pedmax=pedmax1[j];
            		chr[j].min=min1[j];
            		chr[j].tmin=tmin1[j];
            		chr[j].max=max1[j]; 
            		chr[j].tmax=tmax1[j];
            		chr[j].t1=t11; 
            		chr[j].t2=t12; 
		}
	    	float dtp=28000;
            	n++;

        	// --- for negative TRG----//
    	        //htrg->Fill(ped1-min1);
	 	// --- for positive TRG----//
            	//htrg->Fill(max1-ped1);
            	//-------------------------//
            	//if( abs(t11-t21)/5<15 ) {h1->Fill(ped2-min2); ttsimcp1->Fill((float)(t11-t21)/5.0);}
            	if( n==total ) {break;}
            	printf("%d\r",n); fflush(stdout);

		data->Fill();	
   		sum1[nch]={0};
            	tmin1[nch]={0}; min1[nch]={4096}; tmax1[nch]={0}; max1[nch]={0};
	}
    }
    
    //h1->Write(); htrg->Write(); h3->Write(); h5->Write();
    //ttsmcpt->Write(); ttsimcp1->Write(); ttsimcp5->Write();
    data->Write();
	
    fout->Close();
    
}
