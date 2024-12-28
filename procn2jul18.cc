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
// ------- Only for BINARY filedat ------//
    uint32_t BinHeader1[HeaderLength];
    uint32_t BinHeader2[HeaderLength];
    uint32_t BinHeader3[HeaderLength];
    float Data1[MaxDataLength];
    float Data2[MaxDataLength];
    float Data3[MaxDataLength];
    float Data4[MaxDataLength];
    int DataLength1;
    int DataLength2;
    int DataLength3;
    int DataLength4;
    float xpos,ypos,chan=0;
//----Insert 08/01/2018----//
	uint64_t tstamp=0;
	uint64_t k=0;
	uint64_t tstamp0=0;
	uint64_t tstprev=0;
//--- CutOff friquency for Butterworth filter---//
    //double fcutoff=1000e6;
    double fcutoff=500e6;
    //double fcutoff=100e6;
    //double fcutoff=1e6;
    int total,tot;
    TString file1,file2,file3,file4,fout,rfile;
    if( argc == 1 || argc == 2)
    {
        cout<< "Default arguments"<<endl;
        //exit(0);
        file1="wave_5.dat";
        file2="wave_1.dat";
        file3="wave_3.dat";
        file4="wave_7.dat";
        fout="test.root";
        total=-100;
        tot=20;
        xpos=0;
        ypos=0;
    }
    if (argc == 7)
    {
        file1=argv[1];
        file2=argv[2];
        file3=argv[3];
        file4=argv[4];
        fout=argv[5];
//    }
        total=atoi(argv[6]);
        tot=20;
        xpos=0;
        ypos=0;
	chan=0;
        cout<<"total:"<< total<<" ToT:"<<tot << endl;
    }
    //-------- Args parsing for one channel and two trigers tests  -------//
    if (argc == 10 && atoi(argv[9])==2)
    {
//        chan=atoi(argv[9]);
        file1=argv[1];
        file2=argv[2];
            file3=argv[3];
 //       rfile=argv[10];
            fout=argv[4];
        total=atoi(argv[5]);
        tot=atoi(argv[6]);
        xpos=atof(argv[7])/1.;
        ypos=atof(argv[8])/1.;
    }
    else
    {
        if (argc==9)
    {
//      file4=argv[4];
        file1=argv[1];
        file2=argv[2];
            fout=argv[3];
    
//    }
        total=atoi(argv[4]);
        tot=atoi(argv[5]);
        xpos=atof(argv[6])/1.;
        ypos=atof(argv[7])/1.;
        cout<<"total:"<< total<<" ToT:"<<tot <<" Xpos="<<xpos<<"mm"<<" Ypos="<<ypos<<"mm\t Digitizer channel:"<<chan<< endl;
    }
    }
    if (argc == 4 || argc == 6)
    {
        file1=argv[1];
        file2=argv[2];
        fout=argv[3];
        if (argc>4)
        {
            total=atoi(argv[4]);
            tot=atoi(argv[5]);
        }
        else
        {
            total = -100;
            tot = 20;
        }
        xpos=0;
        ypos=0;
	chan=0;
    }
    if(argc!=7&&argc!=8&&argc!=9&&argc>1&&argc>2&&argc!=4&&argc!=6&&argc!=10&&argc!=11)
    {
    //     if (argc == 2) {
    //         total=atoi(argv[1]);
    //         tot=20;
    //     }
    //     else
    //     {
             cout<< "Exit2"<<endl;
             exit(0);
    //    }
    }
    cout<<"total2:"<< total<<" ToT2:"<<tot << endl;
//-------- END: Args parsing for one channel tests -------//
//-------- Creation of utime[] from RunDaq ---------------//
/* -----Commented since 20/07/18
    TFile* fr = new TFile(rfile,"READ");
    TTree* daq;
    if (fr->FindObjectAny("daq")>0) {
        daq=(TTree*)fr->Get("daq");
    }
    else{
        cout << "Wrong name of TTree: not daq" << endl;
        exit(1);
    }
    int ndaqev=daq->GetEntries();
    cout<<"DaqEntries:"<<ndaqev<<endl;
    struct serv{
        UInt_t time;
        UInt_t index;
        uint64_t utime;
    } service;
    uint64_t utime[ndaqev],dut[ndaqev],dtu,dtc,ut;
    daq->SetBranchStatus("*",0);
    daq->SetBranchStatus("service",1);
    
    daq->SetBranchAddress("service",&service.time);
    for (int i=0; i<ndaqev; i++)
    {
        daq->GetEntry(i);
        utime[i]=service.utime;
        if (i==0) {dut[i]=0;}
        else {dut[i]=utime[i]-utime[i-1];}
        //dtu=dut[i];
        printf("%d\t%llu\t%llu\r",i,utime[i],dut[i]); fflush(stdout);
        //getchar();

    }
----Commented since 20/07/18--- */    
//---------- END utime[] creation ------------------------//
    data_t ch1,ch1f,ch5,trgch;
    datar_t ch1r,ch1fr,ch5r,trgchr;
    TFile* fd = new TFile(fout,"RECREATE","exbeam_mcp_data");
    TTree* data=new TTree("exbeamdata","iMCP BINP testbeam data");
    data->Branch("chtrg1",&trgch.numev,"nevent/i:number/i:amp[1024]/F:ti[1024]/i:sum/f");
    data->Branch("ch1",&ch1.numev,"nevent/i:number/i:amp[1024]/F:ti[1024]/i:sum/f");
    data->Branch("ch1f",&ch1f.numev,"nevent/i:number/i:amp[1024]/F:ti[1024]/i:sum/f");
    data->Branch("chtrg1r",&trgchr.numev,"nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
    data->Branch("ch1r",&ch1r.numev,"nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
    data->Branch("ch1rf",&ch1fr.numev,"nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
    data->Branch("chtrg2",&ch5.numev,"nevent/i:number/i:amp[1024]/F:ti[1024]/i:sum/f");
    data->Branch("chtrg2r",&ch5r.numev,"nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
    data->Branch("x",&xpos,"x/f");
    data->Branch("y",&ypos,"y/f");
    data->Branch("ch",&chan,"ch/f");
//----- Insert 08/01/2018-----//
/* ---Commented since 20/07/18---   
data->Branch("ts",&tstamp,"ts/L");
    data->Branch("ut",&ut,"ut/L");
    data->Branch("tslot",&k,"tslot/L");
    data->Branch("dtu",&dtu,"dut/L");
    data->Branch("dtc",&dtc,"dct/L");
-----*/
//----- End Insert 08/01/2018-----//
/*	if (argc!=4 && argc!=6)
	{
        data->Branch("ch3",&ch3.numev,"nevent/i:number/i:amp[1024]/f:ti[1024]/i:sum/f");
        data->Branch("ch5",&ch5.numev,"nevent/i:number/i:amp[1024]/f:ti[1024]/i:sum/f");
        data->Branch("ch3r",&ch3r.numev,"nevent/i:ped/f:min/i:tmin/i:max/i:tmax/i:pedmax/i:pedmin/i:t1/i:t2/i");
	}
*/
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
//-------------------------------------//
    int i=0;
    int n=0;
    int number1=0; int number2=0; int number3=0; int number5=0;
    float sum1=0; float sum2=0; float sumf=0,sum5=0;
    float ped1=0,ped2=0,pedf=0,ped5=0;
    float pedmin1=0,pedmin2=0,pedminf=0,pedmin5=0;
    float pedmax1=0,pedmax2=0,pedmaxf=0,pedmax5=0;
    float min1=4096,min2=4096,minf=4096,min5=4096;
    float max1=0,max2=0,maxf=0,max5=0;
    int tmin1,tmin2,tminf,tmin5;
    int tmax1,tmax2,tmaxf,tmax5;
//-------------------------------------//
    if (file1!=NULL&&file2!=NULL)
    {
        FILE* fd1=fopen(file1,"r");
        FILE* fd2=fopen(file2,"r");
//        FILE* fd3;
//        if (file3!=NULL)
//        {
        FILE*    fd3=fopen(file3,"r");
//        }
cout<<"Files are opened!!!: fd1:"<<file1<<"\t fd2:"<<file2<<"\t fd3:"<<file3<<endl;
//        if (argc!=4 && argc!=6)
//        {
//            FILE* fd3=fopen(file3,"r");
//            FILE* fd4=fopen(file4,"r");
//        }
        int l=0;
        while (!feof(fd1))
        {
//cout<<"Inside while(!feof(fd1))"<<endl;
            if(fread(BinHeader1, sizeof(*BinHeader1), HeaderLength, fd1) != HeaderLength || fread(BinHeader2, sizeof(*BinHeader2), HeaderLength, fd2) != HeaderLength || fread(BinHeader3, sizeof(*BinHeader3), HeaderLength, fd3) != HeaderLength /* ||
                        if(feof(fd1) || feof(fd2) /*|| feof(fd3) || feof(fd4)*/)
            {
                fprintf(stderr,"Some fread() error during header reading!\n");
                break;
            }
            DataLength1=(BinHeader1[0]-HeaderLength*sizeof(*BinHeader1))/sizeof(*Data1);
            DataLength2=(BinHeader2[0]-HeaderLength*sizeof(*BinHeader2))/sizeof(*Data2);
            //cout << "DataLength1="<< DataLength1<< " DataLength2="<< DataLength2<< endl;
            //if (file3!=NULL && fread(BinHeader3, sizeof(*BinHeader3), HeaderLength, fd3) != HeaderLength)
            //{
                DataLength3=(BinHeader3[0]-HeaderLength*sizeof(*BinHeader3))/sizeof(*Data3);
                //DataLength4=(BinHeader2[0]-HeaderLength*sizeof(*BinHeader2))/sizeof(*Data4);
            //}
            if(DataLength1>MaxDataLength || DataLength2>MaxDataLength || DataLength3>MaxDataLength /*|| DataLength1!=MaxDataLength || DataLength2!=MaxDataLength || DataLength4 =MaxDataLength*/)
                {
                    fprintf(stderr,"Record is too long (>%d)! %d/%d/%d\n",MaxDataLength,DataLength1,DataLength2,DataLength3);
                    break;
                }
	//----- INsert 08/01/2018 --------//
            //if (BinHeader1[5]!=BinHeader2[5] || BinHeader1[4]!=BinHeader2[4] || BinHeader1[5]!=BinHeader3[5] || BinHeader1[4]!=BinHeader3[4])
            //{
                //fprintf(stderr,"Timesatmp for event from file1 and file2 is different (%d/%d/%d)!\n",BinHeader1[5],BinHeader2[5],BinHeader3[5]);
                //break;
            //}
            //if (file3!=NULL && (BinHeader1[5]!=BinHeader3[5] || BinHeader1[4]!=BinHeader3[4]))
            //{
             //   fprintf(stderr,"Timesatmp for event from file1 and file3 is different (%d/%d)!\n",BinHeader1[5],BinHeader3[5]);
            //    break;
            //}
            //dtc=BinHeader1[5]-tstprev;
//            if (n==0) {tstamp0=BinHeader1[5]; dtc=0; }
            if (BinHeader1[5]<tstprev)
            {
                k++;
//                cout<<"Curr.t:"<<BinHeader1[5]<<"\t Prev.t:"<<tstprev<<"\t Init.t:"<<tstamp0<<endl;
//                dtc=(BinHeader1[5]+(pow(2,31)-1))-tstprev;
            }
		tstamp=(BinHeader1[5]+k*(pow(2,31)-1))-tstamp0;
            tstprev=BinHeader1[5];
            chan=BinHeader2[3];
            Int_t number=BinHeader2[4];
	//----- End INsert 08/01/2018 --------//
            if(fread(Data1, sizeof(*Data1), DataLength1, fd1) != DataLength1 || fread(Data2, sizeof(*Data2), DataLength2, fd2) != DataLength2 /*|| fread(Data3, sizeof(*Data3), DataLength3, fd3) != DataLength3 || fread(Data4, sizeof(*Data4), DataLength4, fd4) != DataLength4*/)
            {
                if(feof(fd1)||feof(fd2)/*||feof(fd3)||feof(fd4)*/)
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
            else
            {
                if(file3!=NULL && fread(Data3, sizeof(*Data3), DataLength3, fd3) != DataLength3 /*|| fread(Data4, sizeof(*Data4), DataLength4, fd4) != DataLength4*/)
                {
                    if(feof(fd3)/*||feof(fd2)||feof(fd3)||feof(fd4)*/)
                    {
                        fprintf(stderr,"EOF3: Some fread() error during data reading!\n");
                        break;
                    }
                    else
                    {
                        fprintf(stderr,"Some fread3() error during data reading!\n");
                        break;
                    }
                }

//cout<<"Data reading"<<endl;
                ped1=0; ped2=0; pedf=0; ped5=0;
                min1=4096; min2=4096; minf=4096; min5=4096;
                max1=0; max2=0; maxf=0; max5=0;
                int t11=0; int t12=1024;
                int t21=0; int t22=1024;
                int tf1=0; int tf2=1024;
                int t51=0; int t52=1024;
                int flmin1=0, flmin2=0, flminf=0, flmin5=0;
                trgchr.numev=n; ch1r.numev=n; ch1fr.numev=n; ch5r.numev=n;
                trgch.number=number; ch1.number=number; ch1f.number=number; ch5.number=number;
                trgch.numev=n; ch1.numev=n; ch1f.numev=n; ch5.numev=n;
                double din[DataLength1];
                double* dou=new double[DataLength1];
                for(int i=0; i<DataLength1; i++) {din[i]=Data2[i]; }
                int res=Butterworth(din, 2e-10, fcutoff, dou);
                //int res=Butterworth(din, 1e-10, fcutoff, dou);
                //int res=Butterworth(din, 1e-9, fcutoff, dou);
//                cout<<"din/dou:"<<din[100]<<"/"<<dou[0]<<" res="<<res<<endl;
                for(int i=0; i<DataLength1; i++)
                {
                    trgch.amp[i]=Data1[i]; trgch.ti[i]=i; ch1.amp[i]=Data2[i]; ch1.ti[i]=i;
                    ch1f.amp[i]=dou[i]; ch1f.ti[i]=i;
                    if (file3!=NULL)
                    {
                        //ch3.amp[i]=Data3[i]; ch3.ti[i]=i;
                        ch5.amp[i]=Data3[i]; ch5.ti[i]=i;
                    }
                    sum1+=Data1[i]; sum2+=Data2[i]; sumf+=dou[i];
                    if (file3!=NULL)
                    {
                        //sum3+=Data3[i];
                        sum5+=Data3[i];
                    }
                    if (i < 100 && i>10)
                    {
                        ped1+=Data1[i]/89;
                        pedmin1=min1; pedmax1=max1;
                        ped2+=Data2[i]/89;
                        pedmin2=min2; pedmax2=max2;
                        pedf+=dou[i]/89;
                        pedminf=minf; pedmaxf=maxf;
                        if (file3!=NULL)
                        {
//                            ped3+=Data3[i]/89;
//                            pedmin3=min3; pedmax3=max3;
                            ped5+=Data3[i]/89;
                            pedmin5=min5; pedmax5=max5;
                        }
/**/
                    }
                    if (Data1[i]<min1) {min1=Data1[i]; tmin1=i;}
                    if (Data2[i]<min2) {min2=Data2[i]; tmin2=i;}
                    if (dou[i]<minf) {minf=dou[i]; tminf=i;}
                    if (file3!=NULL)
                    {
//                        if (Data3[i]<min3) {min3=Data3[i]; tmin3=i;}
                        if (Data3[i]<min5) {min5=Data3[i]; tmin5=i;}
                    }
/**/
                    if (Data1[i]>max1) {max1=Data1[i]; tmax1=i;}
                    if (Data2[i]>max2) {max2=Data2[i]; tmax2=i;}
                    if (dou[i]>maxf) {maxf=dou[i]; tmaxf=i;}
                    if (file3!=NULL)
                    {
//                        if (Data3[i]>max3) {max3=Data3[i]; tmax3=i;}
                        if (Data3[i]>max5) {max5=Data3[i]; tmax5=i;}
                    }
/**/
//            i++;
// ------ For POSITIVE TRG and SIGNAL Pulses------//
if (tot<0)
{
                    if(i>100 && (Data1[i]-ped1)>tot*(-1) && flmin1==0 && t11==0)   {t11=i; flmin1=1;}
                    if(i>100 && (Data1[i]-ped1)<tot*(-1) && flmin1==1 && t12==1024){t12=i; flmin1=0;}
                    if(i>100 && (Data2[i]-ped2)>tot*(-1) && flmin2==0 && t21==0)   {t21=i; flmin2=1;}
                    if(i>100 && (Data2[i]-ped2)<tot*(-1) && flmin2==1 && t22==1024){t22=i; flmin2=0;}
                    if(i>100 && (dou[i]-pedf)>tot*(-1) && flminf==0 && tf1==0)   {tf1=i; flminf=1;}
                    if(i>100 && (dou[i]-pedf)<tot*(-1) && flminf==1 && tf2==1024){tf2=i; flminf=0;}
}
// ------ For NEGATIVE TRG and SIGNAL Pulses------//
if (tot>0){
            if(i>100 && (ped1-Data1[i])>tot && flmin1==0 && t11==0)   {t11=i; flmin1=1;}
            if(i>100 && (ped1-Data1[i])<tot && flmin1==1 && t12==1024){t12=i; flmin1=0;}
                    if(i>100 && (ped2-Data2[i])>tot && flmin2==0 && t21==0)   {t21=i; flmin2=1;}
                    if(i>100 && (ped2-Data2[i])<tot && flmin2==1 && t22==1024){t22=i; flmin2=0;}
                    if(i>100 && (pedf-dou[i])>tot && flminf==0 && tf1==0)   {tf1=i; flminf=1;}
                    if(i>100 && (pedf-dou[i])<tot && flminf==1 && tf2==1024){tf2=i; flminf=0;}
}
// ------ For NEGATIVE TRG and SIGNAL Pulses------//
                    if (file3!=NULL)
                    {
//                        if(i>100 && (ped3-Data3[i])>tot && flmin3==0 && t31==0)   {t31=i; flmin3=1;}
//                        if(i>100 && (ped3-Data3[i])<tot && flmin3==1 && t32==1024){t32=i; flmin3=0;}
                        if(i>100 && (ped5-Data3[i])>tot && flmin5==0 && t51==0)   {t51=i; flmin5=1;}
                        if(i>100 && (ped5-Data3[i])<tot && flmin5==1 && t52==1024){t52=i; flmin5=0;}
                    }
 /**/
                }
                //Int_t
		float dtp=28000;
                //ut=utime[l]-utime[0];
                //cout<<"Tstamp/Utime:"<<tstamp*(2*8.5e-9)<<"/"<<ut*1e-6<<endl;
                //dtu=dut[l];
                //if (abs(tstamp*(1.706670e-02)-ut)<dtp)
                //{
                //    cout<<"Tstamp==utime"<<endl;
                n++;
                //    l++;
                //}
                //else
                //{
                
//                cout<<"D(Tstamp)!=D(utime): Dtc/Dut:"<<tstamp*(1.706667e-02)<<"/"<<ut<<endl;
                    //l++;
                //}
                trgch.sum=sum1/DataLength1; ch1.sum=sum2/DataLength2; ch1f.sum=sumf/DataLength2;
                trgchr.ped=ped1; ch1r.ped=ped2; ch1fr.ped=pedf;
                trgchr.pedmin=pedmin1; ch1r.pedmin=pedmin2; ch1fr.pedmin=pedminf;
                trgchr.pedmax=pedmax1; ch1r.pedmax=pedmax2; ch1fr.pedmax=pedmaxf;
                trgchr.min=min1; ch1r.min=min2; ch1fr.min=minf;
                trgchr.tmin=tmin1; ch1r.tmin=tmin2; ch1fr.tmin=tminf;
                trgchr.max=max1; ch1r.max=max2; ch1fr.max=maxf;
                trgchr.tmax=tmax1; ch1r.tmax=tmax2; ch1fr.tmax=tmaxf;
                trgchr.t1=t11; ch1r.t1=t21; ch1fr.t1=tf1;
                trgchr.t2=t12; ch1r.t2=t22; ch1fr.t2=tf2;

               if (file3!=NULL)
                {
/*                    ch3.sum=sum3/DataLength3;*/
                    ch5.sum=sum5/DataLength3;
//                    ch3r.ped=ped3;
                    ch5r.ped=ped5;
//                    ch3r.pedmin=pedmin3;
                    ch5r.pedmin=pedmin5;
//                    ch3r.pedmax=pedmax3;
                    ch5r.pedmax=pedmax5;
//                    ch3r.min=min3;
                    ch5r.min=min5;
//                    ch3r.tmin=tmin3;
                    ch5r.tmin=tmin5;
//                    ch3r.max=max3;
                    ch5r.max=max5;
//                    ch3r.tmax=tmax3;
                    ch5r.tmax=tmax5;
//                    ch3r.t1=t31;
                    ch5r.t1=t51;
//                    ch3r.t2=t32;
                    ch5r.t2=t52;
                }
/**/
                //if (abs(tstamp*(1.706670e-02)-ut)<dtp)
                //{
                data->Fill();
                //}
		// --- for negative TRG----//
                htrg->Fill(ped1-min1);
		// --- for positive TRG----//
                //htrg->Fill(max1-ped1);
		//-------------------------//
                if (abs(t11-t21)/5<15) {h1->Fill(ped2-min2); ttsimcp1->Fill((float)(t11-t21)/5.0);}
                if (n==total) {break;}
                printf("%d\r",n); fflush(stdout);
                sum1=0; sum2=0; sumf=0; sum5=0;
                tmin1=0; min1=4096; tmax1=0; max1=0;
                tmin2=0; min2=4096; tmax2=0; max2=0;
                tminf=0; minf=4096; tmaxf=0; maxf=0;
                tmin5=0; min5=4096; tmax5=0; max5=0;
            }
        }
    }
        //else{ fprintf(stderr,"Wrong filename!\n");}
    //}

        cout<<"SUM5:"<<sum1<<":"<<sum2<<":"<<sumf<<endl;
//   TCanvas* mcan;
//    mcan=new TCanvas("mcan","MCP Amps & TTS", 0, 0, 800, 600);
//    mcan->Divide(1,2,0.01,0.01,0);
//    mcan->cd(1);    h1->Draw();
//    ttsimcp1->SetXTitle("TTS, ns");
//    mcan->cd(2);    ttsimcp1->Draw();
    //    mcan->Update();
//    ttsimcp5->SetXTitle("TTS, ns");
//    mcan->cd(6);  ttsimcp5->SetLineColor(2);  ttsimcp5->Draw(); ttsimcp1->SetLineColor(4); ttsimcp1->Draw("same");
//    mcan->Update();
    
    h1->Write(); htrg->Write(); h3->Write(); h5->Write();
    ttsmcpt->Write(); ttsimcp1->Write(); ttsimcp5->Write();
    data->Write();
//    mcan->Write();
//    getchar();
//    fd.Close();
    //theApp->Run();
    fd->Close();
}
