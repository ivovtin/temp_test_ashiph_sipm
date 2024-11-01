#include<iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TMath.h"
#include "TH1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#include <ctime>
#include <chrono>

#define MaxDataLength 1024
//#define HeaderLength 8
#define HeaderLength 6

//const float rate = 5*pow(10,9); // частота дискретизации
const float rate = 1*pow(10,9); // частота дискретизации
const int trgthr = 1800; //порог на срабатывания триггера.
Bool_t b_slow_interferences = kTRUE; // заменить на kFALSE чтобы не пытаться убрать медленную наводку  ??
//ham
//Bool_t b_invert_signal = kFALSE;// инверсия сигнала если выбирать  kTRUE если сигнал отрицательный(чтобы пьедестал был слева)
//ketek
Bool_t b_invert_signal = kFALSE;// инверсия сигнала если выбирать  kTRUE если сигнал отрицательный(чтобы пьедестал был слева)

//ham
//const float hmint = -22E-9, hmaxt = 3E-9, hshift = -90E-9; //область интегрирования , hshift сдвижка для шумов, секунды
//ketek
const float hmint = 50E-9, hmaxt = 150E-9, hshift = -110E-9; //область интегрирования , hshift сдвижка для шумов, секунды


//const float hmint = 75E-9, hmaxt = 115E-9, hshift = -50E-9; //область интегрирования , hshift сдвижка для шумов, секунды
const float gtfitmin = -50E-9, gtfitax = -35E-9; //область фитирования для вычетания подложки

int nbins = 1024;
float mincortime, maxcortime;
float minAmpW, maxAmpW, minAmpT, maxAmpT,minAmpWcor, maxAmpWcor;
int maxoscnum = 1000000; //чтобы не ждать очень долго обработки. Читается только это количество первых осциллограмм
//int maxoscnum = 5000; //чтобы не ждать очень долго обработки. Читается только это количество первых осциллограмм

Int_t npeaks = 30; //ждя вычитания фона

void timecor(vector< vector< float > >&  voscW, vector< vector< float > >&  voscT,   vector< vector< float > >& vt, vector< vector< float > >& vtcor   ) {
    //cout << (voscW.size()) << endl;
    float tcor = 0;
    for ( long unsigned int j = 0; j<voscW.size(); j++ ) {
        vt.push_back(vector<float>());
        vtcor.push_back(vector<float>());

        for (long unsigned int i = 0; i < voscW[0].size(); i++) {
            vt[j].push_back( (float)i / rate );
        }
    //time correction
        for (long unsigned int i = 0; i < voscW[0].size(); i++) {
            if (voscT[j][i] <trgthr) {
                tcor =  vt[j][i];
                break;
            }
        }
        for (long unsigned int i = 0; i < voscW[0].size(); i++) {
            vtcor[j].push_back(   vt[j][i]  - tcor);
            //vtcor[j].push_back(   vt[j][i]);
            if (j == 0 && i == 0 ) {
                mincortime = vt[j][i]  - tcor;
                maxcortime = vt[j][i]  - tcor;
            } else {
                if ( (vt[j][i]  - tcor) < mincortime ) { mincortime = vt[j][i]  - tcor; }
                if ( (vt[j][i]  - tcor) > maxcortime ) { maxcortime = vt[j][i]  - tcor; }
            }
        }
    }
}

void filldata(TString fw, TString ft, vector< vector< float > >&  voscW, vector< vector< float > >&  voscT,  Bool_t invercion = kFALSE ) {

    uint32_t BinHeaderWAVE[HeaderLength];
    uint32_t BinHeaderTRG[HeaderLength];
    float DataWAVE[MaxDataLength];
    float DataTRG[MaxDataLength];
    int DataLengthWAVE;
    int DataLengthTRG;

    FILE* fwaveH=fopen(fw,"r");
    FILE* ftrgH=fopen(ft,"r");
    cout<<"Files are opened: "<<fwaveH<<"\t"<<ftrgH<<endl;
    fprintf(fwaveH, "blah");
    int n=0;
    while (!feof(fwaveH) && n < maxoscnum) {
    
        if(fread(BinHeaderWAVE, sizeof(*BinHeaderWAVE), HeaderLength, fwaveH) != HeaderLength || fread(BinHeaderTRG, sizeof(*BinHeaderTRG), HeaderLength, ftrgH) != HeaderLength /* || if(feof(fd1) || feof(fd2) */)       {
            fprintf(stderr,"Some fread() error during header reading!\n");
            break;
        }
        DataLengthWAVE=(BinHeaderWAVE[0]-HeaderLength*sizeof(*BinHeaderWAVE))/sizeof(*DataWAVE);
        DataLengthTRG=(BinHeaderTRG[0]-HeaderLength*sizeof(*BinHeaderTRG))/sizeof(*DataTRG);
        //if (DataLengthWAVE != DataLengthTRG ) {
        if (DataLengthWAVE != DataLengthTRG || DataLengthTRG != 1024) {
            fprintf(stderr,"Signal event  data length not equal trigger event length!\n");
            cout << " DataLengthWAVE = " << DataLengthWAVE << endl;
            cout << " DataLengthTRG = " << DataLengthTRG << endl;

            break;
        }

        if(DataLengthWAVE>MaxDataLength || DataLengthTRG>MaxDataLength  ) {
            fprintf(stderr,"Record is too long (>%d)! %d/%d\n", MaxDataLength, DataLengthWAVE, DataLengthTRG);
            break;
        }

        if(fread(DataWAVE, sizeof(*DataWAVE), DataLengthWAVE, fwaveH) != DataLengthWAVE || fread(DataTRG, sizeof(*DataTRG), DataLengthTRG, ftrgH) != DataLengthTRG )  {
            if(feof(fwaveH)||feof(ftrgH))    {
                fprintf(stderr,"EOF: Some fread() error during data reading!\n");
                break;
            }
            else    {
                fprintf(stderr,"Some fread() error during data reading!\n");
                break;
            }
        }
        //cout<<"Data reading"<<endl;
        voscW.push_back(vector<float>());
        voscT.push_back(vector<float>());

        for(int i=0; i<DataLengthWAVE; i++)        {
            if (invercion == kTRUE)  {voscW[n].push_back(-DataWAVE[i]);}
            else  {voscW[n].push_back(DataWAVE[i]);}
            voscT[n].push_back(DataTRG[i]);
            //cout << DataWAVE[i] << "\t" <<DataTRG[i] << endl;
        }
        n++;
    }
    cout << "Read "<< n << " waveforms" << endl;
    fclose(fwaveH);
    fclose(ftrgH);

}

void slow_interferences_correction(vector< vector< float > >&  vosc,  vector< vector< float > >&  voscWaveAmpCor ) {
//найдем подложку у осциллограмм, вычтем ее и положим результат в вектор                
    for (long unsigned int i = 0; i < vosc.size(); i++) { 
        voscWaveAmpCor.push_back(vector<float>());
        auto h = new TH1F("h","osc",1023,0,1023);
        for (int j = 0; j < vosc[i].size(); j++ ) {
            h->SetBinContent(j,vosc[i][j]) ;         
        }
        auto s = new TSpectrum();
        auto hb = s->Background(h,40,"same");
        //if (hb) c->Update();
        for (int j = 0; j < vosc[i].size(); j++) {
            voscWaveAmpCor[i].push_back(   h->GetBinContent(j)  - hb->GetBinContent(j));
        }
    }
    
/*
    for (long unsigned int i = 0; i < vosc.size(); i++) { 

        TGraph  gcor(vosc[i].size(), &vosct[i][0], &vosc[i][0]);
        gcor.Fit("pol0","Q","",gtfitmin, gtfitax);
        fitl = gcor.GetFunction("pol0");
        Double_t p0 = fitl->GetParameter(0);
        vampcor.push_back(p0);
    }
    //сдвинем осциллограммы по вертикали для компенсации "медленной" наводки
    for (int i = 0; i < vosc.size(); i++) {
        for (long unsigned int j =0; j < vosc[0].size(); j++) {
            vosc[i][j] = vosc[i][j] - vampcor[i];
        }
    }
    */
}
void integrate(vector< vector< float > >&  vosc, vector< vector< float > >&  vosct, float first,float last, vector< float >&  vint) {
    float integ = 0;
    float timestep = vosct[0][1] - vosct[0][0];
    for (int i = 0; i < vosc.size(); i++ ) {
        integ = 0;
        for (int j = 0; j< vosc[0].size(); j++) {
            if ( vosct[i][j] >= first &&  vosct[i][j] < last) {
                integ = integ + timestep * vosc[i][j];
            }
        }
        vint.push_back(integ);
    }
}

const std::string timeFormat{"%Y-%m-%d_%H-%M-%S"};

int64_t getTimeStamp(const std::string& timeString)
{
    std::istringstream m_istream{timeString};
    std::tm m_tm { 0 };
    std::time_t m_timet { 0 };
    m_istream >> std::get_time(&m_tm, timeFormat.c_str());
    m_timet = std::mktime(&m_tm);
    //m_timet *= 1000; // convert to milliseconds
    return m_timet;
}


void histo(const char* dirname,  int channel = 1, TString outfile = "histogram.root" ) //
{
    TH1::AddDirectory(false);//ownership
    TH2::AddDirectory(false);

//вектора для данных
    std::vector< vector< float > > voscWave;
    std::vector< vector< float > > voscTrg;
    std::vector< vector< float > > vt;
    std::vector< vector< float > > vtcor;
    //std::vector<  float > vampcor;
    std::vector< vector< float > > voscWaveAmpCor;
    void* dir;
    dir = gSystem->OpenDirectory(dirname);
    if( !dir ) {
        cerr<<dirname<<": No such directory"<<endl;
       	return 0;
    }

    //читаем данные и заплняем ветора с данными
    TString filewave=TString(dirname) + "/"  + "wave_" + channel + ".dat";
    TString filetrg=TString(dirname)  + "/"  + "TR_0_0.dat";
    //cout << filewave << endl;
    filldata(filewave, filetrg, voscWave, voscTrg, b_invert_signal) ;


    //создаем канву    
    Double_t w = 1500;
    Double_t h = 750;
    auto c = new TCanvas("c", "c", w, h);
    //auto c1 = new TCanvas("c1", "c1", w, h);
    c->Divide(2,3);

    //заполняем ветора времени вычитая время прихода триггера чтобы компенсировать разброс времен и заполняем гистограммы  
    timecor(voscWave, voscTrg, vt, vtcor   );

    c->cd(1);
    //нарисуем несколько гистограмм с осциллограммами для примера.
    int ngraph=10;//число графиков для отрисовки
    if (voscWave.size() < ngraph) {ngraph =voscWave.size(); }
    for (long unsigned int i = 0; i < ngraph; i++) { 
	    auto h = new TH1F("h","osc",1024,0,1023);
	    for (int j = 0; j < voscWave[i].size(); j++ ) {
		h->SetBinContent(j,voscWave[i][j]) ;         
	    }
	    if (i==0) {h->Draw();}
	    else { h->Draw("same");}
    }
    c->cd(2);
    //Это должно попытаться убрать медленную наводку, чтобы отключить kTRUE заменить на kFALSE
    if (b_slow_interferences==kTRUE) { 
    //if (b_slow_interferences==kFALSE) { 
        slow_interferences_correction(voscWave, voscWaveAmpCor);
    }

    //найдем минимум и максимум среди исходных гистограмм и откоректированнх, если они есть
    for (int i = 0; i < voscWave.size(); i++) {
        for (int j =0; j < voscWave[0].size(); j++) {
            if (i == 0 && j == 0 ) {
                maxAmpW = voscWave[i][j];
                minAmpW = voscWave[i][j];   
                maxAmpT = voscTrg[i][j];   
                minAmpT = voscTrg[i][j];   
                if (b_slow_interferences==kTRUE) {
                    maxAmpWcor = voscWaveAmpCor[i][j];
                    minAmpWcor = voscWaveAmpCor[i][j];   
                }
                else {
                    maxAmpWcor = 0;
                    minAmpWcor = 0;   

                }
            } else {
                if (voscWave[i][j] > maxAmpW ) {maxAmpW = voscWave[i][j]; }
                if (voscWave[i][j] < minAmpW ) {minAmpW = voscWave[i][j]; }
                if (voscTrg[i][j] > maxAmpT ) {maxAmpT = voscTrg[i][j]; }
                if (voscTrg[i][j] < minAmpT ) {minAmpT = voscTrg[i][j]; }
                if (b_slow_interferences==kTRUE) {
                    if (voscWaveAmpCor[i][j] > maxAmpWcor ) {maxAmpWcor = voscWaveAmpCor[i][j]; }
                    if (voscWaveAmpCor[i][j] < minAmpWcor ) {minAmpWcor = voscWaveAmpCor[i][j]; }
                }
                else {
                    maxAmpWcor = 0;
                    minAmpWcor = 0;   
                }
            }
        }
    }  
    //    cout << "minAmpW = " << minAmpW << "\t" << "maxAmpW = " << maxAmpW << endl;
    //    cout << "minAmpT = " << minAmpT << "\t" << "maxAmpT = " << maxAmpT << endl;

    //нарисуем осциллограммы с вычтенной подложкой, если есть
    if (b_slow_interferences==kTRUE) {
        if (voscWaveAmpCor.size() < ngraph) {ngraph =voscWaveAmpCor.size(); }
        for (long unsigned int i = 0; i < ngraph; i++) { 
	        auto h = new TH1F("h","osc",1024,0,1023);
	        for (int j = 0; j < voscWaveAmpCor[i].size(); j++ ) {
		        h->SetBinContent(j,voscWaveAmpCor[i][j]) ;         
	        }
	        if (i==0) {h->Draw();}
	        else { h->Draw("same");}
        }
    
    
    }


    c->cd(3);
    //нарисуем все осциллограммы  на 2д гистограмме
    TH2F *hosc = new TH2F("hosc",";time,s;Amp,channels",(maxcortime-mincortime)*rate-1,mincortime,maxcortime,int(maxAmpW-minAmpW-40),minAmpW+40,maxAmpW); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
    //TProfile *hoscH = new TProfile("hoscH","Hamamatsu;time,s;Amp,channels",(maxcortime-mincortime)*rate,mincortime,maxcortime); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
    for (int i = 0; i < vtcor.size(); i ++  ) {
        for (long unsigned int j = 0; j < vtcor[i].size(); j++ ) {
            hosc->Fill(vtcor[i][j],voscWave[i][j]);
        }
    }
    //hoscH->Draw("colz");
    hosc->Draw("colz");
    //hosc->Draw("box");
    //линии границ для гисторамм
    TLine *clinehs1 = new TLine(hmint,minAmpW+40,hmint,maxAmpW);
    TLine *clinehs2 = new TLine(hmaxt,minAmpW+40,hmaxt,maxAmpW);
    TLine *clinehn1 = new TLine(hmint+hshift,minAmpW+40,hmint+hshift,maxAmpW);
    TLine *clinehn2 = new TLine(hmaxt+hshift,minAmpW+40,hmaxt+hshift,maxAmpW);
    clinehs1->SetLineColor(kRed);
    clinehs1->Draw();
    clinehs2->SetLineColor(kRed);
    clinehs2->Draw();
    clinehn1->SetLineColor(kBlue);
    clinehn1->Draw();
    clinehn2->SetLineColor(kBlue);
    clinehn2->Draw();

    //нарисуем все осциллограммы без фона если они есть  на 2д гистограмме
    if (b_slow_interferences==kTRUE) {
        c->cd(4);
        TH2F *hoscC = new TH2F("hoscC",";time,s;Amp,channels",(maxcortime-mincortime)*rate-1,mincortime,maxcortime,int(maxAmpWcor+40),-40,maxAmpWcor); //
    //TProfile *hoscH = new TProfile("hoscH","Hamamatsu;time,s;Amp,channels",(maxcortime-mincortime)*rate,mincortime,maxcortime); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
        for (int i = 0; i < vtcor.size(); i ++  ) {
            for (long unsigned int j = 0; j < vtcor[i].size(); j++ ) {
                hoscC->Fill(vtcor[i][j],voscWaveAmpCor[i][j]);
            }
        }
        //hoscH->Draw("colz");
        hoscC->Draw("colz");
        //hosC->Draw("box");
        //линии границ для гисторамм
        TLine *clinehs1cor = new TLine(hmint,-40,hmint,maxAmpWcor);
        TLine *clinehs2cor = new TLine(hmaxt,-40,hmaxt,maxAmpWcor);
        TLine *clinehn1cor = new TLine(hmint+hshift,-40,hmint+hshift,maxAmpWcor);
        TLine *clinehn2cor = new TLine(hmaxt+hshift,-40,hmaxt+hshift,maxAmpWcor);
        clinehs1cor->SetLineColor(kRed);
        clinehs1cor->Draw();
        clinehs2cor->SetLineColor(kRed);
        clinehs2cor->Draw();
        clinehn1cor->SetLineColor(kBlue);
        clinehn1cor->Draw();
        clinehn2cor->SetLineColor(kBlue);
        clinehn2cor->Draw();
    }

/*
    c->cd(3);
    //нарисуем все осциллограммы  на 2д гистограмме
    TProfile *hoscP = new TProfile("hoscH","Hamamatsu;time,s;Amp,channels",(maxcortime-mincortime)*rate-1,mincortime,maxcortime); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
    for (int i = 0; i < vtcor.size(); i ++  ) {
        for (long unsigned int j = 0; j < vtcor[i].size(); j++ ) {
            hoscP->Fill(vtcor[i][j],voscWave[i][j]);
        }
    }
    //hoscH->Draw("colz");
    hoscP->Draw();

    c->cd(2);
    //нарисуем все осциллограммы триггера  на 2д гистограмме
    TH2F *hoscT = new TH2F("hoscT",";time,s;Amp,channels",(maxcortime-mincortime)*rate,mincortime,maxcortime,int(maxAmpT-minAmpT),minAmpT,maxAmpT); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
    //TProfile *hoscH = new TProfile("hoscH","Hamamatsu;time,s;Amp,channels",(maxcortime-mincortime)*rate,mincortime,maxcortime); //-40 чтобы обрезать наводку ввиде коротких отрицательных выбросов
    for (int i = 0; i < vtcor.size(); i ++  ) {
        for (long unsigned int j = 0; j < vtcor[i].size(); j++ ) {
            hoscT->Fill(vtcor[i][j],voscTrg[i][j]);
        }
    }
    hoscT->Draw("colz");
*/
    c->cd(5);
    //найдем площадь под осциллограммами 
    std::vector<  float > vSignalCharge; 
    std::vector<  float > vNoiseCharge; 
    
    integrate(voscWave, vtcor, hmint, hmaxt, vSignalCharge) ;
    integrate(voscWave, vtcor, hmint+hshift, hmaxt+hshift, vNoiseCharge) ;
    
    //заполним данными гистограммы
    float minhs = *min_element(vSignalCharge.begin(), vSignalCharge.end());
    float maxhs = *max_element(vSignalCharge.begin(), vSignalCharge.end());
 
    float minhn = *min_element(vNoiseCharge.begin(), vNoiseCharge.end());
    float maxhn = *max_element(vNoiseCharge.begin(), vNoiseCharge.end());

    TH1F *hs = new TH1F("hsig",";",1000,minhn,maxhs);
    for (long unsigned int i = 0; i < vSignalCharge.size(); i++ ) {
        hs->Fill(vSignalCharge[i]);
 
    }
	
    Double_t Smm=0, Srms=0, Nmm=0, Nrms=0;

    /*
    Smm = hs->GetMean();
    Srms= hs->GetRMS();
    cout<<"Signal spectrum:"<<"\t"<<"A="<<Smm<<"+-"<<Srms<<endl;
    */
    hs->SetLineWidth(2);

    TH1F *hn = new TH1F("hnoise",";",1000,minhn,maxhn);
    for (long unsigned int i = 0; i < vNoiseCharge.size(); i++ ) {
        hn->Fill(vNoiseCharge[i]);
    }
    /*
    Nmm = hn->GetMean();
    Nrms= hn->GetRMS();
    cout<<"Noise spectrum:"<<"\t"<<"A="<<Nmm<<"+-"<<Nrms<<endl;
    */
    hn->SetLineWidth(3);
    hn->SetLineColor(4);
    hs->SetLineWidth(2);
    hs->SetLineColor(2);
    
    hs->Draw("");
    hn->Draw("sames");

    //Сохраним гистограммы
    TString postfix = TString::Format("%dch",channel).Data();
    outfile = "histogram_" + postfix + ".root";
    TString outf=TString(dirname)  + "/" + outfile;
    TFile histogram(outf,"RECREATE");
 
    hs->Write();
    hn->Write();

    //найдем площадь под осциллограммами
    
    if (b_slow_interferences==kTRUE) {
        c->cd(6);
        std::vector<  float > vSignalChargeAmpCor; 
        std::vector<  float > vNoiseChargeAmpCor; 
    
        integrate(voscWaveAmpCor, vtcor, hmint, hmaxt, vSignalChargeAmpCor) ;
        integrate(voscWaveAmpCor, vtcor, hmint+hshift, hmaxt+hshift, vNoiseChargeAmpCor) ;
    
        //заполним данными гистограммы
        minhs = *min_element(vSignalChargeAmpCor.begin(), vSignalChargeAmpCor.end());
        maxhs = *max_element(vSignalChargeAmpCor.begin(), vSignalChargeAmpCor.end());
        minhn = *min_element(vNoiseChargeAmpCor.begin(), vNoiseChargeAmpCor.end());
        maxhn = *max_element(vNoiseChargeAmpCor.begin(), vNoiseChargeAmpCor.end());

        TH1F *hsacor = new TH1F("hsigCor",";",1000,minhn,maxhs);
        for (long unsigned int i = 0; i < vSignalChargeAmpCor.size(); i++ ) {
            hsacor->Fill(vSignalChargeAmpCor[i]); 
        }
        hsacor->SetLineWidth(2);
	
    	Smm = hsacor->GetMean();
    	Srms= hsacor->GetRMS();
    	cout<<"Signal spectrum:"<<"\t"<<"A="<<Smm<<"+-"<<Srms<<endl;

        TH1F *hnacor = new TH1F("hnoiseCor",";",1000,minhn,maxhs);
        for (long unsigned int i = 0; i < vNoiseChargeAmpCor.size(); i++ ) {
            hnacor->Fill(vNoiseChargeAmpCor[i]);
        }

	Nmm = hnacor->GetMean();
    	Nrms= hnacor->GetRMS();
    	cout<<"Noise spectrum:"<<"\t"<<"A="<<Nmm<<"+-"<<Nrms<<endl;

        hnacor->SetLineWidth(3);
        hnacor->SetLineColor(4);
        hsacor->SetLineWidth(2);
        hsacor->SetLineColor(2);
        hnacor->Draw("");
        
        hsacor->Draw("sames");
   
        hsacor->Write();
        hnacor->Write();

    }
    
    //cout<<int(maxAmpK-minAmpK) << "\t" << minAmpK << "\t" << maxAmpK << endl;

    voscWaveAmpCor.clear();     
    voscWave.clear(); 
    voscTrg.clear(); 
    vt.clear(); 
    vtcor.clear(); 

    int64_t resTime = getTimeStamp(dirname);
    cout << "resTime:" << resTime << std::endl;

    TString SoutF = TString::Format("channel_%d.dat",channel).Data();
    ofstream SfileOut;
    SfileOut.open(SoutF,ios_base::app | ios_base::ate);
    SfileOut << resTime << "\t" << Smm << "\t" << Srms << "\n";
    SfileOut.close();

    TString NoutF = TString::Format("noise_channel_%d.dat",channel).Data();
    ofstream NfileOut;
    NfileOut.open(NoutF,ios_base::app | ios_base::ate);
    NfileOut << resTime << "\t" << Nmm << "\t" << Nrms << "\n";
    NfileOut.close();

}


