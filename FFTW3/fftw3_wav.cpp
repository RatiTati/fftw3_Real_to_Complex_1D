#include <iostream>
#include <fstream>
#include <Graph.hh>
#include <WaveReader.hh>
#include <fftw3.h>
#include <cmath>

using namespace std;

double cosFi(WaveReader* wr){
    double Preal = 0;
    double Vrms = 0;
    double Irms = 0;
    for(int i=0; i<wr->GetNsamples(); ++i){
        Preal += -1* (wr->GetData(1)[i] * wr->GetData(0)[i]);
        Vrms += (wr->GetData(1)[i] * wr->GetData(1)[i]);
        Irms += (wr->GetData(0)[i] * wr->GetData(0)[i]);
    }
    return Preal/sqrt(Vrms*Irms);
}

int main(){

    //default constructor
    WaveReader myF;
    myF.SetVpp(10.);
    myF.Open("load.wav", 10000);


    double SampleRate = myF.GetSampleRate();
    int N = myF.GetNsamples();
    double df = SampleRate/N;
    double t = N/SampleRate;

    double* inV;
    double* inI;
    fftw_complex *outV;
    fftw_complex *outI;
    fftw_plan pV;
    fftw_plan pI;

    inV = (double*) fftw_malloc(sizeof(double) * N);
    inI = (double*) fftw_malloc(sizeof(double) * N);

    outV = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
    outI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));


    //filling "in" arrays
    for(int i=0; i<N; i++){
        inV[i] = myF.GetData(0)[i];
    }

    for(int i=0; i<N; i++){
        inI[i] = myF.GetData(1)[i];
    }


    //plan
    pV = fftw_plan_dft_r2c_1d((unsigned long)N, inV, outV,  FFTW_ESTIMATE);
    pI = fftw_plan_dft_r2c_1d((unsigned long)N, inI, outI,  FFTW_ESTIMATE);

    //execute
    fftw_execute(pV);
    fftw_execute(pI);


    //filling Graphs
    Graph myGraphVolt("Graph_Voltage");

    for(int i=0; i<N/2+1; ++i){
        myGraphVolt.AddPoint(i * df, sqrt(outV[i][0] * outV[i][0] + outV[i][1] * outV[i][1]));
    }

    Graph myGraphCurr("Graph_Current");

    for(int i=0; i<N/2+1; ++i){
        myGraphCurr.AddPoint(i * df, sqrt(outI[i][0] * outI[i][0] + outI[i][1] * outI[i][1]));
    }

    fftw_destroy_plan(pV);
    fftw_destroy_plan(pI);
    fftw_free(inV);
    fftw_free(outV);
    fftw_free(inI);
    fftw_free(outI);

    myGraphCurr.Print();
    myGraphVolt.Print();
    myGraphVolt.Draw(0);
    myGraphCurr.Draw(0);


    Graph gVol = ("Voltage_Graph");
    for(int i=0; i<myF.GetNsamples(); ++i){
        gVol.AddPoint(myF.GetTimeAxis()[i], myF.GetData(1)[i]);
    }

    Graph gCurr = ("Current_Graph");
    for(int i=0; i<myF.GetNsamples(); ++i){
        gCurr.AddPoint(myF.GetTimeAxis()[i], myF.GetData(0)[i]);
    }
    gVol.Draw(0);
    gCurr.Draw(0);

    cout<<"Cos(fi) = "<<cosFi(&myF)<<endl;
}
