#include <iostream>
#include <RandomGenerator.hh>
#include <SignalGenerator.hh>
#include <Histo.hh>
#include <Graph.hh>
#include <fftw3.h>
#include <cmath>

using namespace std;

int main(){
//     Histo myH("MyHisto", 100, -5, 40);
//     RandomGenerator r;
//     for(int i=0; i<100000; i++){
//         myH.Fill(r.GetGaus2(12, 4));
//     }
//     double mean = myH.GetMean();
//     cout<<"mean of histo: "<<mean<<endl;
//     double sigma = myH.GetRMS();
//     cout<<"mean of sigma: "<<sigma<<endl;
//     myH.Draw(1);


    SignalGenerator sg;
    sg.AddHarmonics(10, 20, 0);
    sg.AddHarmonics(5, 40, 0);
    sg.SetNoise(30);

    const double SampleRate = 500;
    const unsigned long N = 100000;
    const double df = SampleRate/N;
    const double t = N/SampleRate;

    sg.Generate(t, (unsigned long)N);

    double* in;
    fftw_complex *out;
    fftw_plan p;

    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

    for(int i=0; i<N; i++){
        in[i] = sg.getSignal()[i];
    }


    Graph myG("Graph_in_signal");
    for(int i=0; i<N; i++){
        myG.AddPoint(sg.getTime()[i], sg.getSignal()[i]);
    }

    myG.Draw(0);



    p = fftw_plan_dft_r2c_1d((unsigned long)N, in, out,  FFTW_ESTIMATE);

    fftw_execute(p);

    Graph myGraph("Graph_out_fftw");

    for(int i=0; i<N/2+1; ++i){
        myGraph.AddPoint(i * df, sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]));
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    myGraph.Draw(0);
}
