#ifndef train_GP_C
#define train_GP_C

#include "waveform_median_variance.h"
#include "pick_waveform_samples.h"
#include <HCal.hpp>
// APEX peak-search
#include <GP.hpp>
#include <Histo1D.hpp>
// ROOT 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TPad.h>
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
#include <TGraphErrors.h> 
#include <TStyle.h> 
#include <Math/ProbFuncMathCore.h>
#include <TString.h> 
#include <TBox.h>
#include <TAxis.h> 
// stdlib
#include <fstream> 
#include <sstream>
#include <vector> 
#include <string> 
#include <cmath> 



double waveform_variance(const double *v) {
    static constexpr double norm_fact = 1./((double)HCal::Block::n_waveform_samps); 

    double sum{0.}, sum2{0.};
    for (int i=0; i<HCal::Block::n_waveform_samps; i++) { sum += v[i]; sum2 += v[i]*v[i]; }
    sum  *= norm_fact; 
    sum2 *= norm_fact;
    return sum2 - sum*sum;  
}

void train_GP(
    int row, int col,
    std::string path_infile="waveforms.root", 
    std::string path_out_graphic="plots/waveform_samples_row-7_col-6_logrange_-8.0_1.0"
)
{
    using ROOT::VecOps::RVec; 
    
    const int n_waveforms = 40; 

    auto wf = pick_waveform_samples(row,col, -8.0, +1.0, n_waveforms); 

    std::cout << "number of samples: " << wf.size() << "\n"; 

    double stddev, mean; 
    waveform_median_val_stddev(row,col, mean, stddev); 

    RVec<double> pts_err(HCal::Block::n_waveform_samps, stddev); 

    RVec<double> pts_n(HCal::Block::n_waveform_samps); 
    for (int i=0; i<HCal::Block::n_waveform_samps; i++) pts_n[i] = ((double)i); 

    double ymin{-0.03}, ymax{0.10}; 

    auto canv = new TCanvas; 
    gStyle->SetOptStat(0); 

    if (!path_out_graphic.empty())
        path_out_graphic += ".pdf"; 
    
    int page=0; 
    for (int i=0; i<n_waveforms; i++) {

        auto samps = wf[i].data(); 

        auto g = new TGraphErrors(HCal::Block::n_waveform_samps, 
            pts_n.data(), samps, 
            nullptr,      pts_err.data()
        );
    
        double event_variance = waveform_variance(samps); 

        canv->Clear(); 
        canv->SetTopMargin(0.15); 

        g->SetTitle(Form("#splitline{Waveform sample, row-%i, col-%i}{sample %i, variance: %.4e}",row,col,i,event_variance)); 
        g->SetMinimum(ymin); 
        g->SetMaximum(ymax);
        g->SetMarkerStyle(kFullCircle); 
        g->SetMarkerSize(0.40);  
        g->Draw("AP"); 

        canv->Modified(); 
        canv->Update(); 

        //save each 'file' to a pdf
        
        ++page; 
        switch (page) {

            //first page 
            case 1 : 
                canv->Print(std::string{path_out_graphic + "["}.c_str(), "pdf"); break; 
            
            //last page
            case n_waveforms : 
                canv->Print(std::string{path_out_graphic + "]"}.c_str(), "pdf"); break; 

            //any other page
            default : 
                canv->Print(path_out_graphic.c_str(), "pdf"); 
        }     
    }   

    return; 
}



#endif