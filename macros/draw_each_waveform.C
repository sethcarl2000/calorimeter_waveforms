
#include "find_quantile.h"
#include "waveform_median_variance.h"
#include <HCal.hpp>
// ROOT 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <TH1D.h>
#include <TCanvas.h>
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
// analysis_utils
#include <analysis_utils/ROOT.hpp>
// stdlib
#include <string> 
#include <stdexcept> 
#include <array> 
#include <cmath> 
#include <memory> 
#include <limits> 
#include <iostream> 

namespace {
    constexpr int n_samples = 50; 
    constexpr int n_cols = 12; 
    constexpr int n_rows = 24; 

    constexpr double min_amplitude = 1e-25; 
}; 


double max_variance_subrange(const HCal::Block& b, int window_size=8) {
    double norm = 1./((double)window_size); 

    const int n_steps = HCal::Block::n_waveform_samps - window_size; 

    const auto& wf = b.samps; 

    double max_variance = -1.; 
    for (int i=0; i<n_steps; i++) {

        double sum{0.}, sum2{0.};  
        for (int j=i; j<i+window_size; j++) {
            sum  += wf[j]; 
            sum2 += wf[j]*wf[j]; 
        }

        sum *= norm; sum2 *= norm; 

        if (double variance = (sum2 - sum*sum); variance > max_variance) {
            max_variance = variance; 
        }
    }
    return max_variance; 
}

double max_slope(const HCal::Block& b) {
    double max_slope = -1; 
    for (int i=0; i<HCal::Block::n_waveform_samps-1; i++) {
        
        if (double slope = std::fabs( b.samps[i+1] - b.samps[i] ); slope > max_slope) { 
            max_slope = slope; 
        } 
    }
    return max_slope; 
}



template<typename T> using block_rptr = std::array<ROOT::RDF::RResultPtr<T>, HCal::n_blocks>; 

void draw_each_waveform(
    std::string path_infile="waveforms.root", 
    std::string path_out_graphic="plots/all_samps_log-y"
)
{
    using ROOT::VecOps::RVec; 

    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("T", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 
    //now, create the 'hcal' objects

    block_rptr<TH1D> histos; //[HCal::n_blocks]; 
    block_rptr<TH1D> h_variance; //[HCal::n_blocks]; 
    block_rptr<TH1D> h_log_variance; //[HCal::n_blocks]; 
    block_rptr<TH1D> h_subrange_log_variance; //[HCal::n_blocks]; 
    block_rptr<TH1D> h_log_max_slope; //[HCal::n_blocks]; 

    //there are a lot of samples which appear to have null-values near 0. i'm going to cut any waveform which has these. 
    const double minimum_sample_amplitude = 1e-15; 

    block_rptr<double> b_variance2, b_variance; 

    const int subrange_size = 10; 

    path_out_graphic += ".pdf"; 
    std::cout << "booking histograms..." << std::flush; 
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {
        
            double norm_fact = 1./((double)HCal::Block::n_waveform_samps); 

            const int ind = row*HCal::n_cols + col; 
    
            auto df_samps = df

                    //cut out waveforms which have samples below a given magnitude
                    .Filter([row,col,minimum_sample_amplitude](const HCal& d){
                        for (const auto& xi : d.get_block(row,col).samps) { 
                            if ( std::fabs(xi) < minimum_sample_amplitude ) return false; 
                        }
                        return true; 
                    }, {"hcal_data"})

                    .Define("my_samps", [row,col](const HCal& d){ 
                        const auto& block = d.get_block(row,col); 
                        return RVec<double>( block.samps.cbegin(), block.samps.cend() ); 
                    }, {"hcal_data"})

                    .Define("samp_variance", [norm_fact](const RVec<double>& v){
                        double xx{0.}, x{0.}; 
                        for (double xi : v) { xx += xi*xi; x += xi; }
                        xx *= norm_fact; 
                        x  *= norm_fact; 
                        return xx - x*x; 
                    }, {"my_samps"})

                    .Define("subrange_variance", [row,col,subrange_size](const HCal& d)
                    {
                        return max_variance_subrange(d.get_block(row,col), subrange_size); 
                    }, {"hcal_data"})

                    .Define("log_max_slope", [row,col](const HCal& d)
                    {
                        return std::log10( max_slope(d.get_block(row,col)) ); 
                    }, {"hcal_data"})

                    .Define("log_subrange_variance", [](double var){ return std::log10(var); }, {"subrange_variance"})

                    .Define("samp_log_variance", [](double var){ return std::log10(var); }, {"samp_variance"})

                    .Define("samp_variance2", [](double var){ return var*var; }, {"samp_variance"}); 

            b_variance[ind]  = df_samps.Sum<double>("samp_variance");
            b_variance2[ind] = df_samps.Sum<double>("samp_variance2"); 

            histos[ind] = df_samps   
                .Histo1D<RVec<double>>({"h_draw", 
                    Form("#splitline{All waveforms amplitude samples}{row %4i, col %4i};waveform amplitude;",row,col), 
                    200, -0.03, +0.15}, "my_samps"); 
                    
            h_subrange_log_variance[ind] = df_samps
                .Histo1D<double>({"h_log_subrange_variance", 
                    Form("#splitline"
                        "{Max-window Variance of waveform samples: row %4i, col %4i}"
                        "{window size: %i}"
                        ";log_{10}( max[ <x_{i}^{2}> - <x_{i}>^{2} ] );",
                        row,col,subrange_size), 
                    100, -8, 1}, "log_subrange_variance");
            
            h_log_variance[ind] = df_samps
                .Histo1D<double>({"h_log_variance", 
                    Form("Variance of waveform samples: row %4i, col %4i;log_{10}( max[ <x_{i}^{2}> - <x_{i}>^{2} ] );",row,col), 
                    100, -8, 1}, "samp_log_variance");
            
            h_log_max_slope[ind] = df_samps
                .Histo1D<double>({"h_log_max_slope", 
                    Form("Log of max slope: row %4i, col %4i;log_{10}( max[ x_{i+1} - x_{i} ] );",row,col), 
                    100, -8, 1}, "log_max_slope");
            
            h_variance[ind] = df_samps
                .Histo1D<double>({"h_variance", 
                    Form("#splitline"
                        "{Per-Event Variance of waveform samples: row %4i, col %4i}"
                        "{window size: %i}"
                        ";max[ <x_{i}^{2}> - <x_{i}>^{2} ];",
                        row,col,subrange_size), 
                    200, 0., 1.e-4}, "samp_variance");
        }
    }
    std::cout << "done.\n"; 

    auto count = *df.Count(); 

    std::printf("done fetching %.4e events.\ndrawing...\n", (double)count); 

    auto canv = new TCanvas; 
    
    std::vector<double> pts_median(HCal::n_blocks), pts_stddev(HCal::n_blocks); 

    auto tf1_gaus = new TF1("f_gaus", [](double *X, double *par){

        double A = par[0]; 
        double x0 = par[1]; 
        double sigma = par[2]; 

        double arg = (X[0]-x0) / sigma;
        return A * std::exp( -0.5 * arg*arg );

    }, histos[0]->GetXaxis()->GetXmin(), histos[0]->GetXaxis()->GetXmax(), 3); 
    
    tf1_gaus->SetNpx(1000);
    tf1_gaus->SetLineStyle(kDashed); 

    int page=0; 
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {


            const int ind = row*HCal::n_cols + col; 

            auto& h = histos[ind]; 

            auto& h_log = h_log_variance[ind]; 

            //if (!((row==7) && (col==6))) continue; 
                
            canv->Clear(); 
            canv->SetLogy(1); 
            canv->SetTopMargin(0.15); 

            double median_stddev, median_val; 

            waveform_median_val_stddev(row,col, median_val, median_stddev);  

            tf1_gaus->SetParameter(0, h->GetMaximum()); 
            tf1_gaus->SetParameter(1, median_val); 
            tf1_gaus->SetParameter(2, median_stddev); 

            h->DrawCopy(); 
            
            tf1_gaus->Draw("SAME"); 

            canv->Modified(); 
            canv->Update(); 
            /*
            
            //compute the median
            double median = std::pow( 10., find_quantile(&(*h_log), 0.5) ); 

            double var  = *b_variance[ind] / ((double)count); 
            double var2 = *b_variance2[ind] / ((double)count); 

            double stddev = std::sqrt( var2 - (var*var) ); 

            pts_median[ind] = median; 
            pts_stddev[ind] = stddev; 

            std::printf("block id: %3i (row-%02i, col-%02i): median, stddev: (of samp. variance): %.4e %.4e\n",
                ind, row, col, 
                median, stddev
            ); 
            */ 
            
            //save each 'file' to a pdf
            ++page; 
            switch (page) {

                //first page 
                case 1 : 
                    canv->Print(std::string{path_out_graphic + "["}.c_str(), "pdf"); break; 
                
                //last page
                case (HCal::n_blocks) : 
                    canv->Print(std::string{path_out_graphic + "]"}.c_str(), "pdf"); break; 

                //any other page
                default : 
                    canv->Print(path_out_graphic.c_str(), "pdf"); 
            }
        } 
    }

    canv = new TCanvas; 
    canv->SetLogx(1); 
    canv->SetLogy(1); 
    auto g = new TGraph(HCal::n_blocks, pts_median.data(), pts_stddev.data()); 
    g->SetTitle("Median variance of each block vs. Stddev;Median variance of waveforms (per-block);Sttdev variance of waveforms (per-block)");
    g->SetMarkerStyle(kOpenCircle); 
    g->Draw("AP"); 


    return; 

}