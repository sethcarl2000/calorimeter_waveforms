#include "find_quantile.h"
#include <HCal.hpp>
// ROOT
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TH1D.h> 
#include <ROOT/RVec.hxx>
#include <TString.h> 
#include <TLine.h> 
#include <TCanvas.h> 
// stdlib
#include <iostream> 
#include <fstream> 
#include <cmath> 
#include <string> 
#include <array> 

template<typename T> using block_rptr = std::array<ROOT::RDF::RResultPtr<T>, HCal::n_blocks>; 


void find_median_block_stddev(std::string path_infile="waveforms.root", std::string path_outfile="data/csv/block_medians.csv", std::string path_out_graphic="plots/waveform_samps")
{
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("T", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 

    using ROOT::VecOps::RVec; 

    block_rptr<TH1D> h_var; 
    block_rptr<TH1D> histos; 

    //there are a lot of samples which appear to have null-values near 0. i'm going to cut any waveform which has these. 
    const double minimum_sample_amplitude = 1e-15; 


    std::cout << "booking histograms..." << std::flush;     
    
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {

            const int ind = row*HCal::n_cols + col; 

            double norm_fact = 1./((double)HCal::Block::n_waveform_samps); 

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

                .Define("samp_log_variance", [row,col,norm_fact](const HCal& d){
                    double xx{0.}, x{0.}; 
                    for (double xi : d.get_block(row,col).samps) { xx += xi*xi; x += xi; }
                    xx *= norm_fact; 
                    x  *= norm_fact; 
                    return std::log10( xx - x*x );
                }, {"hcal_data"});   

            histos[ind] = df_samps   
                .Histo1D<RVec<double>>({"h_draw", 
                    Form("#splitline{All waveforms amplitude samples}{row %4i, col %4i};waveform amplitude;",row,col), 
                    200, -0.03, +0.03}, "my_samps"); 

            h_var[ind] = df_samps   

                .Histo1D<double>({"h_var", 
                    Form("#splitline{Per-Event Waveform variance}{row-%i, col-%i};log_{10}( <x_{i}^{2}> - <x_{i}>^{2} )",row,col), 
                    200, -8, 1.}, 
                    "samp_log_variance"); 

        }
    }
    std::cout << "done.\n"; 

    auto count = *df.Count(); 


    auto canv = new TCanvas; 

    std::fstream outfile(path_outfile, std::ios::out | std::ios::trunc); 

    std::string path_outgraphic_var = path_out_graphic + "-variance.pdf"; 
    std::string path_outgraphic_val = path_out_graphic + "-median.pdf"; 

    outfile << 
        "\n"
        "// this is a list of median variances computed for all waveforms.\n"
        "// it was computer using '" << path_infile << "',\n"
        "// and the plots (with marked medians and 99%% quantiles) are under '" << path_out_graphic << "'\n\n"; 

    outfile <<  
        "\n"
        "struct samp_median_stddev_t { double med_val, med_stddev; };\n"
        "\n"; 

    outfile << "static const std::vector<samp_median_stddev_t> waveform_median_variances{\n"; 
    
    int page=0; 
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {  

            //save each 'file' to a pdf
            ++page; 
            auto save_canv_as_pdf = [canv,&page](std::string path){
                switch (page) {

                    //first page 
                    case 1 : 
                        canv->Print(std::string{path + "["}.c_str(), "pdf"); break; 
                    
                    //last page
                    case (HCal::n_blocks) : 
                        canv->Print(std::string{path + "]"}.c_str(), "pdf"); break; 

                    //any other page
                    default : 
                        canv->Print(path.c_str(), "pdf"); 
                }
            }; 

            const int ind = row*HCal::n_cols + col; 


            //if (!((row==7) && (col==6))) continue; 
                
            canv->Clear(); 
            canv->SetTopMargin(0.15); 

            auto& hvar = h_var[ind]; 

            hvar->DrawCopy(); 
            
            //compute the median stddev. 
            double median_stddev = std::pow( 10., find_quantile(&(*hvar), 0.5)/2. ); 
            double quant99 = std::pow( 10., find_quantile(&(*hvar), 0.99) ); 

            // draw median
            TLine *line; 
            line = new TLine(std::log10(median_stddev),0., std::log10(median_stddev),hvar->GetMaximum()); 
            line->SetLineStyle(kSolid); 
            line->SetLineColor(kRed); 
            line->Draw(); 
            
            line = new TLine(std::log10(quant99),0., std::log10(quant99),hvar->GetMaximum()); 
            line->SetLineStyle(kDashed); 
            line->SetLineColor(kBlack); 
            line->Draw(); 

            canv->Modified(); 
            canv->Update(); 

            save_canv_as_pdf(path_outgraphic_var); 


            canv->Clear(); 
            canv->SetTopMargin(0.15); 

            auto& hval = histos[ind]; 

            hval->DrawCopy(); 

            //compute the median
            double median_val = find_quantile(&(*hval), 0.50); 
            quant99    = find_quantile(&(*hval), 0.99); 


            // draw median

            line = new TLine(median_val,0., median_val,hval->GetMaximum()); 
            line->SetLineStyle(kSolid); 
            line->SetLineColor(kRed); 
            line->Draw(); 
            
            line = new TLine(quant99,0.,    quant99,hval->GetMaximum()); 
            line->SetLineStyle(kDashed); 
            line->SetLineColor(kBlack); 
            line->Draw(); 


            canv->Modified(); 
            canv->Update(); 
            //if (!((row==7) && (col==6))) continue; 
                
            save_canv_as_pdf(path_outgraphic_val); 
            
            
            char buff[100];
            if (ind < HCal::n_blocks-1) {
                std::sprintf(buff, "    {%+.8e, %+.8e},\n", median_val, median_stddev);
            } else {
                std::sprintf(buff, "    {%+.8e, %+.8e}\n", median_val, median_stddev); 
            }  
            outfile << buff; 
        } 
    }
    outfile.close(); 
    
}