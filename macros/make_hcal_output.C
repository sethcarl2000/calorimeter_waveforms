
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
// analysis_utils
#include <analysis_utils/ROOT.hpp>
// stdlib
#include <string> 
#include <stdexcept> 
#include <array> 
#include <cmath> 
#include <memory> 
#include <utility> 
#include <iostream> 

void make_hcal_output(
    std::string path_infile="data/*", 
    std::string path_outfile="waveforms.root"
)
{
    using ROOT::VecOps::RVec; 

    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("T", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 
    //now, create the 'hcal' objects

    df    
        .Filter([](const RVec<double>& raw_waveforms){ return raw_waveforms.size() == HCal::n_blocks * HCal::Block::n_waveform_samps; }, {"sbs.hcal.samps"})

        .Define("hcal_data", [](const RVec<double>& raw_waveforms){

            if (raw_waveforms.size() != HCal::n_blocks * HCal::Block::n_waveform_samps) {
                throw std::runtime_error("wrong number of waveform samples!"); 
            }

            HCal data; 

            for (int i=0; i<HCal::n_blocks; i++) {

                int row = std::floor(i/HCal::n_cols); 
                int col = i % HCal::n_cols; 

                for (int t=0; t<HCal::Block::n_waveform_samps; t++) {

                    data.get_block(row,col)(t) = raw_waveforms[(HCal::n_cols*row + col)*HCal::Block::n_waveform_samps + t]; 
                }
            }
            
            return data; 

            }, {"sbs.hcal.samps"})
        
            .Snapshot("T", path_outfile, {"hcal_data"}); 

    return; 
}