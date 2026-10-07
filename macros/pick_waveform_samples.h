#ifndef choose_waveform_samples_h
#define choose_waveform_samples_h

#include <HCal.hpp>
// ROOT 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <TString.h> 
// stdlib
#include <vector> 
#include <string> 
#include <cmath> 
#include <iostream> 
#include <array> 

/// @brief select the first 'n' waveform samples 
/// @param row, col 
/// @param min_log_variance 
/// @param max_log_variance 
/// @param n_waveforms 
/// @param overshoot 
/// @param path_infile 
/// @return 
std::vector<std::array<double,HCal::Block::n_waveform_samps>> pick_waveform_samples(
    int row, int col,
    double min_log_variance=-4.40, double max_log_variance=-3.60, 
    std::size_t n_waveforms=40, std::size_t overshoot=3e4, 
    std::string path_infile="waveforms.root"
)
{
    using Samp = std::array<double, HCal::Block::n_waveform_samps>; 

    using ROOT::VecOps::RVec; 
    
    bool is_MT_enabled = ROOT::IsImplicitMTEnabled(); 

    if (is_MT_enabled) ROOT::DisableImplicitMT(); 

    ROOT::RDataFrame df("T", path_infile);  

    const double minval = 1e-25; 

    auto waveform_variance = [](const Samp& v) {
        static constexpr double norm_fact = 1./((double)HCal::Block::n_waveform_samps); 

        double sum{0.}, sum2{0.};
        for (const auto& xi : v) { sum += xi; sum2 += xi*xi; }
        sum  *= norm_fact; 
        sum2 *= norm_fact;
        return sum2 - sum*sum;  
    };

    std::size_t count=0; 

    auto waveforms = df
    
        .Range(0, n_waveforms + overshoot)

        /*.Filter([row,col,minval](const HCal& d){ 
            for (const auto& xi : d.get_block(row,col).samps) {
                if (std::fabs(xi) < minval) return false; 
            }
            return true; 
        }, {"hcal_data"})*/ 

        .Filter([&count,n_waveforms](){ ++count; return (count <= n_waveforms); }, {})

        .Define("samps", [row,col,&count](const HCal& d){

                return d.get_block(row,col).samps;        
            }, {"hcal_data"})

        .Filter([min_log_variance,max_log_variance,&waveform_variance](const Samp& samp){
            double log_var = std::log10( waveform_variance(samp) ); 
            return (log_var >= min_log_variance && log_var < max_log_variance); 
        }, {"samps"})

        .Take<Samp>("samps"); 
    
    std::vector<Samp> pts_kept( 
        waveforms->cbegin(), 
        waveforms->cbegin() + std::min( n_waveforms, waveforms->size() ) 
    ); 

    // re-enable implicit MT, if it was enabled. 
    if (is_MT_enabled) ROOT::EnableImplicitMT(); 

    if (pts_kept.size() < n_waveforms) {
        if (pts_kept.size() < 0) {
            throw std::runtime_error("in <"+std::string{__func__}+">: no relevant waveform samples found!"); 
            return {}; 
        }
        Error(__func__, "Could not find %zi waveform samples as requested, only found %zi", n_waveforms, pts_kept.size());
    }   

    return pts_kept; 
}


#endif