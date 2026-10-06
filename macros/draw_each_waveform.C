
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


double max_variance_subrange(const HCal::Block& b, int window_size=8) 
{
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


double get_median(TH1D* h) {
    double integral = h->Integral(); 
    auto xax = h->GetXaxis(); 
    double sum=0.; 
    for (int ix=1; ix<=xax->GetNbins(); ix++) {
        double val = h->GetBinContent(ix)/integral; 
        if (val + sum >= 0.5) {
            double step_size = val;
            double overshoot = val + sum - 0.5; 

            double x0 = xax->GetBinCenter(ix-1); 
            double dx = xax->GetBinWidth(1); 

            return x0 + dx*((val - overshoot)/val); 
        }
        sum += val; 
    }

    //something went wrong, if we got here
    return std::numeric_limits<double>::quiet_NaN(); 
}

template<typename T> using block_rptr = std::array<ROOT::RDF::RResultPtr<T>, HCal::n_blocks>; 

void draw_each_waveform(
    std::string path_infile="waveforms.root", 
    std::string path_out_graphic="plots/waveform_log_variance"
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

    block_rptr<double> b_variance2, b_variance; 

    const int subrange_size = 4; 

    path_out_graphic += ".pdf"; 
    std::cout << "booking histograms..." << std::flush; 
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {
        
            double norm_fact = 1./((double)HCal::Block::n_waveform_samps); 

            const int ind = row*HCal::n_cols + col; 
    
            auto df_samps = df
                
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

                    .Define("log_subrange_variance", [](double var){ return std::log10(var); }, {"subrange_variance"})

                    .Define("samp_log_variance", [](double var){ return std::log10(var); }, {"samp_variance"})

                    .Define("samp_variance2", [](double var){ return var*var; }, {"samp_variance"}); 

            b_variance[ind] = df_samps.Sum<double>("samp_variance");
            b_variance2[ind] = df_samps.Sum<double>("samp_variance2");  

            histos[ind] = df_samps   
                .Histo1D<RVec<double>>({"h_draw", Form("All waveforms for row %4i, col %4i;waveform amplitude;",row,col), 200, -2.5e-2, +2.5e-2}, "my_samps"); 
                    
            h_subrange_log_variance[ind] = df_samps
                .Histo1D<double>({"h_log_subrange_variance", Form(
                    "Max-window Variance of waveform samples: row %4i, col %4i;log_{10}( max[ <x_{i}^{2}> - <x_{i}>^{2} ] );",row,col), 
                    100, -8, 1}, "log_subrange_variance");
            
            h_log_variance[ind] = df_samps
                .Histo1D<double>({"h_log_variance", 
                    Form("Variance of waveform samples: row %4i, col %4i;log_{10}( max[ <x_{i}^{2}> - <x_{i}>^{2} ] );",row,col), 
                    100, -8, 1}, "samp_log_variance");
            
            h_variance[ind] = df_samps
                .Histo1D<double>({"h_variance", 
                    Form("Max-window Variance of waveform samples: row %4i, col %4i;max[ <x_{i}^{2}> - <x_{i}>^{2} ];",row,col), 
                    100, 0., 5e-5}, "log_subrange_variance");
        }
    }
    std::cout << "done.\n"; 

    auto count = *df.Count(); 

    std::printf("done fetching %.4e events.\ndrawing...\n", (double)count); 

    auto canv = new TCanvas; 
    
    std::vector<double> pts_median(HCal::n_blocks), pts_stddev(HCal::n_blocks); 

    int page=0; 
    for (int row=0; row<HCal::n_rows; row++) {
        for (int col=0; col<HCal::n_cols; col++) {

            const int ind = row*HCal::n_cols + col; 

            auto& h = h_log_variance[ind]; 

            auto& h_log = h_log_variance[ind]; 

            canv->Clear(); 
            canv->SetLogy(1); 

            h->DrawCopy(); 
            
            canv->Modified(); 
            canv->Update(); 

            //compute the median
            double median = std::pow( 10., get_median(&(*h_log)) ); 

            double var  = *b_variance[ind] / ((double)count); 
            double var2 = *b_variance2[ind] / ((double)count); 

            double stddev = std::sqrt( var2 - (var*var) ); 

            pts_median[ind] = median; 
            pts_stddev[ind] = stddev; 

            std::printf("block id: %3i (row-%02i, col-%02i): median, stddev: (of samp. variance): %.4e %.4e\n",
                ind, row, col, 
                median, stddev
            ); 

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

#if 0 
    auto events = *df    

        .Filter([](const RVec<double>& raw_waveforms){ 
            
            return raw_waveforms.size() == n_cols*n_rows*n_samples; 
            
            }, {"sbs.hcal.samps"})

        .Define("HCal_event", [](const RVec<double>& raw_waveforms){

            if (raw_waveforms.size() != n_cols*n_rows*n_samples) {
                throw std::runtime_error("wrong number of waveform samples!"); 
            }

            HCalEvent evt; 

            for (int i=0; i<n_rows*n_cols; i++) {

                int row = std::floor(i/n_cols); 
                int col = i % n_cols; 

                for (int t=0; t<n_samples; t++) {

                    evt.get_block(row,col)(t) = raw_waveforms[(n_cols*row + col)*n_samples + t]; 
                }
            }
            
            return evt; 

            }, {"sbs.hcal.samps"})

        .Take<HCalEvent>("HCal_event"); 

    std::cout << "Number of events: " << events.size() << "\n"; 
    
    using namespace analysis_utils; 

    auto canv = new TCanvas; 

    auto hist = new TH1D("h_draw", "All waveforms for block;waveform amplitude;", 100, -2.5e-2, +2.5e-2); 

    const double dx = hist->GetXaxis()->GetBinWidth(1);     
    const double xmin{hist->GetXaxis()->GetXmin()}, xmax{hist->GetXaxis()->GetXmax()}; 

    path_out_graphic += ".pdf"; 


    // fit the histogram with a gaussian fcn
    using ROOT::Math::normal_cdf; 
    auto gfcn = new TF1("f_gaus", [](const double* X, const double* par){
        double A = par[0]; 
        double sigma = std::fabs(par[2]);

        double arg = (X[0] - par[1])/sigma; 

        return A * std::exp( -0.5*arg*arg ); 
    }, xmin,xmax, 3); 


    auto c_rms = new TCanvas; 
    auto hist_rms = new TH2D("h_rms", "RMS of each waveform;col;row;Amplidute noise RMS", 
        n_cols, -0.5, ((double)n_cols)-0.5,
        n_rows, -0.5, ((double)n_rows)-0.5
    ); 
    hist_rms->GetXaxis()->SetMaxDigits(3); 
    c_rms->cd(); hist_rms->Draw("colz"); 

    auto c_offset = new TCanvas; 
    auto hist_offset = new TH2D("h_offset", "Offset of each waveform;col;row;Amplitude noise offset from Pedistal", 
        n_cols, -0.5, ((double)n_cols)-0.5,
        n_rows, -0.5, ((double)n_rows)-0.5
    ); 
    hist_offset->GetXaxis()->SetMaxDigits(3); 
    c_offset->cd(); hist_offset->Draw("colz"); 

    std::vector<double> pts_rms, pts_offset, pts_rms_err, pts_offset_err; 

    int page=0; 
    for (int row=0; row<n_rows; row++) { 
        for (int col=0; col<n_cols; col++) {

            //
            std::printf("drawing row %4i col %4i\n", row, col); 

            auto h = root::unique_TObject_copy(hist, Form("h_row%i_col%i",row,col)); 
            
            for (const auto& evt : events)
                for (const auto& wf : evt.get_block(row,col).data) { if (std::fabs(wf) > min_amplitude) h->Fill( wf ); }  
        
            canv->cd(); 
            canv->Clear(); 

            h->SetTitle(Form("All Waveforms: row %i, col %i;Amplitude",row,col)); 
            h->Draw(); 

            double sigma = h->GetRMS(); 
            gfcn->SetParameter(1, h->GetMean());
            gfcn->SetParameter(2, sigma);
            gfcn->SetParameter(0, h->GetMaximum()); //h->GetMaximum()/(normal_cdf(dx,sigma)-0.5));

            auto fitresult = h->Fit(gfcn, "L R S Q N"); 

            double rms, offset; 

            if (!fitresult->IsValid()) { 
                Error(__func__, "Something went wrong with gauss fit");

                auto box = new TBox(col-0.5,row-0.5, col+0.5,row+0.5); 
                box->SetFillStyle(3004); 
                box->SetLineColor(kRed); 
                box->SetFillColor(kRed); 

                c_offset->cd(); 
                box->Draw("SAME");
                
                c_rms->cd(); 
                box->Draw("SAME"); 

            } else {

            
                offset = fitresult->Parameter(1); 
                rms    = fitresult->Parameter(2); 

                pts_rms.emplace_back(rms); 
                pts_offset.emplace_back(offset); 

                pts_rms_err.emplace_back(fitresult->ParError(2)); 
                pts_offset_err.emplace_back(fitresult->ParError(1)); 

                std::printf(" row %3i col %i    rms: %+.3e    offset: %+.3e\n", row,col, rms, offset); 

                hist_rms->Fill( col, row, rms ); 
                hist_offset->Fill( col, row, offset ); 
            }
            canv->cd(); 
            gfcn->Draw("SAME"); 



            canv->Modified();
            canv->Update(); 

            ++page; 
            switch (page) {

                //first page 
                case 1 : 
                    canv->Print(std::string{path_out_graphic + "["}.c_str(), "pdf"); break; 
                
                //last page
                case (n_cols*n_rows) : 
                    canv->Print(std::string{path_out_graphic + "]"}.c_str(), "pdf"); break; 

                //any other page
                default : 
                    canv->Print(path_out_graphic.c_str(), "pdf"); 
            }
        }
    }

    gStyle->SetOptStat(0); 
    gStyle->SetPalette(kGreyScale); 
    TColor::InvertPalette(); 

    hist_rms->SetMinimum(0.); 
    c_rms->SetRightMargin(0.15); 
    c_rms->Modified(); 
    c_rms->Update(); 

    hist_offset->SetMinimum(0.); 
    c_offset->SetRightMargin(0.15); 
    c_offset->Modified(); 
    c_offset->Update(); 

    canv = new TCanvas;
    
    canv->SetLeftMargin(0.15);  
    auto g = new TGraphErrors(pts_rms.size(), pts_rms.data(), pts_offset.data(), pts_rms_err.data(), pts_offset_err.data()); 
    g->SetTitle("Waveform noise for all blocks;Waveform noise RMS;Waveform noise Offset");
    g->SetMarkerStyle(kPlus); 
    g->SetMarkerSize(0.75);
    g->GetXaxis()->SetNdivisions(8); 
    g->GetYaxis()->SetNdivisions(8); 
    canv->SetLogx(1); canv->SetLogy(1); 
    g->Draw("APZ");

    canv = new TCanvas; 
    auto h_1d_offset = new TH1D("h_1d_offset", "Waveform noise offset from pedistal;offset x 10^{3};", 50, 0., 2.); 
    for (auto offset : pts_offset) h_1d_offset->Fill( offset*1e3 ); 
    h_1d_offset->GetXaxis()->SetNdivisions(4); 
    h_1d_offset->Draw(); 

    canv = new TCanvas; 
    auto h_1d_rms = new TH1D("h_1d_rms", "Waveform noise rms from pedistal;rms x 10^{3};", 50, 0., 8.); 
    for (auto rms : pts_rms) h_1d_rms->Fill( rms*1e3 ); 
    h_1d_rms->GetXaxis()->SetNdivisions(4); 
    h_1d_rms->Draw(); 


    //output the larget waveform in the first hist as a histogram
    int evt_num =0; 
    auto& evt = events[evt_num]; 

    double highest=-1e30; 
    Waveform *wf_highest=nullptr; 
    int irow,icol; 
    for (int row=0; row<n_rows; row++) {
        for (int col=0; col<n_cols; col++) {
            
            double highest_this_block=-1e30;     
            for (const auto& val : evt.get_block(row,col).data) {
                highest_this_block = std::max(val, highest_this_block); 
            }
            if (highest_this_block > highest) {
                wf_highest = &evt.get_block(row,col); 
                irow=row; 
                icol=col;
                highest = highest_this_block;  
            }
        }
    }

    auto hist_highest = new TH1D("h_wf", Form("Larget-amplitude wavefor for event %i;time index;amplitude",evt_num), n_samples, -0.5, ((double)n_samples)-0.5); 
    for (int t=0; t<n_samples; t++) 
        hist_highest->Fill( t, (*wf_highest)(t) ); 

    new TCanvas; 
    hist_highest->SetMarkerStyle(kPlus); 
    hist_highest->Draw("PL"); 

    std::printf("evt. num highest block: row %3i     col %3i\n", irow,icol); 
    return;
#endif

}