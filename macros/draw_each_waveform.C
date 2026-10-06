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
#include <iostream> 

namespace {
    constexpr int n_samples = 50; 
    constexpr int n_cols = 12; 
    constexpr int n_rows = 24; 

    constexpr double min_amplitude = 1e-25; 
}; 

struct Waveform {
    std::array<double, n_samples> data;
    double& operator()(int i) { return data[i]; } 
    const double& operator()(int i) const { return data[i]; }
}; 

struct HCalEvent {
    std::array<Waveform, n_rows*n_cols> blocks; 

    Waveform& get_block(int row, int col) { return blocks[row*n_cols + col]; }
    const Waveform& get_block(int row, int col) const { return blocks[row*n_cols + col]; }
}; 

void draw_each_waveform(
    std::string path_infile="data/e1209016_replayed_2034_stream0_2_seg0_0_firstevent0_nevent10000.root", 
    std::string path_out_graphic="plots/waveforms_test"
)
{
    using ROOT::VecOps::RVec; 



    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("T", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 
    //now, create the 'hcal' objects

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
}