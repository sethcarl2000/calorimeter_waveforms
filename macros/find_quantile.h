#ifndef find_quantile_h
#define find_quantile_h

#include <TH1D.h> 
#include <TAxis.h>
// stdlib
#include <limits> 

/// @brief Given Histogram, find x-value corresponding to given quantile (from lhs)
/// @param h histogram 
/// @param p p-value to stop at (integrating from the left)
/// @return x-value corresponding to this quantile
double find_quantile(TH1D* h, double p) {
    double integral = h->Integral(); 
    auto xax = h->GetXaxis(); 
    double sum=0.; 
    for (int ix=1; ix<=xax->GetNbins(); ix++) {
        double val = h->GetBinContent(ix)/integral; 
        if (val + sum >= p) {
            double step_size = val;
            double overshoot = val + sum - p; 

            double x0 = xax->GetBinCenter(ix-1); 
            double dx = xax->GetBinWidth(1); 

            return x0 + dx*((val - overshoot)/val); 
        }
        sum += val; 
    }

    //something went wrong, if we got here
    return std::numeric_limits<double>::quiet_NaN(); 
}


#endif