#ifndef HCal_hpp
#define HCal_hpp

#include <TObject.h> 
// stdlib
#include <array> 

struct HCal {

    static constexpr int n_rows = 24; 
    static constexpr int n_cols = 12; 

    static constexpr int n_blocks = n_rows*n_cols; 

    struct Block {

        static constexpr int n_waveform_samps = 50; 
        std::array<double, n_waveform_samps> samps; 

        double& operator()(int t) { return samps[t]; }
        double operator()(int t) const { return samps[t]; }
    }; 

    std::array<Block, n_blocks> blocks; 

    //row-major block indexing 
    Block& get_block(int row, int col) { return blocks[row*n_cols + col]; } 
    const Block& get_block(int row, int col) const { return blocks[row*n_cols + col]; } 
}; 


#endif