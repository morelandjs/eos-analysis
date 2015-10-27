#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "urqmd_reader.h"

using namespace std;

urqmd_reader::urqmd_reader(const char* filename){

    // read from gzipped file
    ifstream file(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istreambuf inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
 
    //Convert streambuf to istream
    std::istream event(&inbuf);
    string line;

    // read particles from the event
    while(getline(event,line)){   
    
        // skip headers
        if(line.length() != 307) continue;
        replace(line.begin(), line.end(), 'D', 'E');

        // properties
        int typ_ = stoi(line.substr(145, 11));
        int chg_ = stoi(line.substr(159, 3));
    
        // four-position
        double t_ = stod(line.substr(180, 16));
        double x_ = stod(line.substr(196, 16));
        double y_ = stod(line.substr(212, 16));
        double z_ = stod(line.substr(228, 16));

        // four momentum
        double e_ = stod(line.substr(244, 16));
        double px_ = stod(line.substr(260, 16));
        double py_ = stod(line.substr(276, 16));
        double pz_ = stod(line.substr(292, 16));
    
        // pseudorapidity
        double p_ = pow(px_*px_ + py_*py_ + pz_*pz_, 0.5);
        double eta_ = 0.5*log((p_ + pz_)/(p_ - pz_));

        // apply particle cuts
        if(typ_ != 101 || abs(eta_) > 1.5 || chg_ == 0) continue;
        
        // store four-position and four-momentum
        vector<double> x = {t_, x_, y_, z_};
        vector<double> p = {e_, px_, py_, pz_};

        // add particle to list
        xlist.push_back(x);
        plist.push_back(p);
    }
}
