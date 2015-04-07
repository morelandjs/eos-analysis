#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "urqmd_reader.h"

using namespace std;

urqmd_reader::urqmd_reader(char filename[]){

    // open the file
    ifstream event(filename, ifstream::in);
    string line;

    // read particles from the event
    while(getline(event,line)){        

        // skip headers
        if(line.length() != 434) continue;
        replace(line.begin(), line.end(), 'D', 'E');

        // properties
        int typ_ = stoi(line.substr(217, 4));
        int chg_ = stoi(line.substr(225, 2));
    
        // four-position
        double t_ = stod(line.substr(243, 23));
        double x_ = stod(line.substr(267, 23));
        double y_ = stod(line.substr(291, 23));
        double z_ = stod(line.substr(315, 23));

        // four momentum
        double e_ = stod(line.substr(339, 23));
        double px_ = stod(line.substr(363, 23));
        double py_ = stod(line.substr(387, 23));
        double pz_ = stod(line.substr(411, 23));
    
        // pseudorapidity
        double p_ = pow(px_*px_ + py_*py_ + pz_*pz_, 0.5);
        double eta_ = 0.5*log((p_ + pz_)/(p_ - pz_));

        // apply particle cuts
        if(typ_ != 101 || abs(eta_) > 2.5) continue;

        // store four-position and four-momentum
        vector<double> x = {t_, x_, y_, z_};
        vector<double> p = {e_, px_, py_, pz_};

        // add particle to list
        xlist.push_back(x);
        plist.push_back(p);
    }
}
