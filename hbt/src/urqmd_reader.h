#ifndef URQMD_READER_H
#define URQMD_READER_H

#include <vector>

using namespace std;

class urqmd_reader
{
    private:
        vector<vector<double>> xlist;
        vector<vector<double>> plist;
    public:
        urqmd_reader(char filename[]);
        vector<vector<double>> get_xlist(){return xlist;}
        vector<vector<double>> get_plist(){return plist;}
};

#endif
