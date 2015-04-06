#include "urqmd_reader.h"
#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
    // read particle data for two events
    cout << "before" << endl;    
    urqmd_reader event1(argv[1]);
    cout << "event 1 read" << endl;
    urqmd_reader event2(argv[2]);
    cout << "event 2 read" << endl;

    // retrieve particle data
    //vector<vector<double>> xlist = event1.get_xlist(), plist = event1.get_plist();
    //vector<vector<double>> xlist_mixed = event1.get_plist(), plist_mixed = event2.get_plist();
    
    return 0;
}
