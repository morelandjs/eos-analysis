#include "urqmd_reader.h"
#include <iostream>
#include <algorithm>
using namespace std;

vector<double> boost_z(vector<double> v_, double beta_z)
{
    double t_ = v_[0], x_ = v_[1], y_ = v_[2], z_ = v_[3];
    double gamma = 1./pow(1.-beta_z*beta_z, 0.5);
    double t = gamma*(t - beta_z*z_);
    double x = x_;
    double y = y_;
    double z = gamma*(z_ - beta_z*t_);
    
    vector<double> v{t, x, y, z};
    
    return v;
}

int main(int argc, char **argv)
{
    // fix q_out, q_side, q_long, kt partitions
    double kt_min = 0.0, kt_max = 1.0, dkt = 0.05; int nkt = (kt_max - kt_min)/dkt;
    double qo_min = -3.0, qo_max = 3.0, dqo = 0.5; int nqo = (qo_max - qo_min)/dqo;
    double qs_min = -0.5, qs_max = 0.5, dqs = 0.1; int nqs = (qs_max - qs_min)/dqs;
    double ql_min = -3.0, ql_max = 3.0, dql = 0.5; int nql = (ql_max - ql_min)/dql;

    // initialize correlation arrays
    double Csame[nkt][nqo][nqs][nql];
    double Cmixed[nkt][nqo][nqs][nql];

    // read particle data for two events
    urqmd_reader event1(argv[1]);
    urqmd_reader event2(argv[2]);

    // retrieve particle data
    vector<vector<double>> x1v = event1.get_xlist(), p1v = event1.get_plist();
    vector<vector<double>> x2v = event2.get_plist(), p2v = event2.get_plist();

    // calculate same event correlations
    for(int i = 0; i < x1v.size(); i++){
        for(int j = i + 1; j < x1v.size(); j++){
	    
	        // select particle pair
	        vector<double> x1 = x1v[i], p1 = p1v[i];
	        vector<double> x2 = x1v[j], p2 = p1v[j];
            
	        // boost into LCM frame
	        double beta_z = (p1[3] + p2[3])/(p1[0] + p2[0]);
	        x1 = boost_z(x1, beta_z);
	        x2 = boost_z(x2, beta_z);
	        p1 = boost_z(p1, beta_z);
	        p2 = boost_z(p2, beta_z);

            // calculate vector dx
            vector<double> dx{x2[1] - x1[1], x2[2] - x1[2], x2[3] - x1[3]}; 

            // calculate q and k vectors
            vector<double> q{p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]};
            vector<double> k{(p1[1] + p2[1])/2, (p1[2] + p2[2])/2, (p1[3] + p2[3])/2};

            // calculate kt
            double kt = pow(k[0]*k[0] + k[1]*k[1], 0.5);
            
            // project q_out, q_side, q_long
            double theta = atan2(q[1], q[0]);
            vector<double> r_hat{cos(theta), sin(theta), 0};
            vector<double> theta_hat{-sin(theta), cos(theta), 0};
            vector<double> z_hat{0, 0, 1};
            double qo = inner_product(q.begin(), q.end(), z_hat.begin(), 0.);
            double qs = inner_product(q.begin(), q.end(), theta_hat.begin(), 0.);
            double ql = inner_product(q.begin(), q.end(), r_hat.begin(), 0.);

            //cout << "kt: " << kt << endl;
            //cout << "qo: " << qo << endl;
            //cout << "qs: " << qs << endl;
            //cout << "ql: " << ql << endl << endl;
    
            // test if q vector is off grid
            if(kt < kt_min || qo < qo_min || qs < qs_min || ql < ql_min) continue;
            if(kt > kt_max || qo > qo_max || qs > qs_max || ql > ql_max) continue;
            
            // calculate ikt, iqs, iqo, iql
            int ikt = (kt - kt_min)/dkt;
            int iqo = (qo - qo_min)/dqo;
            int iqs = (qs - qs_min)/dqs;
            int iql = (ql - ql_min)/dql;

            // add correlations
            cout << ikt << "\t" << iqo << "\t" << iqs << "\t" << iql << endl;
            Csame[ikt][iqo][iqs][iql] += 1.; // + cos(inner_product(q.begin(), q.end(), dx.begin(), 0.));
        }
    }

    return 0;
}
