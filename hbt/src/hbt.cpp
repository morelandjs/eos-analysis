#include "urqmd_reader.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "boost/multi_array.hpp"
#include "H5Cpp.h"
using namespace std;
using namespace H5;

// Lorentz boost along the z-axis (beam direction)
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

// write output using hdf5 writer
void write_correlations(const boost::multi_array<double, 4>& data, char filename[])
{
    // filename and dataset name
    const H5std_string FILE_NAME(filename);
    const H5std_string DATASET_NAME("correlations");

    // create new file using default property lists
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    
    // data set dimensions
    hsize_t dims[4] = {data.shape()[0], data.shape()[1], data.shape()[2], data.shape()[3]};
    DataSpace dataspace(4, dims);

    // create the dataset
    DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace, H5P_DEFAULT);

    // write to file
    dataset.write(data.data(), PredType::NATIVE_DOUBLE);
}

// main function: calculates hbt correlations
int main(int argc, char **argv)
{
    // fix q_out, q_side, q_long, kt partitions
    double kt_min =  0.0, kt_max = 1.0; const int nkt = 20; double dkt = (kt_max - kt_min)/nkt;
    double qo_min = -2.0, qo_max = 2.0; const int nqo = 40; double dqo = (qo_max - qo_min)/nqo;
    double qs_min = -2.0, qs_max = 2.0; const int nqs = 40; double dqs = (qs_max - qs_min)/nqs;
    double ql_min = -2.0, ql_max = 2.0; const int nql = 40; double dql = (ql_max - ql_min)/nql;

    // initialize correlation arrays
    boost::multi_array<double, 4> Csame(boost::extents[nkt][nqo][nqs][nql]);
    boost::multi_array<double, 4> Cmixed(boost::extents[nkt][nqo][nqs][nql]);
    double n_same = 0, n_mixed = 0;

    // read particle data for two events
    urqmd_reader event1(argv[1]);
    urqmd_reader event2(argv[2]);

    // retrieve particle data
    vector<vector<double>> x1v = event1.get_xlist(), p1v = event1.get_plist();
    vector<vector<double>> x2v = event2.get_plist(), p2v = event2.get_plist();

    // loop through particles in given event
    for(int i = 0; i < x1v.size(); i++){

        // loop through partners in the same event
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
            double kt = sqrt(k[0]*k[0] + k[1]*k[1]);
            
            // project q_out, q_side, q_long
            vector<double> rvec{(x1[1] + x2[1])/2, (x1[2] + x2[2])/2, (x1[3] + x2[3])/2};
            double r = sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1]);
            vector<double> r_hat{rvec[0]/r, rvec[1]/r, 0};
            vector<double> theta_hat{-rvec[1]/r, rvec[0]/r, 0};
            vector<double> z_hat{0, 0, 1};
            double qo = inner_product(q.begin(), q.end(), z_hat.begin(), 0.);
            double qs = inner_product(q.begin(), q.end(), theta_hat.begin(), 0.);
            double ql = inner_product(q.begin(), q.end(), r_hat.begin(), 0.);

            // test if q vector is off grid
            if(kt < kt_min || qo < qo_min || qs < qs_min || ql < ql_min) continue;
            if(kt > kt_max || qo > qo_max || qs > qs_max || ql > ql_max) continue;
            
            // calculate ikt, iqs, iqo, iql
            int ikt = (kt - kt_min)/dkt;
            int iqo = (qo - qo_min)/dqo;
            int iqs = (qs - qs_min)/dqs;
            int iql = (ql - ql_min)/dql;

            // add same event correlations
            Csame[ikt][iqo][iqs][iql] += 1. + cos(inner_product(q.begin(), q.end(), dx.begin(), 0.));
            n_same += 1.;
        }

        // loop through partners in a mixed event
        for(int j = 0; j < x2v.size(); j++){

            // select particle pair
	        vector<double> x1 = x1v[i], p1 = p1v[i];
	        vector<double> x2 = x2v[j], p2 = p2v[j];
            
	        // boost into LCM frame
	        double beta_z = (p1[3] + p2[3])/(p1[0] + p2[0]);
	        x1 = boost_z(x1, beta_z);
	        x2 = boost_z(x2, beta_z);
	        p1 = boost_z(p1, beta_z);
	        p2 = boost_z(p2, beta_z);

            // calculate vector dx
            vector<double> xmid{(x1[1] + x2[1])/2, (x1[2] + x2[2])/2, (x1[3] + x2[3])/2};
            vector<double> dx{x2[1] - x1[1], x2[2] - x1[2], x2[3] - x1[3]}; 

            // calculate q and k vectors
            vector<double> q{p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]};
            vector<double> k{(p1[1] + p2[1])/2, (p1[2] + p2[2])/2, (p1[3] + p2[3])/2};

            // calculate kt
            double kt = pow(k[0]*k[0] + k[1]*k[1], 0.5);
            
            // project q_out, q_side, q_long
            vector<double> rvec{(x1[1] + x2[1])/2, (x1[2] + x2[2])/2, (x1[3] + x2[3])/2};
            double r = sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1]);
            vector<double> r_hat{rvec[0]/r, rvec[1]/r, 0};
            vector<double> theta_hat{-rvec[1]/r, rvec[0]/r, 0};
            vector<double> z_hat{0, 0, 1};
            double qo = inner_product(q.begin(), q.end(), z_hat.begin(), 0.);
            double qs = inner_product(q.begin(), q.end(), theta_hat.begin(), 0.);
            double ql = inner_product(q.begin(), q.end(), r_hat.begin(), 0.);
            
            // test if q vector is off grid
            if(kt < kt_min || qo < qo_min || qs < qs_min || ql < ql_min) continue;
            if(kt > kt_max || qo > qo_max || qs > qs_max || ql > ql_max) continue;
            
            // calculate ikt, iqs, iqo, iql
            int ikt = (kt - kt_min)/dkt;
            int iqo = (qo - qo_min)/dqo;
            int iqs = (qs - qs_min)/dqs;
            int iql = (ql - ql_min)/dql;

            // add mixed event correlations
            Cmixed[ikt][iqo][iqs][iql] += 1.;
            n_mixed += 1;
        }
    }

    // normalize by n-pairs
    for(int ikt = 0; ikt < nkt; ikt++){
        for(int iqo = 0; iqo < nqo; iqo++){
           for(int iqs = 0; iqs < nqs; iqs++){
              for(int iql = 0; iql < nql; iql++){
                Csame[ikt][iqo][iqs][iql] /= n_same;
                Cmixed[ikt][iqo][iqs][iql] /= n_mixed;
              }
           }
        }
    }

    // write Csame and Cmixed to hdf5 files 
    write_correlations(Csame, "csame.hd5");
    write_correlations(Cmixed, "cmixed.hd5");
    return 0;
}
