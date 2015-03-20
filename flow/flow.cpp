#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <complex>
using namespace std;

// bool for POI (particles of interest)
bool POI(double itype, double pt, double ptmin, double ptmax){
  // set particle ID for POI
  if(abs(int(itype)) == 101 && ptmin <= pt && pt < ptmax) return true;
  else return false;
}

// bool for RFP (reference frame particles)
bool RFP(double itype){
  // set particle ID for RFP
  if(int(itype) == 101) return true;
  else return false;
}

// flow data needed to characterize each event
struct flowdata
{
  complex<double> two;            // <2 >
  complex<double> two_weight;     // <2 > event weights
  complex<double> four;           // <4 >
  complex<double> four_weight;    // <4 > event weights
  vector<complex<double> > twop;  // <2'>
  vector<double>  twop_weight;    // <2'> event weights
  vector<complex<double> > fourp; // <4'>
  vector<double> fourp_weight;    // <4'> event weights
};

//###########################################################################
int main()
{ 
  // set harmonic to calculate
  double mode = 3.;

  // pt range and step size
  double ptmin=0.0, ptmax = 2.5, dpt=0.25;
  int npt = int((ptmax-ptmin)/dpt);

  // open the flow file
  string line; ifstream infile;   
  infile.open("tmp/s95-PCE-dflow-20-30.dat");

  // initialize vars
  double ityp0, pt0, phi0;
  vector<double> ityp, pt, phi;
  vector<flowdata> batch;

  // loop over events
  int nev = 0;
  do{

    // initialize flowdata event
    flowdata event;

    // clear existing particle containers
    ityp.clear();
    pt.clear();
    phi.clear();

    // read event data until line skip
    do{
      getline(infile,line);

      // store particle data
      if(line.length() != 0){
	istringstream(line) >> ityp0 >> phi0 >> pt0;
	ityp.push_back(ityp0);
	pt.push_back(pt0);
	phi.push_back(phi0);
      }

    // stop on line skip
    }while(line.length()!= 0);

    // print event status to terminal
    cout << "nev: " << nev <<  endl;

    // construct Qn and Q2n vectors over all RFP
    double M = 0.0;
    complex<double> Qn(0.0,0.0), Q2n(0.0,0.0);
    for(int i=0; i<ityp.size(); i++){
      if(RFP(ityp[i])){
	complex<double> dQn(cos(1.*mode*phi[i]),sin(1.*mode*phi[i]));
	complex<double> dQ2n(cos(2.*mode*phi[i]),sin(2.*mode*phi[i]));
	Qn += dQn;
	Q2n += dQ2n;
	M += 1.;
      }
    }

    // construct single event <2> vector
    event.two = (norm(Qn) - M)/(M*(M-1.) + 1e-15);
    event.two_weight = M*(M-1.);

    // construct single event <4> vector
    complex<double> A1 = norm(Qn)*norm(Qn) + norm(Q2n) - 2.*real(Q2n*conj(Qn)*conj(Qn));
    complex<double> A2 = -2.*(2.*(M-2.)*norm(Qn) - M*(M-3.));
    event.four = (A1+A2)/(M*(M-1.)*(M-2.)*(M-3.) + 1e-15);
    event.four_weight = (M*(M-1.)*(M-2.)*(M-3.));

    // loop over momentum
    for(int ipt=0; ipt<npt; ipt++){

      // momentum cuts
      double ptlow = ptmin+ipt*dpt, pthigh = ptmin+(ipt+1)*dpt;

      // initialize variables
      double mp = 0., mq = 0.;
      complex<double> pn(0.0,0.0), qn(0.0,0.0), q2n(0.0,0.0);
      for(int i=0; i<ityp.size(); i++){

	// construct pn vector over POI
	if(POI(ityp[i],pt[i],ptlow,pthigh)){
	  complex<double> dpn(cos(1.*mode*phi[i]),sin(1.*mode*phi[i])); 
	  pn += dpn;
	  mp += 1.;
	}
	// construct qn and q2n vectors over POI intersect RFP
	if(POI(ityp[i],pt[i],ptlow,pthigh) && RFP(ityp[i])){
	  complex<double> dqn(cos(1.*mode*phi[i]),sin(1.*mode*phi[i])); 
	  complex<double> dq2n(cos(2.*mode*phi[i]),sin(2.*mode*phi[i])); 
	  qn += dqn;
	  q2n += dq2n;
	  mq += 1.;
	}
      }

      // record <2'> for each event
      event.twop.push_back((pn*conj(Qn) - mq)/(mp*M - mq + 1e-15));
      event.twop_weight.push_back(mp*M - mq);

      // record <4'> for each event
      event.fourp.push_back((pn*Qn*conj(Qn)*conj(Qn) - q2n*conj(Qn)*conj(Qn)
       - pn*Qn*conj(Q2n) - 2.*M*pn*conj(Qn) - 2.*mq*norm(Qn) + 7.*qn*conj(Qn)
       - Qn*conj(qn) + q2n*conj(Q2n) + 2.*pn*conj(Qn) + 2.*mq*M - 6.*mq)
	 /((mp*M-3.*mq)*(M-1.)*(M-2.) + 1e-15));
      event.fourp_weight.push_back(((mp*M-3.*mq)*(M-1.)*(M-2.)));
    }

    // store the event data to the batch vector
    batch.push_back(event);

    // increment event counter
    nev++;

  // end event loop
  }while(!infile.eof()); 

  // close the read file
  infile.close();

  // calculate <<2>> and <<4>> over all events
  complex<double> num2 = 0.0, denom2 = 0.0, num4 = 0.0, denom4 = 0.0;
  for(int iev=0; iev < nev-1; iev++){
    flowdata event = batch[iev];
    num2 += event.two_weight*event.two;
    num4 += event.four_weight*event.four;
    denom2 += event.two_weight;
    denom4 += event.four_weight;
  }

  // return <<2>> and <<4>>
  complex<double> twoavg = num2/denom2, fouravg = num4/denom4; 

  // calculate <<2'>> and <<4'>> over all events
  vector<complex<double> > num2p(npt,0.0), num4p(npt,0.0);
  vector<double> denom2p(npt,0.0), denom4p(npt,0.0);
  for(int iev=0; iev < nev-1; iev++){
    flowdata event = batch[iev];
    for(int ipt=0; ipt<npt; ipt++){
      num2p[ipt] += event.twop_weight[ipt]*event.twop[ipt];
      num4p[ipt] += event.fourp_weight[ipt]*event.fourp[ipt];
      denom2p[ipt] += event.twop_weight[ipt];
      denom4p[ipt] += event.fourp_weight[ipt];
    }
  } 

  // return <<2'>> and <<4'>>
  vector<complex<double> > twopavg(npt,0.0), fourpavg(npt,0.0);
  for(int ipt=0; ipt<npt; ipt++){
    twopavg[ipt] = num2p[ipt]/denom2p[ipt];
    fourpavg[ipt] = num4p[ipt]/denom4p[ipt];
  }

  // ERROR TERMS

  // weight error terms
  complex<double>  W2=0.0, W2sq=0.0, W4=0.0, W4sq=0.0, W2W4=0.0;
  vector<complex<double> > W2p(npt,0.0), W2psq(npt,0.0), W4p(npt,0.0), W4psq(npt,0.0);
  vector<complex<double> > W2pW4p(npt,0.0), W2W2p(npt,0.0), W4W4p(npt,0.0), W2W4p(npt,0.0), W2pW4(npt,0.0);

  // s error terms
  complex<double> s2sqnum=0.0, s4sqnum=0.0;
  vector<complex<double> > s2psqnum(npt,0.0), s4psqnum(npt,0.0);
  
  // covariance terms
  complex<double> W2two(npt,0.0), W4four(npt,0.0), W2W4twofour(npt,0.0);
  vector<complex<double> > W2ptwop(npt,0.0), W4pfourp(npt,0.0);
  vector<complex<double> > W2W2ptwotwop(npt,0.0), W4W4pfourfourp(npt,0.0); 
  vector<complex<double> > W2W4ptwofourp(npt,0.0), W2pW4twopfour(npt,0.0), W2pW4ptwopfourp(npt,0.0);
  
  for(int iev=0; iev < nev-1; iev++){ 
    flowdata event = batch[iev];

    // weight error terms
    W2 += event.two_weight;
    W2sq += event.two_weight*event.two_weight;
    W4 += event.four_weight;
    W4sq += event.four_weight*event.four_weight;
    W2W4 += event.two_weight*event.four_weight;
    W2two += event.two*event.two_weight;
    W4four += event.four*event.four_weight;
    W2W4twofour += event.two_weight*event.four_weight*event.two*event.four;

    // s terms
    s2sqnum += event.two_weight*norm(event.two-twoavg);
    s4sqnum += event.four_weight*norm(event.four-fouravg);
    
    for(int ipt=0; ipt<npt; ipt++){

      // weight error terms
      W2p[ipt] += event.twop_weight[ipt];
      W2psq[ipt] += event.twop_weight[ipt]*event.twop_weight[ipt];
      W4p[ipt] += event.fourp_weight[ipt];
      W4psq[ipt] += event.fourp_weight[ipt]*event.fourp_weight[ipt];
      W2pW4p[ipt] += event.twop_weight[ipt]*event.fourp_weight[ipt];
      W2W2p[ipt] += event.two_weight*event.twop_weight[ipt];
      W4W4p[ipt] += event.four_weight*event.fourp_weight[ipt];
      W2W4p[ipt] += event.two_weight*event.fourp_weight[ipt];
      W2pW4[ipt] += event.twop_weight[ipt]*event.four_weight;

      // s terms
      s2psqnum[ipt] += event.twop_weight[ipt]*norm(event.twop[ipt]-twopavg[ipt]);
      s4psqnum[ipt] += event.fourp_weight[ipt]*norm(event.fourp[ipt]-fourpavg[ipt]);

      // covariance terms
      W2ptwop[ipt] += event.twop_weight[ipt]*event.twop[ipt];
      W4pfourp[ipt] += event.fourp_weight[ipt]*event.fourp[ipt];
      W2W2ptwotwop[ipt] += event.two_weight*event.twop_weight[ipt]*event.two*event.twop[ipt];
      W4W4pfourfourp[ipt] += event.four_weight*event.fourp_weight[ipt]*event.four*event.fourp[ipt];
      W2W4ptwofourp[ipt] += event.two_weight*event.fourp_weight[ipt]*event.two*event.fourp[ipt];
      W2pW4twopfour[ipt] += event.twop_weight[ipt]*event.four_weight*event.twop[ipt]*event.four;
      W2pW4ptwopfourp[ipt] += event.twop_weight[ipt]*event.fourp_weight[ipt]*event.twop[ipt]*event.fourp[ipt];
    }
  }

  // s terms
  complex<double> s2sq = (s2sqnum/W2)/(1.-W2sq/(W2*W2));
  complex<double> s4sq = (s4sqnum/W4)/(1.-W4sq/(W4*W4));
  vector<complex<double> > s2psq(npt,0.0), s4psq(npt,0.0);

  for(int ipt=0; ipt<npt; ipt++){
    s2psq[ipt] += (s2psqnum[ipt]/W2p[ipt])/(1.-W2psq[ipt]/(W2p[ipt]*W2p[ipt]));
    s4psq[ipt] += (s4psqnum[ipt]/W4p[ipt])/(1.-W4psq[ipt]/(W4p[ipt]*W4p[ipt]));
  }

  // covariance terms
  complex<double> cov24 = (W2W4twofour/W2W4-W2two*W4four/(W2*W4))/(1.-W2W4/(W2*W4)); 
  vector<complex<double> > cov22p(npt,0.0), cov44p(npt,0.0);
  vector<complex<double> > cov24p(npt,0.0), cov2p4(npt,0.0), cov2p4p(npt,0.0);

  for(int ipt=0; ipt<npt; ipt++){
    cov22p[ipt] = (W2W2ptwotwop[ipt]/W2W2p[ipt]-W2two*W2ptwop[ipt]/(W2*W2p[ipt]))/(1.-W2W2p[ipt]/(W2*W2p[ipt])); 
    cov44p[ipt] = (W4W4pfourfourp[ipt]/W4W4p[ipt]-W4four*W4pfourp[ipt]/(W4*W4p[ipt]))/(1.-W4W4p[ipt]/(W4*W4p[ipt])); 
    cov24p[ipt] = (W2W4ptwofourp[ipt]/W2W4p[ipt]-W2two*W4pfourp[ipt]/(W2*W4p[ipt]))/(1.-W2W4p[ipt]/(W2*W4p[ipt])); 
    cov2p4[ipt] = (W2pW4twopfour[ipt]/W2pW4[ipt]-W2ptwop[ipt]*W4four/(W2p[ipt]*W4))/(1.-W2pW4[ipt]/(W2p[ipt]*W4)); 
    cov2p4p[ipt] = (W2pW4ptwopfourp[ipt]/W2pW4p[ipt]-W2ptwop[ipt]*W4pfourp[ipt]/(W2p[ipt]*W4p[ipt]))/(1.-W2pW4p[ipt]/(W2p[ipt]*W4p[ipt]));
  }

  // calc v{2} and v{4}
  complex<double> v2 = sqrt(twoavg);
  complex<double> v4 = pow(2.*twoavg*twoavg-fouravg,1./4.);

  // v{2} and v{4} errors
  complex<double> v2err = sqrt((1./(4.*twoavg))*W2sq/(W2*W2)*s2sq);
  complex<double> v4err = sqrt(pow(2.*twoavg*twoavg-fouravg,-3./2.)*(twoavg*twoavg*(W2sq/(W2*W2))*s2sq
			       + (1./16.)*(W4sq/(W4*W4))*s4sq-(1./2.)*twoavg*(W2W4/(W2*W4))*cov24));

  // calculate v'{2} and v'{4}
  vector<complex<double> > v2p(npt,0.0), v4p(npt,0.0);
  for(int ipt=0; ipt<npt; ipt++){
    v2p[ipt] = twopavg[ipt]/sqrt(twoavg);
    v4p[ipt] = -(fourpavg[ipt]-2.*twopavg[ipt]*twoavg)/pow(2.*twoavg*twoavg-fouravg,3./4.);
  }

  // calculate v'{2} and v'{4} errors
  vector<complex<double> > v2perr(npt,0.0), v4perr(npt,0.0);
  for(int ipt=0; ipt<npt; ipt++){
    v2perr[ipt] = sqrt((1./(4.*pow(twoavg,3.)))*(pow(twopavg[ipt],2.)*W2sq/(W2*W2)*s2sq
		  + 4.*pow(twoavg,2.)*W2psq[ipt]/(W2p[ipt]*W2p[ipt])*s2psq[ipt]
		  - 4.*twoavg*twopavg[ipt]*W2W2p[ipt]/(W2*W2p[ipt])*cov22p[ipt]));

    v4perr[ipt] = sqrt(pow(2.*twoavg*twoavg-fouravg,-7./2.)*(pow(2.*twoavg*twoavg*twopavg[ipt]-3.*twoavg*fourpavg[ipt]
                      +2.*fouravg*twopavg[ipt],2.)*(W2sq/(W2*W2))*s2sq+(9./16.)*pow(2.*twoavg*twopavg[ipt]
                      -fourpavg[ipt],2.)*(W4sq/(W4*W4))*s4sq+4.*twoavg*twoavg*pow(2.*twoavg*twoavg-fouravg,2.)
                      *(W2psq[ipt]/(W2p[ipt]*W2p[ipt]))*s2psq[ipt]+pow(2.*twoavg*twoavg-fouravg,2.)*(W4psq[ipt]/(W4p[ipt]*W4p[ipt]))
                      *s4psq[ipt]-(3./2.)*(2.*twoavg*twopavg[ipt]-fourpavg[ipt])*(2.*twoavg*twoavg*twopavg[ipt]
                      -3.*twoavg*fourpavg[ipt]+2.*fouravg*twopavg[ipt])*(W2W4/(W2*W4))*cov24-4.*twoavg*(2.*twoavg*twoavg-fouravg)
                      *(2.*twoavg*twoavg*twopavg[ipt]-3.*twoavg*fourpavg[ipt]+2.*fouravg*twopavg[ipt])*(W2W2p[ipt]/(W2*W2p[ipt]))
                      *cov22p[ipt]+2.*(2.*twoavg*twoavg-fouravg)*(2.*twoavg*twoavg*twopavg[ipt]-3.*twoavg*fourpavg[ipt]
                      +2.*fouravg*twopavg[ipt])*(W2W4p[ipt]/(W2*W4p[ipt]))*cov24p[ipt]+3.*twoavg*(2.*twoavg*twoavg-fouravg)
                      *(2.*twoavg*twopavg[ipt]-fourpavg[ipt])*(W2pW4[ipt]/(W2p[ipt]*W4))*cov2p4[ipt]-(3./2.)
                      *(2.*twoavg*twoavg-fouravg)*(2.*twoavg*twopavg[ipt]-fourpavg[ipt])*(W4W4p[ipt]/(W4*W4p[ipt]))*cov44p[ipt]
                      -4.*twoavg*pow(2.*twoavg*twoavg-fouravg,2.)*(W2pW4p[ipt]/(W2p[ipt]*W4p[ipt]))*cov2p4p[ipt])); 
  }

  // open results file
  ofstream outfile;
  outfile.open("results/s95-PCE-pions-20-30.dat",std::ios_base::trunc);

  // print results to file
  for(int ipt=0; ipt<npt; ipt++){
    outfile << setprecision(12) << setw(22) << ptmin+ipt*dpt+dpt/2. 
	    << setprecision(12) << setw(22) << sqrt(norm(v2p[ipt]))
            << setprecision(12) << setw(22) << sqrt(norm(v2perr[ipt]))
	    << setprecision(12) << setw(22) << sqrt(norm(v4p[ipt])) 
            << setprecision(12) << setw(22) << sqrt(norm(v4perr[ipt])) << endl;
  }

  // print v{2} and v{4} integrated to terminal with errors
  cout << "v{2}: " << "\t" << sqrt(norm(v2)) << "\t" <<  sqrt(norm(v2err)) << endl;
  cout << "v{4}: " << "\t" << sqrt(norm(v4)) << "\t" <<  sqrt(norm(v4err)) << endl;
  
  // close the print file
  outfile.close();

  return 1;
}