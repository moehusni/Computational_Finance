#include<iostream>
#include<algorithm>
#include<fstream>
#include<vector>
#include <iomanip>
#include<cmath>
using namespace std;


// The cumulative normal function N(x) is given by
double cum_norm(double x){
	const double root = sqrt(0.5);
	return 0.5*(1.0 + erf(x * root));
}
//The symbols S, K, r, q, σ, t0 and T have theit usual meanings.
double d1_func(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double t0
		){
	// Define the variables d1 and d2 as follows:
	double d1 = ((log(S0/K) + (r - q)*(T - t0))/(sigma * sqrt(T - t0))) + 0.5 * sigma *sqrt(T - t0);
	return d1;
}

double d2_func(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double t0
		){
	      return  d1_func(S0, K, r, q, sigma, T, t0 ) - sigma *sqrt(T - t0);
}
// Write functions to calculate the fair value and Delta of
// European call and put options using the Black–Scholes–Merton formula.
double FairValue_Eur_call(
    	double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double t0,
        double d1,
        double d2
		){
		double c = S0 * exp(-1 * q * (T - t0)) * cum_norm(d1) - K * exp(-1 * r * (T - t0)) * cum_norm(d2);
		return c;
}

double FairValue_Eur_put(
    	double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double t0,
        double d1,
        double d2
		){
		double p =  K * exp(-1 * r * (T - t0)) * cum_norm(-1 * d2)- S0 * exp(-1 * q * (T - t0)) * cum_norm(-1* d1);
		return p;
}

double delta_c_func(
        double q,
		double T,
		double t0,
		double d1
		){
		double delta_c = exp(-1 * q * (T - t0)) * cum_norm(d1);
		return delta_c;
}
double delta_p_func(
        double q,
		double T,
		double t0,
		double d1
		){
		double delta_p = -1 * exp(-1 * q * (T - t0)) * cum_norm(-1 * d1);
		return delta_p;
}
//Verify values satisfies put-call party
void verify_put_call(double c_bsm, double p_bsm,
                     double S0, double K, double r, double q, double T, double t0){
	if(abs(abs(c_bsm - p_bsm) - abs(S0*exp(-1 * q * (T - t0)) - K * exp(-1 * r * (T - t0)))) <= 0.0000001){
		cout<<"c_bsm = "<< c_bsm <<" p_bsm = "<<p_bsm<<" satisfy put-call parity "<<endl<<endl;
	}else{
		cout<<"c_bsm = "<< c_bsm <<" p_bsm = "<<p_bsm<<" NOT satisfy put-call parity "<<endl<<endl;
	}
}

//Verify the values satify the relation for Delta
void verify_relation_delta(double delta_c , double delta_p,
						   double q, double T, double t0){
	if(abs(delta_c - delta_p - exp(-1 * q * (T - t0))) <= 0.0000001){
		cout<<"delta_c = "<<delta_c<<" delta_p = "<<delta_p<<" satisfy the relation for Delta"<<endl<<endl;
	}else{
		cout<<"delta_c = "<<delta_c<<" delta_p = "<<delta_p<<" NOT satisfy the relation for Delta"<<endl<<endl;
	}
}

int hw_7_BSM_test(){

	// 8 BlackñScholesñMerton valuations for examples in Lecture 17a
	// • Let us employ the parameter values for the examples in Lecture 17a
	double S0 =100;
	double K = 100;
	double r = 0.1;
	double q = 0;
	double sigma = 0.5;
	double T = 0.3;
	double t0 = 0;

	// Expected values for
	// d1 ≃ 0.246475, d2 ≃ −0.02739
	double d1 = d1_func(S0, K, r, q, sigma, T, t0 );
	double d2 = d2_func(S0, K, r, q, sigma, T, t0 );



	// The Black–Scholes–Merton values for the European call and put are:
	// Expected cBSM ≃ 12.2721 , pBSM ≃ 9.31668
	double c_bsm = FairValue_Eur_call(S0, K, r, q, sigma, T, t0, d1, d2);
	double p_bsm = FairValue_Eur_put(S0, K, r, q, sigma, T, t0, d1, d2);

	// The Black–Scholes–Merton values for the Delta of the European call and put are:
	// Expected ∆c ≃ 0.597343, ∆p ≃ −0.402657
	// Remember that the Delta of a put is negative.
	double delta_c = delta_c_func(q, T, t0, d1);
	double delta_p = delta_p_func(q, T, t0, d1);


	cout<<"d1 = " <<d1 << " ,Expected value d1(0.246475) " << endl << endl;
	cout<<"d2 = " <<d2 << " ,Expected value d2(-0.02739) " << endl << endl;
	cout<<"c_bsm = "<<c_bsm<<" ,Expected value c_bsm(12.2721)" << endl << endl;
	cout<<"p_bsm = "<<p_bsm<<" ,Expected value p_bsm(9.31668)" << endl << endl;
	cout<<"delta_c = "<<delta_c<< " ,Expected value delta_c(0.597343)" << endl << endl ;
	cout<<"delta_p = "<<delta_p<< " ,Expected value delta_p(-0.402657)"<< endl << endl;
	// Verify that the above values satisfy put-call parity
	verify_put_call(c_bsm, p_bsm, S0, K, r, q, T, t0);
	// Verify that the above values satisfy the relation for Delta
	verify_relation_delta(delta_c, delta_p, q, T, t0);


	// 9 BlackñScholesñMerton valuations for examples in Lecture 19a
	// . In the worked examples in Lecture 19a, I calculated American call and put options.
	// . Nevertheless, we can employ the parameter values in Lecture 19a and calculate the fair values
	// of the corresponding European call and put options.
	// The parameter values are:
	// S0 = 100 ,
	//  K = 100 ,
	//  r = 0.1 ,
	//  q = 0.1 ,
	//  σ = 0.5 ,
	//  T = 0.4 ,
	// t0 = 0.0 ,
	q = 0.1;
    T = 0.4;

	d1 = d1_func(S0, K, r, q, sigma, T, t0 );
	d2 = d2_func(S0, K, r, q, sigma, T, t0 );
	c_bsm = FairValue_Eur_call(S0, K, r, q, sigma, T, t0, d1, d2);
	p_bsm = FairValue_Eur_put(S0, K, r, q, sigma, T, t0, d1, d2);
	delta_c = delta_c_func(q, T, t0, d1);
	delta_p = delta_p_func(q, T, t0, d1);

	// The values of d1 and d2 are equal and opposite
	cout<<"d1 = "<<d1<< " ,Expected value d1(0.158114) "<<endl<<endl;
	cout<<"d2 = "<<d2<< " ,Expected value d2(-0.158114)"<< endl<<endl;

	// Scholes–Merton values for the European call and put are equal:
	cout<<"c_bsm = "<<c_bsm<< ",Expected value c_bsm(12.07068)"<<endl<<endl;
	cout<<"p_bsm = "<<p_bsm<< ",Expected value p_bsm(12.07068)"<<endl<<endl;

	// The Black–Scholes–Merton values for the Delta of the European call and put are:
	cout<<"delta_c = "<<delta_c<<" ,Expected value delta_c(0.540748)"<<endl<<endl;
	cout<<"delta_p = "<<delta_p<<" ,Expected value delta_p(−0.420041)"<<endl<<endl;

	// Verify that the above values satisfy put-call parity
	verify_put_call(c_bsm, p_bsm, S0, K, r, q, T, t0);

	// Verify that the above values satisfy the relation for Delta
	verify_relation_delta(delta_c, delta_p, q, T, t0);

	return 0;
}
int main(){
	hw_7_BSM_test();
	return 0;
}

