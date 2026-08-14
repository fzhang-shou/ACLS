// ACLS
// R follows AR1
// year effect of F follows AR1 process + add on age-dependent selectivity
// use an efficient pnorm function
// fit to catch length composition data
// fit to catch age composition data
// fit to CPUE data
// fit to total catch
// add effective sample size to the length/age frequency data

#include <TMB.hpp> 
#include <iostream>

template<class Type>  
Type fast_pnorm(Type x) {  
    const Type b0 = 0.2316419;  // Constants for the approximation  
    const Type b1 = 0.319381530;  
    const Type b2 = -0.356563782;  
    const Type b3 = 1.781477937;  
    const Type b4 = -1.821255978;  
    const Type b5 = 1.330274429;  
    const Type P = 0.398942280401; // 1 / sqrt(2π)  

    Type t = 1.0 / (1.0 + b0 * fabs(x)); // Transitional value  
    Type poly = t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * b5)))); // Polynomial approximation  
    Type normal_pdf = P * exp(-0.5 * x * x); // Normal PDF at x  

    Type result = 1.0 - normal_pdf * poly;  

    // If x is negative, flip the result for the lower tail  
    if (x < 0) result = 1.0 - result;  

    return result;  
}  

// Compute the length-age transition probability matrix  
template<class Type>  
matrix<Type> compute_length_age_transition(int n_ages, vector<Type> bin_border, vector<Type> length_at_age, Type cv_L) {  
    int n_bins = bin_border.size()+1; // Number of bins 
    matrix<Type> P(n_ages, n_bins);     // Probability matrix  
    P.setZero();
	vector<Type> p_border(n_bins - 1);
    for (int i = 0; i < n_ages; i++) {
		p_border.setZero();
        for(int j = 0; j < (n_bins - 1); j++){
			p_border(j) = fast_pnorm((bin_border(j) - length_at_age(i)) / (cv_L*length_at_age(i)));
		}
		P(i,0) = p_border(0);
        for (int j = 1; j < (n_bins-1); j++) {  
            P(i, j) = p_border(j) - p_border(j-1); // Probability of falling in this length-bin  
        }
		P(i, (n_bins-1)) = 1.0 - p_border(n_bins-2);
    } 
	
	for(int i=0; i<n_ages; ++i){
		for(int j=0; j<n_bins; ++j){
			if(P(i,j)<1e-20){P(i,j)=1e-20;};
		}
	}
    return P;  
}  

// Von Bertalanffy Growth: Compute length at age  
template<class Type> 
vector<Type> compute_length_at_age(vector<Type> ages, Type Linf, Type K, Type t0) {  
    vector<Type> lmid = Linf * (1.0 - exp(-K * (ages - t0))); // Vectorized Von Bertalanffy equation  	
    return lmid;
}

// Selectivity function
template<class Type>
vector<Type> compute_selectivity_at_age(vector<Type> ages, Type A50, Type A95){
	int A = ages.size();
	Type slope = log(0.95/0.05)/(A95 - A50);
	Type intercept = -1 * A50 * slope;
	vector<Type> sel(A);
	for(int i=0; i < A; i++){
		Type logit_pt = intercept + slope*ages(i);
		sel(i) = exp(logit_pt)/(1.0 + exp(logit_pt));
	}
	return sel;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // input data;
  DATA_VECTOR(obs_log_C); // annual total catch biomass (by Y)
  DATA_VECTOR(obs_log_I); // annual abundance index (by Y)
  DATA_MATRIX(obs_len_comp); // observed catch length composition (L by Y)
  DATA_MATRIX(obs_age_comp); // observed catch age composition (L by Y)
  DATA_VECTOR(ESS_len); // effective sample size of length data (by Y)
  DATA_VECTOR(ESS_age); // effective sample size of age data (by Y)
  //DATA_VECTOR(sel_len); // catch gear selectivity at length (by L)
  //DATA_VECTOR(sel_age); // catch gear selectivity at age (by A)
  DATA_VECTOR(len_border); // breaking border between length bins (by L-1) 
  DATA_VECTOR(age); // ages (by A)
  DATA_INTEGER(Y); // number of years
  DATA_INTEGER(A); // number of ages
  DATA_INTEGER(L); // number of length bins
  DATA_VECTOR(wgt); // weight at age (by A)
  DATA_VECTOR(mat); // maturing probability at age (by A)
  DATA_SCALAR(M); // assumed constant natural mortality
   
  Type zero = 0.0; 
  Type one = 1.0;

  //define parameters;
  
  // fixed effects
  PARAMETER(log_init_Z); // intital total mortality used to set initial age structure
  PARAMETER(log_std_log_N0); // variability of initial number at age
  
  PARAMETER(mean_log_R); // mean recruitment  
  PARAMETER(log_std_log_R); // recruitment standard error
  PARAMETER(logit_log_R); // define the AR1 coefficient of recruitment variation
  
  PARAMETER(mean_log_F);
  PARAMETER(log_std_log_F);
  PARAMETER(logit_log_F);
  PARAMETER(log_sel_A50); // parameter of selectivity 
  PARAMETER(log_sel_A95); // parameter of selectivity    
  
  PARAMETER(log_vbk); // VB parameter
  PARAMETER(log_Linf); // VB parameter
  PARAMETER(log_t0); // VB parameter
  PARAMETER(log_cv_len); // cv of length-at-age
  
  PARAMETER(log_q); // catchability
  PARAMETER(log_std_index); // measurement error of abundance index (CPUE)
  PARAMETER(log_std_catch); // measurement error of total catch biomass

  
  // random effects
  PARAMETER_VECTOR(dev_log_R); // recruitment deviation
  PARAMETER_VECTOR(dev_log_F); // deviation of fishing moratlity across years
  PARAMETER_VECTOR(dev_log_N0); // deviation of initial number at age (except for the first age), so size = A-1
  
  // derived parameters
  Type init_Z = exp(log_init_Z);
  Type std_log_N0 = exp(log_std_log_N0);
  
  Type std_log_R = exp(log_std_log_R);
  Type phi_log_R = exp(logit_log_R)/(one + exp(logit_log_R));
  
  Type std_log_F = exp(log_std_log_F);
  Type phi_log_F = exp(logit_log_F)/(one + exp(logit_log_F));
  Type sel_A50 = exp(log_sel_A50);
  Type sel_A95 = exp(log_sel_A95);
  
  Type vbk = exp(log_vbk);
  Type Linf = exp(log_Linf);
  Type t0 = exp(log_t0);
  Type cv_len = exp(log_cv_len); 
  
  Type std_index = exp(log_std_index);
  Type std_catch = exp(log_std_catch);  
  
  vector<Type> N(Y); N.setZero();
  vector<Type> F_y(Y); F_y.setZero();
  vector<Type> log_N(Y); log_N.setZero();
  matrix<Type> NA_imm(A,Y); NA_imm.setZero(); // immatured number at age
  matrix<Type> log_NA_imm(A,Y); log_NA_imm.setZero();
  matrix<Type> NA_mat(A,Y); NA_mat.setZero(); // matured number at age
  matrix<Type> log_NA_mat(A,Y); log_NA_mat.setZero();
  matrix<Type> NA_tot(A,Y); NA_tot.setZero();
  matrix<Type> NL(L,Y); NL.setZero();
  matrix<Type> Z(A,Y); Z.setZero();
  matrix<Type> F(A,Y); F.setZero();
  
  
  using namespace density; // call functions in TMB density namespace

  // --- Precomputations ---  
  vector<Type> length_at_age = compute_length_at_age(age, Linf, vbk, t0);       // Length-at-age from growth curve  
  matrix<Type> P = compute_length_age_transition(A, len_border, length_at_age, cv_len); // Length-age transition  
  
  //compute Z,F and M
  vector<Type> sel_age = compute_selectivity_at_age(age, sel_A50, sel_A95);
  for(int i = 0;i < A;++i){
    for(int j = 0;j < Y;++j){
        F(i,j) = exp(mean_log_F + dev_log_F(j)) * sel_age(i) + 1e-20; 
        Z(i,j)=F(i,j)+M;
        F_y(j)= exp(mean_log_F + dev_log_F(j))+ 1e-20;
     }
  }
  
  // get annual recruitment
  vector<Type> log_Rec = dev_log_R + mean_log_R; // may need to plug in SR model here
  vector<Type> Rec = exp(log_Rec); // - half * std_log_R * std_log_R); // bias correction
 
  // The cohort model;
  //initializing first year
  log_NA_imm(0,0) = log_Rec(0);
  NA_imm(0,0) = Rec(0);
  log_NA_mat(0,0) = log(1e-20);
  NA_mat(0,0) = 1e-20;
  //Type cz=zero;
  for(int i = 1;i < A;++i){ 
    log_NA_imm(i,0) = log_NA_imm(i-1,0) -init_Z  + dev_log_N0(i-1) + log(1.0 - mat(i));
	NA_imm(i,0) = exp(log_NA_imm(i,0));
	log_NA_mat(i,0) = (log_NA_imm(i-1,0) -init_Z  + dev_log_N0(i-1)) + log(mat(i));
	NA_mat(i,0) = exp(log_NA_mat(i,0));
  }
  //compute numbers at age
  for(int j = 1;j < Y;++j){
    log_NA_imm(0,j) = log_Rec(j);
	NA_imm(0,j) = Rec(j);
	log_NA_mat(0,j) = log(1e-20);
	NA_mat(0,j) = 1e-20;
    for(int i = 1;i < A;++i){
		log_NA_imm(i,j) = log_NA_imm(i-1,j-1) - Z(i-1,j-1) + log(1.0 - mat(i));
		NA_imm(i,j) = exp(log_NA_imm(i,j));
		log_NA_mat(i,j) = log_NA_imm(i-1,j-1) - Z(i-1,j-1) + log(mat(i));
		NA_mat(i,j) = exp(log_NA_mat(i,j));
	}
  }
  
  for(int i=0; i<A; ++i){
	  for(int j=0; j<Y; ++j){
		  NA_tot(i,j) = NA_imm(i,j) + NA_mat(i,j);
	  }
  }
  
  // compute number at length and total number
  for(int i=0; i<Y; ++i){
	NL.col(i) = P.transpose() * NA_tot.col(i);
	N(i)=NA_tot.col(i).sum();
	log_N(i)=log(N(i));
  }

  // compute SSB and biomass at age
  vector<Type> B(Y); B.setZero();
  vector<Type> SSB(Y); SSB.setZero();
  matrix<Type> BA(A,Y); BA.setZero();
  matrix<Type> SBA(L,Y); SBA.setZero();
  
  for(int i=0; i<A; ++i){
	  for(int j=0; j<Y; ++j){
		  BA(i,j) = NA_tot(i,j) * wgt(i);
		  SBA(i,j) = NA_mat(i,j) * wgt(i);
	  }
  }
  for(int i=0; i<Y; ++i){
	B(i) = BA.col(i).sum();
	SSB(i) = SBA.col(i).sum();
  }  

  //calculate catch at age, and catch age composition
  matrix<Type> CA(A,Y); CA.setZero();
  matrix<Type> CA_imm(A,Y); CA_imm.setZero();
  matrix<Type> CA_mat(A,Y); CA_mat.setZero();
  for(int i=0; i<A; ++i){
    for(int j=0; j<Y; ++j){
    CA_imm(i,j) = NA_imm(i,j) * (one - exp(-Z(i,j))) * (F(i,j) / Z(i,j));
    CA_mat(i,j) = NA_mat(i,j) * (one - exp(-F(i,j)));
    CA(i,j) = CA_imm(i,j) + CA_mat(i,j);

    }
  }

  matrix<Type> pred_age_comp(A,Y);  pred_age_comp.setZero();
  vector<Type> temp_age_comp(A); temp_age_comp.setZero();
  for (int i = 0; i < Y; i++) {
	temp_age_comp.array() = CA.col(i).array();		
    temp_age_comp /= temp_age_comp.sum(); // Normalize age composition
	pred_age_comp.col(i) = temp_age_comp;
  }

  //calculate catch at length, and catch length composition
  matrix<Type> CL(L,Y); CL.setZero();
  for (int i = 0; i < Y; i++) {
	CL.col(i) = P.transpose() * CA.col(i);
  }
  
  matrix<Type> pred_len_comp(L,Y);  pred_len_comp.setZero();
  vector<Type> temp_len_comp(L); temp_len_comp.setZero();
  for (int i = 0; i < Y; i++) {
	temp_len_comp.array() = CL.col(i).array();		
    temp_len_comp /= temp_len_comp.sum(); // Normalize age composition
	pred_len_comp.col(i) = temp_len_comp;
  }

  // calculate total catch
  vector<Type> C(Y);
  vector<Type> log_C(Y);
  matrix<Type> CB(A,Y);
  for(int i=0; i<A; i++){
	  for(int j=0; j<Y; j++){
		  CB(i,j) = CA(i,j) * wgt(i);
	  }
  }
  for(int i=0; i<Y; i++){
	  C(i) = CB.col(i).sum();
	  log_C(i) = log(C(i));
  }
  
  // negative log-likelihoods
  Type nll = zero;
  
  // abundance index, the measurement error
  vector<Type> pred_log_I = log_q + log_N;
  vector<Type> resid_log_I = obs_log_I - pred_log_I;
  nll -= dnorm(resid_log_I,zero,std_index,true).sum();
  
  // total catch
  vector<Type> resid_log_C = obs_log_C - log_C;
  nll -= dnorm(resid_log_C,zero,std_catch,true).sum();
	
  // catch age composition
  for (int i = 0; i < Y; i++) {
	vector<Type> obs_age = obs_age_comp.col(i);
	vector<Type> pred_age = pred_age_comp.col(i);
    nll -= ESS_age(i) * dmultinom(obs_age, pred_age, true);  // weighted by effective sample size
  }  
  
  // catch length composition
  for (int i = 0; i < Y; i++) {
	vector<Type> obs_len = obs_len_comp.col(i);
	vector<Type> pred_len = pred_len_comp.col(i);
    nll -= ESS_len(i) * dmultinom(obs_len, pred_len, true);  // weighted by effective sample size
  }  
  
  // logF
  nll += SCALE(AR1(phi_log_F),std_log_F)(dev_log_F); // AR1 evaluation to generate nll of logF
  
  //recruitment estimates 
  nll += SCALE(AR1(phi_log_R),std_log_R)(dev_log_R); // AR1 evaluation to generate nll for logR
  
  // initial population size
  nll -= dnorm(dev_log_N0,zero,std_log_N0,true).sum(); // iid evaluation of the initial population size
  
 // report results
  REPORT(P);
  REPORT(N)
  REPORT(NL); 
  REPORT(NA_imm);
  REPORT(NA_mat);
  REPORT(NA_tot);
  REPORT(dev_log_N0);
  REPORT(std_log_N0);  
  REPORT(B);
  REPORT(BA);
  REPORT(SSB);
  REPORT(SBA);
  REPORT(CA);
  REPORT(C);
  REPORT(CL);
  REPORT(CA_mat);
  REPORT(CA_imm);
  REPORT(Z);
  REPORT(F);
  REPORT(mean_log_F);
  REPORT(dev_log_F);
  REPORT(phi_log_F);
  REPORT(Rec);
  REPORT(mean_log_R);
  REPORT(dev_log_R);
  REPORT(std_log_R);
  REPORT(phi_log_R);
  REPORT(pred_log_I);  
  REPORT(resid_log_I);
  REPORT(std_index);
  REPORT(std_catch);
  REPORT(pred_age_comp);
  REPORT(pred_len_comp);
  REPORT(F_y);
  REPORT(vbk);
  REPORT(Linf);
  REPORT(t0);
  REPORT(cv_len); 
  
  ADREPORT(Rec);
  ADREPORT(N);
  ADREPORT(NA_tot);
  ADREPORT(NL);
  ADREPORT(B);
  ADREPORT(SSB);
  ADREPORT(F);
 
  return nll;
}
