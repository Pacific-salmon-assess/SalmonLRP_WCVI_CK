#include <TMB.hpp>

// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
  // Scalar version
  template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }


// Start Model

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(P_3);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(stk_predS);
  DATA_IVECTOR(yr);
  DATA_VECTOR(logSurv_3);
  DATA_VECTOR(logSurv_4);
  DATA_INTEGER(N_Stks);
  DATA_VECTOR(muLSurv);
  DATA_SCALAR(Tau_dist);
  DATA_SCALAR(gamma_mean);
  DATA_SCALAR(gamma_sig);
  DATA_INTEGER(Bern_Logistic);
  DATA_VECTOR(LM_S);
  DATA_VECTOR(LM_Agg_Abund);
  DATA_IVECTOR(LM_yr);
  DATA_IVECTOR(LM_stk);
  DATA_VECTOR(Pred_Abund);
  DATA_VECTOR(Pred_Spwn);
  DATA_VECTOR(cap_mean);
  DATA_SCALAR(cap_sig);
  DATA_SCALAR(p);
  DATA_SCALAR(Sgen_sig);
  DATA_SCALAR(B_penalty_mu);
  DATA_SCALAR(B_penalty_sigma);
  
  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(cap);
  PARAMETER_VECTOR(logSigma);
  PARAMETER(gamma);
  PARAMETER_VECTOR(logSgen);
  PARAMETER(B_0);
  PARAMETER(B_1);
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector<Type> LogR_Pred(N_Obs);
  vector<Type> LogR_Pred_3(N_Obs);
  vector<Type> LogR_Pred_4(N_Obs);
  vector <Type> sigma=exp(logSigma);
  vector <Type> SMSY(N_Stks);  
  vector <Type> SRep(N_Stks); 
  vector <Type> LogSMSY(N_Stks);
  vector <Type> Sgen = exp(logSgen);
  vector <Type> logProd(N_Stks);
  vector <Type> A(N_Stks);
  vector <Type> B(N_Stks);

  
  // Ricker likelihood based on Arbeider et al. 2020 Interior Fraser Coho RPA res. doc.
  for(int i=0; i<N_Obs; i++){
    
    B(stk(i)) = (logA(stk(i)) + gamma * muLSurv(stk(i))) / cap(stk(i));
    
    LogR_Pred_3(i) = logA(stk(i)) + gamma*logSurv_3(i) + log(P_3(i)*S(i)) - B(stk(i)) * S(i);
    LogR_Pred_4(i) = logA(stk(i)) + gamma*logSurv_4(i) + log((1-P_3(i))*S(i)) - B(stk(i)) * S(i);
    LogR_Pred(i) = log(exp(LogR_Pred_3(i)) + exp(LogR_Pred_4(i)));
    // get same answer whether put likelihood on log(R/S) or R_pred
    ans += -dnorm(LogR_Pred(i) - log(S(i)), logR(i) - log(S(i)),  sigma(stk(i)), true);
  }
  
  // Calculate SMSY using Lambert's W function
  // Approach from Scheurell 2016
  for(int i=0; i<N_Stks; i++){
	ans += -dgamma(pow(sigma(i),-2), Tau_dist, 1/Tau_dist, true);
    // Jacobian adjustment for prior on sigma (add for TMBstan runs)
	// ans -= log(2) - 2*logSigma(i);
	// Mu-survival adjusted productivity for SMSY, Sgen calcs
	logProd(i) = logA(i) + gamma*muLSurv(i);
    A(i) = exp(logProd(i));
  // Calculate SMSY using Lambert W function
  SMSY(i) = (1 - LambertW(exp(1-logProd[i])) ) / B(i);
  }
  
  // Add priors ====================
  // Gamma prior
  ans += -dnorm(gamma, gamma_mean, gamma_sig, true);
  
  // prior for CU-specific caps
  for (int i = 0; i < N_Stks; i++){
    ans += -dnorm(cap(i), cap_mean(i), cap_sig, true);
  }
  
  // Estimate Sgen =========================
  LogSMSY = logProd + logSgen - B * Sgen;
  vector <Type> Diff = exp(LogSMSY)-SMSY;
  ans += -sum(dnorm(Diff, 0, Sgen_sig, true ));
  
  // Calculate SRep =======================
  SRep = logProd / B;
  
  // Get estimates for plotting SR model fit with CIs =========
  int N_SpwnPreds = Pred_Spwn.size();
  vector<Type> LogRec_Preds(N_SpwnPreds);
  vector<Type> Rec_Preds(N_SpwnPreds);
  
  for(int i=0; i<N_SpwnPreds; i++){
    LogRec_Preds(i) = logA(stk_predS(i)) + gamma*muLSurv(stk_predS(i)) + log(Pred_Spwn(i)) - B(stk_predS(i)) * Pred_Spwn(i);
  }
  Rec_Preds = exp(LogRec_Preds);
  
  
  // Compile "data" for logistic model ====================
  // get number of years, number of obs
  int Logistic_Mod_Yrs = LM_Agg_Abund.size();
  int N_LM_Obs = LM_S.size();
  // create vector with number of pops above their BM (Sgen)
  vector <Type> N_Above_BM(Logistic_Mod_Yrs);
  // Also version with 0/1's if doing bernoulli
  vector <Type> All_Above_BM(Logistic_Mod_Yrs);
  // set both to 0
  All_Above_BM.setZero();
  N_Above_BM.setZero();
  
  // vectors of number of "trials" for logistic mods
  vector<Type> N_bin(Logistic_Mod_Yrs);
  vector<Type> N_bern(Logistic_Mod_Yrs);
  // vector for logistic likelihood
  vector<Type> LogitP(Logistic_Mod_Yrs);
  
  for(int i=0; i<N_LM_Obs; ++i){
    //check if Spawners above Sgen
    if(LM_S(i) > Sgen(LM_stk(i))){
      N_Above_BM(LM_yr(i)) += 1;
    }
  } // end for loop over obs
  
  // create logistic likelihood and
  // Fill Ns with values
  for(int i=0; i<Logistic_Mod_Yrs; ++i){
    LogitP(i) = B_0 + B_1*LM_Agg_Abund(i);
    N_bin(i) = N_Stks;
    N_bern(i) = 1;
    // Also fill in bernoulli version with 1s when all above LRP
    if(N_Above_BM(i) == N_Stks){
      All_Above_BM(i) = 1;
    }
  }
  
  // Now fit logistic model ========================
  if(Bern_Logistic == 1){
    ans += -sum(dbinom_robust(All_Above_BM, N_bern, LogitP, true));
  } else if(Bern_Logistic == 0) {
    ans +=  -sum(dbinom_robust(N_Above_BM, N_bin, LogitP, true));
  }
  
  //Get final BM
  Type Agg_LRP = (log(p/(1-p)) - B_0)/(B_1);
  
  // Add prior penalty on aggregate abundance at low proportion of CUs above benchmark (if Bernoulli == F) 
  // ---- or penalty on aggregate abundance with a 1% probability of all CUs being above benchmark (if Bernoulli = T)
  Type p_min = 0.01; // value p_min = low proportion of CUs greater than benchmark (Sgen)
  ans += -dnorm((log(p_min/(1-p_min)) - B_0) / B_1, B_penalty_mu , B_penalty_sigma, true ); // likelihood penalty on aggregate abundance
  
  // Get estimates for plotting CIs
  int N_Preds = Pred_Abund.size();
  vector<Type> Logit_Preds(N_Preds);
  
  Logit_Preds = B_0 + B_1*Pred_Abund;
  
  REPORT(N_Above_BM);
  REPORT(All_Above_BM);
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  ADREPORT(SRep);
  ADREPORT(Agg_LRP);
  ADREPORT(A);
  ADREPORT(B);
  ADREPORT(Logit_Preds);
  ADREPORT(Rec_Preds);
  REPORT(ans);
  
  return ans;
  
}
