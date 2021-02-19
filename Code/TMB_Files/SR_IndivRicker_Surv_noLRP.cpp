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
  DATA_INTEGER(Bayes);
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
  DATA_VECTOR(Pred_Spwn);
  DATA_SCALAR(Sgen_sig);
  
  
  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);
  PARAMETER(gamma);
  PARAMETER_VECTOR(logSgen);
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector<Type> LogR_Pred(N_Obs);
  vector<Type> LogR_Pred_3(N_Obs);
  vector<Type> LogR_Pred_4(N_Obs);
  vector<Type> R_Resid(N_Obs);
  vector <Type> sigma=exp(logSigma);
  vector <Type> SMSY(N_Stks);  
  vector <Type> SRep(N_Stks); 
  vector <Type> LogSMSY(N_Stks);
  vector <Type> Sgen = exp(logSgen);
  vector <Type> B = exp(logB);
  vector <Type> logProd(N_Stks);
  vector <Type> A(N_Stks);

  
  // Ricker likelihood based on Arbeider et al. 2020 Interior Fraser Coho RPA res. doc.
  for(int i=0; i<N_Obs; i++){
    LogR_Pred_3(i) = logA(stk(i)) + gamma*logSurv_3(i) + log(P_3(i)*S(i)) - exp(logB(stk(i))) * S(i);
    LogR_Pred_4(i) = logA(stk(i)) + gamma*logSurv_4(i) + log((1-P_3(i))*S(i)) - exp(logB(stk(i))) * S(i);
    LogR_Pred(i) = log(exp(LogR_Pred_3(i)) + exp(LogR_Pred_4(i)));
    // get same answer whether put likelihood on log(R/S) or R_pred
    ans += -dnorm(LogR_Pred(i) - log(S(i)), logR(i) - log(S(i)),  sigma(stk(i)), true);
  }
  
  // Calculate SMSY using Lambert's W function
  // Approach from Scheurell 2016
  for(int i=0; i<N_Stks; i++){
	ans += -dgamma(pow(sigma(i),-2), Tau_dist, 1/Tau_dist, true);
  if(Bayes == 1){
    // Jacobian adjustment; only needed for Bayesian runs using TMBstan
    ans -= log(2) - 2*logSigma(i);
  }
  // Mu-survival adjusted productivity for SMSY, Sgen calcs
  logProd(i) = logA(i) + gamma*muLSurv(i);
  A(i) = exp(logProd(i));
  // Calculate SMSY using Lambert W function
  SMSY(i) =  (1 - LambertW(exp(1-logProd(i)))) / B(i) ;
  }
  
  // Add priors ====================
  // Gamma prior
  ans += -dnorm(gamma, gamma_mean, gamma_sig, true);
  
  // Estimate Sgen =========================
  LogSMSY = logProd + logSgen - B * Sgen;
  vector <Type> Diff = exp(LogSMSY)-SMSY;
  ans += -sum(dnorm(Diff, 0, Sgen_sig, true ));
  
  // Calculate SRep =======================
  SRep = logProd / B;
  
  // Get estimates for plotting SR model fit with CIs =========
  int N_SpwnPreds = Pred_Spwn.size();
  vector<Type> LogRec_Preds(N_SpwnPreds); // create vector of predicted log(recruits)
  vector<Type> Rec_Preds(N_SpwnPreds);
  
  for(int i=0; i<N_SpwnPreds; i++){  // predict recruitment based on Ricker equation with stock-specific parameters
    LogRec_Preds(i) = logA(stk_predS(i)) + gamma*muLSurv(stk_predS(i)) + log(Pred_Spwn(i)) - exp(logB(stk_predS(i))) * Pred_Spwn(i);
  }
  Rec_Preds = exp(LogRec_Preds);
  
  R_Resid = exp(logR) - exp(LogR_Pred);
  
  REPORT(R_Resid);
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  ADREPORT(SRep);
  ADREPORT(A);
  ADREPORT(Rec_Preds);
  
  return ans;
  
}
