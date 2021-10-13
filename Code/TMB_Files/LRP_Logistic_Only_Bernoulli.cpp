#include <TMB.hpp>


  
  
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(LM_Agg_Status);
  DATA_VECTOR(LM_Agg_Abund);
  DATA_IVECTOR(LM_yr);
  //DATA_IVECTOR(LM_stk); // FLAG maybe take this out. Redundant?
  DATA_VECTOR(Pred_Abund);
  DATA_SCALAR(p);
 
  PARAMETER(B_0);
  PARAMETER(B_1);
  
  Type ans=0.0;
 
 
  // Compile "data" for logistic model
  // get number of years, number of obs
  int Logistic_Mod_Yrs = LM_Agg_Abund.size();
 

 // vector for logistic likelihood
 vector<Type> LogitP(Logistic_Mod_Yrs);
 vector<Type> N_bern(Logistic_Mod_Yrs);
 
 // create logistic likelihood and
 // Fill Ns with values
 for(int i=0; i<Logistic_Mod_Yrs; ++i){
   LogitP(i) = B_0 + B_1*LM_Agg_Abund(i);
   N_bern(i) = 1;
 }
 
 // Now fit logistic model
  ans += -sum(dbinom_robust(LM_Agg_Status, N_bern, LogitP, true));

 
 //Get final BM
 Type Agg_LRP = (log(p/(1-p)) - B_0)/(B_1);
 
 // Get estimates for plotting CIs
 int N_Preds = Pred_Abund.size();
 vector<Type> Logit_Preds(N_Preds);
 
 Logit_Preds = B_0 + B_1*Pred_Abund;
 
 ADREPORT(Agg_LRP);
 ADREPORT(Logit_Preds);
 REPORT(ans);
 
 return ans;
 
  
}
