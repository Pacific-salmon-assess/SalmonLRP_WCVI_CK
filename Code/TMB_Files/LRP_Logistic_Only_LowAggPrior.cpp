#include <TMB.hpp>


  
  
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_Stks);
  DATA_INTEGER(Bern_Logistic);
  DATA_VECTOR(LM_CU_Status);
  DATA_VECTOR(LM_Agg_Abund);
  DATA_IVECTOR(LM_yr);
  DATA_IVECTOR(LM_stk);
  DATA_VECTOR(Pred_Abund);
  DATA_SCALAR(p);
  DATA_SCALAR(B_penalty_mu); // for likelihood penalty on low aggregate abundance; average value of aggregate abundance at low p
  DATA_SCALAR(B_penalty_sigma); // for likelihood penalty on low aggregate abundnace; sd of value of aggregate abundance at low p

 
  PARAMETER(B_0);
  PARAMETER(B_1);
  
  Type ans=0.0;
 
 
  // Compile "data" for logistic model
  // get number of years, number of obs
  int Logistic_Mod_Yrs = LM_Agg_Abund.size();
  int N_LM_Obs = LM_CU_Status.size();
  // create vector with number of pops above their benchmark (Sgen or percentile)
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
   if(LM_CU_Status(i) == 1){
     N_Above_BM(LM_yr(i)) += 1;
   }
 } // end for loop over obs
 

// for(int i=0; i<N_LM_Obs; ++i){
//   //check if Spawners above Sgen
//   if(LM_S(i) > Sgen(LM_stk(i))){
//     N_Above_BM(LM_yr(i)) += 1;
//   }
// } // end for loop over obs

 
 
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
 
 // Now fit logistic model
 if(Bern_Logistic == 1){
   ans += -sum(dbinom_robust(All_Above_BM, N_bern, LogitP, true));
 } else if(Bern_Logistic == 0) {
   ans +=  -sum(dbinom_robust(N_Above_BM, N_bin, LogitP, true));
 }
 
 //Get final BM
 Type Agg_LRP = (log(p/(1-p)) - B_0)/(B_1);
 
 // Add prior penalty on Aggregate abbundance at low proportion of CUs above benchmark
 Type p_min = 0.01; // value p_min = low proportion of CUs greater than benchmark (Sgen)
 ans += -dnorm(  (log(p_min/(1-p_min)) - B_0) / B_1, B_penalty_mu , B_penalty_sigma, true ); // likelihood penalty on aggregate abundance
 
 // Get estimates for plotting CIs
 int N_Preds = Pred_Abund.size();
 vector<Type> Logit_Preds(N_Preds);
 
 Logit_Preds = B_0 + B_1*Pred_Abund;
 
 REPORT(N_Above_BM);
 REPORT(All_Above_BM);
 ADREPORT(Agg_LRP);
 ADREPORT(Logit_Preds);
 
 return ans;
 
  
}
