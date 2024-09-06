#include <TMB.hpp>


template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_INTEGER(biasCor);
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  //DATA_VECTOR(predS);
  //DATA_IVECTOR(yr);

  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);


  Type ans=0.0;
  int N_Obs = S.size();
  //int N_Preds predS.size();
  vector <Type> LogR_Pred(N_Obs);
  //vector <Type> LogR_Pred2(N_predS);
  vector <Type> sigma = exp(logSigma);
  vector <Type> A = exp(logA);


  // Ricker likelihood
  for (int i = 0; i<N_Obs; i++){
    LogR_Pred(i) = logA(stk(i)) + log(S(i)) - exp(logB(stk(i))) * S(i);
    //ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stk(i)), true);
    if(biasCor == 0) ans += -dnorm(logR(i) - LogR_Pred(i), Type(0),  sigma(stk(i)), true);
    if(biasCor == 1) ans += -dnorm(logR(i) - LogR_Pred(i), - pow(sigma(stk(i)),2)/Type(2),  sigma(stk(i)), true);

  }

  //Ricker predictions over a larger range of S
  //for (int j = 0; j<N_predS; i++){
  //  LogR_Pred2(i) = logA(stk(i)) + log(predS(i)) - exp(logB(stk(i))) * predS(i);
  //  //ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stk(i)), true);
  //  if(biasCor == 0) ans += -dnorm(logR(i) - LogR_Pred2(i), Type(0),  sigma(stk(i)), true);
  //  if(biasCor == 1) ans += -dnorm(logR(i) - LogR_Pred2(i), - pow(sigma(stk(i)),2)/Type(2),  sigma(stk(i)), true);
 // }


  ADREPORT(A);
  ADREPORT(LogR_Pred);
  //ADREPORT(LogR_Pred2);

  return ans;

}

