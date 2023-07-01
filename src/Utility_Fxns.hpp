// Purpose: An assortment of utility functions and likelihoods
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 3/18/23

template <class Type> // Function to call different selectivity parameterizations
// @param age = indexed age within age loop
// @param sel_model = integer of selectivity model 
// == 0, logistic
// @param ln_selpars vector of log selectivity parameters.
Type Get_Selex(Type age, 
               int sel_model, 
               vector<Type> ln_selpars) {
  
  // Create container to return predicted selectivity value
  Type selex = 0;
  
  if(sel_model == 0) { // logistic selectivity
    // Extract out and exponentiate the parameters here
    Type a50 = exp(ln_selpars(0)); // a50
    Type k = exp(ln_selpars(1)); // slope
    selex = Type(1.0) / (Type(1) + exp(-1 * ((age - a50)/k)));
  }
  
  if(sel_model == 1) { // gamma dome-shaped selectivity 
    // Extract out and exponentiate the parameters here
    Type delta = exp(ln_selpars(0)); // slope parameter
    Type amax = exp(ln_selpars(1)); // age at max selex
    
    // Now, calculate/derive power parameter + selex values
    Type p = 0.5 * (sqrt( pow(amax, 2) + (4 * pow(delta, 2))) - amax);
    selex = pow( (age / amax), (amax/p)  ) * exp( (amax - age) / p ); 
    
  }
  
  if(sel_model == 2) { // double logistic selectivity
    // Extract out and exponentiate the parameters here
    Type slope_1 = exp(ln_selpars(0)); // Slope of ascending limb
    Type slope_2 = exp(ln_selpars(1)); // Slope of descending limb
    Type infl_1 = exp(ln_selpars(2)); // Inflection point of ascending limb
    Type infl_2 = exp(ln_selpars(3)); // Inflection point of descending limb

    // Calculate logistic curve 1
    Type logist_1 = 1.0/(1.0 + exp(-slope_1 * (age - infl_1)));
    // Calculate logistic curve 2
    Type logist_2 = 1/(1.0 + exp(-slope_2 * (age - infl_2)));
    // Return selectivity - product of two logistic selectivity curves
    selex = logist_1 * (1 - logist_2);

  }
  
  if(sel_model == 3) { // exponential logistic (Thompson 1994)
    
    // Extract out and exponentiate the parameters here
    Type gamma = exp(ln_selpars(0)); // Controls Doming
    Type beta = exp(ln_selpars(1)); // Controls Peak of Ascending Limb
    Type alpha = exp(ln_selpars(2)); // Controls Ascending Limb

    // Define equations to minimize mistakes
    Type first = (1 / (1 - gamma));
    Type second = pow((1 - gamma) / gamma, gamma);
    Type third = exp( alpha * gamma * (beta - age ) );
    Type fourth = 1 + exp(alpha * (beta - age));

    // Return selectivity
    selex = first * second * (third/fourth);
  }
  
  return selex;
} // end function

template<class Type>
// Function to compute likelihood for dirichlet-multinomial (follows linear parameterization of
// Thorson et al. 2017)
// @param obs = Observed vector (in proportions)
// @param pred = Predicted vector (in proportions)
// @param Input_N = Input sample size
// @param Dir_Param = parameter for DM
// @param give_log = whether or not to compute log of likelihood
Type ddirmult( vector<Type> obs, 
               vector<Type> pred, 
               Type Input_N, 
               Type Dir_Param, 
               int do_log = 1){
  
  // Pre-processing
  int n_a = obs.size(); // Number of age-classes/bins
  Type Ntotal = Input_N; // Input sample size assumed
  vector<Type> p_pred = pred; // Predicted vector
  vector<Type> p_obs = obs; // Observed vector
  Type dirichlet_param = exp(Dir_Param) * Ntotal; // Dirichlet Alpha Parameter

  // NOTE: These calculations are in log space (Thorson et al. 2017)
  // 1st term -- integration constant that could be dropped
  Type logLike = lgamma(Ntotal + 1);
  for(int a = 0; a < n_a; a++) logLike -= lgamma(Ntotal * p_obs(a) + 1);

  // 2nd term in formula
  logLike += lgamma(dirichlet_param) - lgamma(Ntotal + dirichlet_param);

  // Summation in 3rd term
  for(int a = 0; a < n_a; a++){
    logLike += lgamma( (Ntotal * p_obs(a)) + (dirichlet_param * p_pred(a)) );
    logLike -= lgamma(dirichlet_param * p_pred(a));
  } // end a loop
  
  // Type phi = dirichlet_param.sum();
  // Type logLike = lgamma(Ntotal + 1.0) + lgamma(phi) - lgamma(Ntotal + phi);
  // for(int a = 0; a < n_a; a++) {
  //   logLike += -lgamma(p_obs(a) + 1.0) +
  //     lgamma(p_obs(a) + dirichlet_param(a)) -
  //     lgamma(dirichlet_param(a));
  // } // end a loop
  
  if(do_log == 1) return logLike; else return exp(logLike);
} // end function

template <class Type>
// Function to compute likelihood for dirichlet 
// @param p_obs = Observed vector (in proportions)
// @param p_pred = Predicted vector (in proportions)
// @param Ntotal = Input sample size
// @param Dir_Param = parameter for DM
// @param Ntotal = Input sample size
// @param give_log = whether or not to compute log of likelihood
Type ddirichlet(vector<Type> p_obs, 
                vector<Type> p_pred, 
                Type Ntotal, 
                Type Dir_Param,
                int do_log = 1)
{
  // Pre-processing
  Type lambda = Ntotal/pow(exp(Dir_Param), 2) - 1; // Compute Lambda parameters
  vector<Type> alpha = lambda * p_pred;  // Compute alpha dirichlet parameters
  Type phi = alpha.sum(); // Sum to get phi
  int n_ages = p_obs.size(); // Number of ages

  // Define likelihood
  Type logLike = lgamma(phi);
  
  // Loop through to compute likelihood
  for(int a = 0; a < n_ages; a++) {
    logLike +=  -lgamma(alpha(a)) + (alpha(a) - 1.0) * log(p_obs(a));
  } // end a loop
  
  if(do_log == 1) return logLike;
  else return exp(logLike);
}

template<class Type>
// Function to compute the likelihoods of a variety of models (multinomial, dirichlet, and dirichlet-multinomial)
// @param p_obs = Observed vector (in proportions)
// @param p_pred = Predicted vector (in proportions)
// @param Ntotal = Input sample size
// @param Dir_Param = parameter for dirichlet or dirichlet multinomial (log space)
// @param Ntotal = Input sample size
// @param acomp_LL_model = Integer to denote likelihood model (== 0, multinomial, == 1, dirichlet-multinomial, == 2, dirichlet)
// @param give_log = whether or not to compute log of likelihood
Type get_acomp_nLL(vector<Type> p_obs, 
                   vector<Type> p_pred, 
                   Type Ntotal, 
                   Type Dir_Param,
                   int acomp_LL_model, 
                   int do_log) {
  
  Type logLike = 0; // Likelihood storage
  
  // Compute likelihoods
  if(acomp_LL_model == 0) { // Multinomial
    vector<Type> obs_vec = p_obs * Ntotal; // Get Observed vector in numbers
    logLike += dmultinom(obs_vec, p_pred, true); // Compute likelihood
  } // end if
  
  if(acomp_LL_model == 1) { // Dirichlet-Multinomial
    logLike += ddirmult(p_obs, p_pred, Ntotal, Dir_Param, true);  // cCompute likelihood
  } // end if
  
  if(acomp_LL_model == 2) { // Dirichlet
    logLike += ddirichlet(p_obs, p_pred, Ntotal, Dir_Param, true);
  } // end if
  
  if(do_log == 1) return logLike;  else return exp(logLike);
} // end function

template <class Type> 
// Function to compute effective sample sizes for given composition likelihoods
// @param Dir_Param = parameter for dirichlet or dirichlet multinomial (log space)
// @param Ntotal = Input sample size
// @param acomp_LL_model = Integer to denote likelihood model (== 0, multinomial, == 1, dirichlet-multinomial, == 2, dirichlet)
Type get_Neff(Type Ntotal,
              Type Dir_Param,
              int acomp_LL_model) {
  
  Type Neff = 0; // Define for storage
  
  // Compute compositions
  if(acomp_LL_model == 0) { // Multinomial
    Neff = Ntotal;
  } // end if
  
  if(acomp_LL_model == 1) { // Dirichlet-Multinomial (Saturating Version)
    Neff = (1 + (Ntotal * exp(Dir_Param))) / (1 + exp(Dir_Param));
  } // end if
  
  if(acomp_LL_model == 2) { // Dirichlet
    Neff = Ntotal / pow(exp(Dir_Param), 2);
  } // end if
  
  return Neff; 
} // end function

