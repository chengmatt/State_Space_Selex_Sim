// General single species age-structured stock assessment
// that accommodates multiple fishery fleets written in TMB
// Creator: Matthew LH. Cheng (UAF-CFOS)

// TO DO:
// Beverton Holt SR
// Set Selectivities

#include<TMB.hpp>
template<class Type>

Type objective_function<Type>::operator() ()
{
  using namespace density; // Define namespace to use multivariate distributions

  // DATA SECTION ----------------------------------------------
  // Define general model dimensions
  DATA_VECTOR(ages); // Vector of ages
  DATA_VECTOR(years); // Vector of years
  DATA_INTEGER(n_sexes); // Number of sexes
  DATA_INTEGER(n_fleets); // Number of fishery fleets we want to model
  DATA_INTEGER(n_fish_indices); // Number of fishery indices of abundance 
  DATA_INTEGER(n_srv_indices); // Number of survey indices of abundance 
  DATA_INTEGER(n_fish_comps); // Number of fishery comps
  DATA_INTEGER(n_srv_comps); // Number of survey comps
  
  int n_ages = ages.size(); // Determine length of age vector
  int n_years = years.size(); // Determine length of years vector
  
  // Observations ----------------------------------------------
  DATA_MATRIX(obs_catches); // Matrix of catch from each fleet; n_years * n_fleets
  DATA_ARRAY(obs_fish_age_comps); // Array of fishery age comps; n_years * n_ages * n_fleets * n_sexes
  DATA_ARRAY(obs_srv_age_comps); // Array of survey age comps; n_years * n_ages * n_srv_indices * n_sexes
  DATA_MATRIX(obs_fish_age_Neff); // Array of fishery age comps; n_years * n_fleets 
  DATA_MATRIX(obs_srv_age_Neff); // Array of survey age comps; n_years * n_srv_indices 
  DATA_MATRIX(obs_fish_indices); // Matrix of Fishery Indices of Abundance
  DATA_MATRIX(obs_srv_indices); // Matrix of Survey Indices of Abundance
  DATA_VECTOR(fish_cv); // Vector of fishery index CVs
  DATA_VECTOR(srv_cv); // Vector of survey index CVs
  
  // Biological Information Inputs ----------------------------------------------
  DATA_ARRAY(WAA); // Weight-at-age array; n_years * n_ages * n_sexes
  DATA_ARRAY(MatAA); // Maturity-at-age array; n_years * n_ages * n_sexes (But only Sex0 used for calcs)
  DATA_VECTOR(Sex_Ratio); // Sex ratio; n_sexes - Females, Males
  
  // Controls on assessment ----------------------------------------------
  DATA_INTEGER(rec_model); // Recruitment model, == 0 Mean Recruitment, == 1 Beverton Holt 
  DATA_ARRAY(F_Slx_model); // Fishery Selectivity Model, == 0 Logistic n_fleets * n_years * n_sexes
  DATA_ARRAY(S_Slx_model); // Survey Selectivity Model, == 0 Logistic n_fleets * n_years * n_sexes

  // PARAMETER SECTION ----------------------------------------------
  PARAMETER_VECTOR(ln_N1_Devs); // log deviations for initial abundance
  PARAMETER(ln_MeanRec); // log mean recruitment 
  PARAMETER(ln_SigmaRec); // log sigma for recruitment
  PARAMETER_VECTOR(ln_RecDevs); // log recruitment deviations
  PARAMETER(ln_M); // log natural mortality
  // PARAMETER_VECTOR(ln_MeanF);
  PARAMETER_MATRIX(ln_Fy); // Annual Fishing Mortality ; n_years * n_fleets
  PARAMETER_VECTOR(ln_q_fish); // log catchability for fishery; n_fish_indices 
  PARAMETER_VECTOR(ln_q_srv); // log catchability for survey; n_srv_indices 
  PARAMETER(ln_a50_f); // a50 logistic fishery
  PARAMETER(ln_k_f); // k logistic fishery
  PARAMETER(ln_a50_s); // a50 logistic survey
  PARAMETER(ln_k_s); // k logistic survey
  
  // Parameter Transformations ----------------------------------------------
  Type M = exp(ln_M); // Natural Mortality
  Type ln_SigmaRec2 = ln_SigmaRec * ln_SigmaRec; // Variance of recruitment sigma

  // Predicted Quantities
  array<Type> pred_fish_age_comps(obs_fish_age_comps.dim); // Predicted Fishery Age Comps
  array<Type> pred_srv_age_comps(obs_srv_age_comps.dim); // Predicted survey age comps
  matrix<Type> pred_fish_indices(n_years, n_fish_indices); // Predicted fishery indices
  matrix<Type> pred_srv_indices(n_years, n_srv_indices); // Predicted survey indices
  matrix<Type> pred_catches(n_years, n_fleets); // Predicted fishery catches
  
  // Stored Quantities
  array<Type> NAA(n_years + 1, n_ages, n_sexes); // Numbers at age; n_years + 1 (forward projection)
  array<Type> ZAA(n_years, n_ages, n_sexes); // Total Mortality
  array<Type> SAA(n_years, n_ages, n_sexes); // Survival at age
  array<Type> FAA(n_years, n_ages, n_fleets, n_sexes); // Fishing Mortality
  array<Type> CAA(n_years, n_ages, n_fleets, n_sexes); // Catch at Age
  
  matrix<Type> ln_F(n_years, n_fleets); // Get F by fleet in log space
  matrix<Type> F(n_years, n_fleets); // Get F by fleet in normal space
  vector<Type> Total_Fy(n_years); // Total F summed across fleets
  vector<Type> Total_Rec(n_years); // Total Recruitment
  vector<Type> Total_Biom(n_years); // Total Biomass
  vector<Type> SSB(n_years + 1); // Spawning stock biomass; n_years + 1 (forward projection)

  // Selectivities
  array<Type> F_Slx(n_years, n_ages, n_fleets, n_sexes); // Fishery Selectivities
  array<Type> Total_Fishery_Numbers(n_years, n_fleets, n_sexes); // Store Total Fishery Numbers for Proportions
  array<Type> S_Slx(n_years, n_ages, n_srv_indices, n_sexes); // Survey Selectivities
  array<Type> Total_Survey_Numbers(n_years, n_fleets, n_sexes); // Store Total Survey Numbers for Proportions

  // Define objects to store negative log-likelihoods
  vector<Type> catch_nLL(n_fleets); 
  vector<Type> fish_index_nLL(n_fish_indices); 
  vector<Type> srv_index_nLL(n_srv_indices); 
  matrix<Type> fish_comp_nLL(n_fish_comps, n_sexes); 
  matrix<Type> srv_comp_nLL(n_srv_comps, n_sexes);
  vector<Type> rec_nLL(2); // Recruitment likelihood penalty
  Type jnLL = 0; // Joint Negative log Likelihood
  
  // TESTING
  DATA_VECTOR(Init_N); 
  

  // MODEL STRUCTURE ----------------------------------------------
  // y = year, a = age, s = sex, f = fishery fleet, sf = survey fleet
  
  // Selectivity ----------------------------------------------
  
  // Fishery
  for(int f = 0; f < n_fleets; f++) {
    for(int s = 0; s < n_sexes; s++) {
      for(int y = 0; y < n_years; y++) {
        for(int a = 0; a < n_ages; a++) {
          
          if(F_Slx_model(f, y, s) == 0) { // time invariant logistic
            F_Slx(y, a, f, s) = Type(1.0) / (Type(1.0) + exp(-(ages(a) - exp(ln_a50_f) ) /exp(ln_k_f) ));
          } // time invariant logistic
          
        } // a loop
      } // y loop
    } // s loop
  } // f loop
  
  // Survey
  for(int si = 0; si < n_srv_indices; si++) {
    for(int s = 0; s < n_sexes; s++) {
      for(int y = 0; y < n_years; y++) {
        for(int a = 0; a < n_ages; a++) {
          
          if(S_Slx_model(si, y, s) == 0) { // time invariant logistic
            S_Slx(y, a, si, s) = Type(1.0) / (Type(1.0) + exp(-(ages(a) - exp(ln_a50_s) ) / exp(ln_k_s) ));
          } // time invariant logistic
          
        } // a loop
      } // y loop
    } // s loop
  } // f loop

  // Deaths and Removals (Fishing Mortality and Natural Mortality) --------------------------------
  
  // Get total fishing mortality
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fleets; f++) {
      
      ln_F(y, f) = ln_Fy(y, f); // F in log space
      F(y, f) = exp(ln_F(y,f)); // F in normal space
      // Increment F by fleet to get total
      Total_Fy(y) += F(y, f);
      
    for(int s = 0; s < n_sexes; s ++) {
      for(int a = 0; a < n_ages; a++) {
        // Calculate F_at_age
        FAA(y, a, f, s) = F(y, f) * F_Slx(y, a, f, s);
        // Increment to add FAA to ZAA
        ZAA(y, a, s) += FAA(y, a, f, s);
      } // a loop
    } // s loop
    
    } // f loop
  } // y loop

  // Finish calculating Total Mortality and Survival
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s ++) {
        // Increment M to ZAA
        ZAA(y, a, s) += M;
       // Calculate survival
       SAA(y, a, s) = exp(Type(-1.0) * ZAA(y, a, s));
      } // a loop
    } // s loop
  } // y loop

  // Initialization ----------------------------------------------
  
  for(int s = 0; s < n_sexes; s++) {
    for(int a = 1; a < n_ages; a++){
      if(a != n_ages - 1) { // not plus group
        NAA(0, a, s) = exp(ln_MeanRec + ln_N1_Devs(a - 1) -M * Type(a) -(0.5 * ln_SigmaRec2) ) * Sex_Ratio(s);
      } else{
        NAA(0, n_ages - 1, s) = (exp(ln_MeanRec -M * Type( n_ages - 1) ) / (1 - exp(-M)) ) * Sex_Ratio(s);
      }  // Plus group calculation for initializing population
    } // n_ages - 1 a loop
  } // s loop
  
  // TESTING
  // for(int a = 0; a < n_ages; a++) NAA(0, a, 0) = Init_N(a);
  
  // Recruitment ----------------------------------------------
  
  if(rec_model == 0) { // Mean Recruitment 
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        NAA(y, 0, s) = exp( ln_MeanRec + ln_RecDevs(y) -(0.5 * ln_SigmaRec2)) * Sex_Ratio(s); 
      } // s loop
    } // y loop
  } // if for mean recruitment

  // Get total recruitment
  for(int y = 0; y < n_years; y++) for(int s = 0; s < n_sexes; s++) Total_Rec(y) += NAA(y, 0, s);
  
  // Project Numbers At Age Forward ----------------------------------------------
  
  for(int y = 0; y < n_years; y++) {
    for(int s = 0; s < n_sexes; s++) {
      for(int a = 0; a < n_ages; a++) {
        
        // Project ages and years forward
        if(a != (n_ages - 1)) { // Not in Plus Group
          NAA(y + 1, a + 1, s) = NAA(y, a, s) * SAA(y, a, s); 
        } else{ // Add Plus Group
          NAA(y + 1, n_ages - 1, s) = NAA(y + 1, n_ages - 1, s) + (NAA(y, n_ages - 1, s) * 
                                      SAA(y, n_ages - 1, s));
        } // plus group calculation
        
        // Increment Numbers at age to get total biomass
        Total_Biom(y) += NAA(y, a, s) * WAA(y, a, 0);
      } // end ages loop
    } // end sex loop
  } // end year loop
  
  // Sum to get SSB
  for(int y = 0; y < n_years; y++)
    for(int a = 0; a < n_ages; a++) 
      SSB(y) += NAA(y, a, 0) * MatAA(y, a, 0) * WAA(y, a, 0);
      
  // Catch ----------------------------------------------
  
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fleets; f++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Baranov's Catch Equation
          CAA(y, a, f, s) = NAA(y, a, s) * FAA(y, a, f, s) * (Type(1.0) - SAA(y, a, s)) / ZAA(y, a, s);
          // Get Aggregated Catch - Increment catch in biomass
          pred_catches(y, f) += CAA(y, a, f, s) * WAA(y, a, s);
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Indices of Abundance ----------------------------------------------
  // (Assumes that indices are observed at the start of each year)
  
  // Fishery Index of Abundance
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          pred_fish_indices(y, fi) += exp(ln_q_fish(fi)) * NAA(y, a, s) * WAA(y, a, s) * F_Slx(y, a, fi, s);
        } // a loop
      } // s loop
    } // fi loop
  } // y loop
  
  // Survey Index of Abundance
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          pred_srv_indices(y, si) +=  exp(ln_q_srv(si)) * NAA(y, a, s) * S_Slx(y, a, si, s);
        } // a loop
      } // s loop
    } // si loop
  } // y loop
  
  // Compositions ----------------------------------------------  
  
  // Fishery Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int fc = 0; fc < n_fish_comps; fc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Get predicted comps here prior to normalizing w/ catch at age
          pred_fish_age_comps(y, a, fc, s) = CAA(y, a, fc, s);
          // Increment to get total numbers at age for a given fleet
          Total_Fishery_Numbers(y, fc, s) += pred_fish_age_comps(y, a, fc, s);
        } // a loop
      } // s loop
    } // fc loop
  } // y loop
  
  // Normalize to sum to 1
  for(int y = 0; y < n_years; y++) 
    for(int fc = 0; fc < n_fish_comps; fc++) 
      for(int s = 0; s < n_sexes; s++) 
        for(int a = 0; a < n_ages; a++) 
          pred_fish_age_comps(y, a, fc, s) /= Total_Fishery_Numbers(y, fc, s);
  
  // Survey Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int sc = 0; sc < n_srv_comps; sc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Get predicted comps here prior to normalizing
          pred_srv_age_comps(y, a, sc, s) = NAA(y, a, s) * S_Slx(y, a, sc, s);
          // Increment to get total numbers at age for a given fleet
          Total_Survey_Numbers(y, sc, s) += pred_srv_age_comps(y, a, sc, s);
        } // a loop
      } // s loop
    } // sc loop
  } // y loop
  
  // Normalize to sum to 1
  for(int y = 0; y < n_years; y++) 
    for(int sc = 0; sc < n_srv_comps; sc++) 
      for(int s = 0; s < n_sexes; s++) 
        for(int a = 0; a < n_ages; a++) 
          pred_srv_age_comps(y, a, sc, s) /= Total_Survey_Numbers(y, sc, s);
  
  // Likelihoods ----------------------------------------------
  
  // Catch likelihood (Log-normal likelihood) ----------------------------------------------
  
  // Catch observed w/ minimal error
  for(int y = 0; y < n_years; y ++) { 
    for(int f = 0; f < n_fleets; f++) {
      catch_nLL(f) -= dnorm(log(obs_catches(y, f)), log(pred_catches(y, f)), Type(0.01), true);
      
      SIMULATE{ // Simulate catch
        obs_catches(y, f) = exp(rnorm(log(pred_catches(y, f)), Type(0.01) ));
      } // Simulation statement
      
    } // f loop
  } // y loop
  
  // Index likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> fish_sd(n_fish_indices);   // Convert fishery index CV to standard deviation
  vector<Type> srv_sd(n_srv_indices);   // Convert survey index CV to standard deviation
  for(int fi = 0; fi < n_fish_indices; fi++) fish_sd(fi) = sqrt(log( (2*fish_cv(fi)) + 1));
  for(int si = 0; si < n_srv_indices; si++) srv_sd(si) = sqrt(log( (2*srv_cv(si)) + 1));
  
  // Likelihood for fishery index
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      // Likelihood calculations
      fish_index_nLL(fi) -= dnorm(log(obs_fish_indices(y, fi)), log(pred_fish_indices(y, fi)), fish_sd(fi), true);
      
      SIMULATE{ // Simulate Fishery Index
        obs_fish_indices(y, fi) = exp(rnorm(log(pred_fish_indices(y, fi)), fish_sd(fi))); 
      } // Simulation statement
      
    } // fi loop
  } // y loop
  
  // Likelihood for survey index
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      // Likelihood calculations
      srv_index_nLL(si) -= dnorm(log(obs_srv_indices(y, si)), log(pred_srv_indices(y, si)), srv_sd(si), true); 
      
      SIMULATE{ // Simulate Survey Index
        obs_srv_indices(y, si) = exp(rnorm(log(pred_srv_indices(y, si)), srv_sd(si))); 
      } // Simulation statement
      
    } // si loop
  } // y loop
  
  // Composition likelihoods (Multinomial likelihood) ----------------------------------------------
  
  Type c = 1e-03; // Constant to add to multinomial
  
  // Fishery Compositions
  vector<Type> obs_F_vec(n_ages); // Obs fishery vector to hold and pass values to nLL
  vector<Type> pred_F_vec(n_ages); // Pred fishery vector to hold and pass values to nLL
  
  for(int y = 0; y < n_years; y++) {
    for(int fc = 0; fc < n_fish_comps; fc++) {
      for(int s = 0; s < n_sexes; s++) { 
        // Pull out observed age vector
        obs_F_vec = obs_fish_age_comps.col(s).col(fc).transpose().col(y) + c;
        // Pull out predicted age vector
        pred_F_vec = pred_fish_age_comps.col(s).col(fc).transpose().col(y) + c;
        // Evaluate log-likelihood
        fish_comp_nLL(fc, s) -= dmultinom(obs_F_vec.vec(), pred_F_vec.vec(), true);
      } // s loop
  } // fc loop
} // y loop
  
  vector<Type> obs_S_vec(n_ages); // Obs survey vector to hold and pass values to nLL
  vector<Type> pred_S_vec(n_ages); // Pred survey vector to hold and pass values to nLL

  for(int y = 0; y < n_years; y++) {
    for(int sc = 0; sc < n_srv_comps; sc++) {
      for(int s = 0; s < n_sexes; s++) {
        // Pull out observed age vector
        obs_S_vec = obs_srv_age_comps.col(s).col(sc).transpose().col(y) + c;
        // Pull out predicted age vector
        pred_S_vec = pred_srv_age_comps.col(s).col(sc).transpose().col(y) + c;
        // Evaluate log-likelihood
        srv_comp_nLL(sc, s) -=  dmultinom(obs_S_vec, pred_S_vec, true);
      } // s loop
    } // fc loop
  } // y loop
  
  // Recruitment related stuff (likelihoods + derived quantites) ----------------------------------------------
  
  for(int y = 0; y < ln_N1_Devs.size(); y++) { // Mean = log-normal correction
    rec_nLL(0) -= dnorm(ln_N1_Devs(y), Type(0), exp(ln_SigmaRec), true);
  } // Penalty for initial recruitment

  for(int y = 0; y < ln_RecDevs.size(); y++) { // Mean = log-normal correction
    rec_nLL(1) -= dnorm(ln_RecDevs(y), Type(0), exp(ln_SigmaRec), true);
  } // Penalty for recruitment (mean should be 0)
  
  // Add to joint nLL
  jnLL = rec_nLL.sum() + srv_comp_nLL.sum() + fish_comp_nLL.sum() + 
    srv_index_nLL.sum() + fish_index_nLL.sum() + catch_nLL.sum();
  
  // REPORT SECTION ----------------------------------------------
  REPORT(NAA); // Numbers at age
  REPORT(ZAA); // Total Mortality
  REPORT(FAA); // Fishing Mortality
  REPORT(CAA); // Catch at Age
  REPORT(pred_catches); // Aggregate Catch by fleet
  REPORT(SSB); // Spawning Stock Biomass
  REPORT(pred_srv_indices); // Survey Indices
  REPORT(pred_fish_indices); // Fishery Indices
  REPORT(F_Slx); // Fishery Selectivity
  REPORT(S_Slx); // Survey Selectivity
  REPORT(pred_fish_age_comps); // Predicted fishery age comps
  REPORT(pred_srv_age_comps); // Predicted survey age comps
  REPORT(ln_q_fish); // catchability for fishery
  REPORT(ln_q_srv); // catchability for survey
    
  //  Likelihoods
  REPORT(catch_nLL);
  REPORT(srv_index_nLL);
  REPORT(fish_index_nLL);
  REPORT(fish_comp_nLL);
  REPORT(srv_comp_nLL);
  REPORT(rec_nLL);
  REPORT(jnLL);
  
  // ADREPORT 
  ADREPORT(SSB);
  ADREPORT(Total_Fy);
  ADREPORT(Total_Rec);
  ADREPORT(Total_Biom);
  ADREPORT(pred_srv_indices);
  ADREPORT(pred_fish_indices);
  ADREPORT(pred_catches);

  return jnLL;
  
} // end objective function

