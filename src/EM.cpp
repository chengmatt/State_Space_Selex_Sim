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
  DATA_VECTOR(catch_cv); // Vector of catch CVs from each fleet; n_fleets
  DATA_ARRAY(obs_fish_age_comps); // Array of fishery age comps; n_years * n_ages * n_fleets * n_sexes
  DATA_ARRAY(obs_srv_age_comps); // Array of survey age comps; n_years * n_ages * n_srv_indices * n_sexes
  DATA_ARRAY(obs_fish_age_Neff); // Array of fishery age comps; n_years * n_ages * n_fleets * n_sexes
  DATA_ARRAY(obs_srv_age_Neff); // Array of survey age comps; n_years * n_ages * n_srv_indices * n_sexes
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
  DATA_VECTOR(lb_q_fish); // Lower bound for inverse logit fishery catchability
  DATA_VECTOR(ub_q_fish); // Upper bound for inverse logit fishery catchability
  DATA_VECTOR(lb_q_srv); // Lower bound for inverse logit survey catchability
  DATA_VECTOR(ub_q_srv); // Upper bound for inverse logit survey catchability
  
  // PARAMETER SECTION ----------------------------------------------
  PARAMETER_VECTOR(ln_N1_Devs); // log deviations for initial abundance
  PARAMETER(ln_MeanRec); // log mean recruitment 
  PARAMETER(ln_SigmaRec); // log sigma for recruitment
  PARAMETER_VECTOR(ln_RecDevs); // log recruitment deviations
  PARAMETER(ln_M); // log natural mortality
  PARAMETER_MATRIX(ln_Fy); // Annual Fishing Mortality ; n_years * n_fleets
  PARAMETER_VECTOR(logit_q_fish); // logit catchability for fishery; n_fish_indices (bounded between 0 and 1)
  PARAMETER_VECTOR(logit_q_srv); // logit catchability for survey; n_srv_indices (bounded between 0 and 1)
  PARAMETER(ln_a50_f); // a50 logistic fishery
  PARAMETER(ln_k_f); // k logistic fishery
  PARAMETER(ln_a50_s); // a50 logistic survey
  PARAMETER(ln_k_s); // k logistic survey
  
  // Parameter Transformations ----------------------------------------------
  Type M = exp(ln_M); // Natural Mortality

  // Predicted Quantities
  array<Type> pred_fish_age_comps(obs_fish_age_comps.dim); // Predicted Fishery Age Comps
  array<Type> pred_srv_age_comps(obs_srv_age_comps.dim); // Predicted survey age comps
  matrix<Type> pred_fish_indices(n_years, n_fish_indices); // Predicted fishery indices
  matrix<Type> pred_srv_indices(n_years, n_srv_indices); // Predicted survey indices
  matrix<Type> pred_catches(n_years, n_fleets); // Predicted fishery catches
  
  // Stored Quantities
  array<Type> NAA(n_years + 1, n_ages, n_sexes); // Numbers at age; n_years + 1 because year 0 is for numbers at age prior to fishing
  array<Type> ZAA(n_years, n_ages, n_sexes); // Total Mortality
  array<Type> FAA(n_years, n_ages, n_fleets, n_sexes); // Fishing Mortality
  array<Type> CAA(n_years, n_ages, n_fleets, n_sexes); // Catch at Age
  
  array<Type> Fish_F_yas(n_years, n_ages, n_sexes); // Total fishing mortality (sum of fleets)
  vector<Type> Total_Fy(n_years); // Get sum of estimated Fs for all fleets
  vector<Type> Total_Rec(n_years); // Total Recruitment
  vector<Type> SSB(n_years + 1); // Spawning stock biomass; n_years + 1 because year 0 is for numbers at age prior to fishing

  // Selectivities
  array<Type> F_Slx(n_years, n_ages, n_fleets, n_sexes); // Fishery Selectivities
  array<Type> Total_Fishery_Numbers(n_years, n_fleets, n_sexes); // Store Total Fishery Numbers for Proportions
  array<Type> S_Slx(n_years, n_ages, n_srv_indices, n_sexes); // Survey Selectivities
  array<Type> Total_Survey_Numbers(n_years, n_fleets, n_sexes); // Store Total Survey Numbers for Proportions

  // Cactchabilities
  vector<Type> invlogit_q_fish(n_fish_indices);
  vector<Type> invlogit_q_srv(n_srv_indices);
  
  // Define objects to store negative log-likelihoods
  vector<Type> catch_nLL(n_fleets); 
  vector<Type> fish_index_nLL(n_fish_indices); 
  vector<Type> srv_index_nLL(n_srv_indices); 
  matrix<Type> fish_comp_nLL(n_fish_comps, n_sexes); 
  matrix<Type> srv_comp_nLL(n_srv_comps, n_sexes);
  Type rec_nLL = 0; // Recruitment likelihood
  Type jnLL = 0; // Joint Negative log Likelihood


  // MODEL STRUCTURE ----------------------------------------------
  // y = year, a = age, s = sex, f = fishery fleet, sf = survey fleet
  
  // Initialization ----------------------------------------------

    for(int s = 0; s < n_sexes; s++) {
      for(int a = 0; a < n_ages; a++){
        if(a < n_ages - 1) { // not plus group
          NAA(0, a, s) = exp(ln_MeanRec + ln_N1_Devs(a) -M * ages(a) ) * Sex_Ratio(s);
        } else{
          NAA(0, a, s) = (exp(ln_MeanRec -M * (ages(a) ) ) / (1 - exp(-M))) * Sex_Ratio(s);
        } // plus group
      } // a loop
    } // s loop

  // Get B0 / SSB at Time 0 
  for(int a = 0; a < n_ages; a++) SSB(0) += NAA(0,a,0) * MatAA(0,a,0) * WAA(0,a,0);

  
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
  
  // Deaths and Removals (Fishing Mortality and Natural Mortality) ---------------------------------
  
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fleets; f ++) {
      for(int s = 0; s < n_sexes; s ++) {
        for(int a = 0; a < n_ages; a++) {
          // Calculate F_at_age
          FAA(y, a, f, s) = (exp(ln_Fy(y, f) )) * F_Slx(y, a, f, s);
          // Add fishing mortality to get Total F by year age and sex
          Fish_F_yas(y, a, s) += FAA(y, a, f, s); 
          // Get total mortality by year, age, and sex 
          ZAA(y, a, s) +=  Fish_F_yas(y, a, s) + M; 
        } // a loop
      } // s loop
         // Get total fishing mortality here
         Total_Fy(y) += exp( ln_Fy(y, f));
    } // f loop
  } // y loop
  
  
  // Project Numbers At Age Forward ----------------------------------------------
  
  for(int y = 1; y <= n_years; y++) {
    for(int s = 0; s < n_sexes; s++) {
      
      // Recruitment ----------------------------------------------
      if(rec_model == 0) { // Mean Recruitment 
        if(y < n_years) NAA(y, 0, s) = exp( ln_MeanRec + ln_RecDevs(y) ) * Sex_Ratio(s); 
      } 
      
      for(int a = 1; a < n_ages; a++) { // Starting at age 1
        
        if(a != (n_ages - 1)) { // Not in Plus Group
          // Project ages forward in time 
          NAA(y, a, s) = NAA(y - 1, a - 1, s) * (exp(-ZAA(y - 1, a - 1, s)) ); 
        } else{ // Add Plus Group
          NAA(y, a, s)  = NAA(y - 1, a - 1, s) * exp(-ZAA(y - 1, a - 1, s)) +
            NAA(y - 1, a, s) * exp(-ZAA(y - 1, a, s));
        } // plus group calculation
        
        // Calculate SSB when done w/ + group calculation for sex 0 (females)
        if(a == n_ages - 1 && s == 0) { 
          for(int a = 0; a < n_ages; a++) { 
            SSB(y) += NAA(y, a, s) * MatAA(y, a, s) * WAA(y, a, s);
          } // inner age loop
        } // if statement
        
      } // end ages loop
    } // end sex loop
  } // end year loop
  
  
  // Catch ----------------------------------------------
  
  for(int y = 1; y <= n_years; y++) {
    for(int f = 0; f < n_fleets; f++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Baranov's Catch Equation
          CAA(y - 1, a, f, s) = ((FAA(y - 1, a, f, s) / ZAA(y - 1, a, s)) * NAA(y - 1, a, s) * 
                                (Type(1.0) - exp(-ZAA(y - 1, a, s))) * WAA(y - 1, a, s));
          // Get Aggregated Catch - Increment catch in biomass
          pred_catches(y - 1, f) += CAA(y - 1, a, f, s);
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Indices of Abundance ----------------------------------------------
  // (Assumes that indices are observed at the start of each year)
  
  // Fishery Index of Abundance
  for(int fi = 0; fi < n_fish_indices; fi++) { // Inverse logit transform
    invlogit_q_fish(fi) = lb_q_fish(fi) + (ub_q_fish(fi) - lb_q_fish(fi)) * 
                          (1 / (1 + exp(-logit_q_fish(fi)))); 
  } // end fi loop
  
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          pred_fish_indices(y, fi) += NAA(y, a, s) * WAA(y, a, s) * F_Slx(y, a, fi, s);
        } // a loop
      } // s loop
         // Inverse logit q and apply to index after summing
         pred_fish_indices(y, fi) = invlogit_q_fish(fi) * pred_fish_indices(y, fi); 
    } // fi loop
  } // y loop
  
  // Survey Index of Abundance
  for(int si = 0; si < n_srv_indices; si++) { // Inverse logit transform
    invlogit_q_srv(si) = lb_q_srv(si) + (ub_q_srv(si) - lb_q_srv(si)) * 
                         (1 / (1 + exp(-logit_q_srv(si)))); 
  } // end si loop
  
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          pred_srv_indices(y, si) +=  NAA(y, a, s) * S_Slx(y, a, si, s);
        } // a loop
      } // s loop
         // Inverse logit q and apply to index after summing
         pred_srv_indices(y, si) = invlogit_q_srv(si) * pred_srv_indices(y, si);
    } // si loop
  } // y loop

    
  // Compositions ----------------------------------------------  
  
  // Fishery Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int fc = 0; fc < n_fish_comps; fc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Get predicted comps here prior to normalizing
          pred_fish_age_comps(y, a, fc, s) = NAA(y, a, s) * F_Slx(y, a, fc, s);
          // Increment to get total numbers at age for a given fleet
          Total_Fishery_Numbers(y, fc, s) += pred_fish_age_comps(y, a, fc, s);
          
          if(a == n_ages - 1) { 
            for(int a = 0; a < n_ages; a++) { // Normalize to sum to 1
              pred_fish_age_comps(y, a, fc, s) /= Total_Fishery_Numbers(y, fc, s);
            } // inner a loop
          } // if n_ages - 1
          
        } // a loop
      } // s loop
    } // fc loop
  } // y loop
  
  
  // Survey Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int sc = 0; sc < n_srv_comps; sc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Get predicted comps here prior to normalizing
          pred_srv_age_comps(y, a, sc, s) = NAA(y, a, s) * S_Slx(y, a, sc, s);
          // Increment to get total numbers at age for a given fleet
          Total_Survey_Numbers(y, sc, s) += pred_srv_age_comps(y, a, sc, s);
          
          // Normalize to sum to 1
          if(a == n_ages - 1) {
            for(int a = 0; a < n_ages; a++) { 
              pred_srv_age_comps(y, a, sc, s) /= Total_Survey_Numbers(y, sc, s);
            } // inner age loop
          } // if n_ages - 1
          
        } // a loop
      } // s loop
    } // sc loop
  } // y loop
  
  // Likelihoods ----------------------------------------------
  
  // Catch likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> catch_sd(n_fleets);   // Convert catch CV to standard deviation
  for(int f = 0; f < n_fleets; f++) catch_sd(f) = sqrt(log(pow(catch_cv(f), Type(2.0)) + 1));

  for(int y = 0; y < n_years; y ++) { 
    for(int f = 0; f < n_fleets; f++) {
      catch_nLL(f) -= dnorm(log(obs_catches(y, f)), log(pred_catches(y, f)), catch_sd(f), true);
    } // f loop
  } // y loop
  
  SIMULATE{ // Simulate catch
    for(int y = 0; y < n_years; y ++) {
      for(int f = 0; f < n_fleets; f++) {
        obs_catches(y, f) = exp(rnorm(log(pred_catches(y, f)), catch_sd(f) ));
      } // f loop
    } // y loop
  } // end simulate
  
  // Index likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> fish_sd(n_fish_indices);   // Convert fishery index CV to standard deviation
  vector<Type> srv_sd(n_srv_indices);   // Convert survey index CV to standard deviation
  for(int fi = 0; fi < n_fish_indices; fi++) fish_sd(fi) = sqrt(log(pow(fish_cv(fi), Type(2.0)) + 1));
  for(int si = 0; si < n_srv_indices; si++) srv_sd(si) = sqrt(log(pow(srv_cv(si), Type(2.0)) + 1));
  
  // Likelihood for fishery index
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      fish_index_nLL(fi) -= dnorm(log(obs_fish_indices(y, fi)), log(pred_fish_indices(y, fi)),
                            fish_sd(fi), true); 
    } // fi loop
  } // y loop
  
  SIMULATE{ // Simulate Fishery Index
    for(int y = 0; y < n_years; y++) {
      for(int fi = 0; fi < n_fish_indices; fi++) {
        obs_fish_indices(y, fi) = exp(rnorm(log(pred_fish_indices(y, fi)), fish_sd(fi))); 
      } // fi loop
    } // y loop
  } // end simulate for fishery index
  
  // Likelihood for survey index
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      srv_index_nLL(si) -= dnorm(log(obs_srv_indices(y, si)), log(pred_srv_indices(y, si)), srv_sd(si), true); 
    } // si loop
  } // y loop
  
  SIMULATE{ // Simulate Survey Index
    for(int y = 0; y < n_years; y++) {
      for(int si = 0; si < n_srv_indices; si++) {
        obs_srv_indices(y, si) = exp(rnorm(log(pred_srv_indices(y, si)), srv_sd(si))); 
      } // si loop
    } // y loop
  } // end simulate for survey index
  
  // Composition likelihoods (Multinomial likelihood) ----------------------------------------------
  
  // Fishery Compositions
  vector<Type> obs_F_vec(n_ages); // Obs fishery vector to hold and pass values to nLL
  vector<Type> pred_F_vec(n_ages); // Pred fishery vector to hold and pass values to nLL
  
  for(int y = 0; y < n_years; y++) {
    for(int fc = 0; fc < n_fish_comps; fc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          obs_F_vec(a) = obs_fish_age_comps(y, a, fc, s);
          pred_F_vec(a) = pred_fish_age_comps(y, a, fc, s);
      } // a loop
        // true = gives log probability density - calculate nLL here
        fish_comp_nLL(fc, s) -= dmultinom(obs_F_vec, pred_F_vec, true);
    } // s loop
  } // fc loop
} // y loop
  
  vector<Type> obs_S_vec(n_ages); // Obs survey vector to hold and pass values to nLL
  vector<Type> pred_S_vec(n_ages); // Pred survey vector to hold and pass values to nLL
  
  for(int y = 0; y < n_years; y++) {
    for(int sc = 0; sc < n_srv_comps; sc++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          obs_S_vec(a) = obs_srv_age_comps(y, a, sc, s);
          pred_S_vec(a) = pred_srv_age_comps(y, a, sc, s);
        } // a loop
        // true = gives log probability density - calculate nLL here
        srv_comp_nLL(sc, s) -= dmultinom(obs_S_vec, pred_S_vec, true);
      } // s loop
    } // fc loop
  } // y loop
  
  // Recruitment related stuff (likelihoods + derived quantites) ----------------------------------------------
  
  // Get total recruitment
  for(int y = 0; y < n_years; y++) for(int s = 0; s < n_sexes; s++) Total_Rec(y) += NAA(y, 0, s);

  for(int y = 0; y < ln_RecDevs.size(); y++) {
    rec_nLL -= dnorm(ln_RecDevs(y), Type(0), ln_SigmaRec, true);
  } // Penalty for recruitment
  
  for(int y = 0; y < ln_N1_Devs.size(); y++) {
    rec_nLL -= dnorm(ln_N1_Devs(y), Type(0), ln_SigmaRec, true);
  } // Penalty for initial recruitment
  
  // Sum Joint negative log-likelihood
  for(int f = 0; f < n_fleets; f++) jnLL += catch_nLL(f); // Catch nLL
  for(int fi = 0; fi < n_fish_indices; fi++) jnLL += fish_index_nLL(fi); // Fish Index nLL
  for(int si = 0; si < n_srv_indices; si++) jnLL += srv_index_nLL(si); // Survey Index nLL
  for(int fc = 0; fc < n_fish_comps; fc++) for(int s = 0; s < n_sexes; s++) jnLL += fish_comp_nLL(fc, s); // Fishery Comps nLL
  for(int sc = 0; sc < n_srv_comps; sc++) for(int s = 0; s < n_sexes; s++) jnLL += srv_comp_nLL(sc, s); // Survey Comps nLL
  jnLL += rec_nLL; // Add recruitment penalty
  
  // REPORT SECTION ----------------------------------------------
  REPORT(NAA); // Numbers at age
  REPORT(ZAA); // Total Mortality
  REPORT(FAA); // Fishing Mortality
  REPORT(CAA); // Catch at Age
  REPORT(pred_catches); // Aggregate Catch by fleet
  REPORT(Fish_F_yas); // Total Fishing Mortality summed across fleets
  REPORT(SSB); // Spawning Stock Biomass
  REPORT(pred_srv_indices); // Survey Indices
  REPORT(pred_fish_indices); // Fishery Indices
  REPORT(F_Slx); // Fishery Selectivity
  REPORT(S_Slx); // Survey Selectivity
  REPORT(pred_fish_age_comps); // Predicted fishery age comps
  REPORT(pred_srv_age_comps); // Predicted survey age comps
  
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
  ADREPORT(pred_srv_indices);
  ADREPORT(pred_fish_indices);
  ADREPORT(pred_catches);

  return jnLL;
  
} // end objective function

