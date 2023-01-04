// General single species age-structured stock assessment
// that accommodates multiple fishery fleets written in TMB
// Creator: Matthew LH. Cheng (UAF-CFOS)

// TO DO:
// Beverton Holt SR
// Set Selectivities

#include<TMB.hpp>
#include "Get_Selex.hpp"

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
  DATA_VECTOR(catch_cv); // Vector of catch CVs
  DATA_VECTOR(fish_cv); // Vector of fishery index CVs
  DATA_VECTOR(srv_cv); // Vector of survey index CVs
  
  // Biological Information Inputs ----------------------------------------------
  DATA_ARRAY(WAA); // Weight-at-age array; n_years * n_ages * n_sexes
  DATA_ARRAY(MatAA); // Maturity-at-age array; n_years * n_ages * n_sexes (But only Sex0 used for calcs)
  DATA_VECTOR(Sex_Ratio); // Sex ratio; n_sexes - Females, Males
  
  // Controls on assessment ----------------------------------------------
  DATA_INTEGER(rec_model); // Recruitment model, == 0 Mean Recruitment, == 1 Beverton Holt 
  DATA_IVECTOR(F_Slx_model); // Fishery Selectivity Model, == 0 Logistic n_fleets
  DATA_IMATRIX(F_Slx_Blocks); // Fishery Selectivity Time Blocks, n_years * n_fish_fleets; 
  // this is set up such that the selectivity within a fleet and across sexes is constant 
  DATA_IVECTOR(S_Slx_model); // Survey Selectivity Model, == 0 Logistic n_fleets 
  DATA_IMATRIX(S_Slx_Blocks); // Survey Selectivity Time Blocks, n_years * n_srv_fleets; 
  
  // PARAMETER SECTION ----------------------------------------------
 
  // Biological parameters
  PARAMETER(ln_M); // log natural mortality
  PARAMETER(ln_MeanRec); // log mean recruitment 
  PARAMETER(ln_SigmaRec); // log sigma for recruitment
  PARAMETER_VECTOR(ln_N1_Devs); // log deviations for initial abundance
  PARAMETER_VECTOR(ln_RecDevs); // log recruitment deviations
  
  // Indices of abundance
  PARAMETER_VECTOR(ln_q_fish); // log catchability for fishery; n_fish_indices 
  PARAMETER_VECTOR(ln_q_srv); // log catchability for survey; n_srv_indices 
  
  // Fishery and selectivity parameters
  PARAMETER_MATRIX(ln_Fy); // Annual Fishing Mortality ; n_years * n_fleets
  PARAMETER_ARRAY(ln_fish_selpars); // Fishery Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  PARAMETER_ARRAY(ln_srv_selpars); // Survey Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  
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
  vector<Type> SSB(n_years); // Spawning stock biomass; n_years + 1 (forward projection)
  vector<Type> Depletion(n_years); // Depletion SSB / SSB(0)

  // Selectivities
  array<Type> F_Slx(n_years, n_ages, n_fish_comps, n_sexes); // Fishery Selectivities
  array<Type> Total_Fishery_Numbers(n_years, n_fish_comps, n_sexes); // Store Total Fishery Numbers for Proportions
  array<Type> S_Slx(n_years, n_ages, n_srv_comps, n_sexes); // Survey Selectivities
  array<Type> Total_Survey_Numbers(n_years, n_srv_comps, n_sexes); // Store Total Survey Numbers for Proportions

  // Define objects to store negative log-likelihoods
  vector<Type> catch_nLL(n_fleets); 
  vector<Type> fish_index_nLL(n_fish_indices); 
  vector<Type> srv_index_nLL(n_srv_indices); 
  array<Type> fish_comp_nLL(n_years, n_fish_comps, n_sexes); 
  array<Type> srv_comp_nLL(n_years, n_srv_comps, n_sexes);
  vector<Type> rec_nLL(2); // Recruitment likelihood penalty 
  // (index 1 = penalty for initial recruitment, 2 = penalty for recruitment deviations)
  Type jnLL = 0; // Joint Negative log Likelihood
  
  // TESTING
  DATA_VECTOR(Init_N); 
  

  // MODEL STRUCTURE ----------------------------------------------
  // y = year, a = age, s = sex, f = fishery fleet, sf = survey fleet
  
  // Selectivity ----------------------------------------------
  
  // Fishery
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fish_comps; f++) {

      // Index fishery blocks here
      int b = F_Slx_Blocks(y, f);

      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // a + 1 because TMB indexes starting at 0
          F_Slx(y,a,f,s) = Get_Selex(a + 1, F_Slx_model(f), 
                ln_fish_selpars.transpose().col(f).col(s).col(b).vec()); 
          
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_srv_comps; f++) {
      
      // Index fishery blocks here
      int b = S_Slx_Blocks(y, f);
      
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // a + 1 because TMB indexes starting at 0
          S_Slx(y,a,f,s) = Get_Selex(a + 1, S_Slx_model(f), 
                ln_srv_selpars.transpose().col(f).col(s).col(b).vec()); 
          // transposing selectivity array to coerce to vector
             
        } // a loop       
      } // s loop
    } // f loop
  } // y loop

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
        //(incrementing rather than transposing and summing is faster)
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
  
  for(int y = 0; y < n_years; y++) {
    for(int s = 0; s < n_sexes; s++) {
      
    if(rec_model == 0) { // Mean Recruitment 
    NAA(y, 0, s) = exp( ln_MeanRec + ln_RecDevs(y) -(0.5 * ln_SigmaRec2)) * Sex_Ratio(s); 
    } // if for rec_model == 0
  
    } // s loop
  } // y loop

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
        
        if(s == 0) { // SSB calculation here if sex = 0 (females)
          SSB(y) += NAA(y, a, 0) * MatAA(y, a, 0) * WAA(y, a, 0);
        } // if statement for ssb calculations
        
        // Increment Numbers at age to get total biomass
        Total_Biom(y) += NAA(y, a, s) * WAA(y, a, 0);
      } // end ages loop
      
      // Increment total recruitment here
      Total_Rec(y) += NAA(y, 0, s);
      
    } // end sex loop
    
    // Calculate depletion rates once we are done 
    Depletion(y) = SSB(y) / SSB(0); 
    
  } // end year loop
  
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
          pred_fish_indices(y, fi) += NAA(y, a, s) * WAA(y, a, s) * F_Slx(y, a, fi, s);
        } // a loop
      } // s loop    
      pred_fish_indices(y, fi) = exp(ln_q_fish(fi)) * pred_fish_indices(y, fi);
    } // fi loop
  } // y loop
  
  // Survey Index of Abundance
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          pred_srv_indices(y, si) += NAA(y, a, s) * S_Slx(y, a, si, s);
        } // a loop
      } // s loop
      pred_srv_indices(y, si) = exp(ln_q_srv(si)) * pred_srv_indices(y, si);
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
  vector<Type> catch_sd(n_fleets); // Empty container to store sd calculation
  for(int f = 0; f < n_fleets; f++) catch_sd(f) = sqrt(log( (catch_cv(f)*catch_cv(f)) + 1));  // Calculate catch cv to sd

  
  // Catch observed w/ minimal error
  for(int y = 0; y < n_years; y ++) { 
    for(int f = 0; f < n_fleets; f++) {
      catch_nLL(f) -= dnorm(log(obs_catches(y, f)), 
                            log(pred_catches(y, f)) - (0.5 * catch_sd(f) * catch_sd(f)), // bias correction
                            catch_sd(f), true);
      
      SIMULATE{ // Simulate catch
        obs_catches(y, f) = exp(rnorm(log(pred_catches(y, f)) - (0.5 * catch_sd(f) * catch_sd(f)), 
                                catch_sd(f) ));
      } // Simulation statement
      
    } // f loop
  } // y loop
  
  
  // Index likelihood (Log-normal likelihood) ----------------------------------------------
  
  vector<Type> fish_sd(n_fish_indices);   // Convert fishery index CV to standard deviation
  vector<Type> srv_sd(n_srv_indices);   // Convert survey index CV to standard deviation
  for(int fi = 0; fi < n_fish_indices; fi++) fish_sd(fi) = sqrt(log( (fish_cv(fi)*fish_cv(fi)) + 1));
  for(int si = 0; si < n_srv_indices; si++) srv_sd(si) = sqrt(log( (srv_cv(si)*srv_cv(si)) + 1));
  
  // Likelihood for fishery index
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      // Likelihood calculations
      fish_index_nLL(fi) -= dnorm(log(obs_fish_indices(y, fi)), 
                            log(pred_fish_indices(y, fi)) - (0.5 * fish_sd(fi) * fish_sd(fi)), // bias correction
                            fish_sd(fi), true);
      
      SIMULATE{ // Simulate Fishery Index
        obs_fish_indices(y, fi) = exp(rnorm(log(pred_fish_indices(y, fi)) - (0.5 * fish_sd(fi) * fish_sd(fi)), 
                                  fish_sd(fi))); 
      } // Simulation statement
      
    } // fi loop
  } // y loop
  
  // Likelihood for survey index
  for(int y = 0; y < n_years; y++) {
    for(int si = 0; si < n_srv_indices; si++) {
      // Likelihood calculations
      srv_index_nLL(si) -= dnorm(log(obs_srv_indices(y, si)), 
                           log(pred_srv_indices(y, si)) - (0.5 * srv_sd(si) * srv_sd(si)), // bias correction
                           srv_sd(si), true); 
      
      SIMULATE{ // Simulate Survey Index
        obs_srv_indices(y, si) = exp(rnorm(log(pred_srv_indices(y, si)) - (0.5 * srv_sd(si) * srv_sd(si)), 
                                 srv_sd(si))); 
      } // Simulation statement
      
    } // si loop
  } // y loop
  
  // Composition likelihoods (Multinomial likelihood) ----------------------------------------------
  
  Type c = 1e-10; // Constant to add to multinomial
  
  // Fishery Compositions
  vector<Type> obs_fish_age_vec(n_ages); // Obs fishery vector to hold and pass values to nLL
  vector<Type> pred_fish_age_vec(n_ages); // Pred fishery vector to hold and pass values to nLL
  
  for(int y = 0 ; y < n_years; y++) {
    for(int fc = 0; fc < n_fish_comps; fc++) {
      for(int s = 0; s < n_sexes; s++) { 
        // Pull out observed age vector and multiply by the effective sample size
        obs_fish_age_vec = obs_fish_age_comps.col(s).col(fc).transpose().col(y) * obs_fish_age_Neff(y, fc);
        // Pull out predicted age vector
        pred_fish_age_vec = (pred_fish_age_comps.col(s).col(fc).transpose().col(y) + c);
        // Evaluate log-likelihood
        fish_comp_nLL(y, fc, s) -= dmultinom(obs_fish_age_vec.vec(), pred_fish_age_vec.vec(), true);
      } // s loop
  } // fc loop
} // y loop
  
  vector<Type> obs_srv_age_vec(n_ages); // Obs survey vector to hold and pass values to nLL
  vector<Type> pred_srv_age_vec(n_ages); // Pred survey vector to hold and pass values to nLL

  for(int y = 0; y < n_years; y++) {
    for(int sc = 0; sc < n_srv_comps; sc++) {
      for(int s = 0; s < n_sexes; s++) {
        // Pull out observed age vector and multiply by the effective sample size
        obs_srv_age_vec = obs_srv_age_comps.col(s).col(sc).transpose().col(y) * obs_srv_age_Neff(y, sc);
        // Pull out predicted age vector
        pred_srv_age_vec = (pred_srv_age_comps.col(s).col(sc).transpose().col(y) + c);
        // Evaluate log-likelihood
        srv_comp_nLL(y, sc, s) -=  dmultinom(obs_srv_age_vec, pred_srv_age_vec, true);
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
  ADREPORT(Depletion); 
  ADREPORT(Total_Fy);
  ADREPORT(Total_Rec);
  ADREPORT(Total_Biom);
  ADREPORT(pred_srv_indices);
  ADREPORT(pred_fish_indices);
  ADREPORT(pred_catches);

  return jnLL;
  
} // end objective function

