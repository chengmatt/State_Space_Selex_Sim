// General single species age-and sex-structured stock assessment
// that accommodates multiple fishery fleets written in TMB
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date updated: 2/8/23

#include<TMB.hpp>
#include "Get_Selex.hpp"

template<class Type>

Type objective_function<Type>::operator() ()
{
  using namespace density; // Define namespace to use SEPARABLE, AR1, SCALE
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
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
  
  // OBSERVATIONS ----------------------------------------------
  
  DATA_MATRIX(obs_catches); // Matrix of catch from each fleet; n_years * n_fleets
  DATA_ARRAY(obs_fish_age_comps); // Array of fishery age comps; n_years * n_ages * n_fleets * n_sexes
  DATA_ARRAY(obs_srv_age_comps); // Array of survey age comps; n_years * n_ages * n_srv_indices * n_sexes
  DATA_ARRAY(obs_fish_age_Neff); // Array of fishery age comps; n_years * n_fleets * n_sexes 
  DATA_ARRAY(obs_srv_age_Neff); // Array of survey age comps; n_years * n_srv_indices * n_sexes 
  DATA_MATRIX(obs_fish_indices); // Matrix of Fishery Indices of Abundance
  DATA_MATRIX(obs_srv_indices); // Matrix of Survey Indices of Abundance
  DATA_VECTOR(catch_cv); // Vector of catch CVs
  DATA_VECTOR(fish_cv); // Vector of fishery index CVs
  DATA_VECTOR(srv_cv); // Vector of survey index CVs
  
  // BIOLOGICAL INPUTS  ----------------------------------------------
  
  DATA_ARRAY(WAA); // Weight-at-age array; n_years * n_ages * n_sexes
  DATA_ARRAY(MatAA); // Maturity-at-age array; n_years * n_ages * n_sexes (But only Sex0 used for calcs)
  DATA_VECTOR(Sex_Ratio); // Sex ratio; n_sexes - Females, Males
  
  // ASSESSMENT CONTROLS ----------------------------------------------
  
  DATA_INTEGER(rec_model); // Recruitment model, == 0 Mean Recruitment
  DATA_IVECTOR(F_Slx_model); // Fishery Selectivity Model, == 0 Logistic n_fleets
  DATA_IMATRIX(F_Slx_Blocks); // Fishery Selectivity Time Blocks, n_years * n_fish_comps; 
  // this is set up such that the selectivity within a fleet and across sexes is constant
  DATA_IVECTOR(S_Slx_model); // Survey Selectivity Model, == 0 Logistic n_fleets 
  DATA_IMATRIX(S_Slx_Blocks); // Survey Selectivity Time Blocks, n_years * n_srv_fleets; 
  
  // Fishery Random Effects Selectivity ------------------------------
  
  DATA_IMATRIX(F_Slx_re_model); // Fishery Selectivity Random Effects Model, 
  // n_fish_comps * n_sexes, == 0 AR(1y) == 1 2DAR(1)
  
  // DATA INDICATORS  ----------------------------------------------
  
  // Indicator 0 == not fitting, 1 == fit
  DATA_IMATRIX(use_catch); // Whether or not to use catch data; n_years x n_fleets
  DATA_IMATRIX(use_fish_index); // Whether or not to use fishery index; n_years x n_fish_indices
  DATA_IMATRIX(use_srv_index); // Whether or not to use survey index; n_years x n_srv_indices
  DATA_IARRAY(use_fish_comps); // Whether or not to use fishery comps; n_years x n_fish_comps x n_sexes
  DATA_IARRAY(use_srv_comps); // Wheter or not to use survey comps; n_years x n_srv_comps x n_sexes
  
  // PARAMETER SECTION ----------------------------------------------
  
  // Biological parameters
  PARAMETER(ln_M); // log natural mortality
  PARAMETER_VECTOR(ln_RecPars); // Vector of recruitment parameters
  PARAMETER(ln_SigmaRec); // log sigma for recruitment
  PARAMETER_VECTOR(ln_N1_Devs); // log deviations for initial abundance
  PARAMETER_VECTOR(ln_RecDevs); // log recruitment deviations
  
  // Indices of abundance
  PARAMETER_VECTOR(logit_q_fish); // log catchability for fishery; n_fish_indices 
  PARAMETER_VECTOR(logit_q_srv); // log catchability for survey; n_srv_indices 
  
  // Fishery and selectivity parameters
  PARAMETER_MATRIX(ln_Fy); // Annual Fishing Mortality; n_years * n_fleets
  PARAMETER_ARRAY(ln_fish_selpars); // Fishery Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  PARAMETER_ARRAY(ln_srv_selpars); // Survey Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  
  // Random effects for fishery selectivity
  PARAMETER_ARRAY(ln_fish_selpars_re); // Fishery Selectivity Parameters, n_years * n_ages * n_fish_comps * n_sexes
  PARAMETER_ARRAY(fixed_sel_re_fish); // Correlations and sigma random effects
  // n_fixed_re_pars * n_fish_comps * n_sexes
  
  // Fishery Selectivity Random Effects Array Dimensions
  vector<int> fishsel_re_dim = ln_fish_selpars_re.dim; // Get dimensions of random effects array
  int n_re_years = fishsel_re_dim(0); // Rows in random effects array (Years)
  int n_re_pars = fishsel_re_dim(1); // Columns in random effects array (Ages or Parameters)
  
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
  array<Type> sum_FAA(n_years, n_ages, n_sexes); // Fishing mortality summed across fleets
  array<Type> CAA(n_years, n_ages, n_fleets, n_sexes); // Catch at Age
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
  matrix<Type> catch_nLL(n_years, n_fleets); 
  matrix<Type> fish_index_nLL(n_years,n_fish_indices); 
  matrix<Type> srv_index_nLL(n_years, n_srv_indices); 
  array<Type> fish_comp_nLL(n_years, n_fish_comps, n_sexes); 
  array<Type> srv_comp_nLL(n_years, n_srv_comps, n_sexes);
  Type rec_nLL = 0; // Recruitment likelihood penalty 
  Type fish_sel_re_nLL = 0; // Fishery selectivity random effects
  Type jnLL = 0; // Joint Negative log Likelihood
  
  // Set nLL components to zeros
  catch_nLL.setZero();
  fish_index_nLL.setZero();
  srv_index_nLL.setZero();
  fish_comp_nLL.setZero();
  srv_comp_nLL.setZero();

  // MODEL STRUCTURE ----------------------------------------------
  // y = year, a = age, s = sex, f = fishery fleet, sf = survey fleet
  
  // Selectivity ----------------------------------------------
  
  // Selectivity random effects  ------------------------------------------------------
  for(int f = 0; f < n_fish_comps; f++) {
    for(int s = 0; s < n_sexes; s++) {
      
      if(F_Slx_re_model(f, s) == 0) { // Random Walk Model
        
        for(int p = 0; p < n_re_pars; p++) {
          for(int y = 0; y < n_re_years; y++) {
            // penalize deviations
            fish_sel_re_nLL -= dnorm(ln_fish_selpars_re(y, p, f, s), Type(0.0), 
                                     fixed_sel_re_fish(p, f, s), true);
          } // y loop
        } // p loop
        
      } // end random walk if statement
      
      if(F_Slx_re_model(f, s) == 1) { // AR(1) by year 
        
        for(int p = 0; p < n_re_pars; p++) {
          // Create container to fill in with empty array
          array<Type> tmp_F_selpars_re(n_re_years);
          for(int y = 0; y < n_re_years; y++){
            tmp_F_selpars_re(y) = ln_fish_selpars_re(y, p, f, s);
          } // y loop
          
          // Sigma for Variance
          Type sigma_fish = exp( fixed_sel_re_fish(0, p, f, s) ); 
          // Correlation by year
          Type rho_y = Type(2)/(Type(1) + exp(-Type(2) * fixed_sel_re_fish(1, p, f, s) )) - Type(1); 
          // Variance of the AR process
          Type sigma_sel =  pow(pow(sigma_fish,2) / (1-pow(rho_y,2)),0.5); 
          // Evaluate likelihood here
          fish_sel_re_nLL += SCALE(AR1(rho_y), sigma_sel)(tmp_F_selpars_re); 
          
        } // p loop
      } // if AR(1) by year
      
    } // s loop
  } // f loop
  
  // Fishery Selectivity -------------------------------------------------------
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fish_comps; f++) {  
      
      // Index fishery blocks here
      int b = F_Slx_Blocks(y, f);
      
      for(int s = 0; s < n_sexes; s++) {
        
        // Define and extract selectivity parameters to feed into our 
        // Get_Selex function (our fixed selex parameters)
        vector<Type> tmp_ln_fish_selpars = ln_fish_selpars.transpose().col(f).col(s).col(b).vec();
        
        if(F_Slx_re_model(f, s) == 0) {  // random walk on parameters here
          
          // Temporary container object to store random walk objects deviations
          array<Type> tmp_re_devs_vec(n_re_years, n_re_pars);
          
          for(int p = 0; p < n_re_pars; p++) {
            if(y == 0) { // Year = 0
              tmp_re_devs_vec(y, p) = ln_fish_selpars_re(0, p, f, s);
            } else{
              tmp_re_devs_vec(y, p) = tmp_re_devs_vec(y - 1, p) + ln_fish_selpars_re(y, p, f, s);
            } // else = adding in random walk deviations
            
            // Add random walk deviates to parameters
            tmp_ln_fish_selpars(p) += tmp_re_devs_vec(y, p); // Exponentiated within Get_Selex
            
          } // end parameter (p) loop
        } // end if statement for random walk
        
        if(F_Slx_re_model(f, s) == 1) { // AR1 process by year
          
          for(int p = 0; p < n_re_pars; p++) {
            tmp_ln_fish_selpars(p) += ln_fish_selpars_re(y, p, f, s);
          } // p loop
        } // end if statement for AR1_y
        
        for(int a = 0; a < n_ages; a++) {
          // a + 1 because TMB indexes starting at 0
          F_Slx(y,a,f,s) = Get_Selex(a + 1, F_Slx_model(f), tmp_ln_fish_selpars);
          
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Survey Selectivity --------------------------------------------------------
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_srv_comps; f++) {
      
      // Index fishery blocks here
      int b = S_Slx_Blocks(y, f);
      
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Coerce selectivity parameters to vector
          vector<Type> tmp_ln_srv_selpars = ln_srv_selpars.transpose().col(f).col(s).col(b).vec();
          
          // a + 1 because TMB indexes starting at 0
          S_Slx(y,a,f,s) = Get_Selex(a + 1, S_Slx_model(f), tmp_ln_srv_selpars);
          
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Deaths and Removals (Fishing Mortality and Natural Mortality) --------------------------------
  
  // Get total fishing mortality
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fleets; f++) {
      
      // F in normal space
      F(y, f) = exp(ln_Fy(y,f)); 
      // Increment F by fleet to get total
      Total_Fy(y) += F(y, f);
      
      for(int a = 0; a < n_ages; a++) {
        for(int s = 0; s < n_sexes; s ++) {
          // Calculate F_at_age
          FAA(y, a, f, s) = F(y, f) * F_Slx(y, a, f, s);
          // Increment to add FAA to ZAA 
          sum_FAA(y, a, s) += FAA(y, a, f, s);
        } // s loop
      } // a loop
      
    } // f loop
  } // y loop
  
  // Finish calculating Total Mortality and Survival
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s ++) {
        // Sum FAA and M to get ZAA
        ZAA(y, a, s) = sum_FAA(y, a, s) + M;
        // Calculate survival fraction
        SAA(y, a, s) = exp(Type(-1.0) * ZAA(y, a, s));
      } // a loop
    } // s loop
  } // y loop
  
  // Initialization ----------------------------------------------
  Type SBPR0 = 0; // Spawning Biomass Per recruit container
  for(int a = 0; a < n_ages; a++) SBPR0 += exp(-M * Type(a)) * WAA(0, a, 0) * MatAA(0, a, 0); // Calculate SBPR0 here
  Type ssb0 = SBPR0 * exp(ln_RecPars(0)); // SSB0
  
  
  for(int s = 0; s < n_sexes; s++) {
    for(int a = 0; a < n_ages; a++){
      
      // Define initial recruitment parameter
      Type ln_RecInit = ln_RecPars(0); 
      
      // Fill in initial age-structure
      if(a < n_ages - 1) {
        NAA(0, a, s) = exp(ln_RecInit + ln_N1_Devs(a) -(0.5 * ln_SigmaRec2) -M * Type(a) ) * Sex_Ratio(s);
      } else{
        NAA(0, a, s) = exp(ln_RecInit -M * Type(a)) / (Type(1) - exp(-M)) * Sex_Ratio(s);
      }
      
      // Calculate SSB, Depletion and SBPR0 at first time point here
      if(s == 0) {
        SSB(0) += NAA(0, a, 0) * WAA(0, a, 0) * MatAA(0, a, 0);
        Depletion(0) = Type(1); 
      } // if for sex == 0 (females)
      
    } //  a loop
  } // s loop
  
  // Population Dynamics Equations ----------------------------------------------
  for(int y = 0; y < n_years; y++) {
    for(int s = 0; s < n_sexes; s++) {
      
      // Recruitment ----------------------------------------------
      if(y >= 1) { 
        
        if(rec_model == 0) { // Mean Recruitment
          Type ln_MeanRec = ln_RecPars(0); // Mean Recruitment parameter
          NAA(y, 0, s) = exp( ln_MeanRec + ln_RecDevs(y - 1) 
                                - Type(0.5) * ln_SigmaRec2) * Sex_Ratio(s);
        } // if for rec_model == 0
        
        if(rec_model == 1) { // Beverton Holt Recruitment
          
          // Define parameters 
          Type R0 = exp(ln_RecPars(0)); // Virgin Recruitment
          Type h = exp(ln_RecPars(1)); // Steepness

          // Get determininstic BH rec
          Type ln_det_BH_rec = log( (Type(4) * h * R0 * SSB(y - 1))  / 
                            ( ssb0*(Type(1)-h) + SSB(y - 1) * (Type(5)*h-Type(1)) ));
          
          // Get recruitment with process error here
          NAA(y, 0, s) =  exp( ln_det_BH_rec + ln_RecDevs(y-1) - Type(0.5) * ln_SigmaRec2 )  * Sex_Ratio(s);
          
        } // if Beverton Holt Recruitment
      } // only estimate recruitment if y >= 1
      
      // Increment total recruitment here
      Total_Rec(y) += NAA(y, 0, s);
      
      // Project Numbers At Age Forward ----------------------------------------
      
      for(int a = 0; a < n_ages; a++) {
        
        // Project ages and years forward
        if(a < (n_ages - 1)) { // Not in Plus Group
          NAA(y + 1, a + 1, s) = NAA(y, a, s) * SAA(y, a, s); 
        } 
        
        if(a == n_ages - 1){ // Increment previous year's plus group back in 
          NAA(y + 1, n_ages - 1, s) += (NAA(y, n_ages - 1, s) * SAA(y, n_ages - 1, s));
        } 
        
        // Calculate SSB
        if(y >= 1 && s == 0) {
          SSB(y) += NAA(y, a, 0) * WAA(y, a, 0) * MatAA(y, a, 0); 
          if(a == n_ages - 1) { // Get depletion
            Depletion(y) = SSB(y) / SSB(0); 
          } // if statement
        } // y >= 1
        
        // Increment Numbers at age to get total biomass
        Total_Biom(y) += NAA(y, a, s) * WAA(y, a, s);
        
      } // end ages loop
    } // end sex loop
  } // end year loop
  
  // Catch ----------------------------------------------
  pred_catches.setZero(); // set zero
  
  for(int f = 0; f < n_fleets; f++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Baranov's Catch Equation
          CAA(y, a, f, s) = NAA(y, a, s) * FAA(y, a, f, s) *  
            (Type(1.0) - SAA(y, a, s))  / ZAA(y, a, s);
          
          // Get Aggregated Catch - Increment catch in biomass
          pred_catches(y, f) += ( CAA(y, a, f, s) * WAA(y, a, s) );
          
        } // a loop
      } // s loop
      
    } // y loop
  } // f loop
  
  // Indices of Abundance ----------------------------------------------
  // (Assumes that indices are observed at the start of each year)
  pred_fish_indices.setZero(); // set zero
  pred_srv_indices.setZero(); // set zero
  
  // Fishery Index of Abundance
  for(int fi = 0; fi < n_fish_indices; fi++) {
    for(int y = 0; y < n_years; y++) { 
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Increment to get total index, conditional on selectivity
          pred_fish_indices(y, fi) += NAA(y, a, s) * WAA(y, a, s) *  F_Slx(y, a, fi, s);
        } // a loop
      } // s loop
      
      // Inverse logit transform
      Type tmp_q_fish = Type(0) + (Type(1) - Type(0))/(1 + exp(-logit_q_fish(fi))); 
      // Scale index by catchability here
      pred_fish_indices(y, fi) = tmp_q_fish  * pred_fish_indices(y, fi); 
      
    } // y loop
  } // fi loop
  
  // Survey Index of Abundance
  for(int si = 0; si < n_srv_indices; si++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Increment to get index conditional on selectivity
          pred_srv_indices(y, si) +=  NAA(y, a, s) * S_Slx(y, a, si, s);
        } // a loop
      } // s loop
      
      // Inverse logit transform here
      Type tmp_q_srv = Type(0) + (Type(1) - Type(0))/(1 + exp(-logit_q_srv(si))); 
      // Scale index by catchability
      pred_srv_indices(y, si) = tmp_q_srv * pred_srv_indices(y, si); 
      
    } // y loop
  } // si loop  
  
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
          
          // Divide by the total when done w/ predicting fish age comps (normalize to sum to 1)
          if(a == n_ages - 1) for(int a = 0; a < n_ages; a++) 
            pred_fish_age_comps(y, a, fc, s) /= Total_Fishery_Numbers(y, fc, s);
          
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
          
          // Divide by the total when done w/ predicting fish age comps (normalize to sum to 1)
          if(a == n_ages - 1) for(int a = 0; a < n_ages; a++) 
            pred_srv_age_comps(y, a, sc, s) /= Total_Survey_Numbers(y, sc, s);
          
        } // a loop
      } // s loop
    } // sc loop
  } // y loop
  
  
  // Likelihoods ----------------------------------------------
  
  // Catch likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> catch_sd(n_fleets); // Empty container to store sd calculation
  for(int f = 0; f < n_fleets; f++) catch_sd(f) = sqrt(log( (catch_cv(f)*catch_cv(f)) + 1));
  
  // Catch observed w/ minimal error
  for(int y = 0; y < n_years; y ++) {
    for(int f = 0; f < n_fleets; f++) {
      
      // Get likelihood here
      catch_nLL(y, f) -= use_catch(y, f) * dnorm(log(obs_catches(y, f)),
                log(pred_catches(y, f)) -(Type(0.5)* exp(2 * log(catch_sd(f))) ), 
                catch_sd(f), true);
      
      SIMULATE{ // Simulate catch
        obs_catches(y, f) = exp(rnorm(log(pred_catches(y, f)),  catch_sd(f) ));
      } // Simulation statement
      
    } // f loop
  } // y loop
  
  
  // Index likelihood (Log-normal likelihood) ----------------------------------------------
  
  vector<Type> fish_sd(n_fish_indices);   // Convert fishery index CV to standard deviation
  vector<Type> srv_sd(n_srv_indices);   // Convert survey index CV to standard deviation
  for(int si = 0; si < n_srv_indices; si ++) srv_sd(si) = sqrt(log( (srv_cv(si)*srv_cv(si)) + 1));
  for(int fi = 0; fi < n_fish_indices; fi++)  fish_sd(fi) = sqrt(log( (fish_cv(fi)*fish_cv(fi)) + 1));
  
  // Likelihood for fishery index
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      
      // Likelihood calculations
      fish_index_nLL(y, fi) -= use_fish_index(y, fi) * dnorm(log(obs_fish_indices(y, fi)), 
                     log(pred_fish_indices(y, fi))- (Type(0.5) * exp(2 * log(fish_sd(fi)) )),
                     fish_sd(fi), true);
      
      SIMULATE{ // Simulate Fishery Index
        obs_fish_indices(y, fi) = exp(rnorm(log(pred_fish_indices(y, fi)), fish_sd(fi))); 
      } // Simulation statement
      
    } // fi loop
  } // y loop
  
  // Likelihood for survey index
  for(int si = 0; si < n_srv_indices; si++) {
    for(int y = 0; y < n_years; y++) {
      
      // Likelihood calculations
      srv_index_nLL(y, si) -= use_srv_index(y, si) * dnorm(log(obs_srv_indices(y, si)), 
                    log(pred_srv_indices(y, si))- (Type(0.5) * exp(2*log(srv_sd(si)) )), 
                    srv_sd(si), true); 
      
      SIMULATE{ // Simulate Survey Index
        obs_srv_indices(y, si) = exp(rnorm(log(pred_srv_indices(y, si)), srv_sd(si))); 
      } // Simulation statement
      
    } // y loop
  } // si loop
  
  // Composition likelihoods (Multinomial likelihood) ----------------------------------------------
  
  Type c = 1e-10; // Constant to add to multinomial
  
  // Fishery Compositions
  vector<Type> obs_fish_age_vec(n_ages); // Obs fishery vector to hold and pass values to nLL
  vector<Type> pred_fish_age_vec(n_ages); // Pred fishery vector to hold and pass values to nLL
  
  for(int fc = 0; fc < n_fish_comps; fc++) {
    for(int y = 0 ; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) { 
        
        // Pull out observed age vector and multiply by the effective sample size
        obs_fish_age_vec = obs_fish_age_comps.col(s).col(fc).transpose().col(y) * obs_fish_age_Neff(y, fc, s);
        
        // Pull out predicted age vector
        pred_fish_age_vec = (pred_fish_age_comps.col(s).col(fc).transpose().col(y) + c);
        
        // Evaluate log-likelihood
        fish_comp_nLL(y, fc, s) -= use_fish_comps(y, fc, s) * dmultinom(obs_fish_age_vec.vec(), 
                      pred_fish_age_vec.vec(), true);
      } // s loop
    } // y loop
  } // fc loop
  
  vector<Type> obs_srv_age_vec(n_ages); // Obs survey vector to hold and pass values to nLL
  vector<Type> pred_srv_age_vec(n_ages); // Pred survey vector to hold and pass values to nLL
  
  for(int sc = 0; sc < n_srv_comps; sc++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Pull out observed age vector and multiply by the effective sample size
        obs_srv_age_vec = obs_srv_age_comps.col(s).col(sc).transpose().col(y) * obs_srv_age_Neff(y, sc, s);
        
        // Pull out predicted age vector
        pred_srv_age_vec = (pred_srv_age_comps.col(s).col(sc).transpose().col(y) + c);
        
        // Evaluate log-likelihood
        srv_comp_nLL(y, sc, s) -=  use_srv_comps(y, sc, s) * dmultinom(obs_srv_age_vec.vec(), 
                     pred_srv_age_vec.vec(), true);
        
      } // s loop
    } // y loop
  } // sc loop
  
  // Recruitment related stuff (likelihoods + derived quantites) ----------------------------------------------
  
  for(int y = 0; y < ln_N1_Devs.size(); y++) { // Mean = log-normal correction
    rec_nLL -= dnorm(ln_N1_Devs(y), Type(0), exp(ln_SigmaRec), true);
  } // Penalty for initial recruitment
  
  for(int y = 0; y < ln_RecDevs.size(); y++) { // Mean = log-normal correction
    rec_nLL -= dnorm(ln_RecDevs(y), Type(0), exp(ln_SigmaRec), true);
  } // Penalty for recruitment (mean should be 0)
  
  // Add to joint nLL   
  jnLL = rec_nLL + sum(srv_comp_nLL)+ sum(fish_comp_nLL) + 
    sum(srv_index_nLL) + sum(fish_index_nLL) + sum(catch_nLL) + fish_sel_re_nLL;
  
  // REPORT SECTION ----------------------------------------------
  REPORT(NAA); // Numbers at age
  REPORT(ZAA); // Total Mortality
  REPORT(FAA); // Fishing Mortality
  REPORT(CAA); // Catch at Age  
  REPORT(pred_catches); // Aggregate Catch by fleet
  REPORT(pred_srv_indices); // Survey Indices
  REPORT(pred_fish_indices); // Fishery Indices
  REPORT(F_Slx); // Fishery Selectivity
  REPORT(S_Slx); // Survey Selectivity
  REPORT(pred_fish_age_comps); // Predicted fishery age comps
  REPORT(pred_srv_age_comps); // Predicted survey age comps
  REPORT(ln_fish_selpars_re); // Selectivity random effects
  REPORT(SBPR0); // Spawning biomass per recruit
  REPORT(ssb0); // SSB0
  
  //  Likelihoods
  REPORT(catch_nLL);
  REPORT(srv_index_nLL);
  REPORT(fish_index_nLL);
  REPORT(fish_comp_nLL);
  REPORT(srv_comp_nLL);
  REPORT(rec_nLL);
  REPORT(fish_sel_re_nLL);
  REPORT(jnLL);
  
  // ADREPORT 
  ADREPORT(SSB);
  ADREPORT(Depletion); 
  ADREPORT(Total_Fy);
  ADREPORT(Total_Rec);
  ADREPORT(Total_Biom);
  
  return jnLL;
  
} // end objective function
