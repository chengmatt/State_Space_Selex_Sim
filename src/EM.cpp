// General single species age-and sex-structured stock assessment
// that accommodates multiple fishery fleets written in TMB
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date updated: 3/11/23

#include<TMB.hpp>
#include "Utility_Fxns.hpp"

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
  DATA_ARRAY(obs_fish_age_comps); // Array of fishery age comps; n_years * n_ages * n_fleets * n_sexes (in numbers)
  DATA_ARRAY(obs_srv_age_comps); // Array of survey age comps; n_years * n_ages * n_srv_indices * n_sexes (in numbers)
  DATA_ARRAY(obs_fish_age_Input_N); // Array of fishery age comps; n_years * n_fleets * n_sexes 
  DATA_ARRAY(obs_srv_age_Input_N); // Array of survey age comps; n_years * n_srv_indices * n_sexes 
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
  DATA_IMATRIX(F_Slx_model); // Fishery Selectivity Model, dim = year, fleet
  DATA_IMATRIX(F_Slx_Blocks); // Fishery Selectivity Time Blocks, n_years * n_fish_comps; 
  // this is set up such that the selectivity within a fleet and across sexes is constant
  DATA_IVECTOR(S_Slx_model); // Survey Selectivity Model, == 0 Logistic n_fleets 
  DATA_IMATRIX(S_Slx_Blocks); // Survey Selectivity Time Blocks, n_years * n_srv_fleets; 
  DATA_IMATRIX(fish_comp_likelihoods);  // Vector of integers denoting likelihoods for fish comps (n_fish_fleets, n_sexes)
  DATA_IMATRIX(srv_comp_likelihoods); // Vector of integers denoting likelihoods for survey comps (n_srv_fleets, n_sexes)
  
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
  PARAMETER_VECTOR(ln_N1Devs); // log recruitment deviations for inital age structure
  PARAMETER_VECTOR(ln_RecDevs); // log recruitment deviations
  
  // Indices of abundance
  PARAMETER_VECTOR(ln_q_fish); // log catchability for fishery; n_fish_indices 
  PARAMETER_VECTOR(ln_q_srv); // log catchability for survey; n_srv_indices 
  
  // Fishery and selectivity parameters
  PARAMETER_MATRIX(ln_Fy); // Annual Fishing Mortality; n_years * n_fleets
  PARAMETER_ARRAY(ln_fish_selpars); // Fishery Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  PARAMETER_ARRAY(ln_srv_selpars); // Survey Selectivity Parameters, n_comps * n_sexes * n_blocks * n_pars
  
  // Compositional likelihood parameters
  PARAMETER_MATRIX(ln_DM_Fish_Param); // Vector of theta parameters for Dirichlet-Multinomial for fishery (n_fleets, n_sexes)
  PARAMETER_MATRIX(ln_DM_Srv_Param); // Vector of theta parameters for Dirichlet-Multinomial for survey (n_srv_fleets, n_sexes)
  
  // Random effects for fishery selectivity
  PARAMETER_ARRAY(ln_fish_selpars_re); // Fishery Selectivity Parameters, n_years * n_ages * n_fish_comps * n_sexes
  PARAMETER_ARRAY(ln_fixed_sel_re_fish); // Correlations and sigma random effects
  // n_fixed_re_pars * n_fish_comps * n_sexes
  
  // Fishery Selectivity Random Effects Array Dimensions
  vector<int> fishsel_re_dim = ln_fish_selpars_re.dim; // Get dimensions of random effects array
  int n_re_years = fishsel_re_dim(0); // Rows in random effects array (Years)
  int n_re_pars = fishsel_re_dim(1); // Columns in random effects array (Ages or Parameters)
  
  // Parameter Transformations ----------------------------------------------
  Type M = exp(ln_M); // Natural Mortality
  Type SigmaRec = exp(ln_SigmaRec); // Recruitment Sigma in normal space
  Type SigmaRec2 = pow(SigmaRec, 2); // Variance of recruitment sigma
  
  // Predicted Quantities
  array<Type> pred_fish_age_comps(obs_fish_age_comps.dim); // Predicted Fishery Age Comps
  array<Type> pred_srv_age_comps(obs_srv_age_comps.dim); // Predicted survey age comps
  matrix<Type> pred_fish_indices(n_years, n_fish_indices); // Predicted fishery indices
  matrix<Type> pred_srv_indices(n_years, n_srv_indices); // Predicted survey indices
  matrix<Type> pred_catches(n_years, n_fleets); // Predicted fishery catches
  
  // Stored Quantities
  array<Type> NAA(n_years + 1, n_ages, n_sexes); // Numbers at age; n_years + 1 
  array<Type> ZAA(n_years, n_ages, n_sexes); // Total Mortality
  array<Type> SAA(n_years, n_ages, n_sexes); // Survival at age
  array<Type> FAA(n_years, n_ages, n_fleets, n_sexes); // Fishing Mortality
  array<Type> sum_FAA(n_years, n_ages, n_sexes); // Fishing mortality summed across fleets
  matrix<Type> Fy(n_years, n_fleets); // Fishing Mortality by year and fleet
  matrix<Type> Exploit_Biom(n_years, n_fleets); // Exploitable Biomass for a given fleet
  matrix<Type> Harvest_Rate(n_years, n_fleets); // Harvest rate calculations for a given fleet
  array<Type> CAA(n_years, n_ages, n_fleets, n_sexes); // Catch at Age
  vector<Type> Total_Fy(n_years); // Total F summed across fleets
  vector<Type> Total_Harvest_Rate(n_years); // Total harvest rate
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
  array<Type> Fish_Neff(n_years, n_fish_comps, n_sexes); 
  array<Type> Srv_Neff(n_years, n_srv_comps, n_sexes);
  Type rec_nLL = 0; // Recruitment likelihood penalty 
  Type Fpen_nLL = 0; // Penalty for Fs
  Type fish_sel_re_nLL = 0; // Fishery selectivity random effects
  Type jnLL = 0; // Joint Negative log Likelihood
  
  // Set nLL components to zeros
  catch_nLL.setZero();
  fish_index_nLL.setZero();
  srv_index_nLL.setZero();
  fish_comp_nLL.setZero();
  srv_comp_nLL.setZero();
  
  // MODEL STRUCTURE ----------------------------------------------
  // y = year, a = age, s = sex; in general, f = fishery fleet, sf = survey fleet
  
  // Selectivity ----------------------------------------------
  // Selectivity random effects  ------------------------------------------------------
  for(int f = 0; f < n_fish_comps; f++) {
    for(int s = 0; s < n_sexes; s++) {
      
      if(F_Slx_re_model(f, s) == 0) { // Random Walk Model
        
        for(int p = 0; p < n_re_pars; p++) {
          for(int y = 0; y < n_re_years; y++) {
            // penalize deviations
            fish_sel_re_nLL -= dnorm(ln_fish_selpars_re(y, p, f, s), 
                                     Type(0.0), exp(ln_fixed_sel_re_fish(p, f, s)), true);
          } // y loop
        } // p loop
      } // end random walk if statement
      
      if(F_Slx_re_model(f, s) == 1) { // Semi-parametric deviations (year x age)
        
        for(int a = 0; a < n_ages; a++) {
          for(int y = 0; y < n_re_years; y++) {
            // penalize deviations
            fish_sel_re_nLL -= dnorm(ln_fish_selpars_re(y, a, f, s), 
                                     Type(0.0), exp(ln_fixed_sel_re_fish(f, s)), true);
          } // y loop
        } // p loop
      } // end random walk if statement
      
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
        vector<Type> tmp_ln_fish_selpars = ln_fish_selpars.transpose().col(f).col(s).col(b);
        
        if(F_Slx_re_model(f, s) == 0) {  // random walk on parameters here
          
          // Temporary container object to store selectivity deviations
          array<Type> tmp_seldevs_vec(n_years, n_re_pars);
          
          for(int p = 0; p < n_re_pars; p++) {
            if(y == 0) { // Year = 0 (Initial conditions) (Initial Condition)
              tmp_seldevs_vec(y, p) = tmp_ln_fish_selpars(p); // Fixed effect for first year 
            } else{ // Random walk deviations after first year
              tmp_seldevs_vec(y, p) = tmp_seldevs_vec(y - 1, p) + ln_fish_selpars_re(y - 1, p, f, s);
            } // else = adding in random walk deviations
            
            // Redefine tmp_ln_fish_pars with updated RW deviates
            tmp_ln_fish_selpars(p) = tmp_seldevs_vec(y, p); // Exponentiated within Get_Selex
            
          } // end parameter (p) loop
        } // end if statement for random walk
    
        for(int a = 0; a < n_ages; a++) {
          
          // Get age index
          Type fish_age_idx = ages(a);
          
          // Get selex value
          if(F_Slx_re_model(f, s) == 1) { // if semi-parametric
            if(y == 0) {
              F_Slx(y,a,f,s) = Get_Selex(fish_age_idx, F_Slx_model(y,f), tmp_ln_fish_selpars); // estimate parametric form
            } // end first year
            if(y > 0) {
              F_Slx(y,a,f,s) = F_Slx(y-1 ,a,f,s) * exp(ln_fish_selpars_re(y-1, a, f, s));
            }
          } else{ // not semi-parametric (either random walk or time-invariant)
            F_Slx(y,a,f,s) = Get_Selex(fish_age_idx, F_Slx_model(y,f), tmp_ln_fish_selpars);
          } // end if else 
             
        } // a loop
        if(F_Slx_re_model(f, s) == 1) { // semi-parametric curvature penalty
          // Extract out fishery selectivity
          if(y > 0) {
            vector<Type> selex_vec = F_Slx.col(s).col(f).transpose().col(y);
            // Set up curvature pentalty - squared second differences
            // Get first difference
            int n = selex_vec.size();
            vector<Type> ans(n - 1);
            for (int i = 0; i < n - 1; i++) ans(i) = selex_vec(i+1) - selex_vec(i);
            // Get second difference
            int n1 = ans.size();
            vector<Type> ans1(n1 - 1);
            for (int i1 = 0; i1 < n1 - 1; i1++) ans1(i1) = ans(i1+1) - ans(i1);
            fish_sel_re_nLL += pow(sum(ans1),2); // add into likelihood penalty
          } // end if y > 0
        } // end curvature penalty
      } // s loop
    } // f loop
  } // y loop
  
  // Survey Selectivity --------------------------------------------------------
  for(int y = 0; y < n_years; y++) {
    for(int sf = 0; sf < n_srv_comps; sf++) {
      
      // Index fishery blocks here
      int b = S_Slx_Blocks(y, sf);
      
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Get age index
          Type srv_age_idx = ages(a);
          
          // Coerce selectivity parameters to vector
          vector<Type> tmp_ln_srv_selpars = ln_srv_selpars.transpose().col(sf).col(s).col(b);
          
          // Get selex value here
          S_Slx(y,a,sf,s) = Get_Selex(srv_age_idx, S_Slx_model(sf), tmp_ln_srv_selpars);
          
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Removals (Fishing Mortality) --------------------------------
  // Get FAA and fishing mortality
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fleets; f++) {
      for(int a = 0; a < n_ages; a++) {
        for(int s = 0; s < n_sexes; s++) {
          // Calculate F_at_age
          FAA(y, a, f, s) =  exp(ln_Fy(y, f)) * F_Slx(y, a, f, s);
          // Increment to add FAA to ZAA 
          sum_FAA(y, a, s) += FAA(y, a, f, s);
        } // s loop
      } // a loop
      // Increment F
      Total_Fy(y) +=  exp(ln_Fy(y, f));
    } // f loop 
  } // y loop
  
  // Get Removals and Deaths here (Fishing and Natural Mortality)
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s++) {
        ZAA(y, a, s) = sum_FAA(y, a, s) + M; // Sum FAA and M to get ZAA
        SAA(y, a, s) = exp(Type(-1.0) * ZAA(y, a, s)); // Calculate survival fraction
      } // s loop
    } // a loop
  } // y loop
  
  // Initialization and SBPR calculations ----------------------------------------------
  vector<Type> SBPR_N(n_ages); // Numbers Per recruit container
  vector<Type> SBPR_SSB0(n_ages); // Spawning Biomass Per recruit container
   
  // Loop through spawning biomass per recruit calculations
  for(int a = 0; a < n_ages; a++) {
    if(a == 0) SBPR_N(a) = Type(1);
    if(a > 0 && a < n_ages - 1) SBPR_N(a) = SBPR_N(a - 1) * exp(-M); 
    if(a == n_ages - 1) SBPR_N(a) = ((SBPR_N(a - 1) * exp(-M)) / (1 - exp(-M)));
    SBPR_SSB0(a) = SBPR_N(a) * MatAA(0, a, 0) * WAA(0, a, 0); 
  } // a loop
  
  // Get B0 here
  Type ssb0 = SBPR_SSB0.sum() * exp(ln_RecPars(0)); // SSB0
  
  // Initialize age-structure 
  for(int a = 0; a < n_ages; a++) {
    for(int s = 0; s < n_sexes; s++){
      // TESTING
      // NAA(0, a, s) = exp(ln_N1Devs(a)) * Sex_Ratio(s);
      if(a < n_ages - 1) { // not plus group
        NAA(0, a, s) = exp(ln_N1Devs(a)) * exp(ln_RecPars(0)) * exp(-M * Type(a)) * Sex_Ratio(s);
      } else{
        NAA(0, a, s) = exp(ln_RecPars(0)) * exp(-M * Type(a)) / (1 - exp(-M)) * Sex_Ratio(s);
      } // plus group
    } // end s loop 
  } // end a loop
   
   // Population Dynamics Equations ----------------------------------------------
   for(int y = 0; y < n_years; y++) {
     for(int s = 0; s < n_sexes; s++) {
       
       // Recruitment ----------------------------------------------
       if(y >= 1) { 
         
         if(rec_model == 0) { // Mean Recruitment
           Type ln_MeanRec = ln_RecPars(0); // Mean Recruitment parameter
           NAA(y, 0, s) = exp( ln_MeanRec + ln_RecDevs(y-1) )  * Sex_Ratio(s);
         } // if for rec_model == 0
          
         if(rec_model == 1) { // Beverton Holt Recruitment
           
           // Define parameters 
           Type R0 = exp(ln_RecPars(0)); // Virgin Recruitment
           Type h = exp(ln_RecPars(1)); // Steepness in normal space
           
           // Define BH components
           Type BH_first_part = Type(4) * h * R0 * SSB(y - 1);
           Type BH_sec_part = (ssb0 * (Type(1) - h)) + SSB(y - 1) * ((Type(5)*h) - 1);
           Type ln_BH_tmp_rec = log(BH_first_part/BH_sec_part); // Deterministic BH in log space
           
           // Get recruitment with process error here
           NAA(y, 0, s) =  exp(ln_BH_tmp_rec + ln_RecDevs(y - 1))  * Sex_Ratio(s);
           
         } // if BH Recruitment
       } // only estimate recruitment if y >= 1
       
       // Project Numbers At Age Forward ----------------------------------------
       for(int a = 0; a < n_ages; a++) {
         
         // Project ages and years forward
         if(a < (n_ages - 1)) { // Not in Plus Group
           NAA(y + 1, a + 1, s) = NAA(y, a, s) * SAA(y, a, s); 
         } // if a < n_ages - 1
         
         if(a == n_ages - 1){ // Increment previous year's plus group back in 
           NAA(y + 1, n_ages - 1, s) += (NAA(y, n_ages - 1, s) * SAA(y, n_ages - 1, s));
         } // if a == n_ages - 1
         
         // Increment Numbers at age to get total biomass
         Total_Biom(y) += NAA(y, a, s) * WAA(y, a, s);
         
       } // end ages loop
       // Increment total recruitment here
       Total_Rec(y) += NAA(y, 0, s);
     } // end sex loop
     
     // Calculate SSB here (Indexing 0 for females)
     SSB(y) = sum(NAA.col(0).transpose().col(y) * 
       MatAA.col(0).transpose().col(y) *
       WAA.col(0).transpose().col(y));
     
   } // end year loop
  
  // Catch and Related Quantities ----------------------------------------------
  pred_catches.setZero(); // set zero
  for(int f = 0; f < n_fleets; f++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Baranov's Catch Equation
          CAA(y, a, f, s) = NAA(y, a, s) * (Type(1.0) - exp(-ZAA(y, a, s))) *
                           (FAA(y, a, f, s) / ZAA(y, a, s));
          
          // Compute Exploitable Biomass; N * Selex * WAA
          Exploit_Biom(y, f) += NAA(y, a, s) * F_Slx(y, a, f, s) * WAA(y, a, s); 
            
          // Get Aggregated Catch - Increment catch in biomass
          pred_catches(y, f) += ( CAA(y, a, f, s) * WAA(y, a, s) );
          
        } // a loop
      } // s loop
      
      // Get Harvest Rates calculations here
      Harvest_Rate(y, f) = pred_catches(y, f) / Exploit_Biom(y, f);
      // Get Total Harvest rate here
      Total_Harvest_Rate(y) += Harvest_Rate(y, f);
      
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
          pred_fish_indices(y, fi) += exp(ln_q_fish(fi)) * NAA(y, a, s) * 
                                      WAA(y, a, s) *  F_Slx(y, a, fi, s);
        } // a loop
      } // s loop
    } // y loop
  } // fi loop 
  
  // Survey Index of Abundance
  for(int si = 0; si < n_srv_indices; si++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) { 
        for(int a = 0; a < n_ages; a++) {
          // Increment to get scale index conditional on selectivity
          pred_srv_indices(y, si) += exp(ln_q_srv(si)) * NAA(y, a, s) *
                                     S_Slx(y, a, si, s);
        } // a loop
      } // s loop
    } // y loop 
  } // si loop  
  
  // Calculate Compositions and Related Quantities ----------------------------------------------  
  
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
          if(a == n_ages - 1) {
            for(int a = 0; a < n_ages; a++) {
              pred_fish_age_comps(y, a, fc, s) /= Total_Fishery_Numbers(y, fc, s);
            } // a loop
          } // if a = plus group

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
          if(a == n_ages - 1) {
            for(int a = 0; a < n_ages; a++) {
              pred_srv_age_comps(y, a, sc, s) /= Total_Survey_Numbers(y, sc, s);
            } // a loop
          } // if a = plus group
          
        } // a loop
      } // s loop
    } // sc loop
  } // y loop
  
  // Catch likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> catch_sd(n_fleets); // Empty container to store sd calculation
  for(int f = 0; f < n_fleets; f++) catch_sd(f) = sqrt(log( (pow(catch_cv(f),2) + 1)));
  
  // Catch observed w/ minimal error
  for(int y = 0; y < n_years; y ++) {
    for(int f = 0; f < n_fleets; f++) {
      
      // Get likelihood here
      catch_nLL(y, f) -= use_catch(y, f) * dnorm(log(obs_catches(y, f)),
                         log(pred_catches(y, f)) - pow(catch_sd(f), 2)/2, catch_sd(f), true);
      
      SIMULATE{ // Simulate catch
        obs_catches(y, f) = exp(rnorm(log(pred_catches(y, f)) - 
                            pow(catch_sd(f), 2)/2, catch_sd(f) ));
      } // Simulation statement
      
    } // f loop
  } // y loop
  
   
  // Index likelihood (Log-normal likelihood) ----------------------------------------------
  vector<Type> fish_sd(n_fish_indices); // Convert fishery index CV to standard deviation
  vector<Type> srv_sd(n_srv_indices); // Convert survey index CV to standard deviation
  for(int si = 0; si < n_srv_indices; si ++) srv_sd(si) = sqrt(log( pow(srv_cv(si), 2) + 1));
  for(int fi = 0; fi < n_fish_indices; fi++) fish_sd(fi) = sqrt(log( pow(fish_cv(fi), 2) + 1));
  
  // Likelihood for fishery index
  for(int y = 0; y < n_years; y++) {
    for(int fi = 0; fi < n_fish_indices; fi++) {
      
      // Likelihood calculations
      fish_index_nLL(y, fi) -= use_fish_index(y, fi) * dnorm(log(obs_fish_indices(y, fi)), 
                               log(pred_fish_indices(y, fi)) - pow(fish_sd(fi),2)/2, 
                               fish_sd(fi), true);
      
      SIMULATE{ // Simulate Fishery Index
        obs_fish_indices(y, fi) = exp(rnorm(log(pred_fish_indices(y, fi)) 
                                    - pow(fish_sd(fi),2)/2, fish_sd(fi))); 
      } // Simulation statement
      
    } // fi loop
  } // y loop
  
  // Likelihood for survey index
  for(int si = 0; si < n_srv_indices; si++) {
    for(int y = 0; y < n_years; y++) {
      
      // Likelihood calculations
      srv_index_nLL(y, si) -= use_srv_index(y, si) * dnorm(log(obs_srv_indices(y, si)), 
                              log(pred_srv_indices(y, si)) - pow(srv_sd(si), 2)/2, 
                              srv_sd(si), true); 
                
      SIMULATE{ // Simulate Survey Index
        obs_srv_indices(y, si) = exp(rnorm(log(pred_srv_indices(y, si)) 
                                    -pow(srv_sd(si), 2)/2, srv_sd(si))); 
      } // Simulation statement
      
    } // y loop
  } // si loop
  
  // Composition likelihoods  ----------------------------------------------
  // Fishery Compositions
  vector<Type> obs_fish_age_vec(n_ages); // Obs fishery vector to hold and pass values to nLL
  vector<Type> pred_fish_age_vec(n_ages); // Pred fishery vector to hold and pass values to nLL
  matrix<Type> DM1_FI(n_years, n_ages);
  matrix<Type> DM2_FI(n_years, n_ages);
  
  for(int fc = 0; fc < n_fish_comps; fc++) {
    for(int y = 0 ; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Pre-processing - extract out quantities
        vector<Type> obs_fish_age = obs_fish_age_comps.col(s).col(fc).transpose().col(y); // Pull out observed vector
        vector<Type> pred_fish_age = pred_fish_age_comps.col(s).col(fc).transpose().col(y); // Pull out predicted vector 
        Type Fish_Input_N = obs_fish_age_Input_N(y, fc, s); // Input Sample Size
        Type ln_fish_theta = ln_DM_Fish_Param(fc, s); // Dispersion parameter if needed
        
        // Get likelihood here
        fish_comp_nLL(y, fc, s) -= use_fish_comps(y, fc, s) * 
                                   get_acomp_nLL(obs_fish_age, pred_fish_age, Fish_Input_N,
                                                 ln_fish_theta, fish_comp_likelihoods(fc, s), true);
         
        // Get effective sample sizes
        Fish_Neff(y, fc, s) = get_Neff(Fish_Input_N, ln_fish_theta, fish_comp_likelihoods(fc, s));
 
      } // s loop
    } // y loop
  } // fc loop
  
  vector<Type> obs_srv_age_vec(n_ages); // Obs survey vector to hold and pass values to nLL
  vector<Type> pred_srv_age_vec(n_ages); // Pred survey vector to hold and pass values to nLL
  
  for(int sc = 0; sc < n_srv_comps; sc++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Pre-processing - extract out quantities
        vector<Type> obs_srv_age= obs_srv_age_comps.col(s).col(sc).transpose().col(y); // Pull out observed vector
        vector<Type> pred_srv_age = pred_srv_age_comps.col(s).col(sc).transpose().col(y); // Pull out predicted vector
        Type Srv_Input_N = obs_srv_age_Input_N(y, sc, s); // Get Input sample size
        Type ln_srv_theta = ln_DM_Srv_Param(sc, s); // Dispersion parameter if needed
        
        // Get likelihood here
        srv_comp_nLL(y, sc, s) -= use_srv_comps(y, sc, s) * 
                                  get_acomp_nLL(obs_srv_age, pred_srv_age, Srv_Input_N,
                                                ln_srv_theta, srv_comp_likelihoods(sc, s), true);
        
        // Get effective sample sizes
        Srv_Neff(y, sc, s) = get_Neff(Srv_Input_N, ln_srv_theta, srv_comp_likelihoods(sc, s));
        
      } // s loop
    } // y loop 
  } // sc loop
  
  
  // Recruitment and derived quantities ---------------------------------------------
  for(int y = 0; y < ln_RecDevs.size(); y++) { 
    rec_nLL -= dnorm(ln_RecDevs(y), -(SigmaRec2/Type(2)), SigmaRec, true);
  } // Penalty for all recruitment deviations
  
  for(int y = 0; y < ln_N1Devs.size(); y++) {
    rec_nLL -= dnorm(ln_N1Devs(y), -(SigmaRec2/Type(2)), SigmaRec, true);
  } // Penalty for all recruitment deviations
  
  // Calculate Depletion at time 0
  Depletion = SSB/ssb0;    
  
  // Add everything back to joint nLL   
  jnLL = rec_nLL + srv_comp_nLL.sum()+ fish_comp_nLL.sum() + 
    srv_index_nLL.sum() + fish_index_nLL.sum() + catch_nLL.sum() + 
    fish_sel_re_nLL;
  
  // REPORT SECTION ----------------------------------------------
  REPORT(NAA); // Numbers at age
  REPORT(ZAA); // Total Mortality
  REPORT(FAA); // Fishing Mortality
  REPORT(CAA); // Catch at Age  
  REPORT(SAA); // Survival at Age
  REPORT(pred_catches); // Aggregate Catch by fleet
  REPORT(pred_srv_indices); // Survey Indices
  REPORT(pred_fish_indices); // Fishery Indices
  REPORT(Exploit_Biom); // Exploitable Biomass
  REPORT(Harvest_Rate); // Harvest Rates
  REPORT(F_Slx); // Fishery Selectivity
  REPORT(S_Slx); // Survey Selectivity
  REPORT(pred_fish_age_comps); // Predicted fishery age comps
  REPORT(pred_srv_age_comps); // Predicted survey age comps
  REPORT(ln_fish_selpars_re); // Selectivity random effects
  REPORT(SBPR_SSB0); // Spawning biomass per recruit
  
  //  Likelihoods
  REPORT(catch_nLL);
  REPORT(srv_index_nLL);
  REPORT(fish_index_nLL);
  REPORT(fish_comp_nLL);
  REPORT(srv_comp_nLL);
  REPORT(rec_nLL);
  REPORT(Fpen_nLL);
  REPORT(fish_sel_re_nLL);
  REPORT(jnLL);
  
  // ADREPORT 
  ADREPORT(SSB);
  ADREPORT(Depletion); 
  ADREPORT(Total_Fy);
  ADREPORT(Total_Rec);
  ADREPORT(Total_Harvest_Rate);
  ADREPORT(Exploit_Biom); 
  ADREPORT(Total_Biom);
  ADREPORT(ssb0);
  ADREPORT(Fish_Neff);
  ADREPORT(Srv_Neff);
  
  return jnLL;
  
} // end objective function
