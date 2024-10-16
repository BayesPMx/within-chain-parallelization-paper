// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, and Ka (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// Matrix-exponential solution using Torsten (the matrix-exponential seems to be
//   faster than the analytical solution for this model)
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0

functions{

  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y; 

  }
  
  vector depot_1cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke = cl/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];           // depot
    dydt[2] = ka*y[1] - ke*y[2];  // central
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  real<lower = 0> location_tvka;  // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
  
  int<lower = 1, upper = 3> solver; // 1 = analytical, 2 = matrix exponential, 3 = rk45
 
}
transformed data{ 

  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 3;                    // Number of random effects
  int n_cmt = 2;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_ka}; 
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] KE;
  
  vector[n_obs] ipred;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    vector[n_total] dv_ipred;
    matrix[n_total, 2] x_ipred;                      
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_ka = col(eta, 3);
    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    KE = CL ./ VC;
    
    for(j in 1:n_subjects){
    
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], KA[j]})';
        
      }else if(solver == 2){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -CL[j]/VC[j];
      
        x_ipred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
                         
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j]}, bioav, tlag)';
        
      }
                           
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }
    
    ipred = dv_ipred[i_obs];
    
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[TVCL/TVVC, ];

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    for(i in 1:n_obs){
      real sigma_tmp = ipred[i]*sigma_p;
      if(bloq_obs[i] == 1){
        target += log_diff_exp(normal_lcdf(lloq_obs[i] | ipred[i], sigma_tmp),
                               normal_lcdf(0.0 | ipred[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred[i], sigma_tmp); 
      }else{
        target += normal_lpdf(dv_obs[i] | ipred[i], sigma_tmp) -
                  normal_lccdf(0.0 | ipred[i], sigma_tmp);
      }
    }
  }
}
// generated quantities{
//   
//   real<lower = 0> sigma_sq_p = square(sigma_p);
// 
//   real<lower = 0> omega_cl = omega[1];
//   real<lower = 0> omega_vc = omega[2];
//   real<lower = 0> omega_ka = omega[3];
// 
//   real<lower = 0> omega_sq_cl = square(omega_cl);
//   real<lower = 0> omega_sq_vc = square(omega_vc);
//   real<lower = 0> omega_sq_ka = square(omega_ka);
// 
//   real cor_cl_vc;
//   real cor_cl_ka;
//   real cor_vc_ka;
//   real omega_cl_vc;
//   real omega_cl_ka;
//   real omega_vc_ka;
// 
//   vector[n_obs] pred;
//   vector[n_obs] dv_ppc;
//   vector[n_obs] log_lik;
//   vector[n_obs] res;
//   vector[n_obs] wres;
//   vector[n_obs] ires;
//   vector[n_obs] iwres;
//  
//   {
// 
//     matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
//     matrix[n_random, n_random] Omega = quad_form_diag(R, omega);
// 
//     vector[n_total] dv_pred;
//     matrix[n_total, 2] x_pred;
// 
//     cor_cl_vc = R[1, 2];
//     cor_cl_ka = R[1, 3];
//     cor_vc_ka = R[2, 3];
// 
//     omega_cl_vc = Omega[1, 2];
//     omega_cl_ka = Omega[1, 3];
//     omega_vc_ka = Omega[2, 3];
// 
//     for(j in 1:n_subjects){
//       
//       if(solver == 1){
//       
//         x_pred[subj_start[j]:subj_end[j],] =
//           pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
//                            amt[subj_start[j]:subj_end[j]],
//                            rate[subj_start[j]:subj_end[j]],
//                            ii[subj_start[j]:subj_end[j]],
//                            evid[subj_start[j]:subj_end[j]],
//                            cmt[subj_start[j]:subj_end[j]],
//                            addl[subj_start[j]:subj_end[j]],
//                            ss[subj_start[j]:subj_end[j]],
//                            {TVCL, TVVC, TVKA})';
//                          
//       }else if(solver == 2){
//         
//         matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
// 
//         K_tv[1, 1] = -TVKA;
//         K_tv[2, 1] = TVKA;
//         K_tv[2, 2] = -TVCL/TVVC;
// 
//         x_pred[subj_start[j]:subj_end[j],] =
//           pmx_solve_linode(time[subj_start[j]:subj_end[j]],
//                            amt[subj_start[j]:subj_end[j]],
//                            rate[subj_start[j]:subj_end[j]],
//                            ii[subj_start[j]:subj_end[j]],
//                            evid[subj_start[j]:subj_end[j]],
//                            cmt[subj_start[j]:subj_end[j]],
//                            addl[subj_start[j]:subj_end[j]],
//                            ss[subj_start[j]:subj_end[j]],
//                            K_tv, bioav, tlag)';
//         
//       }else{
//         
//         x_pred[subj_start[j]:subj_end[j],] =
//           pmx_solve_rk45(depot_1cmt_ode,
//                          n_cmt,
//                          time[subj_start[j]:subj_end[j]],
//                          amt[subj_start[j]:subj_end[j]],
//                          rate[subj_start[j]:subj_end[j]],
//                          ii[subj_start[j]:subj_end[j]],
//                          evid[subj_start[j]:subj_end[j]],
//                          cmt[subj_start[j]:subj_end[j]],
//                          addl[subj_start[j]:subj_end[j]],
//                          ss[subj_start[j]:subj_end[j]],
//                          {TVCL, TVVC, TVKA}, bioav, tlag)';
//         
//       }
// 
//       dv_pred[subj_start[j]:subj_end[j]] =
//         x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
// 
//     }
// 
//     pred = dv_pred[i_obs];
// 
//   }
// 
//   res = dv_obs - pred;
//   ires = dv_obs - ipred;
// 
//   for(i in 1:n_obs){
//     real ipred_tmp = ipred[i];
//     real sigma_tmp = ipred_tmp*sigma_p;
//     dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
//     if(bloq_obs[i] == 1){
//       log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
//                                 normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
//                    normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
//     }else{
//       log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) -
//                    normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
//     }
//     wres[i] = res[i]/sigma_tmp;
//     iwres[i] = ires[i]/sigma_tmp;
//   }
//   
// }
