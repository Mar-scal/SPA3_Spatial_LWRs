// ------------------------------------------------------------------------
// modeling of weight~height+age+... : spatial slope
// ------------------------------------------------------------------------

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() (){
    
    using namespace density;

    
    DATA_VECTOR(weight);
    DATA_MATRIX(varmat);
    DATA_VECTOR(depth);
    DATA_VECTOR(temperature);
    DATA_FACTOR(ind_loc); // length of obs
    DATA_FACTOR(ind_year); // length of obs


    // SPDE objects
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);

    
    PARAMETER_VECTOR(beta);
    PARAMETER_MATRIX(beta_s);
    PARAMETER_MATRIX(beta_t);
    PARAMETER_ARRAY(beta_st);
    PARAMETER_MATRIX(beta_depth);
    PARAMETER_MATRIX(beta_temperature);
    PARAMETER_VECTOR(log_kappa_s);
    PARAMETER_VECTOR(log_tau_s);
    PARAMETER_MATRIX(log_kappa_st);
    PARAMETER_MATRIX(log_tau_st);
    PARAMETER(log_phi);


    int n_loc = beta_s.rows(); 
    int n_year = beta_t.rows();
    int n_var = varmat.cols();
    int n_obs = varmat.rows();


    // ---------------------------
    // Joint negative log-likelihood
    Type jnll = Type(0.0);
    
    
    // Temporal effect
    // not implemented and left as a fixed effect


    // Spatial effect: Random intercept and Random slope
    
    for(int i_var = 0; i_var < n_var; i_var++){
        Eigen::SparseMatrix<Type> Q = exp(4*log_kappa_s(i_var))*M0 + Type(2.0)*exp(2*log_kappa_s(i_var))*M1 + M2;
        vector<Type> tmp(n_loc);
        for (int i_loc = 0; i_loc < n_loc; i_loc++){
            tmp(i_loc) = beta_s(i_loc, i_var);
        }
        jnll += SCALE(GMRF(Q), 1/exp(log_tau_s(i_var)))(tmp);
    }


    // Spatiotemporal effect: Random intercept and Random slope for each year
    
    for(int i_year = 0; i_year < n_year; i_year++){
        for(int i_var = 0; i_var < n_var; i_var++){
            Eigen::SparseMatrix<Type> Q = exp(4*log_kappa_st(i_var))*M0 + Type(2.0)*exp(2*log_kappa_st(i_var))*M1 + M2;
            vector<Type> tmp(n_loc);
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
                tmp(i_loc) = beta_st(i_loc, i_year, i_var);
            }
            jnll += SCALE(GMRF(Q), 1/exp(log_tau_st(i_var)))(tmp);
        }
    }



    
    // Observation likelihood
    int nobs = weight.size();
    vector<Type> mu(nobs); mu.setZero(); 
    vector<Type> mu_fixed(nobs); mu_fixed.setZero(); 
    vector<Type> mu_s(nobs); mu_s.setZero(); 
    vector<Type> mu_t(nobs); mu_t.setZero(); 
    vector<Type> mu_st(nobs); mu_st.setZero(); 
    vector<Type> mu_depth(nobs); mu_depth.setZero(); 
    vector<Type> mu_temperature(nobs); mu_temperature.setZero(); 
    for (int i = 0; i < nobs; i++){
        for(int j = 0; j < varmat.cols(); j++){
            mu_fixed(i) += varmat(i, j) * beta(j);
            mu_s(i) += varmat(i, j) * beta_s(ind_loc(i), j);
            mu_t(i) += varmat(i, j) * beta_t(ind_year(i), j);
            mu_st(i) += varmat(i, j) * beta_st(ind_loc(i), ind_year(i), j);
            mu_depth(i) += varmat(i, j) * beta_depth(ind_year(i), j) * depth(i);
            mu_temperature(i) += varmat(i, j) * beta_temperature(ind_year(i), j) * temperature(i);
        }
        mu(i) = mu_fixed(i) + mu_s(i) + mu_t(i) + mu_st(i) + mu_depth(i) + mu_temperature(i); 
        jnll -= dnorm(weight(i), exp(mu(i)), exp(log_phi), true);
        // jnll -= dgamma(weight(i), exp(mu(i))*exp(mu(i))/exp(log_phi), exp(log_phi)/exp(mu(i)), true);
    }
    
    
    
    // Report
    REPORT(beta);
    REPORT(beta_s);
    REPORT(beta_t);
    REPORT(beta_st);
    REPORT(beta_depth);
    REPORT(beta_temperature);
    REPORT(mu);
    REPORT(mu_fixed);
    REPORT(mu_s);
    REPORT(mu_t);
    REPORT(mu_st);
    REPORT(mu_depth);
    REPORT(mu_temperature);
    ADREPORT(mu);
    ADREPORT(beta);
    ADREPORT(beta_s);
    ADREPORT(beta_t);
    ADREPORT(beta_st);
    ADREPORT(beta_depth);
    ADREPORT(beta_temperature);

    
    return(jnll);
    
}

