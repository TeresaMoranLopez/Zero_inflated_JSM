
functions { 
  //Kronecker product
    matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:cols(A)) { 
    for (j in 1:rows(A)) { 
    kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
    } 
    } 
    return kron; 
    } 
    } 

////////////////////////////////////////
data { 
    int<lower=1> N;                 // number of sites
    int<lower=1> J;                 // Number of species
    real<lower=0, upper = 1> y[N,J]; // matrix of species abundances
    int<lower=0, upper=1> y_zero[N,J] ; 
    int<lower=1> Kpres;                 //Covariates (+intercept) for presence
    int<lower=1> Kcov;                 //Covariates (+intercept) for presence
    int<lower=1> L_Jpres;               // Traits affecting response as presence
    int<lower=1> L_Jcov;               // Traits affecting response as cover
    matrix[N, Kpres] Xpres;       // Design matrix for covariates of presence
    matrix[N, Kcov] Xcov ;       //Design matrix for covariates of cover
    matrix[J, L_Jpres] TTpres;            // Species traits for presence
    matrix[J, L_Jcov] TTcov;            // Species traits for cover
    matrix[J, J] C;             // phylogenetic correlation matrix
    vector[J] ones;               // vector on 1s
    }
///////////////////////////////////////////
 parameters {
   //Presence of species
    corr_matrix[Kpres] Omega_p;     // correlation matrix for var-covar of betas (species environmental responses for presence)
    vector<lower=0>[Kpres] tau_p;   // scales for the variance covariance of betas
    vector[L_Jpres * Kpres] z_p; // coeffs for traits for presence and abundance
    vector[J * Kpres] betas_p; //coeffs for x for presence and abundance
    real<lower=0,upper=1> rho_p;  // correlation between phylogeny and betas
    real<lower=0> phi; // dispersion parameter beta
   //Species cover
    corr_matrix[Kcov] Omega_c;     // correlation matrix for var-covar of betas (now for cover)
    vector<lower=0>[Kcov] tau_c;   // scales for the variance covariance of betas
    vector[L_Jpres * Kcov] z_c; // coeffs for traits for presence and abundance
    vector[J * Kcov] betas_c; //coeffs for x for presence and abundance
    real<lower=0,upper=1> rho_c;  // correlation between phylogeny and betas
    }
///////////////////////////////////////////////
  transformed parameters { 
    ///Related to presence
    matrix[Kpres, Kpres] Sigma_p = quad_form_diag(Omega_p, tau_p);
    matrix[J*Kpres, J*Kpres] S_p = kronecker(Sigma_p, rho_p * C + (1-rho_p) *    diag_matrix(ones)); // 
    matrix[L_Jcov, Kpres] Z_p = to_matrix(z_p, L_Jpres, Kpres);    
    vector[J * Kpres] m_p = to_vector(TTpres * Z_p);        // coeffs  pres
    matrix[J, Kpres] b_p = to_matrix(betas_p, J, Kpres);  // 
    ////Now I estimate alpha (per species an site). Alpha = probability of presence
    matrix[N,J] alpha = inv_logit(Xpres * b_p');
    ///Related to cover
    matrix[Kcov, Kcov] Sigma_c = quad_form_diag(Omega_c, tau_c);
    matrix[J*Kcov, J*Kcov] S_c = kronecker(Sigma_c, rho_c * C + (1-rho_c) *    diag_matrix(ones)); // 
    matrix[L_Jcov, Kcov] Z_c = to_matrix(z_c, L_Jcov, Kcov);    
    vector[J * Kcov] m_c = to_vector(TTcov * Z_c);        // coeffs  pres
    matrix[J, Kcov] b_c = to_matrix(betas_c, J, Kcov);  
    //Now I estiamte mu (expected cover per species and site) 
    //And calculate the parameters of the beta distribution according to mu and phi
    matrix[N,J]  mu = inv_logit(Xcov*b_c');
    matrix[N,J] p = mu*phi;
    matrix[N,J] q = phi - (mu*phi);
  }
  
  model {
    // priors for presence
    Omega_p ~ lkj_corr(2);
    tau_p ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
    rho_p ~ beta(2,2);
    betas_p ~ multi_normal(m_p, S_p);
    z_p ~ normal(0,1);// you may want to use a wider prior with sd = 10
    //priors for cover
    // priors for presence
    Omega_c ~ lkj_corr(2);
    tau_c ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
    rho_p ~ beta(2,2);
    betas_c ~ multi_normal(m_c, S_c);
    z_c ~ normal(0,1); // you may want to use a wider prior with sd = 10
    phi ~ inv_gamma(.001, .001);//phi ~ cauchy(0, 5) it could also be something like this
  
   //likelihood
   for(n in 1:N)
   {
     for(j in 1:J)
     {
       y_zero[n,j] ~ bernoulli(1-alpha[n,j]);#Probability of absence
       if(y[n,j]  > 0)
       {
         y[n,j] ~beta(p[n,j], q[n,j]);#If present, cover
       }
     }
   }
  }

