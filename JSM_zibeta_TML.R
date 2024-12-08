# Data simulation ---------------------------------------------------------
rm(list = ls())
library(ape)
library(mvtnorm)
library(emulator)
library(matrixcalc)
set.seed(1234)

ns <-  50  # number of species
np <- 1   # number of covariates (without intercept)
nt <- 2   # number of traits (without intercept)
spp = paste("sp", 1:ns, sep = "")

# Create species ----------------------------------------------------------

# Sample species traits
dummy <- rbinom(ns, 1, 0.6)    # Dummy variable
quant <- rnorm(ns, 0, 1)  # Quantitative trait
#Create the trait matrix, remember that theu first column will be the intercept
TT <- as.matrix(cbind(rep(1, ns), dummy, quant))
rownames(TT) = spp
# Simulate phylogeny and estimate correlation
tree <- rtree(n = ns)    
C <- vcv(tree, corr = TRUE) # species correlation based on phylogeny
sum(diag(C)) == ns
isSymmetric(C)
colnames(C) = spp
rownames(C) = spp


# Simulate occupancy and cover --------------------------------------------

#Matrix Z [nt,np], traits per covariates, when it multiplies the trait matrix
#you will have the expected value of species responses to a certain covariate
#as a function of their traits.

Z_pres = array(NA, c((nt+1), (np+1)))#In both cases number of traits and parameters 
#we will sum a column of ones, the intercept

Z_pres[,1] = runif((nt+1), min = qlogis(0.05), max = qlogis(0.20))#Basal occupancies not too large
Z_pres[,2]= rnorm((nt+1), 1, 1)

# Expected species-specific responses to covariates related to presence
M_pres <- TT %*% Z_pres
#Variance-covariance matrix for presence

cor = rnorm(1, 0, 0.25)#Correlation among parameters only one because I have intercept-covariate1 only
omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= (np+1), nrow= (np+1))#Correlation matrix
taus = abs(runif((np+1), 0.1, 0.5))#variance of intercept and covariate
Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= (np+1), nrow= (np+1))
Sigma = round(Sigma, 4)#Avoid problems with floating
if(is.positive.definite(Sigma) == FALSE)#Needs to be positive definite because it is a variance-covar matrix
{
  while(is.positive.definite(Sigma) == FALSE)
  {
    cor = rnorm(1, 0, 0.25)#Correlation among parameters only one because I have intercept-covariate1 only
    omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= (np+1), nrow= (np+1))#Correlation matrix
    taus = abs(runif((np+1), 0.1, 0.5))#variance of intercept and covariate
    Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= (np+1), nrow= (np+1))
    Sigma = round(Sigma, 4)#Avoid problems with floating
  }
  
}

rho_pres = runif(1, 0.25, 0.75)  # correlation with phylogeny
thetas_pres <- rmvnorm(1, mean = as.vector(M_pres), kronecker(Sigma, rho_pres*C + (1-rho_pres) * diag(ns)))
Theta_pres <- matrix(thetas_pres[1,], ns, np+1)
hist(plogis(Theta_pres[,1]))#Check that you have variability intercept of presences 


#Define parameters for cover (proportion) and sample thetas (species environmental responses in terms of cover)
Z_cov = array(NA, c((nt+1), (np+1)))#

Z_cov[,1] = runif((nt+1), min = qlogis(0.05), max = qlogis(0.20))#Basal covers not too large
Z_cov[,2]= rnorm((nt+1), 1, 1)
M_cov <- TT %*% Z_cov

#Create sigma for cover
cor2 = rnorm(1, 0, 0.25)#Correlation among parameters only one because I have intercept-covariate1 only
omega2 = matrix(c(1, cor2, cor2, 1), byrow = T, ncol= (np+1), nrow= (np+1))#Correlation matrix
taus2 = abs(runif((np+1), 0.1, 0.5))#variance of intercept and covariate
Sigma2 = matrix(as.numeric(quad.form(omega2, diag(taus2))),  byrow = T, ncol= (np+1), nrow= (np+1))
Sigma2 = round(Sigma2, 4)#Avoid problems with floating
if(is.positive.definite(Sigma2) == FALSE)#Needs to be positive definite because it is a variance-covar matrix
{
  while(is.positive.definite(Sigma2) == FALSE)
  {
    cor2 = rnorm(1, 0, 0.25)#Correlation among parameters only one because I have intercept-covariate1 only
    omega2 = matrix(c(1, cor2, cor2, 1), byrow = T, ncol= (np+1), nrow= (np+1))#Correlation matrix
    taus2 = abs(runif((np+1), 0.1, 0.5))#variance of intercept and covariate
    Sigma2 = matrix(as.numeric(quad.form(omega2, diag(taus2))),  byrow = T, ncol= (np+1), nrow= (np+1))
    Sigma2 = round(Sigma2, 4)#Avoid problems with floating
  }
  
}

rho_cov = runif(1, 0.25, 0.75)  # correlation with phylogeny
thetas_cov <- rmvnorm(1, mean = as.vector(M_cov), kronecker(Sigma2, rho_cov*C + (1-rho_cov) * diag(ns)))
Theta_cov <- matrix(thetas_cov[1,], ns, np+1)
hist(plogis(Theta_cov[,1]))#Ok mean cover has variability (mean because I will scale the X variables
#and hence, when X = mean, it is zero and y = intercept (first column)



# Simulate communities based on species responses and environment-------------------------

nsu <- 75# sampling units
# sample level predictor (environmental covariate)
x = array(NA, c(nsu, (np+1)))
x[,1] = 1 #Intercept
x[,2] = scale(rnorm(nsu, mean = 0, sd = 1))
#If you want to practice with dummy variables (e.g., removal/control) You cound use instead
#x[,2]<- rbinom(nsu, 1, 0.5) #Dummy variabe (e.g., removal/Control)

#Presence
logit_p <- x %*% t(Theta_pres)#To make it faster I have vectorized this part
logit_p_v = c(logit_p)
#When we create a vector from a matrix it will join is by column
#In this case this means sp1site1, sp1site2, sp1siteN...spJsite1...spJsiteN
#Now we simulate presence based on these probabilities
z <- rbinom(nsu*ns, size = 1, prob = plogis(logit_p_v))#



#Cover
mu = plogis(x %*% t(Theta_cov))#exponencial porque el link de la poisson es un log
phi  = 10#Dispersion parameter of the beta distribution. You could also sample it for instance from an inverse gamma distribution
# Calculate parameters for beta distribution
p <- mu * phi
q <- phi - mu * phi
##Simulate responses. y = simulated cover*presence
cover = rbeta(nsu*ns, p, q)
y <- cover*z 

#Beta distribution is bounded in 0 and 1. In here to avoid an extra
#parameterization for values of 1 (which would be another bernoulli)
#I will substitute values of 1 by values of 0.99

tmp = which(y >=0.99)
if(length(tmp) >0)
{
  y[tmp] = 0.99
}

#Now I will transform y data into a matrix (input data for the model)

y_m = matrix(y, byrow = F, nrow = nsu, ncol = ns)

#Also I need to create a matrix that identifies zeroes

y_zero = array(NA, c(nsu, ns))
for(i in 1:nsu)
{
  y_zero[i,] = ifelse(y_m[i,] == 0, 1, 0)
}


#Now I save simulated data to fit the model

save(ns, nsu, C, TT, y, y_m, y_zero, x,
     Z_pres, Theta_pres, rho_pres,
     Z_cov, Theta_cov, rho_cov, file = "simulated_data.RData")



# Stan model --------------------------------------------------------------

cat(file = "Joint_zi_beta.stan", "
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

")

# Adjust model ------------------------------------------------------------
rm(list = ls())
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
load("simulated_data.RData")

stan_dat <- list(
  N = nrow(x),
  J = nrow(TT),
  y = y_m,
  y_zero = y_zero,
  Kpres = ncol(x),
  Kcov = ncol(x),
  L_Jpres = ncol(TT), #In our simulations we had the same traits for presence and cover, but they may differ
  L_Jcov = ncol(TT),
  Xpres = x, #In our simulations we had the same covariates for presence and cover, but htey can vary
  Xcov = x,
  TTpres = TT,
  TTcov = TT,
  C = C,
  ones = rep(1, nrow(TT))
 )

pars <- c("Omega_p", "tau_p", "rho_p", "betas_p", "z_p",
          "Omega_c", "tau_c", "rho_c", "betas_c", "z_c", "phi")

fit <- stan(file = 'Joint_zi_beta.stan', 
            data = stan_dat,
            pars = pars, 
            iter = 10000, #Usually you should put more, but this is just to see if we capture model parameters
            chains = 1)#You will also need more than 1 chain
#control = list(adapt_delta = 1.2)#Read https://mc-stan.org/rstanarm/reference/adapt_delta.html
#Usually with these JSM models we have divergent transitions, in my experience this is not problematic
#it may reflect some "weird" topographies in the likelihood. It happens when I capture model parameters and the
#model to fit exactly reflect the process of simulated data. Anyways, you can change the acceptance probability
#to try to fix it. 

save(fit, file = "model_fit.RData")

# Evaluate parameterization ----------------------------------------------
library(ggplot2)
fit_summary <- summary(fit)$summary
op <- par(mfrow = c(1,2))
hist(fit_summary[,10], main = "R-hat")
hist(fit_summary[,9], main = "n-eff" )
par(op)
View(fit_summary)#Do not worry about Omega being NaN. It will always be 1, because it is a correlation of a variable with itself
max(na.omit(fit_summary[,10]))#Rhat < 1.1, or 1.01
min(na.omit(fit_summary[,9]))#Minimum Neff

# Evaluate z (effect of traits on responses) ------------------------------


# On presence -------------------------------------------------------------

z_pres_f <- fit_summary[grepl("z_p", rownames(fit_summary)),]
#In vector z trait effects is organized in p1tr1, p1tr2... p1trt---
#pptr1...pptrt. Since we have
df <- data.frame(x = paste("parameter_", 1:nrow(Z_pres),sep = ""),
                 tz = c(Z_pres),
                 fz = z_pres_f[,1],
                 L = z_pres_f[,4],
                 U = z_pres_f[,8])

ggplot(df, aes(x = x, y = tz)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ylab("Trait effects (for presence)")+
  xlab ("")
#Here and in the following we are comparing true parameters (simulated, in red)
#vs model estimates (black points and lines).


# On cover ----------------------------------------------------------------

z_cov_f <- fit_summary[grepl("z_c", rownames(fit_summary)),]
#In vector z trait effects is organized in p1tr1, p1tr2... p1trt---
#pptr1...pptrt. In R and stan, when you vectorize a matrix you 
#do it linking column-vectors from the first to the last column
df <- data.frame(x = paste("parameter_", 1:nrow(Z_cov),sep = ""),
                 tz = c(Z_cov),
                 fz = z_cov_f[,1],
                 L = z_cov_f[,4],
                 U = z_cov_f[,8])

ggplot(df, aes(x = x, y = tz)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ylab("Trait effects (for cover)")+
  xlab ("")


# Species-specific responses ----------------------------------------------


# Presence ----------------------------------------------------------------

b_pres_f <- fit_summary[grepl("betas_p", rownames(fit_summary)),]
#Since we are vectorizing by columns it is b1sp1...b1spJ, b2sp1...b2spJ

#Intercept (first 50 rows) and first column of Theta_pres
df <- data.frame(x = paste("sp", 1:nrow(Theta_pres),sep = ""),
                 tz = Theta_pres[,1],
                 fz = b_pres_f[1:ns,1],
                 L = b_pres_f[1:ns,4],
                 U = b_pres_f[1:ns,8])
df$x = factor(df$x, levels = c(paste("sp", 1:nrow(Theta_pres),sep = "")))
ggplot(df, aes(x = x, y = tz)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ylab("Intercept (for presence)")+
  xlab ("")+
  theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust=1))
#If you want to see mean occupancies you may convert these values with plogis
#Like in the following
#You can also look at it with a scaterplot

ggplot(df, aes(x=plogis(tz), y=plogis(fz)) )+ geom_point(size = 3)+
  theme_classic()+
  ylab("Fitted mean occupancies")+
  xlab ("Simulated mean occupancies")+
  geom_abline(intercept =0 , slope = 1, color = "blue")
 #The geometric line is the perfect match 


#Estimates of response to the covariate, next 50 rows and second column
id_row = (ns+1):(ns*2)
df <- data.frame(x = paste("sp", 1:ns,sep = ""),
                 tz = Theta_pres[,2],
                 fz = b_pres_f[id_row,1],
                 L = b_pres_f[id_row,4],
                 U = b_pres_f[id_row,8])
df$x = factor(df$x, levels = c(paste("sp", 1:ns,sep = "")))
ggplot(df, aes(x = x, y = tz)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ylab("Environmental responses (for presence)")+
  xlab ("")+
  theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust=1))

ggplot(df, aes(x=tz, y=fz)) + geom_point(size = 3)+
  theme_classic()+
  ylab("Fitted betas")+
  xlab ("Simulated betas")+
  geom_abline(intercept =0 , slope = 1, color = "blue")



# Rhos --------------------------------------------------------------------

#You may also want to look at phylogeny effects

draws <- extract(fit)#In here I am extracting posterior distributions
plot(density(draws$rho_p), main = "Rho pres")
abline(v=rho_pres, col= "red")

draws <- extract(fit)#In here I am extracting posterior distributions
plot(density(draws$rho_c), main = "Rho cov")
abline(v=rho_pres, col= "red")

#You may want to look at other parameters, like those 
#Related to the variance-covariance of betas to construct Sigma together with phylogeny
#(i.e., variance covariance of the mutivariate normal from which we sample species responses)
#The phi parameter (dispersion parameter of beta)
#etc.
