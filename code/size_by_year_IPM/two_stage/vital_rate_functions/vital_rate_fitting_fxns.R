##  This script contains a set of functions for fitting growth, survival, and
##  and recruiment models to demographic data. There are separate models for
##  juveniles (one-year-old plants) and adults (greater than one year old).
##  There are also two distinct model types: one that contains year-by-size 
##  interactions, and one that does not. All models include random year effects,
##  however.


####
####  GROWTH REGRESSIONS ----
####

#' Estimate growth regression coefficients using lmer.
#' 
#' @author Andrew Tredennick
#' @param df Time series dataframe of log genet sizes at t0 (logarea.t0) and 
#'           t1 (logarea.t1), effects of crowding (W), year as a factor (fyear),
#'           age (age), and group ID (Group).
#' @param juvenile  Logical value for fitting juvenile model (TRUE) or not (FALSE).
#' @param size_by_year Logical value for fitting size by year interaction (TRUE) or not (FALSE).
#' @return A lmer model object for adults or fitted distribution for juveniles.
fit_growth <- function(df, juvenile, size_by_year, seedling_size=NULL) {
  # Require lme4 for GLMMs
  require(lme4)
  
  # Check to make sure all needed columns are in the dataframe
  checknames <- c("logarea.t0", "logarea.t1", "W", "Group", "fyear", "age")
  checkits <- which(checknames %in% names(df))
  if(length(checkits) != length(checknames)) {
    stop(paste0("Missing columns in dataframe: ", checknames[-checkits]))
  }
  
  # If not a juvenile model, fit GLMM
  if(juvenile == FALSE) {
    
    # If sizeXyear interaction requested
    if(size_by_year == TRUE) {
      fit <- lmer(logarea.t1 ~ logarea.t0 * W + (1|Group) + (logarea.t0|fyear), 
                  data = df)
      
      # Fit variance
      x <- predict(fit)
      y <- resid(fit)  
      outVar <- nls((y^2) ~ a*exp(b*x), start=list(a=1,b=0))
      out_var <- c(sigma.a=coef(outVar)[1], sigma.b=coef(outVar)[2])
      # con <- gamlss.control(n.cyc = 1000)
      # tmp_data <- data.frame(x = x, y = y)
      # fit.t <- gamlss(y~1, 
      #                 sigma.fo=~x,
      #                 nu.fo=~1, 
      #                 family=TF(sigma.link="log"), 
      #                 data=tmp_data, 
      #                 control=con)
      # 
      # sigma.a <- fit.t$sigma.coefficients[1]
      # sigma.b <- fit.t$sigma.coefficients[2]
      # nu      <- fit.t$nu.coefficients[1]
      # out_var <- c(sigma.a=sigma.a, sigma.b=sigma.b, sigma.nu=nu)
    }
    
    # If sizeXyear interaction not requested
    if(size_by_year == FALSE) {
      fit <- lmer(logarea.t1 ~ logarea.t0 * W + (1|Group) + (1|fyear), 
                  data = df)
      # Fit variance
      x <- predict(fit)
      y <- resid(fit)
      outVar <- nls((y^2) ~ a*exp(b*x), start=list(a=1,b=0))
      out_var <- c(sigma.a=coef(outVar)[1], sigma.b=coef(outVar)[2])
      # con <- gamlss.control(n.cyc = 1000)
      # tmp_data <- data.frame(x = x, y = y)
      # fit.t <- gamlss(y~1, 
      #                 sigma.fo=~x,
      #                 nu.fo=~1, 
      #                 family=TF(sigma.link="log"), 
      #                 data=tmp_data, 
      #                 control=con)
      # 
      # sigma.a <- fit.t$sigma.coefficients[1]
      # sigma.b <- fit.t$sigma.coefficients[2]
      # nu      <- fit.t$nu.coefficients[1]
      # out_var <- c(sigma.a=sigma.a, sigma.b=sigma.b, sigma.nu=nu)
    }
    
    # Return the lmer model object
    return(list(fit=fit, out_var=out_var))
  }
  
  # If a juvenile model, fit distribution
  if(juvenile == TRUE) {
    nochangers <- which(df$logarea.t1-df$logarea.t0==0)
    df$grow_yn <- 0
    df$grow_yn[-nochangers] <- 1
    
    # print(paste0("Warning: This growth analysis is being done on ", nrow(df)," individuals."))
    # print(paste0("And, only ", nrow(df)-length(nochangers)," individuals actually change size!"))
    # print("But, I'll estimate a probability density anyway")
    
    sizes <- na.omit(df$logarea.t1)
    small <-  which(sizes < log(0.26))
    frac_small <-  length(small)/length(sizes)
    if(length(small) > 0){
      bigsizes <- sizes[-small]
      NegLogLik2 <- function(p,mu2,sigma2) { 
        z <- p*dnorm(sizes, mean=log(0.25), sd=0.2) + (1-p)*dnorm(sizes, mean=mu2, sd=sigma2)
        return(-sum(log(z)))
      }
      MLE_fit <- stats4::mle(minuslogl=NegLogLik2, start= list(p=frac_small,mu2=mean(bigsizes), sigma2=sd(bigsizes)), 
                     method="Nelder-Mead")
    }
    if(length(small) == 0){
      bigsizes <- sizes
      NegLogLik2 <- function(mu2,sigma2) { 
        z <- dnorm(sizes, mean=mu2, sd=sigma2)
        return(-sum(log(z)))
      }
      MLE_fit <- mle(minuslogl=NegLogLik2, start= list(mu2=mean(bigsizes), sigma2=sd(bigsizes)), 
                     method="Nelder-Mead")
    }
  
    return(coef(MLE_fit))
  }
  
} # end grow_fit function



####
####  SURVIVAL REGRESSIONS ----
####

#' Estimate survival regression coefficients using lmer.
#' 
#' @author Andrew Tredennick
#' @param df Time series dataframe of log genet sizes at t0 (logarea) effects of
#'           crowding (W), year as a factor (fyear), age (age), and 
#'           group ID (Group).
#' @param juvenile  Logical value for fitting juvenile model (TRUE) or not (FALSE).
#' @param size_by_year Logical value for fitting size by year interaction (TRUE) or not (FALSE).
#' @return A lmer model object for adults or fitted distribution for juveniles.
fit_survival <- function(df, juvenile, size_by_year) {
  # Require lme4 for GLMMs
  require(lme4)
  
  # Check to make sure all needed columns are in the dataframe
  checknames <- c("logarea", "W", "Group", "fyear", "age")
  checkits <- which(checknames %in% names(df))
  if(length(checkits) != length(checknames)) {
    stop(paste0("Missing columns in dataframe: ", checknames[-checkits]))
  }
  
  # If not a juvenile model, fit GLMM
  if(juvenile == FALSE) {
    
    # If sizeXyear interaction requested
    if(size_by_year == TRUE) {
      fit <- glmer(survives ~ logarea * W + (1|Group) + (logarea|fyear), 
                  data = df, family = "binomial")
    }
    
    # If sizeXyear interaction not requested
    if(size_by_year == FALSE) {
      fit <- glmer(survives ~ logarea * W + (1|Group) + (1|fyear), 
                   data = df, family = "binomial")
    }
    
    # Return the lmer model object
    return(fit)
  }
  
  # If a juvenile model, fit simpler model
  if(juvenile == TRUE) {
    fit <- glmer(survives ~ W + (1|Group) + (1|fyear), 
                 data=df, family="binomial")
    
    return(fit)
  }
  
} # end fit_survival function



####
####  RECRUITMENT REGRESSIONS ----
####
#' Estimate recruitment regression coefficients using JAGS.
#' 
#' @author Andrew Tredennick
#' @param df Time series dataframe of log genet sizes at t0 (logarea) effects of
#'           crowding (W), year as a factor (fyear), age (age), and 
#'           group ID (Group).
#' @param n_init Number of iterations for adaptation phase.
#' @param n_update Number of iterations for update phase.
#' @param n_iters Number of iterations to sample from the posterior distribution.
#' @param n_thin Number of iterations by which to thin the posterior distribution.
#' @return A matrix of statistical results by fitted parameter.
fit_recruitment <- function(df, 
                            n_adapt=5000,
                            n_update=10000, 
                            n_iters=20000, 
                            n_thin=50, 
                            model_file) {
  require(coda)
  require(rjags)
  require(parallel)
  
  ##  Extract objects for JAGS data list
  y        <- df$recruits
  parents1 <- df$cover/100
  parents2 <- df$Gcov/100
  year     <- as.numeric(as.factor(df$year))
  Nyrs     <- length(unique(year))
  N        <- nrow(df)
  Group    <- as.numeric(as.factor(df$Group))
  #if(min(Group) == 0) Group <- Group+1
  Ngroups  <- length(unique(Group))

  
  ##  Fit as negative binomial with random effects in JAGS
  data_list = list(N = N, 
                   y = y, 
                   parents1 = parents1, 
                   parents2 = parents2, 
                   year = year, 
                   Nyrs = Nyrs,
                   Ngroups = Ngroups, 
                   Group = Group)
  
  params <- c("intcpt.yr",
              "intcpt.mu",
              "intcpt.tau",
              "intcpt.gr",
              "g.tau",
              "dd",
              "theta",
              "u",
              "pvalue_mean",
              "pvalue_sd")

  init_func <- function(Nyrs, Ngroups){
    return(list(
      intcpt.yr = runif(Nyrs,-2,2),
      intcpt.mu = runif(1,-2,2),
      intcpt.tau = runif(1,0.01,5),
      intcpt.gr = runif(Ngroups,-2,2),
      g.tau = runif(1,0.01,5),
      dd = runif(1,-3,1),
      theta = runif(1,1,5)))}
  
  cl <- makeCluster(3)
  clusterExport(cl, c("data_list",
                      "init_func",
                      "n_adapt",
                      "n_update",
                      "n_iters",
                      "n_thin",
                      "model_file",
                      "Nyrs",
                      "Ngroups",
                      "params"),
                envir=environment())

  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <-  jags.model(model_file, data = data_list, inits = init_func(Nyrs,Ngroups),
                    n.chains = 1, n.adapt = n_adapt)
    update(jm, n.iter = n_update)
    zmCore = coda.samples(jm, variable.names = params,
                          n.iter = n_iters, thin = n_thin)
    return(as.mcmc(zmCore))
  })

  out_mcmc <- mcmc.list(out)
  stopCluster(cl)
  out_stats <- summary(out_mcmc)$stat
  return(list(out_stats,out_mcmc))
  
} # end function
