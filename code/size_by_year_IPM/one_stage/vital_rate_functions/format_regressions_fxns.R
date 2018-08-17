#' Formates growth regression parameters.
#' 
#' @author Andrew Tredennick
#' @param model_fit Lmer model object.
#' @return A dataframe with coefficients for each year.
format_growth_params <- function(model_fit) {
  require(tidyverse)
  require(dplyr)
  model_params <- coef(model_fit[["fit"]])$fyear %>%
    mutate(sigma.a = model_fit$out_var[1],
           sigma.b = model_fit$out_var[2],
           sigma.nu = model_fit$out_var[3],
           year_id = rownames(coef(model_fit[["fit"]])$fyear),
           year = 1900 + as.numeric(rownames(coef(model_fit[["fit"]])$fyear))) %>%
    dplyr::rename(intercept = `(Intercept)`,
                  logarea.W = `logarea.t0:W`)
  
  return(model_params)
}



#' Formates survival regression parameters.
#' 
#' @author Andrew Tredennick
#' @param model_fit Lmer model object.
#' @param all_years A vector of year ids for all the years present in the data.
#' @return A dataframe with coefficients for each year.
format_surv_params <- function(model_fit, all_years) {
  require(tidyverse)
  require(dplyr)
  
  model_params <- coef(model_fit)$fyear %>%
    mutate(year_id = rownames(coef(model_fit)$fyear),
           year = 1900 + as.numeric(rownames(coef(model_fit)$fyear))) %>%
    dplyr::rename(intercept = `(Intercept)`,
                  logarea.W = `logarea:W`)
  
  return(model_params)
}



#' Formates recruitment regression parameters.
#' 
#' @author Andrew Tredennick
#' @param model_fit Model fit object. 
#' @param years Vector of years from the data, sorted.
#' @return A dataframe with coefficients for each year.
format_rec_params <- function(model_fit, years) {
  require(tidyverse)
  require(dplyr)
  
  year_effects <- as.data.frame(model_fit) %>%
    dplyr::select(Mean) %>%
    mutate(parameter = row.names(model_fit)) %>%
    filter(grepl("intcpt.yr",parameter)) %>%
    mutate(year = years) %>%
    dplyr::select(-parameter) %>%
    dplyr::rename(intcpt.yr = Mean)
  
  fixed_effects <- as.data.frame(model_fit) %>%
    dplyr::select(Mean) %>%
    mutate(parameter = row.names(model_fit)) %>%
    filter(!grepl("intcpt.gr",parameter)) %>%
    filter(!grepl("intcpt.yr",parameter))
  long_names <- expand.grid(fixed_effects$parameter, years)
  long_effects <- expand.grid(fixed_effects$Mean, years)
  long_effects$parameter <- long_names$Var1
  colnames(long_effects) <- c("mean", "year", "parameter")
  
  long_effects <- long_effects %>%
    spread(parameter, mean)
  
  model_params <- year_effects %>%
    left_join(long_effects, by = "year") %>%
    dplyr::select(-year, year)
}


