#' Calculate the Meng and Rubin's (1992) D-statistic (denoted D3) for a specified replication, data condition, number of imputations
#'
#'
#' @param pooled_tracker_df (data.frame) Tracker with both free and constrained (fixed) loglikelihood values
#' @param rep_D3 (integer) Replication to calculate D3
#' @param z_D3 (integer) Data condition number to calculate D3
#' @param Mmax_D3 (integer) Number of imputations to use. Note that a D3 may be calculated using fewer imputed data sets if not all data sets converged.
#'
#' @return (data.frame) A data frame with K-1 rows the following columns:
#' rep - (integer) Replication
#' z - (integer) data condition
#' M_condition - (integer) Number of requested to us imputations to use.
#' H0 -  Reduced model number of components.
#' H1 - Full model number of components.
#' dQ - Difference in number of parameters between full and reduced model.
#' M_D_3 - Number of imputed data sets used to calculate D_3 and ARIV_2 statistics. May be less than M_condition (and Mmax_D3) specified if the model failed converge for some of the data sets.
#' ARIV_2 - Average relative increase in variance statistic calculated a proposed by Meng and Rubin (1992)
#' v_4 - The denominator degrees of freedom for the F reference distribution to which the sampling distriubtion of D_3 is supposed to follow (see Meng & Rubin, 1992)
#' D_3 - The Meng and Rubin (1992) D-statistic. Should be nonnegative.
#' @export
#'
#' @examples
#' Meng_Rubin_D3(pooled_tracker_df = lpa_pool_tracker_df2, rep_D3 = unique(lpa_pool_tracker_df2$rep),z_D3 = 10, Mmax_D3 = 100)

Meng_Rubin_D3<-function(pooled_tracker_df, rep_D3,z_D3,Mmax_D3){

  require(data.table)
  require(plyr); require(tidyverse)


  ####### Preliminary checks #####
  # Check that all inputs are of the right type.
    if(!is.data.frame(pooled_tracker_df)){stop("pooled_tracker_df is not a data.frame.")}
    if(length(rep_D3)>1){stop("rep_D3 must be single-valued integer")}
    if(length(z_D3)>1){stop("z_D3 must be single-valued integer.")}
    if(length(Mmax_D3)>1){stop("Mmax_D3 must be single-valued integer.")}

  # Check that the integers are in the the tracker
    if( !(rep_D3 %in% pooled_tracker_df$rep) ){stop(" rep_D3 value not found in pooled_tracker_df$rep. ") }
    if( !(z_D3 %in% pooled_tracker_df$z) ){stop(" z_D3 value not found in pooled_tracker_df$z.") }


  ###### Lower-level function to calculate D3 ########
    parse_and_calc_D3<-function(the_tracker){
      # This function parses the tracker into relelvant data sets for the sequential nested model tests.
      # After parsing, the D3 statistics and relevant ancillaries are calculated.

      # Inputs: the_tracker - A data.frame at a given replication, z, and M_condition value.
      #                       The data.frame must be in long format in that the free and constrained (i.e., fixed)
      #                       log-likeihood values must be side-by-side
      # Notes: That when calculating Meng and Rubin's D statistic, only the imputed datastes that converged for BOTH the FULL
      #        and the constrained modelS are used to in the calculation.

      # Create a list that does the parsing and calculating
      parsed_list <-
        lapply(X = seq(min(the_tracker$kfit), max(the_tracker$kfit)-1,by = 1),
               FUN = function(x){
                 tmp_tracker = the_tracker %>% filter(kfit == x | kfit==x+1)
                 tmp_df = ddply(tmp_tracker,.(kfit),summarize,
                                MM_max = max(mm))
                 MM_minmax = min(tmp_df$MM_max)
                 tmp_aggbar = tmp_tracker %>% filter(mm<=MM_minmax) %>%
                   transform(H0 = x, H1 = x+1) %>%
                   ddply(., .(H0,kfit), summarize,
                         H1 = mean(H1), M = max(mm),Parameters = mean(Parameters),
                         LL_bar = mean(LL),
                         LL_c_bar = mean(LL_c))

                 tmp_D3 = tmp_aggbar[-1, c("Parameters","LL_bar","LL_c_bar")] - tmp_aggbar[-2, c("Parameters","LL_bar","LL_c_bar")]
                 tmp_D3 = tmp_D3 %>%  transform(LR_bar = 2*LL_bar,  LR_c_bar= 2*LL_c_bar, dQ = Parameters,
                                                H0 = tmp_aggbar$H0[1], H1 = tmp_aggbar$H1[1], M = tmp_aggbar$M[1] )
                 tmp_D3 = tmp_D3 %>% transform(ARIV_2 = (M+1)/(dQ*(M-1))*(LR_bar-LR_c_bar) )
                 tmp_D3 = tmp_D3 %>% transform(D_3  = LR_c_bar / (dQ*(1+ARIV_2)), # Enders (2010) Eqn 8.32. Note substitution k -> dQ
                                               v_4 = ifelse(dQ*M-dQ>4, #see Enders (2010) pg. 241
                                                            4+(dQ*M-dQ-4)*(  1+(1-(2/(dQ*M-dQ)))*(1/ARIV_2) )^2, # Enders (2010) Eqn 8.34
                                                            0.5*(dQ*M-dQ)*(1+1/dQ)*(1+1/ARIV_2)^2) # Enders (2010) Eqn 8.35
                                      )
                 return(tmp_D3 %>% select(H0,H1,dQ,M,ARIV_2,v_4,D_3))
               }#end FUN = function(x)
        ) #end lapply
      return(parsed_list %>% rbindlist())
    } #end parse_and_calc_D3 function

  ###### Step 1: Wrangle the tracker to put free/fix in wide format ########
    # Freely estimated likelihoods
    D3_tracker_free = pooled_tracker_df %>%
      filter(rep == rep_D3, z==z_D3, mm<=Mmax_D3, fix == FALSE, converged == TRUE) %>%
      select(rep,z,converged,fix,kfit,mm, Parameters,LL) %>%
      arrange(kfit,mm)
    # Likelihood estimates constrained [fixed] to the RR pooled values
    D3_tracker_fix = pooled_tracker_df %>%
      filter(rep == rep_D3, z==z_D3, M_condition==Mmax_D3, mm<=Mmax_D3, fix == TRUE, converged == TRUE) %>%
      select(rep,z,converged,fix,kfit,mm, Parameters,LL) %>%
      arrange(kfit,mm)
    names(D3_tracker_fix)[which(names(D3_tracker_fix)=="LL")] = "LL_c"
    # Put free/fix in wide format
    D3_tracker = D3_tracker_free %>%
      left_join(D3_tracker_fix %>% select(rep,z,kfit,mm,LL_c),
                by = c("rep","z","kfit","mm")) %>%
      select(rep,z,converged,kfit,mm,Parameters,LL,LL_c)

  ###### Step 2: Calculate the D3 statistic ########
    D3_df = parse_and_calc_D3(the_tracker = D3_tracker) %>% data.frame()

  ##### Step 3: Clean and return #################
    names(D3_df)[names(D3_df)=="M"]="M_D_3"
    D3_df = D3_df %>% transform(rep = rep_D3, z = z_D3, M_condition = Mmax_D3) %>%
      select(rep,z,M_condition,H0,H1,dQ,M_D_3,ARIV_2,v_4,D_3)

  return(D3_df)
}
