#' Resolve label switching with population values known.
#'
#' @param out_ftc (list) An object from fit_til_convergence. 
#' @param name (character) Defaults to "glpk"
#' @param t_max (integer) Maximum time for estimation using "glpk". Defaults to 30 sec.
#' @param approximate (integer) Allow for approximate solution. Defauts to 0 (i.e., FALSE)
#' @param round_cplex (integer) 
#' @param trace_cplex (integer)
#' @param ... Other arguments to send to create_naiveMplus_inpfile()
#' @return (list) A list with the following things:
#'           (a) class_switched - (logical) identifies if class switching had occurred.
#'           (b) map_from - (vector) always just 1:K showing class definitions from sample
#'           (c) map_to - (NULL or vector) same as map_from if no class switching.
#'                      Otherwise it is the corresponding population-based classes from sample classes
#'           (d) out_Mplus- (list) readModels() fitting model with starting values that resolve class switching 
#' @export
#' @examples
#'
#'
#'
#'
resolve_label_switch<-function(out_ftc,
                               name = "glpk", t_max = 30, approximate = 0, round_cplex = 0, trace_cplex = 0, ...){

      require(gaussDiff); require(designmatch); require(plyr); require(MplusAutomation)
      parameters_df = out_ftc$out_readModels$parameters$unstandardized
      
      # Get the distance matrix
      agg_cprob = ddply(out_ftc$out_readModels$savedata[,c("SUBPOP","CPROB1","CPROB2","CPROB3")],
                        .(SUBPOP), summarize, CPROB1 = mean(CPROB1), CPROB2 = mean(CPROB2), CPROB3 = mean(CPROB3))
      DMAT = 1-as.matrix(agg_cprob[,-1])
  
      # Get number of class indicators and number of classes based on data condtion z
      out_get_FMM = out_ftc$out_get_FMM
      K_z = out_ftc$out_get_FMM$K_z
      
      #Identify if class switching occured
      class_switched = FALSE;
      if(prod(apply(DMAT,2,"which.min") == c(1:K_z))==0){class_switched = TRUE}

      # If class switching occured, relable classes according to population labels
      map_from = 1:K_z
      map_to = 1:K_z
      if (class_switched){
         obj.bmatch = bmatch(t_ind = c(rep(1,K_z), rep(0,K_z)),dist_mat = t(DMAT),
                             solver = list(name = name, t_max = t_max, approximate = approximate, round_cplex = round_cplex,   trace_cplex = trace_cplex),
                             total_groups = K_z)
         map_to =  obj.bmatch$c_id-K_z
      }
      
      # If class switching, refit the model
      if (class_switched){
        # Extract model text
        model_txt = out_ftc$model_txt
        
        # Create a properly switched model text
          list_model_txt = list(NULL)
          for(k in 1:K_z){
            i_start = which(endsWith(model_txt, paste0("%C#",k,"%")))+1
            i_end = ifelse(k==K_z,length(model_txt),which(endsWith(model_txt, paste0("%C#",k+1,"%"))))-1
            k_to = map_to[k]
            list_model_txt[[k_to]] = c(paste0("%C#",k_to,"%"), model_txt[i_start:i_end])
          }
          pi_vec = gamma2pi(parameters_df$est[parameters_df$LatentClass=="Categorical.Latent.Variables"])
          pi_vec = pi_vec[map_to]
          gamma_vec = log(pi_vec[-K_z]/pi_vec[K_z])
          model_txt_to = c("MODEL:","\n",
                           "%OVERALL%",paste0("[C#",seq(1,K_z-1),"*",round(gamma_vec,5),"];"), "",
                           unlist(list_model_txt))
          
        # Fit the model specified by tne new model text
          inp_txt = do.call(what = "create_naiveMplus_inpfile",
                            args = list(z = out_ftc$z, 
                                             out_get_FMM = out_ftc$out_get_FMM,
                                             dffolderfiles = out_ftc$dffolderfiles, 
                                             temp_wd_p = out_ftc$temp_wd_p,
                                             starts_txt = "0;", type_imputation = FALSE, 
                                             Model_txt = model_txt_to, ...))
          runModels(target = out_ftc$target_wd, 
                    filefilter = list.files(path = out_ftc$target_wd, pattern = out_ftc$pattern_inp),
                    logFile = NULL)
          
        # Read in the model
          out_file= list.files(path = out_ftc$target_wd, pattern = out_ftc$pattern_out, full.names = FALSE)
          out_readModels = readModels(target = out_ftc$target_wd, filefilter = gsub(".out", "", out_file))
          
      } else {
        
        out_readModels = out_ftc$out_readModels
        
      }# end class_switched

      # Create and return out list
      out_list = list(class_switched = class_switched,
                      map_from = map_from, map_to = map_to, 
                      out_readModels = out_readModels)

      return(out_list)

}
