# 2.0.0. Set working directory to current path -----------------------------------
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(this.dir)


# 2.1.0. Library -----------------------------------------------------------------
library(DEoptim)  # Library for global optimization by differential evolution via the Differential Evolution algorithm   
library(deSolve)  # Library for solving differential equations
library(ggplot2) 
library(reshape2) 
library(xlsx)

# 2.2.0. Load data ---------------------------------------------------------------
Dataset <- read.csv("chamber_processed_data.csv", sep = ";", dec = ".", header = TRUE)   # Output data from Script_1


# 2.3.0. Protein turnover model --------------------------------------------------
Model1 <- function(t, NX, parms) {
  with(as.list(parms),{
    
    # Equation M5
    Synthesis <- 2*S_max/(1 + exp(t*td))
    
    # Equation M6
    Degradation <- Dr*NX
    
    # Equation M4  
    dNX <- Synthesis - Degradation
    
    list(dNX)
  })} 


# 2.4.0. Function for fitting ----------------------------------------------------
fitting_function <- function(data = Dataset,
                             list = list(Dataset$VarietyID, Dataset$LightID, Dataset$NitrogenID),     
                             Time_step = 0.1,            # Time step (oC) for simulation
                             Duration_min = 0,           # Starting time (oC) for simulation
                             Duration_max = 600,         # Total duration (oC) for simulation
                             NX0 = 0,                    # Initial value of nitrogen pool
                            
                             # Lower and upper boundaries for optimization 
                             lower = list(NV = c(S_max = 0.01, td = 0.0001, Dr = 0.001),
                                          NJ = c(S_max = 0.01, td = 0.0001, Dr = 0.001),
                                          NC = c(S_max = 0.01, td = 0.0001, Dr = 0.001)), 
                             upper = list(NV = c(S_max = 2, td = 0.01, Dr = 0.03),
                                          NJ = c(S_max = 2, td = 0.01, Dr = 0.03),
                                          NC = c(S_max = 2, td = 0.01, Dr = 0.03)),       
                             
                             NP = 30,                    # Set to 10 times the number of parameters
                             itermax = 20,               # Maximum number of procedure iterations
                             
                             # Set name of the output file
                             filename = paste0("fitting_", Sys.Date(),".csv")) {
  
  
  x2 <- do.call('rbind', lapply(split(data, list, drop = TRUE), function(x){            
    ## x <- split(data, list, drop = TRUE)[[1]]                                   
    
    V <- as.character(x$VarietyID[1])
    LLV <- round(mean(x$LightLevel_mol_m2_d), 1) 
    NLV <- round(mean(x$NitrogenLevel_Mm), 1)    
    
   x1 <- do.call('rbind', lapply(list("NV", "NJ", "NC"), function(nx){        
     ## nx = "NC"                                      
      
     # Print values
     cat("V =", V, "L =", round(LLV,1), "N =", NLV, nx, " \n")       
     
     # Select leaf age and measured photosynthetic nitrogen data
     df <- na.omit(subset(x, select = c("LeafAge_oCd", paste0(nx, "_mmol_m2"))))
      
     # Function to calculate Root-Mean-Square Deviation (RMSD) between predicted and measured data
      ssq <- function(parms){
        ## parms <- upper[nx][[1]]  
        
        # Time steps for simulation, should include the time points where measured data are available
        t <- sort(unique(c(seq(Duration_min, Duration_max, Time_step), df$LeafAge_oCd)))
        # Solve differential equations for a given set of parameters
        out <- data.frame(lsoda(y = c(NX = NX0), times = t, func = Model1, 
                                parms = list(S_max = parms[1], td = parms[2], Dr = parms[3])))    
        
        # Select predicted data at the time points where measured data are available
        outdf <- out[out$time %in% df$LeafAge_oCd,]
        preddf <- outdf[order(outdf$time),]
        
        # Measured data 
        measdf <- melt(df, id.var = "LeafAge_oCd", variable.name = "param", value.name = "NX")
        measdf <- measdf[order(measdf$LeafAge_oCd),]
        
        # Combine predicted and measured data
        out_merge <- merge(preddf, measdf, by.x = "time", by.y = "LeafAge_oCd")
        
        # Calculate RMSD between predicted and measured data
        rmsd <- sqrt(mean((out_merge$NX.x - out_merge$NX.y)^2, na.rm = TRUE))
        
        return(rmsd)
      }
      
      # Start optimization
      opti <- DEoptim(ssq, lower = lower[nx][[1]], upper = upper[nx][[1]], DEoptim.control(VTR = 0, NP = NP, itermax = itermax)) 
      # View optimization result
      summary(opti)
      # Output optimization result
      parms <- c(opti$optim$bestmem)
      Mean <- round(mean(df[,2]), 2)
      RMSD <- round(opti$optim$bestval, 2)
      Accuracy <- 1 - RMSD/Mean
      Iteration <- opti$optim$iter 
      
      # Plot predicted vs measured data
      t <- seq(Duration_min, Duration_max, Time_step)
      out <- data.frame(lsoda(y = c(NX = NX0), times = t, func = Model1, parms = parms))
      pred <- melt(out, id.var = c("time"), variable.name = "par", value.name = nx)
      meas <- melt(df, id.var = c("LeafAge_oCd"), variable.name = "par", value.name = nx)
      value_max <- max(c(unlist(meas[nx]), unlist(pred[nx])))*1.5
      
      p <- ggplot(meas, aes(x = LeafAge_oCd, y = unlist(meas[nx]))) +
        geom_point(color = "red", size = 5) +
        geom_line(data = pred, aes(x = time, y = unlist(pred[nx])), size = 1.5) +
        xlim(0, 600) +
        ylim(0, value_max) +
        ggtitle(paste0("VAR = ", V, ", mean LLV = ", round(LLV, 1), ", mean NLV = ", NLV, ", ", nx)) +
        labs(x = "Leaf age (degree days)", y = paste0(nx, " (mmol/m2)")) +
        theme_bw() +
        theme(
          aspect.ratio = 1/1,
          plot.title = element_text(color = "#993333", size = 14, face = "bold.italic"),
          axis.title.x = element_text(color = "black", size = 20, face = "bold"),
          axis.title.y = element_text(color = "black", size = 20, face = "bold"),
          axis.text.x = element_text(face = "plain", color = "black", size = 20),
          axis.text.y = element_text(face = "plain", color = "black", size = 20)) +
        annotate("text", x = 430, y = value_max*0.9, label = paste0("Mean = ", Mean,  " \n", "RMSD = ", RMSD), size = 8)
      print(p)
      
      # Output values 
      x1 <- data.frame(VarietyID = V, LLV = round(LLV, 2), NLV = round(NLV, 2), NX = nx, 
                        S_max = parms["S_max"], td = parms["td"], Dr = parms["Dr"], 
                        Mean = Mean, RMSD = RMSD, Accuracy = Accuracy, 
                        NP = NP, Iteration = Iteration, Time_stamp = Sys.time())

    })) 
   x1}))
  
  write.table(x2, file = filename, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, quote = TRUE) 
  }



###----- IN SILICO TEST simulation with various fd -----###
####-----------------------------------------------------####
####----- Arguments CHANGABLE for parameterization ------####
####-----------------------------------------------------####
# 2.4.1. Parameterization step 1 fitting -----------------------------------------

# Set seed
set.seed(345)

fitting_function(data = Dataset,                                              
                 list = list(Dataset$VarietyID),                              
                 Time_step = 0.1,                                
                 Duration_min = 0,                                
                 Duration_max = 600,                             
                 NX0 = 0,                                        
         
                 lower = list(NV = c(S_max = 0.01, td = 0.0001, Dr = 0.001),
                              NJ = c(S_max = 0.01, td = 0.0001, Dr = 0.001),
                              NC = c(S_max = 0.01, td = 0.0001, Dr = 0.001)), 
                 upper = list(NV = c(S_max = 2, td = 0.01, Dr = 0.0195),
                              NJ = c(S_max = 2, td = 0.01, Dr = 0.0195),
                              NC = c(S_max = 2, td = 0.01, Dr = 0.0195)),      
         
                 NP = 30,                                         
                 itermax = 10,                                    
                 filename = "parameterization_step_1.csv")


# 2.4.2. Parameterization step 2 fitting -----------------------------------------

# Input result form step 1
output_1 <- read.csv("parameterization_step_1.csv", sep = ",", dec = ".")

# Set seed
set.seed(345)

fitting_function(data = Dataset,
                 list = list(Dataset$VarietyID, Dataset$LightID, Dataset$NitrogenID),    
                 Time_step = 0.1,                                
                 Duration_min = 0,                                
                 Duration_max = 600,                             
                 NX0 = 0,                                        
                 
                 lower = list(NV = c(S_max = 0.01, 
                                     td = unique(output_1$td[output_1$NX == "NV"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NV"])),
                              NJ = c(S_max = 0.01, 
                                     td = unique(output_1$td[output_1$NX == "NJ"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NJ"])),
                              NC = c(S_max = 0.01, 
                                     td = unique(output_1$td[output_1$NX == "NC"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NC"]))), 
                 upper = list(NV = c(S_max = 2, 
                                     td = unique(output_1$td[output_1$NX == "NV"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NV"])),
                              NJ = c(S_max = 2, 
                                     td = unique(output_1$td[output_1$NX == "NJ"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NJ"])),
                              NC = c(S_max = 2, 
                                     td = unique(output_1$td[output_1$NX == "NC"]), 
                                     Dr = unique(output_1$Dr[output_1$NX == "NC"]))),      
                 
                 NP = 30,                                         
                 itermax = 10,                                    
                 filename = "parameterization_step_2.csv")

# 2.4.3. Parameterization step 3 fitting -----------------------------------------

# Input result form step 2
output_2 <- read.csv("parameterization_step_2.csv", sep = ",", dec = ".")

# Set seed
set.seed(345)

# Parameterize potential maximum synthesis rate (Smm), light (kI) and nitrogen (kN) effect
parameterize_result <- do.call('rbind', lapply(split(output_2, list(output_2$VarietyID, output_2$NX), drop = TRUE), function(x){
  ## x <- split(output_2, list(output_2$VarietyID, output_2$NX), drop=TRUE)[[1]]

  fit_res <- nls(S_max ~ Smm*kI*LLV/(Smm+kI*LLV)*NLV/(kN + NLV), 
                 start = c(Smm = max(x$S_max)*1.2, kI = 0.1, kN = 0.5), 
                 data = x)
  ## summary(fit_res)
  
  x$Smm <- round(summary(fit_res)$coefficients[,1]["Smm"], 4)
  x$kI <- round(summary(fit_res)$coefficients[,1]["kI"], 4)
  x$kN <- round(summary(fit_res)$coefficients[,1]["kN"], 4)
  
  x$Smm_se <- round(summary(fit_res)$coefficients[,2]["Smm"], 5)
  x$kI_se <- round(summary(fit_res)$coefficients[,2]["kI"], 5)
  x$kN_se <- round(summary(fit_res)$coefficients[,2]["kN"], 5)
  
  x$Smm_pv <- round(summary(fit_res)$coefficients[,4]["Smm"], 4)
  x$kI_pv <- round(summary(fit_res)$coefficients[,4]["kI"], 4)
  x$kN_pv <- round(summary(fit_res)$coefficients[,4]["kN"], 4)

  x
  }))

parameterize_result_output <- unique(subset(parameterize_result, select = c("VarietyID", "NX", "td", "Dr", "Smm", "kI", "kN", "Smm_se", "kI_se", "kN_se", "Smm_pv", "kI_pv", "kN_pv")))


# 2.5.0. Output parameterization result to file ------------------------------------
write.table(parameterize_result_output, file = "parameterize_result_output.csv", sep = ";", dec = ".", row.names = FALSE, col.names = TRUE, quote = TRUE) 






