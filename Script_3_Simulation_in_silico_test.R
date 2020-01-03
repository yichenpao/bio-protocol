# 3.0.0. Set working directory to current path -----------------------------------
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(this.dir)


# 3.1.0. Library -----------------------------------------------------------------
library(DEoptim)  # Library for global optimization by differential evolution via the Differential Evolution algorithm
library(deSolve)  # Library for solving differential equations
library(dplyr)    # group_by() function from 'dplyr' does not work if 'plyr' is loaded!
## detach("package:plyr", unload=TRUE)
library(ggplot2)
library(magrittr) 
library(xlsx)


####-----------------------------------------------------####
####---- CHANGABLE Constants and empirical functions ----####
####-----------------------------------------------------####
# 3.2.0. Constants for simulation process ------------------------------------------------

# Initial value of nitrogen pool
NX0 = 0

# Size of time step for simulation
Time_step_N = 0.1                  # Photosynthetic functional nitrogen simulation, oCd
Time_step_photo = 0.1              # Daily photosynthesis simulation, h


# 3.2.1. General constants --------------------------------------------------------

# Environemntal constant
Ca = 400             # atmospheric CO2 concentration, umol CO2/mol
D = 1.2              # leaf-to-air vapor pressure deficit, kPa
DL = 14.4            # day length, h 

# Plant constant
Temp_base = 10       # base temperature for cucumber, oC
leafx = 3                          # Leaf rank number appeared around the time of transplanting
PD = 1.33                          # plant/m2
k = 0.695                          # Light extinction coefficient, unitless

# Photosynthetic nitrogen and function capacity 
ChiV = 4.49          # carboxylation capacity per unit Rubisco nitrogen, umol CO2/mmol N/s 
ChiJ = 9.48          # electron transport capacity per unit electron transport nitrogen, umol e-/mmol N/s
ChiC = 0.03384       # conversion coefficients for chlorophyll per light harvesting component nitrogen, mmol Chl/mmol N
ChiCJ = 4.64*10^-4   # conversion coefficients for chlorophyll per electron transport nitrogen, mmol Chl/mmol N

# Photosynthesis and stomatal conductance 
beta = 0.5           # partitioning fraction of photons between photosystem II and I, unitless
theta = 0.7          # constant convexity factor describing the response of J to PPFD, unitless
phi = 0.425          # conversion efficiency of photons to J, umol e-/umol photon
gama_str = 43.02     # CO2 compensation point in the absence of mitochondrial respiration, umol CO2/mol
Kc = 404             # Michaelis-Menten constants of Rubisco for CO2, umol CO2/mol
Ko = 278             # Michaelis-Menten constants of Rubisco for O2, mmmol O2/mol
Oc = 210             # mole fraction of O2 at the site of carboxylation, mmol O2/mol
Km <- Kc*(1 + Oc/Ko) # Equation P8
g0 = 0.009           # species-specific coefficient for gs, mol CO2/m2/s
g1 = 3.51            # species-specific coefficient for gs, unitless
p1 = Ca - gama_str
p2 = Ca + Km
k1 = Ca - gama_str
k2 = Ca + 2*gama_str
gsc_b = (1 + g1/(D^0.5))/Ca


# 3.3.0. Empirical functions of LAI and EA, Equations G4 and G5 ------------------

# Estimate LAI using leaf age, Equation G4   
lai_age <- function(x, PD = PD){
  x <- x[order(-x$LeafNo),]                           # Arrange leaf rank from young to old
  x$LAI_cum_plant <- cumsum(x$LA_cm2_mean)*10^-4      # Cumulative LAI, m2 LA per plant
  x$LAI_cum_ground <- x$LAI_cum_plant*PD              # Cumulative LAI, m2 LA m-2 groud
  x <- x[order(x$LeafNo),]                            # Arrange leaf rank from old to young
  
  fit_LAI <- nls(LAI_cum_ground ~ SSlogis(LeafAge_oCd, lai_max, lai_t0, lai_scal), data = x)
  x$lai_max <- summary(fit_LAI)$coefficients[,"Estimate"]["lai_max"]
  x$lai_t0 <- summary(fit_LAI)$coefficients[,"Estimate"]["lai_t0"]
  x$lai_scal <- summary(fit_LAI)$coefficients[,"Estimate"]["lai_scal"]
  x$LAI <- x$lai_max/(1 + exp((x$lai_t0 - x$LeafAge_oCd)/x$lai_scal))   
  return(x)
}

# Fit EA using leaf age, Equation G5   
ea_age <- function(x){
  fit_EA <- nls(EA_degree_mean ~ (90 - ea_min*exp(-0.5*(log(LeafAge_oCd/ea_t0)/ea_scal)^2)), start = c(ea_min = 150, ea_t0 = x$LeafAge_oCd[x$EA_degree_mean == min(x$EA_degree_mean)], ea_scal = 1), data = x)
  x$ea_min <- summary(fit_EA)$coefficients[,"Estimate"]["ea_min"]
  x$ea_t0 <- summary(fit_EA)$coefficients[,"Estimate"]["ea_t0"]
  x$ea_scal <- summary(fit_EA)$coefficients[,"Estimate"]["ea_scal"]
  x$EA <- 90 - x$ea_min*exp(-0.5*(log(x$LeafAge_oCd/x$ea_t0)/x$ea_scal)^2)
  return(x)
}


# 3.3.1. Empirical functions of Rd and gm, Equations P9 and P10 -----------------

# Respiration rate (Rd), Equation P9
Rd_empirical <- function(data, environment_data){
  
  r_max = 0.308
  r_g = 4.16*10^-4
  r_m = 1.88*10^-4
  
  r_max*data$DPI_4d*exp(-r_g*data$DPI_4d*max(environment_data$GDD_leaf)) + r_m*data$DPI_4d*max(environment_data$GDD_leaf)
  }

# Mesophyll conductance (gm), Equation P10
gm_empirical <- function(data){
  
  gm_m = 0.00164 
  gm_m0 = 0.14 
  gm_t0 = 121  
  gm_scal = 0.86
  
  (gm_m*data$Nph_mmol_m2_adj + gm_m0)*exp(-0.5*(log(data$LeafAge_oCd/gm_t0)/gm_scal)^2)
  }
 



####-----------------------------------------------------####
####------------------ Data preparation -----------------####
####-----------------------------------------------------####
# 3.4.0. Load data ---------------------------------------------------------------

# Structural data of each leaf in the plant obtained by digitization 
df_Structure <- read.csv("example_greenhouse_structure_data.csv", sep = ";", dec = ".", header = TRUE)

# Select scenario to be simulated
df_Structure <- subset(df_Structure, MeasureDate == "10.05.2017")


# Environmental input
df_Environment <- read.table("example_greenhouse_environment_data.csv", sep = ";", dec = ".", header = TRUE)

# Parameters
df_para <- read.table("parameterize_result_output.csv", sep = ";", dec = ".", header = TRUE)  # Output data from Script_2


# 3.4.1. Data preparation --------------------------------------------------------

# Discard the first two leaves because they appeared before experiment started 
# and they have usually small area and little contribution to canopy photosynthesis at later stage
df_Structure <- subset(df_Structure, LeafNo != 1 & LeafNo != 2)
# Set date format
df_Structure$MeasureDate <- as.Date(df_Structure$MeasureDate, "%d.%m.%Y")

# Set date format
df_Environment$Date <- as.Date(df_Environment$Date,"%d.%m.%Y")

# Calculate mean nitrogen level and growing degree days (GDD) from environmental input data
df_Environment <- df_Environment[order(df_Environment$Date),]
colname_df <- names(df_Environment)
Supply_N_prefix <- "Supply_N_"
Substrate_N_prefix <- "Substrate_N_"
Temperature_prefix <- "Tmean_L_"

suffix1 <- sub(Supply_N_prefix, '', colname_df[which(grepl(Supply_N_prefix, colname_df))])
suffix2 <- sub(Temperature_prefix, '', colname_df[which(grepl(Temperature_prefix, colname_df))])

for (i in suffix1 ){
  df_Environment[,paste0('NLV_', i)] <-rowMeans(df_Environment[,c(paste0(Supply_N_prefix,i), paste0(Substrate_N_prefix,i))], na.rm = TRUE)
}

for (j in suffix2){
  df_Environment[,paste0("GDD_L_", j)] <- cumsum(df_Environment[,paste0(Temperature_prefix, j)] - Temp_base)
}

# Calculate mean and standard error of leaf area from the digitized data
df_Structure_mean <- as.data.frame(df_Structure %>%
                            group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo) %>%
                            summarise(LA_cm2_mean = mean(LA_cm2),
                                      LA_cm2_sd = sd(LA_cm2),
                                      EA_degree_mean = mean(EA_degree),
                                      EA_degree_sd = sd(EA_degree)))


####-----------------------------------------------------####
####----------------- Defined functions -----------------####
####-----------------------------------------------------####
# 3.5.0. Functions for simulation process ----------------------------------------

# Create sequence of time points for a given day length (DL) with given interval (Time_step_photo, h)
timepoints <- function(DL, Time_step_photo){ seq(12 - DL/2, 12 + DL/2, Time_step_photo) }  

# Create matrix of a parameter (para) in a data frame (data) with row number = time points and column number = leaf number
matrix.para <- function(data, para, para_length = para){
  matrix.x <- matrix(para, byrow = TRUE, nrow = length(timepoints(DL, Time_step_photo)), ncol = length(para_length))
  colnames(matrix.x) <- paste(data$MeasureDate, data$VarietyID, data$LightID, data$NitrogenID, data$LeafNo)
  rownames(matrix.x) <- timepoints(DL, Time_step_photo)
  return(matrix.x)
  }


# 3.5.1. Functions for phyllochron, leaf age and light extinction, Equations G1, G2, G3 and P11 ----

# Calculate phyllochron and leaf age, Equations G1, G2, G3
leaf_age <- function(x, df_Env_Exp){
  x$GDDcanopy <- subset(df_Env_Exp, Date == MeasureDate)[,paste0("GDD_L_", L)]
  x$Phyllochron <- x$GDDcanopy/(max(x$LeafNo) - (leafx - 1))
  x$LeafAge_oCd <- x$GDDcanopy - (x$LeafNo - leafx)*x$Phyllochron
  x$App_GDD <- x$GDDcanopy - x$LeafAge_oCd     # Calculate when the leaf appeared
  return(x)
}

# Light extinction, Equation P11
light_extinction <- function(k, t, lai_max, lai_t0, lai_scal, ea_min, ea_t0, ea_scal){
  as.numeric(exp(-k*(lai_max/(1 + exp((lai_t0 - t)/lai_scal))))*cos((90 - ea_min*exp(-0.5*(log(t/ea_t0)/ea_scal)^2))*pi/180))
  }


# 3.6.0. Function for simulating daily canopy carbon assimilation (DCA) ----------
simDCA <- function(data, list, DPI_multiplier_list){

out <- do.call('rbind', lapply(split(data, list, drop = TRUE), function(x){
                                                      
  # x <- split(Simulation_fd2, list(Simulation_fd2$ExpID, Simulation_fd2$MeasureDate, Simulation_fd2$VarietyID, Simulation_fd2$LightID, Simulation_fd2$NitrogenID, Simulation_fd2$fd), drop=TRUE)[[1]]
  
  # x <- split(data, list(data$ExpID, data$MeasureDate, data$VarietyID, data$LightID, data$NitrogenID, data$fd), drop=TRUE)[[1]]
  
  E <<- as.character(x$ExpID[1])                                                      
  MeasureDate <<- as.character(x$MeasureDate[1])
  V <<- as.character(x$VarietyID[1])
  L <<- as.character(x$LightID[1])
  N <<- as.character(x$NitrogenID[1])   
  
  x1 <- do.call('rbind', lapply(DPI_multiplier_list, function(DPI_multiplier){
   ## DPI_multiplier <- DPI_multiplier_list[[1]]
    
    x$DPI_multiplier <- DPI_multiplier
    
    # Daily light integral above canopy (DPIabovecanopy)
    DPIabovecanopy_matrix <- matrix.para(data = x, para = x$MeanDPI_plant*DPI_multiplier, para_length = x$LeafNo)

    # Matric for diurnal PPFD above the canopy (PPFDabovecanopy), Equation P12
    DiurnalPPFD_matrix <- (DPIabovecanopy_matrix*pi/(2*DL))*((10^6)/3600)*cos(pi*(timepoints(DL, Time_step_photo) - 12)/DL)   
    # Check simulated light integral
    x$DPI_actual <- colSums(DiurnalPPFD_matrix*Time_step_photo*3600*10^-6)   
 
    # Calculate capacities by photo nitrogen, Equations M1a, M1b, M1c
    x$Vcmax <- x$NV_mmol_m2_adj*ChiV
    x$Jmax <- x$NJ_mmol_m2_adj*ChiJ
    x$Chl <- x$NC_mmol_m2_adj*ChiC + x$NJ_mmol_m2_adj*ChiCJ
    x$alpha <- x$Chl/(x$Chl + 0.076)
  
    # Calculate mesophyll conductance (gm) using empirical relationship
    x$gm <- gm_empirical(data = x)
    x$g0m <- g0 + x$gm
   
    # Prepare matrices of each leaf at each time point for photosynthesis simulation
    Vcmax_matrix <- matrix.para(data = x, para = x$Vcmax)
    Jmax_matrix <- matrix.para(data = x, para = x$Jmax)
    leaf_absorption_matrix <- matrix.para(data = x, para = x$alpha)
    Rd_matrix <- matrix.para(data = x, para = x$Rd)
    gm_matrix <- matrix.para(data = x, para = x$gm)
    g0m_matrix <- matrix.para(data = x, para = x$g0m)
    Light_extinction_matrix <- matrix.para(data = x, para = x$Light_extinction_last_d)
    LA_m2_matrix <- matrix.para(data = x, para = x$LA_cm2_mean*10^-4)
    Light_intensity_leaf_matrix <- Light_extinction_matrix*DiurnalPPFD_matrix
    
    # Dependency of carboxylation rate on PPFD, Equation P13
    R_activation_rate <- ((31 + 69/(1 + exp(-0.009*(Light_intensity_leaf_matrix - 500))))*10^-2) 
    Vc_matrix <- Vcmax_matrix*R_activation_rate   
   
    # Dependency of electron transport rate on PPFD, Equation P4
    J_matrix <- (phi*leaf_absorption_matrix*Light_intensity_leaf_matrix + Jmax_matrix
                 - ((phi*leaf_absorption_matrix*Light_intensity_leaf_matrix + Jmax_matrix)^2 
                    - 4*theta*Jmax_matrix*phi*leaf_absorption_matrix*Light_intensity_leaf_matrix)^0.5)/(2*theta)   

    # leaf Photosynthesis analytical solution with daily light distribution, designate spaces for big matrices
    M00 <- Vcmax_matrix*0    
    Av_matrix1 <- M00        
    Aj_matrix1 <- M00       
    An_matrix1 <- M00       
    Av_matrix <- M00        
    Aj_matrix <- M00        
    An_matrix <- M00         
    
    # First calculation of leaf photosynthesis
    for (i in 1:length(timepoints(DL, Time_step_photo))) {
      for (j in 1:length(x$Vcmax)) {
      
        ## i=1
        ## j=2
      
        Rdij <- Rd_matrix[i,j]
        Vcij <- Vc_matrix[i,j]
        Jij <- J_matrix[i,j]
        g0mij <- g0m_matrix[i,j]
        gmij <- gm_matrix[i,j]  
        
        Acij <- c(Re(c(polyroot(c(g0*gmij*p2*Rdij - g0*gmij*p1*Vcij,
                                  g0*gmij*p2 + gsc_b*gmij*p2*Rdij - Rdij*g0mij - gmij*p1*Vcij*gsc_b + g0mij*Vcij,
                                  gsc_b*gmij*p2 - g0mij - Rdij*gsc_b + gsc_b*Vcij,
                                  -gsc_b)))[2]))
        
        Ajij <- c(Re(c(polyroot(c(4*g0*gmij*k2*Rdij - g0*gmij*k1*Jij,
                                  4*g0*gmij*k2 + 4*gsc_b*gmij*k2*Rdij - 4*Rdij*g0mij - gmij*k1*Jij*gsc_b + g0mij*Jij,
                                  4*gsc_b*gmij*k2 - 4*g0mij - 4*Rdij*gsc_b + gsc_b*Jij,
                                 -4*gsc_b)))[2]))
      
        Anij <- c(min(Acij, Ajij))
        Av_matrix1[i,j] <- Acij
        Aj_matrix1[i,j] <- Ajij 
        An_matrix1[i,j] <- Anij
      }
    }
  
    # Calculate and correct stomatal conductance (gs) and chloroplastic CO2 (Cc)
    G_matrix <- g0 + gsc_b*An_matrix1
    G_matrix[G_matrix <= 0] <- g0        # Replace gs <= 0 with g0
    Cc_matrix <- Ca - An_matrix1*((G_matrix + gm_matrix)/(G_matrix*gm_matrix))
    
    # Re-calculate photosynthesis using the corrected gs and Cc
    for (i in 1:length(timepoints(DL, Time_step_photo))) {
      for (j in 1:length(x$Vcmax)) {
        
        ## i=1
        ## j=2
        
        Rdij <- Rd_matrix[i,j]
        Vcij <- Vc_matrix[i,j]
        Jij <- J_matrix[i,j]
        Ccij <- Cc_matrix[i,j]  
        
        Acij <- Vcij*(Ccij - gama_str)/(Ccij + Km) - Rdij
        Ajij <- Jij*(Ccij - gama_str)/(4*Ccij + 8*gama_str) - Rdij
        Anij <- c(min(Acij,Ajij))
        
        Av_matrix[i,j] <- Acij
        Aj_matrix[i,j] <- Ajij 
        An_matrix[i,j] <- Anij
      }
    }
  
    # Calculate daily leaf carbon assimilation (DLA) and canopy carbon assimilation (DCA)
    An_leaf_mol_step <- (An_matrix*LA_m2_matrix)*3600*Time_step_photo*10^-6   # mol per time step per leaf
    x$DLA <- colSums(An_leaf_mol_step)                                        # daily leaf assimilation, mol d-1 per leaf
    x$DCA <- sum(x$DLA)
    
    # Calculate daily leaf light interception (DLI) and canopy leaf light interception (DCI)
    x$DLI_m2 <- x$MeanDPI_plant*x$Light_extinction_last_d                 # mol d-1 m-2
    x$DLI_leaf <- x$DLI_m2*(x$LA_cm2_mean*10^-4)                          # mol d-1 per leaf
    x$DCI <- sum(x$DLI_leaf)                                              # mol d-1 per plant
    x$LI_m2 <- x$DCI*PD                                                   # mol d-1 m-2 ground
    
    # Calculate daily instaneouse use effeciencies
    x$DLUEi <- x$DCA/x$DCI                                   # instantaneous light use efficiency
    x$PNUEi <- x$DCA/(x$Ncanopy_mmol_adj*10^-3)              # instantaneous photo nitrogen use efficiency
    
    x}))
  x1}))

  rownames(out) <- NULL
  return(out) 
}

####-----------------------------------------------------####
####-- Defined function to simulate Nph with various fd -####
####-----------------------------------------------------####
# 3.7.0. Function for simulating photosynthetic nitrogen (Nph) with various fd ----
simN_fd <- function(data, fd_list) {

  out <- do.call('rbind', lapply(split(data, 
                                       list(data$ExpID,
                                            data$MeasureDate,
                                            data$VarietyID,
                                            data$LightID,
                                            data$NitrogenID), 
                                       drop = TRUE), function(x){
                                                  
  ## x <- split(df_Structure_mean, list(df_Structure_mean$ExpID, df_Structure_mean$MeasureDate, df_Structure_mean$VarietyID, df_Structure_mean$LightID, df_Structure_mean$NitrogenID), drop = TRUE)[[1]]
                                             
  ## x <- split(data, list(data$ExpID, data$MeasureDate, data$VarietyID, data$LightID, data$NitrogenID), drop = TRUE)[[1]]
  
  E <<- as.character(x$ExpID[1])                
  MeasureDate <<- as.character(x$MeasureDate[1])
  V <<- as.character(x$VarietyID[1])
  L <<- as.character(x$LightID[1])
  N <<- as.character(x$NitrogenID[1])   
  
  # Select environmental input for a experiment
  df_Env_Exp <- subset(df_Environment, ExpID == E)

  # Calculate phyllochron and leaf age, Equations G1, G2, G3
  x <- leaf_age(x, df_Env_Exp)
  
  # Fit LAI to leaf age, Equation G4  
  x <- lai_age(x, PD)
  
  # Fit EA to leaf age, Equation G5  
  x <- ea_age(x)
  
  x <- do.call('rbind', lapply(split(x, x$LeafNo, drop = TRUE), function(x1){
    ## x1 <- split(x, x$LeafNo, drop=TRUE)[[2]]
    
    cat("[", MeasureDate, E, V, L, N, "]", "Leaf No:", x1$LeafNo[1], "\n")
    
    LID <- paste0("DPI_L_", L)
    NID <- paste0("NLV_", N)
    GDDID <- paste0("GDD_L_", L)
    
    # Calculate date of appearance for each leaf
    x1$App_Date <- min(subset(df_Env_Exp, df_Env_Exp[GDDID] >= x1$App_GDD)$Date)
    
    df_Env_Exp$DPI <- df_Env_Exp[,LID]
    df_Env_Exp$NS <- df_Env_Exp[,NID]
   
    # Select environment of a plant from transplating to measurement
    Env_input_plant <- subset(df_Env_Exp, Date <= MeasureDate)
    
    # Select environment of a leaf from its appearance to measurement
    Env_input_leaf <- subset(Env_input_plant, Date >= x1$App_Date)
    Env_input_leaf$GDD_leaf <- Env_input_leaf[,GDDID] - x1$App_GDD
    
    x1$MeanDPI_exp <- mean(df_Env_Exp$DPI)           # Mean incoming light integral during the whole experimental period
    x1$MeanDPI_plant <- mean(Env_input_plant$DPI)    # Mean incoming light integral during plant growth
    x1$MeanDPI_leaf <- mean(Env_input_leaf$DPI)      # Mean incoming light integral during leaf growth
    x1$MeanNLV_exp <- mean(df_Env_Exp$NS)            # Mean N availablility during the whole experimental period
    x1$MeanNLV_plant <- mean(Env_input_plant$NS)     # Mean N availablility during plant growth
    x1$MeanNLV_leaf <- mean(Env_input_leaf$NS)       # Mean N availablility during leaf growth 
    
    # Simulate light extinction using beer's law and elevation angle correction
    Env_input_leaf$Light_extinction <- light_extinction(k = k, t = Env_input_leaf$GDD_leaf, 
                                                        lai_max = x1$lai_max, lai_t0 = x1$lai_t0, lai_scal = x1$lai_scal,
                                                        ea_min = x1$ea_min, ea_t0 = x1$ea_t0, ea_scal = x1$ea_scal)
    
    Env_input_leaf$DPI_leaf <- Env_input_leaf$DPI*Env_input_leaf$Light_extinction
    
    # Respiration rate, simulated using the last 4 days of light at the leaf
    x1$DPI_4d <- mean(tail(Env_input_leaf$DPI_leaf, n = 4))  
    x1$Rd <- Rd_empirical(data = x1, environment_data = Env_input_leaf)
    
    # Light extinction on the last day for later daily canopy assimilation simulation
    x1$Light_extinction_last_d <- exp(-k*x1$LAI)*cos(x1$EA*pi/180)
    
    # Protein turnover model
    Model2 <- function(t, N_poolX, parms) {
      with(as.list(parms),{
        
        DPI <- head(Env_input_leaf$DPI_leaf[Env_input_leaf$GDD_leaf >= t], n = 1)
        NS <- head(Env_input_leaf$NS[Env_input_leaf$GDD_leaf >= t], n = 1)
        
        # Equation M7
        S_max <- (Smm*kI*DPI)/(Smm + kI*DPI)*(NS/(kN + NS))   
        
        # Equation S1
        Synthesis <- 2*S_max/(1 + exp(t*td*fd))
        
        # Equation M6
        Degradation <- Dr*N_poolX
        
        # Equation M4  
        dN_poolX <- Synthesis - Degradation
        
        list(dN_poolX)
        # return(c(N_poolX = list(dN_poolX), DPI = list(DPI), NS = list(NS), dN_poolX = list(dN_poolX)))
      }) 
    } 
    
    
    x2 <- do.call('rbind', lapply(fd_list, function(fd){
      ## fd <- fd_list[[1]]
      
    x2 <- do.call('rbind', lapply(list("NV", "NJ", "NC"), function(nx){        
      ## nx = "NV" 

      # Time points for simulation
      t = seq(0, x1$LeafAge_oCd, Time_step_N)   #simulation step: 0.1 oCd, until the growing oCd (TS) of the leaf
      # start simulating functional nitrogen over time until TS
      sim <- as.data.frame(lsoda(y = c(N_poolX = NX0), times = t, func = Model2, 
                                    parms = c(Dr = subset(df_para, NX == nx & VarietyID == V)$Dr,
                                              td = subset(df_para, NX == nx & VarietyID == V)$td,
                                              Smm = subset(df_para, NX == nx & VarietyID == V)$Smm,
                                              kI = subset(df_para, NX == nx & VarietyID == V)$kI,
                                              kN = subset(df_para, NX == nx & VarietyID == V)$kN,
                                              fd = fd)))
      
      x1$fd <- fd
      x1$NX <- nx
      x1$NX_mmol_m2 <- round(mean(tail(sim[,"N_poolX"], n = 12/Time_step_N)), 2) 
      x1$NX_mmol_leaf = x1$NX_mmol_m2*(x1$LA_cm2_mean*10^-4)
 
      cat("[", MeasureDate, E, V, L, N, "]", "fd:", fd, nx, x1$NX_mmol_m2, "\n")
      
      x1}))
    x2}))
    x2}))
    # x2}))
  x}))

rownames(out) <- NULL

# Adjust Nph in canopy to control value
out <- out %>%
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, fd, LeafNo) %>%
  mutate(Nph_mmol_m2 = sum(NX_mmol_m2),
         pX = NX_mmol_m2/Nph_mmol_m2) %>%
  
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, fd) %>%
  mutate(Ncanopy_mmol = sum(NX_mmol_leaf)) %>%
  
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID) %>%
  mutate(Ncanopy_mmol_control = unique(Ncanopy_mmol[fd == 1])) %>%
  
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo, NX) %>%
  mutate(
    # Equation S2
    Nph_mmol_m2_adj = Nph_mmol_m2*(Ncanopy_mmol_control/Ncanopy_mmol),
    # Equation S3
    pX_control = NX_mmol_m2[fd == 1]/Nph_mmol_m2[fd == 1],
    NX_mmol_m2_adj = Nph_mmol_m2_adj*pX_control) %>%
  
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, fd) %>%
  mutate(Ncanopy_mmol_adj = sum(NX_mmol_m2_adj*LA_cm2_mean*10^-4)) %>%
  
  group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, fd, LeafNo) %>%
  mutate(
    NV_mmol_m2 = NX_mmol_m2[NX == "NV"],
    NJ_mmol_m2 = NX_mmol_m2[NX == "NJ"],
    NC_mmol_m2 = NX_mmol_m2[NX == "NC"],
    NV_mmol_m2_adj = NX_mmol_m2_adj[NX == "NV"],
    NJ_mmol_m2_adj = NX_mmol_m2_adj[NX == "NJ"],
    NC_mmol_m2_adj = NX_mmol_m2_adj[NX == "NC"],
    pV = pX[NX == "NV"],
    pJ = pX[NX == "NJ"],
    pC = pX[NX == "NC"],
    pV_control = pX_control[NX == "NV"],
    pJ_control = pX_control[NX == "NJ"],
    pC_control = pX_control[NX == "NC"])

out <- as.data.frame(unique(subset(out, select = -c(NX, NX_mmol_m2, NX_mmol_m2_adj, NX_mmol_leaf, pX, pX_control))))

return(out)
}


####-----------------------------------------------------####
####----- IN SILICO TEST simulation with various fd -----####
####-----------------------------------------------------####
# 3.7.1. Simulate Nph with various fd --------------------------------------------
Simulation_fd1 <- simN_fd(data = df_Structure_mean, 
                          fd_list = list(0.5, 1, 2, 3, 4, 5))

# 3.7.2. Simulate DCA with various fd --------------------------------------------

Simulation_result <- simDCA(data = Simulation_fd1, 
                         list = list(Simulation_fd1$ExpID,
                                     Simulation_fd1$MeasureDate,
                                     Simulation_fd1$VarietyID,
                                     Simulation_fd1$LightID,
                                     Simulation_fd1$NitrogenID,
                                     Simulation_fd1$fd),
                         DPI_multiplier_list = c(0.5, 1, 2))


# 3.7.3. Output DCA comparison with various fd to file ----------------------------

# Calculate DCA change between fd scenarios
Simulation_result1 <- do.call('rbind', lapply(split(Simulation_result, 
                                                    list(Simulation_result$ExpID,
                                                         Simulation_result$MeasureDate,
                                                         Simulation_result$VarietyID,
                                                         Simulation_result$LightID,
                                                         Simulation_result$NitrogenID,
                                                         Simulation_result$DPI_multiplier), 
                                                    drop = TRUE), function(x){
                                                      
  # x <- split(Simulation_result, list(Simulation_result$ExpID, Simulation_result$MeasureDate, Simulation_result$VarietyID, Simulation_result$LightID, Simulation_result$NitrogenID, Simulation_result$DPI_multiplier), drop=TRUE)[[1]]
                     
  x$DCA_perct_change <- round((x$DCA-x$DCA[x$fd == 1][1])/abs(x$DCA[x$fd == 1][1])*100, 2)
  x$PNUEi_perct_change <- round((x$PNUEi-x$PNUEi[x$fd == 1][1])/abs(x$PNUEi[x$fd == 1][1])*100, 2)
  rownames(x) <- NULL
  x}))

Test_fd_result <- unique(subset(Simulation_result1, select = c(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, DPI_multiplier, MeanDPI_exp, MeanDPI_plant, MeanNLV_exp, MeanNLV_plant, fd, Ncanopy_mmol, Ncanopy_mmol_control, Ncanopy_mmol_adj, DCI, LI_m2, DCA, DCA_perct_change, PNUEi_perct_change)))

# Output result to file
write.xlsx(Test_fd_result, "Test_fd_result.xlsx", row.names = FALSE, col.names = TRUE)


# 3.7.4. Plot DCA comparison between various fd ----------------------------------

## Test_fd_result <- read.xlsx("Test_fd_result.xlsx", sheetName = "Sheet1")

do.call('rbind', lapply(split(Test_fd_result, list(Test_fd_result$ExpID,
                                                   Test_fd_result$MeasureDate,
                                                   Test_fd_result$VarietyID,
                                                   Test_fd_result$DPI_multiplier), 
                              drop = TRUE), function(x){
                                                      
  # x <- split(Test_fd_result, list(Test_fd_result$ExpID, Test_fd_result$MeasureDate, Test_fd_result$VarietyID, Test_fd_result$DPI_multiplier), drop=TRUE)[[1]]
 
  p <- ggplot(x, aes(x = fd, y = DCA_perct_change, color = LightID, shape = NitrogenID, group = interaction(LightID, NitrogenID))) +
       ylim(-20, 20) +
       geom_point(size = 5) +
       geom_line(size = 0.5) +
       geom_hline(yintercept = 0, linetype = "dashed") +
       ggtitle(paste(x$ExpID[1], x$VarietyID[1], x$MeasureDate[1], "DPI multiplier =", x$DPI_multiplier[1])) +
               labs(x = "fd", y = "DCA change by fd (%)") +
                    theme_bw() +
                    theme(
                      aspect.ratio = 1/1,
                      plot.title = element_text(color = "#993333", size = 14, face = "bold.italic"),
                      axis.title.x = element_text(color = "black", size = 20, face = "bold"),
                      axis.title.y = element_text(color = "black", size = 20, face = "bold"),
                      axis.text.x = element_text(face = "plain", color = "black", size = 20),
                      axis.text.y = element_text(face = "plain", color = "black", size = 20))
  print(p)          
}))


# 3.7.5. Output detailed results with various fd to file ----------------------------
write.xlsx(Simulation_result, "Test_fd_result_detailed.xlsx", row.names = FALSE, col.names = TRUE)

# 3.7.6. Output control nitrogen partitioning in the canopy ----------------------
N_partitioning_control <- subset(Simulation_fd1, fd == 1)

write.table(N_partitioning_control, file = "N_partitioning_control.csv", sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, quote = TRUE)

DCA_control <- subset(subset(Test_fd_result, select = c("ExpID", "MeasureDate","VarietyID", "LightID", "NitrogenID", "DPI_multiplier", "fd", "DCA")), fd == 1)

write.table(DCA_control, file = "DCA_control.csv", sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, quote = TRUE)

####-----------------------------------------------------####
####- Defined functions to simulate Nph with various fp -####
####-----------------------------------------------------####
# 3.8.0. Function for simulating photosynthetic nitrogen (Nph) with various fp ----

simN_fp <- function(data,
                    list,
                    fp_list,
                    Input_N_partitioning_control) {
  
  out <- do.call('rbind', lapply(split(data, list, drop = TRUE), function(x){
                                         
  ## x <- split(df_Structure_mean, list(df_Structure_mean$ExpID, df_Structure_mean$MeasureDate, df_Structure_mean$VarietyID, df_Structure_mean$LightID, df_Structure_mean$NitrogenID), drop = TRUE)[[2]]
                                         
  ## x <- split(data, list(data$ExpID, data$MeasureDate, data$VarietyID, data$LightID, data$NitrogenID), drop = TRUE)[[1]]
                                         
  E <<- as.character(x$ExpID[1])                
  MeasureDate <<- as.character(x$MeasureDate[1])
  V <<- as.character(x$VarietyID[1])
  L <<- as.character(x$LightID[1])
  N <<- as.character(x$NitrogenID[1])   
          
  # Input control Ncanopy
  Input_N_partitioning_control_x <- subset(Input_N_partitioning_control, 
                                           ExpID == E &
                                           VarietyID == V &
                                           LightID == L &
                                           NitrogenID == N)
  Input_N_partitioning_control_x <- Input_N_partitioning_control_x[Input_N_partitioning_control_x$MeasureDate == MeasureDate,]
  x$Ncanopy_mmol_control <- Input_N_partitioning_control_x$Ncanopy_mmol_control[1]  
                                         
  # Select environmental input for a experiment
  df_Env_Exp <- subset(df_Environment, ExpID == E)
                                         
  # Calculate phyllochron and leaf age, Equations G1, G2, G3
  x <- leaf_age(x, df_Env_Exp)
                   
  # Fit LAI to leaf age, Equation G4  
  x <- lai_age(x, PD)
           
  # Fit EA to leaf age, Equation G5  
  x <- ea_age(x)
                                    
  x <- do.call('rbind', lapply(split(x, x$LeafNo, drop = TRUE), function(x1){
    ## x1 <- split(x, x$LeafNo, drop=TRUE)[[2]]
                                           
    ## cat("[", MeasureDate, E, V, L, N, "]", "Leaf No:", x1$LeafNo[1], "\n")
                                           
    LID <- paste0("DPI_L_", L)
    NID <- paste0("NLV_", N)
    GDDID <- paste0("GDD_L_", L)
    
    # Calculate date of appearance for each leaf
    x1$App_Date <- min(subset(df_Env_Exp, df_Env_Exp[GDDID] >= x1$App_GDD)$Date)
    
      df_Env_Exp$DPI <- df_Env_Exp[,LID]
      df_Env_Exp$NS <- df_Env_Exp[,NID]
      
      # Select environment of a plant from transplating to measurement
      Env_input_plant <- subset(df_Env_Exp, Date <= MeasureDate)
      
      # Select environment of a leaf from its appearance to measurement
      Env_input_leaf <- subset(Env_input_plant, Date >= x1$App_Date)
      Env_input_leaf$GDD_leaf <- Env_input_leaf[,GDDID] - x1$App_GDD
      
      x1$MeanDPI_exp <- mean(df_Env_Exp$DPI)           # Mean incoming light integral during the whole experimental period
      x1$MeanDPI_plant <- mean(Env_input_plant$DPI)    # Mean incoming light integral during plant growth
      x1$MeanDPI_leaf <- mean(Env_input_leaf$DPI)      # Mean incoming light integral during leaf growth
      x1$MeanNLV_exp <- mean(df_Env_Exp$NS)            # Mean N availablility during the whole experimental period
      x1$MeanNLV_plant <- mean(Env_input_plant$NS)     # Mean N availablility during plant growth
      x1$MeanNLV_leaf <- mean(Env_input_leaf$NS)       # Mean N availablility during leaf growth 
      
      # Simulate light extinction using beer's law and elevation angle correction
      Env_input_leaf$Light_extinction <- light_extinction(k = k, t = Env_input_leaf$GDD_leaf,
                                                          lai_max = x1$lai_max, lai_t0 = x1$lai_t0, lai_scal = x1$lai_scal,
                                                          ea_min = x1$ea_min, ea_t0 = x1$ea_t0, ea_scal = x1$ea_scal)
      
      Env_input_leaf$DPI_leaf <- Env_input_leaf$DPI*Env_input_leaf$Light_extinction
      
      # Respiration rate, simulated using the last 4 days of light at the leaf
      x1$DPI_4d <- mean(tail(Env_input_leaf$DPI_leaf, n = 4))  
      x1$Rd <- Rd_empirical(data = x1, environment_data = Env_input_leaf)
      
      # Light extinction on the last day for later daily canopy assimilation simulation
      x1$Light_extinction_last_d <- exp(-k*x1$LAI)*cos(x1$EA*pi/180)
      
      # Protein turnover model
      Model2 <- function(t, N_poolX, parms) {
        with(as.list(parms), {
          DPI <- head(Env_input_leaf$DPI_leaf[Env_input_leaf$GDD_leaf >= t], n = 1)
          NS <- head(Env_input_leaf$NS[Env_input_leaf$GDD_leaf >= t], n = 1)
          
          # Equation M7
          S_max <- (Smm*fp*kI*DPI)/(Smm*fp + kI*DPI)*(NS/(kN + NS)) 
          
          # Equation S1
          Synthesis <- 2*S_max/(1 + exp(t*td))
          
          # Equation M6
          Degradation <- Dr*N_poolX
          
          # Equation M4  
          dN_poolX <- Synthesis - Degradation
          
          list(dN_poolX)
          #return(c(NX = list(dNX), DPI = list(DPI), NS = list(NS), dNX = list(dNX)))
        }) 
        } 
      
      x2 <- do.call('rbind', lapply(list("NV", "NJ", "NC"), function(nx){   
        ## nx = "NV" 
        
        # Time points for simulation
        t = seq(0, x1$LeafAge_oCd, Time_step_N)   #simulation step: 0.1 oCd, until the growing oCd (TS) of the leaf
        # start simulating functional nitrogen over time until TS
        sim <- as.data.frame(lsoda(y = c(N_poolX = NX0), times = t, func = Model2, 
                                   parms = c(Dr = df_para$Dr[df_para$NX == nx & df_para$VarietyID == V],
                                             td = df_para$td[df_para$NX == nx & df_para$VarietyID == V],
                                             Smm = df_para$Smm[df_para$NX == nx & df_para$VarietyID == V],
                                             kI = df_para$kI[df_para$NX == nx & df_para$VarietyID == V],
                                             kN = df_para$kN[df_para$NX == nx & df_para$VarietyID == V],
                                             fp = fp_list[[nx]])))
        
        x1$fp <- fp_list[[nx]]
        x1$NX <- nx
        x1$NX_mmol_m2 <- round(mean(tail(sim[,"N_poolX"], n = 12/Time_step_N)), 2) 
        x1$NX_mmol_leaf = x1$NX_mmol_m2*(x1$LA_cm2_mean*10^-4)
        ## cat("[", MeasureDate, E, V, L, N, "]", "fp:", fp_list[[nx]], nx, x1$NX_mmol_m2, "\n")
        
        x1}))
      x2}))
  x}))
  
  rownames(out) <- NULL
  
  # Adjust Nph in canopy to control value
  out <- out %>%
    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo) %>%
    mutate(Nph_mmol_m2 = sum(NX_mmol_m2),
           pX = NX_mmol_m2/Nph_mmol_m2) %>%
  
    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID) %>%
    mutate(Ncanopy_mmol = sum(NX_mmol_leaf),
           Nph_mmol_m2_adj = Nph_mmol_m2*(Ncanopy_mmol_control/Ncanopy_mmol),
           NX_mmol_m2_adj = pX*Nph_mmol_m2_adj,
           Ncanopy_mmol_adj = sum(NX_mmol_m2_adj*LA_cm2_mean*10^-4)) %>%

    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo) %>%
    mutate(
      fpV = fp[NX == "NV"],
      fpJ = fp[NX == "NJ"],
      fpC = fp[NX == "NC"],
      NV_mmol_m2 = NX_mmol_m2[NX == "NV"],
      NJ_mmol_m2 = NX_mmol_m2[NX == "NJ"],
      NC_mmol_m2 = NX_mmol_m2[NX == "NC"],
      NV_mmol_m2_adj = NX_mmol_m2_adj[NX == "NV"],
      NJ_mmol_m2_adj = NX_mmol_m2_adj[NX == "NJ"],
      NC_mmol_m2_adj = NX_mmol_m2_adj[NX == "NC"],
      pV = pX[NX == "NV"],
      pJ = pX[NX == "NJ"],
      pC = pX[NX == "NC"])
  
  out <- as.data.frame(unique(subset(out, select = -c(NX, NX_mmol_m2, NX_mmol_m2_adj, NX_mmol_leaf, pX, fp))))
  
  return(out)
}


# 3.8.1. Function for maximizing DCA with given fp -------------------------------
opt_fp <- function(p){
  
  fp_list <- list(NV = p[1], NJ = p[2], NC = p[3])
  
  out <- simN_fp(data = x,
                 list = x$ExpID,
                 fp_list,
                 Input_N_partitioning_control = Input_N_partitioning_control) 
  
  DCA <- simDCA(out, list = list(out$ExpID), DPI_multiplier_list = list(DPI_multiplier_test))$DCA[1]
    
  cat("[", MeasureDate, E, V, L, N, "]", "DCA:", round(DCA, 4), "DPI multiplier:", DPI_multiplier_test, " fpV:", round(p[1],3), " fpJ:", round(p[2],3), " fpC:", round(p[3],3), "\n")
    
    return(-DCA)
    }


# 3.8.2. Function for optimizing Nph partitioning with various fp ----------------
opt_fp_test <- function(data, lower, upper, NP, itermax, list){
  
  out3 <- do.call('rbind', lapply(split(data, list, drop = TRUE), function(x){
    
    ## x <- split(df_Structure_mean, list(df_Structure_mean$ExpID, df_Structure_mean$MeasureDate, df_Structure_mean$VarietyID, df_Structure_mean$LightID, df_Structure_mean$NitrogenID), drop = TRUE)[[1]]
    
    ## x <- split(data, list(data$ExpID, data$MeasureDate, data$VarietyID, data$LightID, data$NitrogenID), drop = TRUE)[[1]]
    
    x <<- x
    
    E <<- as.character(x$ExpID[1])                
    MeasureDate <<- as.character(x$MeasureDate[1])
    V <<- as.character(x$VarietyID[1])
    L <<- as.character(x$LightID[1])
    N <<- as.character(x$NitrogenID[1])   
      
      ## cat("[", MeasureDate, E, V, L, N, "]", "\n")
      
      # Start optimization procedure
      out1 <- DEoptim(opt_fp, lower = lower, upper = upper, control = DEoptim.control(NP = NP, itermax = itermax))
      
      out2 <- data.frame(ExpID = E, MeasureDate = MeasureDate, VarietyID = V, 
                LightID = L, NitrogenID = N, 
                DPI_multiplier = DPI_multiplier_test,
                fpV = out1$optim$bestmem[1],
                fpJ = out1$optim$bestmem[2],
                fpC = out1$optim$bestmem[3],
                DCA_opt = -out1$optim$bestval, 
                iteration = out1$optim$iter)
      out2}))

  return(out3)}


####-----------------------------------------------------####
####----- IN SILICO TEST simulation with various fp -----####
####-----------------------------------------------------####
# 3.8.3. Optimize Nph partitioning with various fp -------------------------------

# Input control Nph partitioning 
Input_N_partitioning_control <- read.csv("N_partitioning_control.csv", sep = ",", dec = ".", header = TRUE)

# Set seed
set.seed(345)

# User defined DPI multiplier
DPI_multiplier_list = c(0.5, 1, 2)

# Optimize Nph partitioning with various fp
Test_fp_result <- do.call('rbind', lapply(DPI_multiplier_list, function(DPI_multiplier_test){        
  ## DPI_multiplier_test = 1 
  
  DPI_multiplier_test <<- DPI_multiplier_test
  
  out <- opt_fp_test(data = df_Structure_mean,
                    lower = c(0.2, 0.2, 0.2), 
                    upper = c(2, 2, 2),
                    NP = 30,
                    itermax = 20, 
                    list = list(df_Structure_mean$ExpID, 
                                df_Structure_mean$MeasureDate, 
                                df_Structure_mean$VarietyID, 
                                df_Structure_mean$LightID, 
                                df_Structure_mean$NitrogenID))

  out}))


# 3.8.4. Output DCA comparison with various fp to file ----------------------------

# Input control DCA
DCA_control <- read.csv("DCA_control.csv", sep = ",", dec = ".", header = TRUE)

Test_fp_result1 <- merge(Test_fp_result, DCA_control, 
                         by = c("ExpID", "MeasureDate","VarietyID", "LightID", "NitrogenID", "DPI_multiplier"))

# Calculate DCA change by optimal paratitioning
Test_fp_result1$DCA_perct_change <- round((Test_fp_result1$DCA_opt-Test_fp_result1$DCA)/abs(Test_fp_result1$DCA)*100, 2)

# Output results to file
write.xlsx(Test_fp_result1, "Test_fp_result.xlsx", row.names = FALSE, col.names = TRUE)


# 3.8.5. Plot DCA comparison with optimal fp -------------------------------------
do.call('rbind', lapply(split(Test_fp_result1, list(Test_fp_result1$ExpID,
                                                    Test_fp_result1$MeasureDate,
                                                    Test_fp_result1$VarietyID), 
                              drop = TRUE), function(x){
                                
  # x <- split(Test_fp_result1, list(Test_fp_result1$ExpID, Test_fp_result1$MeasureDate, Test_fp_result1$VarietyID), drop=TRUE)[[1]]
                                
  p <- ggplot(x, aes(x = DPI_multiplier, y = DCA_perct_change, color = LightID, shape = NitrogenID, group = interaction(LightID, NitrogenID))) +
       #ylim(-20, 20) +
       geom_point(size = 5) +
       geom_line(size = 0.5) +
       geom_hline(yintercept = 0, linetype = "dashed") +
       ggtitle(paste(x$ExpID[1], x$VarietyID[1], x$MeasureDate[1])) +
       labs(x = "DPI multiplier", y = "DCA change by optimal fp (%)") +
       theme_bw() +
       theme(
       aspect.ratio = 1/1,
       plot.title = element_text(color = "#993333", size = 14, face = "bold.italic"),
       axis.title.x = element_text(color = "black", size = 20, face = "bold"),
       axis.title.y = element_text(color = "black", size = 20, face = "bold"),
       axis.text.x = element_text(face = "plain", color = "black", size = 20),
       axis.text.y = element_text(face = "plain", color = "black", size = 20))
  print(p)
  
}))
# 3.8.6. Output detailed results with optimal fp to file ----------------

# Test_fp_result <- read.xlsx("Test_fp_result.xlsx", sheetName = "Sheet1")

# User defined DPI multiplier
DPI_multiplier_list = c(0.5, 1, 2)

# Input control Nph partitioning 
Input_N_partitioning_control <- read.csv("N_partitioning_control.csv", sep = ",", dec = ".", header = TRUE)
Input_N_partitioning_control_pX <- subset(Input_N_partitioning_control, select = c("ExpID", "MeasureDate", "VarietyID", "LightID", "NitrogenID", "LeafNo", "pV_control", "pJ_control", "pC_control"))

# Simulate detailed result with optimal fp
Test_fp_result_detailed <- do.call('rbind', lapply(split(df_Structure_mean, list(df_Structure_mean$ExpID, df_Structure_mean$MeasureDate, df_Structure_mean$VarietyID, df_Structure_mean$LightID, df_Structure_mean$NitrogenID), drop = TRUE), function(x){
  
  # x <- split(df_Structure_mean, list(df_Structure_mean$ExpID, df_Structure_mean$MeasureDate, df_Structure_mean$VarietyID, df_Structure_mean$LightID, df_Structure_mean$NitrogenID), drop=TRUE)[[1]]
  
  E <<- as.character(x$ExpID[1])                                                      
  MeasureDate <<- as.character(x$MeasureDate[1])
  V <<- as.character(x$VarietyID[1])
  L <<- as.character(x$LightID[1])
  N <<- as.character(x$NitrogenID[1])
  
  out <- do.call('rbind', lapply(DPI_multiplier_list, function(DPI_multiplier_test){        
    ## DPI_multiplier_test = DPI_multiplier_list[[1]] 
    
    x$DPI_multiplier <- DPI_multiplier_test
    
    # Select optimal fp
    fp_list <- subset(Test_fp_result, ExpID == E & MeasureDate == MeasureDate & 
                        VarietyID == V & LightID == L & 
                        NitrogenID == N & DPI_multiplier == DPI_multiplier_test)
    
    fp_list <- list(NV = fp_list$fpV, NJ = fp_list$fpJ, NC = fp_list$fpC)
    
    # Select environmental input for a experiment
    df_Env_Exp <- subset(df_Environment, ExpID == E)
    
    # Calculate phyllochron and leaf age, Equations G1, G2, G3
    x <- leaf_age(x, df_Env_Exp)
    
    # Fit LAI to leaf age, Equation G4  
    x <- lai_age(x, PD)
    
    # Fit EA to leaf age, Equation G5  
    x <- ea_age(x)
    
    out1 <- do.call('rbind', lapply(split(x, x$LeafNo, drop = TRUE), function(x1){
      ## x1 <- split(x, x$LeafNo, drop=TRUE)[[2]]
      
      LID <- paste0("DPI_L_", L)
      NID <- paste0("NLV_", N)
      GDDID <- paste0("GDD_L_", L)
      
      # Calculate day of appearance for each leaf
      x1$App_Date <- min(subset(df_Env_Exp, df_Env_Exp[GDDID] >= x1$App_GDD)$Date)
      
      df_Env_Exp$DPI <- df_Env_Exp[,LID]
      df_Env_Exp$NS <- df_Env_Exp[,NID]
      
      # Select environment of a plant from transplating to measurement
      Env_input_plant <- subset(df_Env_Exp, Date <= MeasureDate)
      # Select environment of a leaf from its appearance to measurement
      Env_input_leaf <- subset(Env_input_plant, Date >= x1$App_Date)
      Env_input_leaf$GDD_leaf <- Env_input_leaf[,GDDID] - x1$App_GDD
      
      x1$MeanDPI_exp <- mean(df_Env_Exp$DPI)                 # Mean incoming light integral during the whole experimental period
      x1$MeanDPI_plant <- mean(Env_input_plant$DPI)          # Mean incoming light integral during plant growth
      x1$MeanDPI_leaf <- mean(Env_input_leaf$DPI)            # Mean incoming light integral during leaf growth
      x1$MeanNLV_exp <- mean(df_Env_Exp$NS)                  # Mean N availablility during the whole experimental period
      x1$MeanNLV_plant <- mean(Env_input_plant$NS)           # Mean N availablility during plant growth
      x1$MeanNLV_leaf <- mean(Env_input_leaf$NS)             # Mean N availablility during leaf growth 
      
      # Simulate light extinction using beer's law and elevation angle correction
      Env_input_leaf$Light_extinction <- light_extinction(k = k, t = Env_input_leaf$GDD_leaf, 
                                                          lai_max = x1$lai_max, lai_t0 = x1$lai_t0, lai_scal = x1$lai_scal,
                                                          ea_min = x1$ea_min, ea_t0 = x1$ea_t0, ea_scal = x1$ea_scal)
      
      Env_input_leaf$DPI_leaf <- Env_input_leaf$DPI*Env_input_leaf$Light_extinction
      
      # Respiration rate, simulated using the last 4 days of light at the leaf
      x1$DPI_4d <- mean(tail(Env_input_leaf$DPI_leaf, n = 4))                           
      x1$Rd <- Rd_empirical(data = x1, environment_data = Env_input_leaf)
      
      # Light extinction on the last day for later daily canopy assimilation simulation
      x1$Light_extinction_last_d <- exp(-k*x1$LAI)*cos(x1$EA*pi/180)
      
      # Protein turnover model
      Model2 <- function(t, N_poolX, parms) {
        with(as.list(parms),{
          
          DPI <- head(Env_input_leaf$DPI_leaf[Env_input_leaf$GDD_leaf >= t], n = 1)
          NS <- head(Env_input_leaf$NS[Env_input_leaf$GDD_leaf >= t], n = 1)
          
          # Equation S4
          S_max <- (Smm*fp*kI*DPI)/(Smm*fp + kI*DPI)*(NS/(kN + NS))   
          
          # Equation M5
          Synthesis <- 2*S_max/(1 + exp(t*td))
          
          # Equation M6
          Degradation <- Dr*N_poolX
          
          # Equation M4  
          dN_poolX <- Synthesis - Degradation
          
          list(dN_poolX)
          # return(c(NX = list(dNX), DPI = list(DPI), NS = list(NS), dNX = list(dNX)))
        }) 
      } 
      
      x2 <- do.call('rbind', lapply(list("NV", "NJ", "NC"), function(nx){        
        ## nx = "NV" 
        
        # Time points for simulation
        t = seq(0, x1$LeafAge_oCd, Time_step_N)   #simulation step: 0.1 oCd, until the growing oCd (TS) of the leaf
        # start simulating functional nitrogen over time until TS
        sim <- as.data.frame(lsoda(y = c(N_poolX = NX0), times = t, func = Model2, 
                                   parms = c(Dr = subset(df_para, NX == nx & VarietyID == V)$Dr,
                                             td = subset(df_para, NX == nx & VarietyID == V)$td,
                                             Smm = subset(df_para, NX == nx & VarietyID == V)$Smm,
                                             kI = subset(df_para, NX == nx & VarietyID == V)$kI,
                                             kN = subset(df_para, NX == nx & VarietyID == V)$kN,
                                             fp = fp_list[[nx]])))
        x1$NX <- nx
        x1$NX_mmol_m2 <- mean(tail(sim[,"N_poolX"], n = 12/Time_step_N)) 
        x1$NX_mmol_leaf = x1$NX_mmol_m2*(x1$LA_cm2_mean*10^-4)
        
        x1})) 
      x2}))
    out1}))
  
  # Input control Ncanopy
  Input_N_partitioning_control_x <- subset(Input_N_partitioning_control, 
                                           ExpID == E & VarietyID == V & LightID == L & NitrogenID == N)
  Input_N_partitioning_control_x <- Input_N_partitioning_control_x[Input_N_partitioning_control_x$MeasureDate == MeasureDate,]
  out$Ncanopy_mmol_control <- Input_N_partitioning_control_x$Ncanopy_mmol_control[1]  
  
  # Adjust Nph in canopy to control value
  out <- out %>%
    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo, DPI_multiplier) %>%
    mutate(Nph_mmol_m2 = sum(NX_mmol_m2),
           pX = NX_mmol_m2/Nph_mmol_m2)
  
  out <- out %>%
    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, DPI_multiplier) %>%
    mutate(Ncanopy_mmol = sum(NX_mmol_leaf),
           Nph_mmol_m2_adj = Nph_mmol_m2*(Ncanopy_mmol_control/Ncanopy_mmol),
           NX_mmol_m2_adj = pX*Nph_mmol_m2_adj,
           Ncanopy_mmol_adj = sum(NX_mmol_m2_adj*LA_cm2_mean*10^-4)
    )
  
  out <- out %>%
    group_by(ExpID, MeasureDate, VarietyID, LightID, NitrogenID, LeafNo, DPI_multiplier) %>%
    mutate(
      NV_mmol_m2 = NX_mmol_m2[NX == "NV"],
      NJ_mmol_m2 = NX_mmol_m2[NX == "NJ"],
      NC_mmol_m2 = NX_mmol_m2[NX == "NC"],
      NV_mmol_m2_adj = NX_mmol_m2_adj[NX == "NV"],
      NJ_mmol_m2_adj = NX_mmol_m2_adj[NX == "NJ"],
      NC_mmol_m2_adj = NX_mmol_m2_adj[NX == "NC"],
      pV_opt = pX[NX == "NV"],
      pJ_opt = pX[NX == "NJ"],
      pC_opt = pX[NX == "NC"])
  
  out <- as.data.frame(unique(subset(out, select = -c(NX, NX_mmol_m2, NX_mmol_m2_adj, NX_mmol_leaf, pX))))
}))

# Merge control values of pX to result
Test_fp_result_detailed1 <- do.call('rbind', lapply(DPI_multiplier_list, function(DPI_multiplier){
  # DPI_multiplier = DPI_multiplier_list[[2]]
  out <- Test_fp_result_detailed[Test_fp_result_detailed$DPI_multiplier == DPI_multiplier,]
  out <- merge(out, Input_N_partitioning_control_pX, 
               by = c("ExpID", "MeasureDate","VarietyID", "LightID", "NitrogenID", "LeafNo"))
  }))

# Output results to file
write.xlsx(Test_fp_result_detailed1, "Test_fp_result_detailed.xlsx", row.names = FALSE, col.names = TRUE)