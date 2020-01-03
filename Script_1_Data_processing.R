# 1.0.0. Set working directory to current path -----------------------------------
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(this.dir)


# 1.1.0. General constants --------------------------------------------------------

# Photosynthesis
beta = 0.5           # partitioning fraction of photons between photosystem II and I, unitless
theta = 0.7          # constant convexity factor describing the response of J to PPFD, unitless
phi = 0.425          # conversion efficiency of photons to J, umol e-/umol photon
gama_str = 43.02     # CO2 compensation point in the absence of mitochondrial respiration, umol CO2/mol
Kc = 404             # Michaelis-Menten constants of Rubisco for CO2, umol CO2/mol
Ko = 278             # Michaelis-Menten constants of Rubisco for O2, mmmol O2/mol
Oc = 210             # mole fraction of O2 at the site of carboxylation, mmol O2/mol
Km <- Kc*(1 + Oc/Ko) # Equation P8

# Photosynthetic nitrogen and function capacity 
ChiV = 4.49          # carboxylation capacity per unit Rubisco nitrogen, umol CO2/mmol N/s 
ChiJ = 9.48          # electron transport capacity per unit electron transport nitrogen, umol e-/mmol N/s
ChiC = 0.03384       # conversion coefficients for chlorophyll per light harvesting component nitrogen, mmol Chl/mmol N
ChiCJ = 4.64*10^-4   # conversion coefficients for chlorophyll per electron transport nitrogen, mmol Chl/mmol N

# Molar mass
MolWt_N = 14         # molar mass of nitrogen, g N/mol
MolWt_Chla = 893.51  # molar mass of chlorophyll a, g Chl a/mol
MolWt_Chlb = 907.49  # molar mass of chlorophyll b, g Chl b/mol

# Plant constant
Temp_base = 10       # base temperature for cucumber, oC


# 1.2.0. Load data ---------------------------------------------------------------
df_harvest <- read.csv("example_chamber_harvest_data.csv", sep = ";", dec = ".", header = TRUE)
df_gas_exchange <- read.csv("example_chamber_gas_exchange_data.csv", sep = ";", dec = ".", header = TRUE)    


# 1.2.1. Merge harvest and gas exchange data -------------------------------------
df <- merge(df_harvest, df_gas_exchange, by = "LeafID")
df$AppearanceDate <- as.Date(df$AppearanceDate, "%d.%m.%Y")
df$HarvestDate <- as.Date(df$HarvestDate, "%d.%m.%Y")


# 1.3.0. Data processing ---------------------------------------------------------

# Leaf mass per area (LMA), g/m2
df$LeafArea_m2 <-  df$LeafArea_cm2*10^-4                            
df$LMA_g_m2 <- df$DryMass_g/df$LeafArea_m2
  
# Convert unit of chlorophyll into mmol/m2
df$Chl_a_b_mmol_m2 <- ((df$Chl_a_mg_g/893.51) + (df$Chl_b_mg_g/907.49))*df$LMA_g_m2 

# Equation P2
df$abs <- (df$Chl_a_b_mmol_m2)/((df$Chl_a_b_mmol_m2) + 0.076)
  
# Equation P3
df$J <- df$abs*beta*df$PPFD*df$PhiPS2        

# Estimate maximum electron transport rate (Jmax) and daytime respiration rate (Rd)
df1 <- do.call('rbind', lapply(split(df, df$LeafID, drop = TRUE), function(x){
  ## x <- split(df, df$LeafID, drop = TRUE)[[1]]

  # Estimate Jmax
   if (any(is.na(x$J)) == TRUE) {
     x$Jmax <- NA
     x$Jmax_se <- NA
     x$Jmax_pv <- NA
   }else{
     # Equation P4
     Jmax_fit <- nls(J ~ (Jmax+phi*PPFD-((Jmax+phi*PPFD)^2-4*theta*Jmax*phi*PPFD)^0.5)/(2*theta), 
                  data = x, 
                  start = list(Jmax = max(x$J)))
     x$Jmax <- summary(Jmax_fit)$parameters[,1]     # Fitted value of Jmax
     x$Jmax_se <- summary(Jmax_fit)$parameters[,2]  # Standard error
     x$Jmax_pv <- summary(Jmax_fit)$parameters[,4]  # p value
   }
  
   # Estimate Rd
   x2 <- subset(x, x$PPFD >= 40 & x$PPFD <= 100)
   x$Rd <- -as.numeric(lm(x2$An ~ x2$PPFD)[[1]][1])   
   x$Rd <- ifelse(x$Rd >= 0, yes = x$Rd, no = 0)
  
  x})) 

# Equation P5
df1$gm_mol_m2_s <- df1$An/(df1$Ci - (gama_str*(df1$J + 8*df1$An + 8*df1$Rd)/(df1$J - 4*df1$An - 4*df1$Rd)))

# Equation P6
df1$Cc_umol_mol <- df1$Ci - df1$An/df1$gm_mol_m2_s
 
# Equation P7 
df1$Vc <- (df1$An + df1$Rd)*((df1$Cc_umol_mol + Km)/(df1$Cc_umol_mol - gama_str))


# 1.3.1. Finalizing data processing and estimate photosynthetic nitrogen ---------
Dataset <- subset(df1, PPFD >= 1295)  
Dataset$Vcmax <- Dataset$Vc

# Equation M1a, M1b and M1c
Dataset$NV_mmol_m2 = Dataset$Vcmax/ChiV                            
Dataset$NJ_mmol_m2 = Dataset$Jmax/ChiJ                            
Dataset$NC_mmol_m2 = (Dataset$Chl_a_b_mmol_m2 - Dataset$NJ_mmol_m2*ChiCJ)/ChiC  

# Equation M2
Dataset$Nph_mmol_m2 <- Dataset$NV_mmol_m2 + Dataset$NJ_mmol_m2 + Dataset$NC_mmol_m2

# Equation M3
Dataset$pV <- Dataset$NV_mmol_m2/Dataset$Nph_mmol_m2
Dataset$pJ <- Dataset$NJ_mmol_m2/Dataset$Nph_mmol_m2
Dataset$pC <- Dataset$NC_mmol_m2/Dataset$Nph_mmol_m2

# Calculate leaf age (oCd), Equation G4
Dataset$LeafAge_oCd <- as.numeric((Dataset$HarvestDate - Dataset$AppearanceDate))*(Dataset$MeanTemp_oC - Temp_base)


# 1.4.0. Output processed data to file -------------------------------------------
write.table(Dataset, "chamber_processed_data.csv", row.names = FALSE, sep = ";", dec = ".")
