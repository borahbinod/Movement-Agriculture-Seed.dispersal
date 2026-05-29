###############################################
## FUNCTION 1: Euclidean distance calculator ##
###############################################

# Calculates straight-line distance between two points/vectors
# Formula:
# sqrt((x1-x2)^2)

euc.dist <- function(x1, x2)
  
  # Difference between coordinates
  # Square them
  # Sum them
  # Take square root
  
  sqrt(sum((x1 - x2) ^ 2))



#########################################
## FUNCTION 2: Bootstrap mean sampler ##
#########################################

boots <- function(x){
  
  # Repeat 1000 bootstrap replicates
  for (i in 1:1000){
    
    # Randomly sample 3 observations WITH replacement
    # from vector x
    
    s1 <- sample(x,
                 size = 3,
                 replace = TRUE)
    
    # Calculate mean of sampled values
    ms <- mean(s1)
    
    # Store bootstrap mean
    s[i] <- ms
  }
  
  # These lines were probably used earlier
  # to summarize bootstrap distribution
  
  #m1 <- mean(s)   # mean of bootstrap means
  #k1 <- sd(s)     # sd of bootstrap means
  #r1 <- c(m1,k1)
  
  # Return vector of 1000 bootstrap means
  return(s)
}



##################################################
## FUNCTION 3: Estimate beta distribution params ##
##################################################

# Converts mean and variance into
# alpha and beta parameters of a beta distribution

# Useful because:
# rbeta() needs alpha and beta
# but empirical data usually provide mean and variance

estBetaParams <- function(mu, var) {
  
  # Estimate alpha parameter
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  
  # Estimate beta parameter
  beta <- alpha * (1 / mu - 1)
  
  # Store parameters together
  beta.params <- c(alpha, beta)
  
  return(beta.params)
}



###################################################
## FUNCTION 4: Estimate gamma distribution params ##
###################################################

# Converts mean and variance into
# shape and scale parameters of gamma distribution

# Used because rgamma() requires shape and scale

estGammaParams <- function(mu,var){
  
  # shape parameter
  shp <- mu^2 / var^2
  
  # scale parameter
  scl <- var^2 / mu
  
  gamma.params <- c(shp,scl)
  
  return(gamma.params)
}



###########################################
## LOAD LANDSCAPE / HABITAT RASTER DATA ##
###########################################

# Load artificial fragmented landscape raster

lscp <- raster("Data_derived/frag_landscape_0.1_50.tif")

# Plot raster
plot(lscp)



#########################################
## CONVERT RASTER TO DATAFRAME FORMAT ##
#########################################

# Convert raster cells to dataframe
# xy = TRUE keeps coordinates

df <- as.data.frame(lscp, xy = TRUE)

# Rename columns
# x = x coordinate
# y = y coordinate
# habitat = habitat type

colnames(df) <- c("x","y","habitat")

# Add unique ID to each raster cell
df$id <- 1:nrow(df)



##############################################
## LOAD STARTING LOCATIONS OF INDIVIDUALS ##
##############################################

# Starting locations for simulated birds

start_pops <- read.csv("Data_derived/start_pops.csv",
                       header = TRUE)

# Raster cell IDs corresponding to start locations
xid <- start_pops$id

# Force those cells to habitat = 1
# ensures birds start in habitat

df[xid,"habitat"] <- 1



#########################################
## CONVERT DATAFRAME BACK TO RASTER ##
#########################################

spg <- df

# Convert to spatial object
coordinates(spg) <- ~x+y

# Define as gridded spatial object
gridded(spg) <- TRUE

# Convert back to raster
lscp <- raster(spg)

# Plot modified raster
plot(lscp)

# Add starting locations on top
points(start_pops[,c("x","y")], pch = 19)



###################################################
## LOAD EMPIRICAL MOVEMENT / BEHAVIOR ESTIMATES ##
###################################################

# Input behavioral estimates

est.data <- read.csv("Data_derived/Estimates.csv",
                     header = TRUE)

# Keep only specialist species
est <- est.data %>%
  filter(Habitat %in% c("Specialist"))



######################################
## GUT RETENTION TIME (GRT) DATA ##
######################################

# GRT = time seeds remain inside bird gut

grt_data <- read.csv("Data_derived/grt.csv",
                     header = TRUE)



###################################################
## INITIALIZE OBJECTS FOR MOVEMENT SIMULATIONS ##
###################################################

tree <- start_pops

# Stores all simulations
sim_list <- list()

# Lists for movement paths
mp01 <- list()
mp02 <- list()
mp03 <- list()



##############################################
## BOOTSTRAP PARAMETER UNCERTAINTY ESTIMATES ##
##############################################

# Step lengths in encamping mode
sl_en <- boots(est$sl.encamp)

# Step lengths in travel mode
sl_tvl <- boots(est$sl.travel)

# Turning angle concentration in encamping
ta_en <- boots(est$ta.encamp)

# Turning angle concentration in traveling
ta_tvl <- boots(est$ta.travel)

# Habitat probabilities in encamping state
pr_en_hab <- boots(est$habitat.encamp)

d3 <- data.frame(pr.en = pr_en_hab,
                 hab = "frag",
                 grp = "spec")

# Matrix probabilities in encamping state
pr_en_mat <- boots(est$habitat.encamp)

d4 <- data.frame(pr.en = pr_en_mat,
                 hab = "matrix",
                 grp = "spec")

# Habitat probabilities in traveling state
pr_tvl_hab <- boots(est$habitat.travel)

d2 <- data.frame(pr.en = pr_tvl_hab,
                 hab = "frag",
                 grp = "gen")

# Matrix probabilities in traveling state
pr_tvl_mat <- boots(est$matrix.travel)

d4 <- data.frame(pr.en = pr_tvl_mat,
                 hab = "matrix",
                 grp = "gen")

# Habitat selection coefficients
sel_hab <- boots(log(est$exp.coef.))

d2 <- data.frame(rss = sel_hab,
                 grp = "spec")



##############################################
## FIT MIXED MODEL TO HABITAT SELECTION ##
##############################################

dx1 <- rbind(d1,d2)

# Gaussian mixed model
obj <- glmer(rss ~ grp + (1|grp),
             data = dx1,
             family = "gaussian")

summary(obj)

# Confidence intervals
confint(obj)

# Comment:
# gen matrix travel: -0.05 [-0.06,-0.05]



#################################################
## CONVERT BOOTSTRAP MEAN/SD TO DISTRIBUTIONS ##
#################################################

# STEP LENGTH DISTRIBUTIONS
# Gamma distributions for movement distances

encamp.params <- estGammaParams(sl_en[1],
                                sl_en[2])

travel.params <- estGammaParams(sl_tvl[1],
                                sl_tvl[2])



#########################################
## TURNING ANGLE CONCENTRATION ##
#########################################

# Beta distributions for directional persistence

encamp.conc <- estBetaParams(ta_en[1],
                             ta_en[2])

travel.conc <- estBetaParams(ta_tvl[1],
                             ta_tvl[2])



#########################################
## STATE PROBABILITIES IN HABITAT ##
#########################################

hab.encamp <- estBetaParams(pr_en_hab[1],
                            pr_en_hab[2])

hab.travel <- estBetaParams(pr_tvl_hab[1],
                            pr_tvl_hab[2])



#########################################
## STATE PROBABILITIES IN MATRIX ##
#########################################

mat.encamp <- estBetaParams(pr_en_mat[1],
                            pr_en_mat[2])

mat.travel <- estBetaParams(pr_tvl_mat[1],
                            pr_tvl_mat[2])



###########################################################
## MAIN SIMULATION LOOP ##
## m = simulation replicate
###########################################################

for (m in 1:100){
  
  trajectory_list <- list()
  
  
  #########################################################
  ## LOOP THROUGH EACH STARTING INDIVIDUAL ##
  #########################################################
  
  for(t in 1:nrow(tree)){
    
    
    #########################################################
    ## DRAW RANDOM MOVEMENT PARAMETERS ##
    #########################################################
    
    # Parameters when bird is in habitat
    
    df.h <- data.frame(
      
      # step lengths
      scl = c(
        rgamma(1,
               shape = encamp.params[1],
               scale = encamp.params[2]),
        
        rgamma(1,
               shape = travel.params[1],
               scale = travel.params[2])
      ) / 10,
      
      # directional persistence
      conc = c(
        rbeta(1,
              encamp.conc[1],
              encamp.conc[2]),
        
        rbeta(1,
              travel.conc[1],
              travel.conc[2])
      ),
      
      # probabilities of state choice
      pr = c(
        rbeta(1,
              hab.encamp[1],
              hab.encamp[2]),
        
        rbeta(1,
              hab.travel[1],
              hab.travel[2])
      )
    )
    
    
    #########################################################
    ## PARAMETERS WHEN BIRD IS IN MATRIX ##
    #########################################################
    
    df.m <- data.frame(
      scl = c(
        rgamma(1,
               shape = encamp.params[1],
               scale = encamp.params[2]),
        
        rgamma(1,
               shape = travel.params[1],
               scale = travel.params[2])
      ) / 10,
      
      conc = c(
        rbeta(1,
              encamp.conc[1],
              encamp.conc[2]),
        
        rbeta(1,
              travel.conc[1],
              travel.conc[2])
      ),
      
      pr = c(
        rbeta(1,
              mat.encamp[1],
              mat.encamp[2]),
        
        rbeta(1,
              mat.travel[1],
              mat.travel[2])
      )
    )
    
    
    
    #########################################################
    ## RANDOM BODY MASS ##
    #########################################################
    
    BodyMass.Value <- runif(
      1,
      min = min(est$Mass),
      max = max(est$Mass)
    )
    
    
    #########################################################
    ## GUT RETENTION TIME ##
    #########################################################
    
    # allometric scaling relationship
    
    grt <- 4.5 * (BodyMass.Value/1000)^0.5 * 60
    
    # stochastic GRT
    value.grt <- rgamma(
      n = 1,
      shape = grt^2/(10)^2,
      scale = (10)^2/grt
    )
    
    
    
    #########################################################
    ## INITIALIZE MOVEMENT TRACK ##
    #########################################################
    
    start_point <- as.data.frame(
      matrix(
        NA,
        nrow = round(value.grt/10),
        ncol = 2
      )
    )
    
    colnames(start_point) <- c("x","y")
    
    # first location = starting tree
    start_point[1,] <- tree[t,c("x","y")]
    
    
    #########################################################
    ## NUMBER OF MOVEMENT STEPS ##
    #########################################################
    
    iter <- round(value.grt/10)
    
    
    
    #########################################################
    ## TIME STEP LOOP ##
    #########################################################
    
    for (i in 1:(iter-1)){
      
      
      #########################################################
      ## DETERMINE CURRENT HABITAT ##
      #########################################################
      
      source.attr <- cbind(
        raster::extract(
          lscp,
          start_point[i,1:2],
          start_point = TRUE
        ),
        start_point[i,1:2]
      )
      
      colnames(source.attr) <- c("habitat","x","y")
      
      
      #########################################################
      ## IF BIRD IS CURRENTLY IN HABITAT ##
      #########################################################
      
      if (source.attr$habitat == 1){
        
        
        #########################################################
        ## CHOOSE BEHAVIORAL STATE ##
        #########################################################
        
        # choose encamp or travel
        c <- df.h[
          sample(
            c(1,2),
            size = 1,
            prob = df.h$pr
          ),
        ]
        
        
        #########################################################
        ## SIMULATE CORRELATED RANDOM WALK ##
        #########################################################
        
        traj1 <- simm.crw(
          date = 1:10,
          h = c$scl,
          r = c$conc,
          x0 = as.numeric(start_point[i,c("x","y")]),
          id = "A1",
          typeII = TRUE,
          proj4string = CRS()
        )
        
        
        #########################################################
        ## EXTRACT TRAJECTORY ##
        #########################################################
        
        df <- as.data.frame(traj1[[1]])
        
        df <- df %>%
          dplyr::select(x=x,y=y)
        
        
        #########################################################
        ## CHECK DESTINATION HABITAT ##
        #########################################################
        
        dest.attr <- cbind(
          raster::extract(
            lscp,
            df[10,1:2],
            df = TRUE
          ),
          df[10,1:2]
        )
        
        colnames(dest.attr) <- c("ID","habitat","x","y")
        
        
        #########################################################
        ## IF DESTINATION IS MATRIX ##
        #########################################################
        
        if (dest.attr$habitat == 0){
          
          # Bird may avoid matrix
          # habitat selection coefficient controls this
          
          s <- sample(
            c(1,10),
            size = 1,
            prob = c(
              exp(rnorm(
                1,
                mean = sel_hab[1],
                sd = sel_hab[2]
              )),
              1
            )
          )
          
          # either stay near origin or move fully
          start_point[i+1,] <- df[s,]
          
        } else {
          
          #########################################################
          ## DESTINATION IS HABITAT ##
          #########################################################
          
          # equal chance of using either point
          
          s <- sample(
            c(1,10),
            size = 1,
            prob = c(0.5,0.5)
          )
          
          start_point[i+1,] <- df[s,]
        }
        
        
      } else {
        
        
        #########################################################
        ## IF BIRD IS CURRENTLY IN MATRIX ##
        #########################################################
        
        e <- df.m[
          sample(
            c(1,2),
            size = 1,
            prob = df.m$pr
          ),
        ]
        
        
        #########################################################
        ## MATRIX MOVEMENT ##
        #########################################################
        
        traj2 <- simm.crw(
          date = 1:10,
          h = e$scl,
          r = e$conc,
          x0 = as.numeric(start_point[i,c("x","y")]),
          id = "A1",
          typeII = TRUE,
          proj4string = CRS()
        )
        
        
        #########################################################
        ## EXTRACT MOVEMENT ##
        #########################################################
        
        df <- as.data.frame(traj2[[1]])
        
        df <- df %>%
          dplyr::select(x=x,y=y)
        
        
        #########################################################
        ## CHECK DESTINATION ##
        #########################################################
        
        dest.attr <- cbind(
          raster::extract(
            lscp,
            df[10,1:2],
            df = TRUE
          ),
          df[10,1:2]
        )
        
        colnames(dest.attr) <- c("ID","habitat","x","y")
        
        
        #########################################################
        ## DESTINATION IS HABITAT ##
        #########################################################
        
        if (dest.attr$habitat == 1){
          
          # bird prefers habitat
          
          s <- sample(
            c(1,10),
            size = 1,
            prob = c(
              1,
              rnorm(
                1,
                mean = sel_hab[1],
                sd = sel_hab[2]
              )
            )
          )
          
          start_point[i+1,] <- df[s,]
          
        } else {
          
          #########################################################
          ## DESTINATION IS MATRIX ##
          #########################################################
          
          s <- sample(
            c(1,10),
            size = 1,
            prob = c(0.5,0.5)
          )
          
          start_point[i+1,] <- df[s,]
        }
      }
      
      
      #########################################################
      ## STORE TRAJECTORY ##
      #########################################################
      
      mp01[[i]] <- df
      
    } ## END i LOOP
    
    
    
    #########################################################
    ## COMBINE TRAJECTORIES ##
    #########################################################
    
    mpdf <- do.call("rbind",mp01)
    
    mpdf$id <- t
    
    start_point$id <- t
    
    trajectory_list[[t]] <- start_point
    
    mp02[[t]] <- mpdf
    
  } ## END t LOOP
  
  
  
  #########################################################
  ## COMBINE ALL INDIVIDUALS ##
  #########################################################
  
  tl <- do.call("rbind",trajectory_list)
  
  tl$sim <- m
  
  sim_list[[m]] <- tl
  
  
  
  #########################################################
  ## COMBINE MOVEMENT PATHS ##
  #########################################################
  
  mpdf01 <- do.call("rbind",mp02)
  
  mpdf01$sim <- m
  
  mp03[[m]] <- mpdf01
  
} ## END m LOOP



#########################################################
## FINAL COMBINED OUTPUT ##
#########################################################

ul <- do.call("rbind",sim_list)

mp_df <- do.call("rbind",mp03)



#########################################################
## ADD METADATA ##
#########################################################

ul$br <- 0.8
mp_df$br <- 0.8

ul$frag <- 0.9
mp_df$frag <- 0.9

ul$habitat <- 0.75
mp_df$habitat <- 0.75



#########################################################
## APPEND TO EXISTING DATASETS ##
#########################################################

al <- read.csv(
  "Data_derived/new/ldd_fragments_05.csv",
  header = TRUE
)

gl <- rbind(al,ul)

write.csv(
  gl,
  "Data_derived/new/ldd_fragments_05.csv",
  row.names = FALSE
)



#########################################################
## SAVE MOVEMENT PATHS ##
#########################################################

mp_al <- read.csv(
  "Data_derived/new/mp_fragments_05.csv",
  header = TRUE
)

mp_gl <- rbind(mp_df,mp_al)

write.csv(
  mp_gl,
  "Data_derived/new/mp_fragments_05.csv",
  row.names = FALSE
)
