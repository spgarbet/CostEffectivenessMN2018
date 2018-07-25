#########              Sick-Sicker Markov model                 #####################

#####################################################################################
## This code forms the basis for the cohort model of the article:                  ## 
## 'Microsimulation modeling for health decision sciences using R: a tutorial'     ##
## Authors: Eline Krijkamp, Fernando Alarid-Escudero, Eva Enns,                    ##
##          Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou                  ##
## citation below                                                                  ##
#####################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	
# Petros Pechlivanoglou, PhD (5) 

# In collaboration of: 		
# 1 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 2 Erasmus MC, Rotterdam, The Netherlands
# 3 Harvard T.H. Chan School of Public Health, Boston, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 5 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Please cite our publications when using this code
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property 
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without written permission.
#####################################################################################

#####################################################################################

rm(list = ls())      # clear memory (removes all the variables from the workspace)

#### 01 Load packages ####
library(ggplot2)
library(dampack)
library(dplyr)
library(scales)
library(ellipse)

#### 02 Load Functions ####
source("../Functions.R")

#### 03 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Treatment", "Treatment")  
## Number of strategies
n.str <- length(v.names.str)
## Markov model parameters
age     <- 25                                 # age at baseline
max.age <- 55                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles
v.n  <- c("H", "S1", "S2", "D")               # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s <- length(v.n)                            # number of health states 

## Transition probabilities (per cycle) and hazard ratios
p.HD    <- 0.005           # probability to die when healthy
p.HS1   <- 0.15          	 # probability to become sick when healthy
p.S1H   <- 0.5           	 # probability to become healthy when sick
p.S1S2  <- 0.105         	 # probability to become sicker when sick
hr.S1   <- 3             	 # hazard ratio of death in sick vs healthy
hr.S2   <- 10            	 # hazard ratio of death in sicker vs healthy 
r.HD    <- - log(1 - p.HD) # rate of death in healthy
r.S1D   <- hr.S1 * r.HD  	 # rate of death in sick
r.S2D   <- hr.S2 * r.HD  	 # rate of death in sicker
p.S1D   <- 1 - exp(-r.S1D) # probability to die in sick
p.S2D   <- 1 - exp(-r.S2D) # probability to die in sicker

## Cost and utility inputs 
c.H     <- 2000            # cost of remaining one cycle in the healthy state
c.S1    <- 4000            # cost of remaining one cycle in the sick state
c.S2    <- 15000           # cost of remaining one cycle in the sicker state
c.Trt   <- 12000           # cost of treatment(per cycle)
c.D     <- 0               # cost of being in the death state
u.H     <- 1               # utility when healthy
u.S1    <- 0.75            # utility when sick
u.S2    <- 0.5             # utility when sicker
u.D     <- 0               # utility when dead
u.Trt   <- 0.95            # utility when being treated

# Discounting factor
d.r <- 0.03                                   # equal discount of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r

#### 04 Define and initialize matrices and vectors ####
#### 04.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_trt <- m.M_trt <- matrix(NA, 
                                nrow = n.t + 1, ncol = n.s,
                                dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_trt) # show first 6 rows of the matrix 

# The cohort starts as healthy
m.M_no_trt[1, ] <- m.M_trt[1, ] <- c(1, 0, 0, 0) # initiate first cycle of cohort trace 

#### 04.2 Transition probability MATRIX ####
# create the transition probability matrix for NO treatment
m.P_notrt  <- matrix(0,
                     nrow = n.s,
                     ncol = n.s,
                     dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_notrt

# fill in the transition probability array
### From Healthy
m.P_notrt["H", "H"]  <- 1 - (p.HS1 + p.HD)
m.P_notrt["H", "S1"] <- p.HS1
m.P_notrt["H", "D"]  <- p.HD
### From Sick
m.P_notrt["S1", "H"]  <- p.S1H
m.P_notrt["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
m.P_notrt["S1", "S2"] <- p.S1S2
m.P_notrt["S1", "D"]  <- p.S1D
### From Sicker
m.P_notrt["S2", "S2"] <- 1 - p.S2D
m.P_notrt["S2", "D"]  <- p.S2D
### From Dead
m.P_notrt["D", "D"] <- 1

# check rows add up to 1
rowSums(m.P_notrt)

# create transition probability matrix for treatment same as NO treatment
m.P_trt <- m.P_notrt

#### 05 Run Markov model ####
for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_trt[t + 1, ] <- t(m.M_no_trt[t, ]) %*% m.P_notrt # estimate the Markov trace for cycle the next cycle (t + 1)
     m.M_trt[t + 1, ] <- t(m.M_trt[t, ])    %*% m.P_trt   # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_trt)  # show the first 6 lines of the matrix


#### 06 Compute and Plot Epidemiological Outcomes ####
#### 06.1 Cohort trace #####
matplot(m.M_no_trt, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

#### 06.2 Overall Survival (OS) #####
#### 06.2 Overall Survival (OS) #####
v.os_no_trt <- 1 - m.M_no_trt[, "D"]       # calculate the overall survival (OS) probability for no treatment
v.os_no_trt <- rowSums(m.M_no_trt[, 1:3])  # alternative way of calculating the OS probability   

plot(0:n.t, v.os_no_trt, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

#### 06.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os_no_trt)                       # summing probablity of OS over time  (i.e. life expectancy)

#### 06.3 Disease prevalence #####
v.prev <- rowSums(m.M_no_trt[, c("S1", "S2")])/v.os_no_trt
plot(v.prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 06.4 Proportion of sick in S1 state #####
v.prop.S1 <- m.M_no_trt[, "S1"] / v.prev
plot(0:n.t, v.prop.S1,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in S1 state", 
     col = "black", type = "l")

#### 07 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_trt <- c(u.H, u.S1, u.S2, u.D)
v.u_trt    <- c(u.H, u.Trt, u.S2, u.D)

v.c_no_trt <- c(c.H, c.S1, c.S2, c.D)
v.c_trt    <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, c.D)

#### 07.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_trt <- m.M_no_trt %*% v.u_no_trt
v.tu_trt    <- m.M_trt    %*% v.u_trt

v.tc_no_trt <- m.M_no_trt %*% v.c_no_trt
v.tc_trt    <- m.M_trt    %*% v.c_trt

#### 07.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_trt <- t(v.tu_no_trt) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
tu.d_trt    <- t(v.tu_trt)    %*% v.dwe

tc.d_no_trt <- t(v.tc_no_trt) %*% v.dwc
tc.d_trt    <- t(v.tc_trt)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_trt, tc.d_trt)
v.tu.d <- c(tu.d_no_trt, tu.d_trt)

# Matrix with discounted costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.tc.d,
                   Effect   = v.tu.d)

#### 08 Compute ICERs of Decision Tree ####
m.cea <- calculate_icers(m.ce)
m.cea

#### 09 Plot frontier of Decision Tree ####
front.ce <- getFrontier(m.ce, plot = F)
plot.frontier(CEmat = m.ce, frontier = front.ce)
ggsave("figs/Markov-SickSicker-CEA-Frontier.png", width = 8, height = 6)
