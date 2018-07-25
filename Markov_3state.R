#####################################################################################
##########             Simple 3-state Markov model in R         #####################
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

rm(list = ls())      # clear memory (removes all the variables from the workspace)

#### 01 Load packages ####
# no packages required
# library(ggplot2)
# library(dampack)
# library(dplyr)
# library(scales)
# library(ellipse)

#### 02 Load Functions ####
source("../Functions.R")

#### 03 Input Model Parameters ####
## Strategy names
v.names.str <- c("Base Case")  
## Number of strategies
n.str <- length(v.names.str)
## Markov model parameters
v.n  <- c("Healthy", "Sick", "Dead")    # state names
n.s  <- length(v.n)                     # number of states
n.t  <- 60                              # number of cycles

p.HD <- 0.02                    # probability to die when healthy
p.HS <- 0.05                    # probability to become sick when healthy
p.SD <- 0.1                     # probability to die when sick

# Costs and utilities  
c.H  <- 400                     # cost of remaining one cycle healthy
c.S  <- 100                     # cost of remaining one cycle sick
c.D  <- 0                       # cost of remaining one cycle dead
u.H  <- 0.8                     # utility when healthy 
u.S  <- 0.5                     # utility when sick
u.D  <- 0                       # utility when dead
d.r  <- 0.03                    # discount rate per cycle, same for costs and effectiveness
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


#### 04 Define and initialize matrices and vectors ####
#### 04.1 Cohort trace ####
# create the cohort trace
m.M <- matrix(NA, 
              nrow = n.t + 1 ,             # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
              ncol = n.s, 
              dimnames = list(0:n.t, v.n))

m.M[1, ] <- c(1, 0, 0)                     # initialize first cycle of Markov trace

#### 04.2 Transition probability MATRIX ####
# create the transition probability matrix
m.P  <- matrix(0,
               nrow = n.s,
               ncol = n.s,
               dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P

# fill in the transition probability matrix
### From Healthy
m.P["Healthy", "Healthy"] <- 1 - p.HD - p.HS
m.P["Healthy", "Sick"]    <- p.HS
m.P["Healthy", "Dead"]    <- p.HD

### From Sick
m.P["Sick", "Sick"] <- 1 - p.SD
m.P["Sick", "Dead"] <- p.SD

### From Dead
m.P["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P)

#### 05 Run Markov model ####
for (t in 1:n.t){                         # loop through the number of cycles
  m.M[t + 1, ] <- m.M[t, ] %*% m.P        # estimate the state vector for the next cycle (t + 1)
}

#### 06 Compute and Plot Epidemiological Outcomes ####
#### 06.1 Cohort trace #####
matplot(m.M, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("right", v.n, col = c("black", "red", "green"), lty = 1:3, bty = "n")  # add a legend to the graph

#### 06.2 Overall Survival (OS) #####
v.os <- 1 - m.M[, "Dead"]                  # calculate the overall survival (OS) probability
v.os <- rowSums(m.M[, 1:2])                # alternative way of calculating the OS probability   

plot(v.os, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

#### 06.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os)                             # summing probablity of OS over time  (i.e. life expectancy)

#### 06.3 Disease prevalence #####
v.prev <- m.M[, "Sick"]/v.os
plot(v.prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 07 Compute Cost-Effectiveness Outcomes ####
#### 07.1 Mean Costs and QALYs ####
# per cycle
v.tc <- m.M %*% c(c.H, c.S, c.D)  # calculate expected costs by multiplying m.M with the cost vector for the different health states   
v.tu <- m.M %*% c(u.H, u.S, u.D)  # calculate expected QALYs by multiplying m.M with the utilities for the different health states   

#### 07.2 Discounted Mean Costs and QALYs ####
v.tc.d <-  t(v.tc) %*% v.dwc   # Discount costs  by multiplying the cost vector with discount weights (v.dw) 
v.te.d <-  t(v.tu) %*% v.dwe   # Discount QALYS  by multiplying the QALYs vector with discount weights (v.dw)

results <- data.frame( "Total Discounted Cost" = v.tc.d, 
                       "Life Expectancy" = v.le, 
                       "Total Discounted QALYs" = v.te.d, 
                       check.names = F)
results
