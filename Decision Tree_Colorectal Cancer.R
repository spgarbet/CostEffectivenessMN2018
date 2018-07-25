######  Follow up after Colorectal Cancer - A Decision Tree in R  ############

##############################################################################
## Credits for the R Code: Petros Pechlivanoglou, Fernando Alarid- Escudero ##
##                         Hawre Jalal and Mohammad Kaviul Kahn             ##
## Credits for the example: A. Gray et al 2011                              ##
## Applied Methods of Cost-effectiveness Analysis in Healthcare             ##
## Cost effectiveness of follow-up strategies after colon cancer surgery    ##
## Strategies: PC(Primary Care), HC (Hospital Care), RP (Routine Practice)  ##
##############################################################################

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
library(ggplot2)
library(dampack)
library(dplyr)
library(scales)
library(ellipse)

#### 02 Load Functions ####
source("../Functions.R")

#### 03 Input Model Parameters ####
## Strategy names
v.names.str <- c("Routine Practice", "Primary Care", "Hospital Care")  
## Number of strategies
n.str <- length(v.names.str)
## Decision tree parameters
p.PCed <- 0.40   # Probability of early detection PC
p.HCed <- 0.45   # Probability of early detection HC
p.RPed <- 0.35   # Probability of early detection RP 

le.ed  <- 7      # Life expectancy after early detection
le.ld  <- 1      # Life expectancy after  late detection 

c.PCed <-  3900  # Total costs after early detection with PC
c.HCed <-  6200  # Total costs after early detection with HC
c.RPed <-  3030  # Total costs after early detection with RP
c.PCld <- 12800  # Total costs after  late detection with PC
c.HCld <- 14400  # Total costs after  late detection with HC
c.RPld <- 12020  # Total costs after  late detection with RP
wtp    <- 10000  # Define willingness to pay threshold


#### 04 Run Decision Tree ####

# the solution of the tree is the sum of the weights (probabilities of each path) times the reward of each path

# for costs
c.RP <- p.RPed * c.RPed + (1 - p.RPed) * c.RPld # Routine Practice cost
c.PC <- p.PCed * c.PCed + (1 - p.PCed) * c.PCld # Primary Care cost
c.HC <- p.HCed * c.HCed + (1 - p.HCed) * c.HCld # Hospital Care cost

# ...and for effects
e.RP <- p.RPed * le.ed + (1 - p.RPed) * le.ld   # Routine Practice life expectancy
e.PC <- p.PCed * le.ed + (1 - p.PCed) * le.ld   # Primary Care life expectancy
e.HC <- p.HCed * le.ed + (1 - p.HCed) * le.ld   # Hospital Care life expectancy``

v.le <- c(e.RP, e.PC, e.HC)                       # vector of life expectancies
v.c  <- c(c.RP, c.PC, c.HC)                       # vector of total costs
names(v.le) <- names(v.c) <- v.names.str          # names for the elements of the two vectors

# Matrix with costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.c,
                   Effect   = v.le)

#### 05 Compute ICERs of Decision Tree ####
m.cea <- calculate_icers(m.ce)
m.cea

#### 06 Plot frontier of Decision Tree ####
front.ce <- getFrontier(m.ce, plot = F)
plot.frontier(CEmat = m.ce, frontier = front.ce)
ggsave("figs/DecTree-CEA-Frontier.png", width = 8, height = 6)
