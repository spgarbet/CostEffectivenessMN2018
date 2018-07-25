#---------------------------------------------------#
#### R functions to compute and plot CE Frontier ####
#---------------------------------------------------#
#####################################################################################
getFrontier <- function(CEmat, maxWTP = Inf, plot = TRUE){
  # Name: getFrontier.R
  # Goal: Find the CEA frontier, up to a given WTP level, by 
  #       identifying strategies with the highest NMB
  # Originally written by: Sze Suen on Feb 25, 2015
  # Citation: Suen S-C, Goldhaber-Fiebert JD. An Efficient, 
  #           Noniterative Method of Identifying the Cost-Effectiveness Frontier. 
  #           Med Decis Making. 2016
  #           https://www.ncbi.nlm.nih.gov/pubmed/25926282
  # Modified by: Fernando Alarid-Escudero on July 20, 2015 
  # Notes: 
  #    ~Frontier strategies are displaced on the R output screen and 
  #      plotted in red on the scatter plot.  
  #
  #    ~User needs to provide a csv file of costs and QALYs 
  #	    (CostQalyInputFile_online_supp.csv) inside the folder specified    
  #	    below (inputFolder). The CSV should have three columns (labeled 
  #     in first row) in this order: 
  #      Strategy number, costs, and QALYs.
  #
  #    ~User can specify the maximum willingness-to-pay level to 
  #      consider (maxWTP).  Can be Inf for infinity.
  #
  #    ~QALY-reducing strategies will be on the frontier if they save
  #      enough money (script assumes maximum willingness to save money 
  #      per QALY lost is equivalent to maximum willingness to pay per QALY
  #      gained). If the user does not wish to consider such policies as
  #      being on the frontier, do not include strategies with negative 
  #      QALYs in the input csv file.
  #
  #    ~Script does not use the one-line code cited in the text
  #      as the max function is slow. This implementation is
  #      faster and methodologically does the same thing.
  #
  #    ~May take a few minutes if thousands of strategies and 
  #       processing resources are low.  Please be patient.
  #
  #    Please cite article if this code is used.
  #
  # USER INPUTS:
  #inputFolder <- "CostEffectivenessFrontier_MDM/"
  #maxWTP <- Inf        # any positive value or Inf
  
  ## Clean everythng from workspace
  #rm(list=ls())
  ####################################################################
  ####################################################################
  
  # check for duplicated strategies
  dups <- CEmat[c(duplicated(CEmat[, 2:3]) | duplicated(CEmat[, 2:3], fromLast = TRUE)), 1]
  
  # initialize some variables
  costsCol <- 2; qalyCol <- 3
  numStrat <- nrow(CEmat)
  
  # find WTP levels to test so that all strategies on frontier will be captured
  # this means testing on either side of all NMB intersections, which are just all the pairwise ICERs
  ICERmat <- matrix(1, numStrat, numStrat)
  suppressWarnings(
  for (i in 1:numStrat) {
    indexStrat <- matrix(1, numStrat, 3)
    indexStrat[, costsCol] <- indexStrat[, costsCol] * CEmat[i, costsCol]
    indexStrat[, qalyCol] <- indexStrat[, qalyCol] * CEmat[i, qalyCol]
    delCostQalys <- CEmat - indexStrat
    ICERmat[, i] <- delCostQalys[, costsCol] / delCostQalys[, qalyCol]
  }  
  )
  intersections <- sort(unique(c(ICERmat)))
  intersections <- intersections[is.finite(intersections)]
  WTPtestPoints <- c(0, intersections [intersections >= 0 & intersections <= maxWTP], maxWTP)
  
  # Find the strategy with the max NMB at each of the WTP test points
  indiciesOfMax <- vector()
  NMBmat <- matrix(0, numStrat, length(WTPtestPoints))
  for (i in 1:length(WTPtestPoints) ) {
    NMBmat[, i] <- (WTPtestPoints[i] * CEmat[, qalyCol]) - CEmat[, costsCol]
  }
  if (is.infinite(maxWTP)) {
    #WTP of infinity means costs are not considered
    NMBmat[, length(WTPtestPoints)] = CEmat[, qalyCol] - (0 * CEmat[, costsCol]); 
  }
  maxVals <- apply(NMBmat, 2, max)  #find strategy that maximizes NMB at each WTP
  for (i in 1:length(WTPtestPoints) ) {  #find all strategies that match max at each WTP
    indiciesOfMax <- c(indiciesOfMax, which( NMBmat[, i] == maxVals[i]))
  }
  v.frontier <- unique(indiciesOfMax)  #find strategy that maximizes NMB at each WTP
  
  if (plot == TRUE){
    # display out: make plot and print to output screen
    plot(CEmat[v.frontier, qalyCol], CEmat[v.frontier, costsCol], col = 'red', pch = 16, 
         xlab = "Effectiveness", ylab = "Cost")
    points(CEmat[, qalyCol], CEmat[, costsCol])
    lines(CEmat[v.frontier, qalyCol], CEmat[v.frontier, costsCol])
    if (length(dups) > 0){
      warning("Strategies have the same costs and benefits (displayed above)")
      print(dups)
    }
  }
  sprintf("Frontier is formed by strategies: %s", paste( sort(CEmat[v.frontier, 1]), collapse=" "))
  
  class(v.frontier) <- "frontier"
  return(v.frontier)
}

plot.frontier <- function(CEmat, frontier, 
                         ncol = 1,
                         coord.flip = F,
                         txtsize = 16)
{
  # A function to plot CE frontier
  # USER INPUTS:
  #   CEmat: A CE matrix arranged as: Col1: Strategy; Col2: Cost; Col3: Effectiveness
  # Create a dataframe from matrix
  CEmat.df <- data.frame(CEmat)
  colnames(CEmat.df)[3] <- "Effectiveness"
  n.strategies <- nrow(CEmat.df)
  # Make Strategies as factor
  CEmat.df$Strategy <- as.factor(CEmat.df$Strategy)
  #
  if (coord.flip == T) {
    ggplot(CEmat.df, aes(Effectiveness, Cost)) +
      geom_point(aes(color = Strategy, shape = Strategy), size = 4) + 
      coord_flip() +
      scale_y_continuous("Cost ($)", labels = comma) +
      ggtitle("Cost-Effectiveness Frontier") +
      geom_point(data = CEmat.df[frontier, ], 
                 aes(Effectiveness, Cost, shape = Strategy, color = Strategy), size = 4) +
      geom_line(data = CEmat.df[frontier, ], aes(Effectiveness, Cost)) +
      scale_shape_manual(values = 0:(n.strategies - 1)) +
      guides(shape = guide_legend(ncol = ncol)) +
      theme_bw(base_size = txtsize) +
      theme(legend.position = "bottom")
  } else {
    ggplot(CEmat.df, aes(Effectiveness, Cost)) +
      geom_point(aes(color = Strategy, shape = Strategy), size = 4) + 
      scale_y_continuous("Cost ($)", labels = comma) +
      ggtitle("Cost-Effectiveness Frontier") +
      geom_point(data = CEmat.df[frontier, ], 
                 aes(Effectiveness, Cost, shape = Strategy, color = Strategy), size = 4) +
      geom_line(data = CEmat.df[frontier, ], aes(Effectiveness, Cost)) +
      scale_shape_manual(values = 0:(n.strategies - 1)) +
      guides(shape = guide_legend(ncol = ncol)) +
      theme_bw(base_size = txtsize) +
      theme(legend.position = "bottom")
  }
}

#-------------------------------------------------#
#### R functions to compute ICERs of ANY model ####
#-------------------------------------------------#
# Source: https://miqdad.freeasinspeech.org.uk/icer_calculator/
compute_icers = function(non_dominated){
  icers = non_dominated %>% 
    arrange(Cost,desc(Effect))
  
  if(nrow(non_dominated)>1){  
    icers[1,"ICER"] = NA
    for(i in 2:nrow(non_dominated)){
      icers[i,"ICER"] = (icers[i,"Cost"] - icers[i-1,"Cost"])/(icers[i,"Effect"] - icers[i-1,"Effect"])
    }
  }
  return(icers)
}

calculate_icers = function(data){
  # check data is in correct format
  if(nrow(data)<2 | ncol(data)!=3){
    return(NULL)
  }
  colnames(data) = c("Strategy","Cost","Effect")
  # remove dominated strategies
  data = data %>% arrange(Cost,desc(Effect))
  dominated = data[FALSE,]
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      if(data[j,"Effect"]<=data[i,"Effect"]){
        dominated = union(dominated,data[j,])
      }
    }
  }
  non_dominated = setdiff(data, dominated)
  
  # remove extendedly dominated strategies
  # extended dominance means ICER > that of both previous and next strategies in order of cost
  ext_dominated = non_dominated[FALSE,]
  if(nrow(non_dominated)>3){
    ext_dominated_count = nrow(ext_dominated)
    prev_ext_dominated_count = -1 # run atleast once
    while(ext_dominated_count>prev_ext_dominated_count & (nrow(non_dominated)-ext_dominated_count)>3){
      non_dominated_set = compute_icers(setdiff(non_dominated,ext_dominated))
      prev_ext_dominated_count = ext_dominated_count
      for(i in 3:(nrow(non_dominated_set)-1)){
        if(non_dominated_set[i,"ICER"]>non_dominated_set[i-1,"ICER"] & non_dominated_set[i,"ICER"]>non_dominated_set[i+1,"ICER"]){
          ext_dominated = union(ext_dominated,non_dominated_set[i,c("Strategy","Cost","Effect")])
        }
      }
      ext_dominated_count = nrow(ext_dominated)
    }
  }
  
  # calculate ICERs for those strategies not dominated
  # nor etendedly dominated
  non_ext_dominated = setdiff(non_dominated,ext_dominated) %>%
    compute_icers() %>% 
    mutate(Status = "ND")
  
  dominated = dominated %>%
    mutate(ICER = NA, Status = "D")
  
  ext_dominated = ext_dominated %>%
    mutate(ICER = NA, Status = "ED")
  
  # recombine all results to produce final output
  results = bind_rows(non_ext_dominated,dominated,ext_dominated)
  
  return(results)
}

#----------------------------------------------------------#
#### R functions for deterministic sensitivity analyses ####
#----------------------------------------------------------#

### ====================================================
###     Function for plotting One-way SA Diagrams
### ====================================================
owsa.plot.det <- function(param, outcomes,
                          paramName,
                          strategyNames, 
                          outcomeName = "Outcome"){
  ## Load required packages
  library(ggplot2)
  library(reshape2)
  library(scales)
  
  owsa.df <- data.frame(outcomes)
  colnames(owsa.df) <- strategyNames
  owsa.df$param <- param
  
  owsa.gg <- ggplot(data = melt(owsa.df, id.vars = "param", 
                                variable.name = "Strategy"), 
                    aes(x = param, y = value, color = Strategy)) +
    geom_line() +
    xlab(paramName) +
    ggtitle("One-way sensitivity analysis", subtitle = outcomeName)+
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  return(owsa.gg)
}

### ====================================================
###     Function for plotting Two-way SA Diagrams
### ====================================================
twsa.plot.det <- function(params, outcomes, strategyNames, outcomeName = "Outcome", mx = T){
  ## Load required packages
  library(ggplot2)
  library(reshape2)
  library(scales)
  
  parm1 <- colnames(params)[1]
  parm2 <- colnames(params)[2]
  if (mx == T){
    strategy <- factor(max.col(outcomes), labels = strategyNames)
  } else {
    strategy <- factor(max.col(-outcomes), labels = strategyNames)
  }
  
  twsa.df <- data.frame(strategy)
  #A simple trick to define my variables in my functions environment
  twsa.df$parm1 <- params[, parm1]
  twsa.df$parm2 <- params[, parm2]
  
  twsa.gg <- ggplot(twsa.df, aes(x = parm1, y = parm2)) + 
    geom_tile(aes(fill = strategy)) +
    theme_bw(base_size = 14) +
    xlab(parm1) +
    ylab(parm2) +
    ggtitle("Two-way sensitivity analysis", subtitle = outcomeName) +
    scale_fill_discrete("Strategy: ", l = 50) +
    theme(legend.position = "bottom")
  
  return(twsa.gg)
}

#### =========================
### Tornado Plots
## Test code:
# paramNames <-  c( "Param 1 [Low/High CI]",
#                   "Param 2 [Low/High CI]",
#                   "Param 3 [-/+ 15%]",
#                   "Param 4 [-/+ 15%]"
# )
# 
# # data structure: ymean	ymin	ymax
# data <- matrix(c(100,	80,	120,
#                  100,	25,	150,
#                  100, 95,	120,
#                  100, 75, 160), nrow = 4, ncol = 3, byrow = TRUE)
# 
# data
# TornadoPlot(Parms = paramNames, Outcomes = data, titleName = "Tornado Plot")
# Parms = paramNames
# Outcomes = data
# titleName = "Tornado Plot"

### ====================================================
###     Function for plotting Tornado Diagrams
### ====================================================
TornadoPlot <- function(Parms, Outcomes, titleName, outcomeName, ylab = "$"){
  # Parm:        vector with parameter names  
  # Outcomes:    matrix including parameter specific outcomes (number of Parm x 3)
  # titleName:   title of the plot (e.g Tornado Plot)
  # outcomeName: name of the outcome shown in the Tornado plot
  
  library(ggplot2)
  library(reshape2)
  library(scales)
  
  # Grouped Bar Plot
  # Determine the overall optimal strategy
  paramNames2 <- Parms
  
  # Combine the parameter list with the data
  ymean <- Outcomes[1, 1]
  
  yMin <- Outcomes[, 2] - ymean
  yMax <- Outcomes[, 3] - ymean
  ySize <- abs(yMax - yMin)  # High value - Low value
  
  rankY<- order(ySize)
  nParams <- length(paramNames2)
  
  Tor <- data.frame(
    Parameter = c(paramNames2[rankY], paramNames2[rankY]),  
    Level = c(rep("Low", nParams), rep("High", nParams)),
    value = ymean + c(yMin[rankY], yMax[rankY]),
    sort = seq(1, nParams)
  )
  
  #re-order the levels in the order of appearance in the data.frame
  Tor$Parameter2 <- ordered(Tor$Parameter, Tor$Parameter[1:(length(Tor$Parameter) / 2)])
  # Tor$Parameter2 <- factor(Tor$Parameter, as.character(Tor$Parameter))
  #Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
  offset_trans <- function(offset = 0) {
    trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x + offset)
  }
  #Plot the Tornado diagram.
  txtsize <- 12
  tor.gg <- ggplot(Tor[Tor$Level == "Low", ], aes(x = Parameter2, y = value, fill = level)) +
    geom_bar(stat = "identity", fill = "blue") +
    ggtitle("Tornado Plot", subtitle = outcomeName) +
    scale_fill_discrete("Parameter Level: ", l = 50) +
    scale_y_continuous(name = ylab, trans = offset_trans(offset = ymean)) +
    scale_x_discrete(name = "Parameter") +
    geom_bar(data = Tor[Tor$Level == "High", ], aes(x = Parameter2, y = value, fill = level), stat = "identity", fill = "red", alpha = 0.5) +
    geom_hline(yintercept = ymean, linetype = "dotted", size = 0.5) +
    theme_bw(base_size = 14) +
    coord_flip() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = txtsize, angle = 0, hjust = 1),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face = "bold", size = txtsize),
          axis.title.y = element_text(face = "bold", size = txtsize),
          axis.text.y = element_text(size = txtsize),
          axis.text.x = element_text(size = txtsize),
          axis.ticks.y = element_blank())
  
  return(tor.gg)
  # ggsave(paste("results/", titleName,".png"))
}

#----------------------------#
#### Scatter CEA frontier ####
#----------------------------#
ScatterCE <- function(strategies, m.c, m.e, 
                      ncol = 1, txtsize = 16) {
  if(!is.data.frame(m.c))
    m.c <- data.frame(m.c)
  if(!is.data.frame(m.e))
    m.e <- data.frame(m.e)
  
  df.c  <- reshape2::melt(m.c, variable.name = "Strategy", id.vars = NULL)
  levels(df.c$Strategy) <- strategies
  df.e <- reshape2::melt(m.e, variable.name = "Strategy", id.vars = NULL)
  levels(df.e$Strategy) <- strategies
  
  ## Create matrix with costs and effectiveness
  m.ce <- cbind(df.c, df.e[, 2])
  colnames(m.ce) <- c("Strategy", "Cost", "Effectiveness")
  
  # Ellipses code
  df.ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(m.ce$Strategy)) {
    # create 100 points per variable around the mean of each group
    df.ell <- rbind(df.ell, 
                    cbind(as.data.frame(with(m.ce[m.ce$Strategy == g,], 
                                             ellipse(cor(Effectiveness, Cost), 
                                                     scale = c(sd(Effectiveness), sd(Cost)), 
                                                     centre = c(mean(Effectiveness), mean(Cost)))
                    )), group = g))
  }
  
  Means <- m.ce %>% 
    dplyr::group_by(Strategy) %>%
    dplyr::summarise(N = length(Cost),
                     Cost.mean = mean(Cost),
                     Eff.mean = mean(Effectiveness)) 
  
  #Define ggplot object
  scatter.ce.gg <- ggplot(m.ce, aes(x = Effectiveness, y = Cost, color = Strategy)) + 
    geom_point(size = 0.7, alpha = 0.2, shape = 21) +
    geom_point(data = Means, aes(x = Eff.mean, y = Cost.mean, shape = Strategy),
               size = 8, fill = "white") +
    geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label = 1:length(strategies)), 
              size = 5, colour = "gray", alpha = 1) +
    geom_path(data = df.ell, aes(x = x, y = y, colour = group), 
              size = 1, linetype = 2, alpha = 1) + # draw ellipse lines
    ggtitle("Cost-Effectiveness Scatterplot") +
    scale_colour_discrete(l = 50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar) +
    scale_x_continuous(breaks =number_ticks(6)) +
    guides(shape = guide_legend(ncol = ncol),
           color = guide_legend(ncol = ncol)) +
    theme_bw(base_size = txtsize) +
    theme(legend.position = "bottom")
  return(scatter.ce.gg)
}

#----------------------------#
#### Expected Loss Curves ####
#----------------------------#
elc <- function(v.wtp, strategies, m.e, m.c) {
  library(reshape2)
  library(ggplot2)
  library(scales)
  n.sim <- nrow(m.e)
  n.str  <- ncol(m.e)
  m.exp.loss <- matrix(0, nrow = length(v.wtp), ncol = n.str)
  for(l in 1:length(v.wtp)) {
    m.nmb <- m.e * v.wtp[l] - m.c # Effectiveness minus Costs, with vector indexing
    max.str <- max.col(m.nmb)
    m.loss <- m.nmb - m.nmb[cbind(1:n.sim, max.str)]
    m.exp.loss[l, ] <- colMeans(m.loss)
  }
  # Optimal strategy based on lowest expected loss
  optimal.str <- max.col(m.exp.loss)
  # Expected loss of optimal strategy
  optimal.el <- m.exp.loss[cbind(1:length(v.wtp), optimal.str)]
  # Format expected loss for plotting
  df.exp.loss <- data.frame(cbind(v.wtp, m.exp.loss, optimal.el))
  colnames(df.exp.loss) <- c("WTP", strategies, "Frontier & EVPI")
  df.exp.loss.plot <- melt(df.exp.loss, 
                           id.vars = "WTP", 
                           variable.name = "Strategy")
  ## Plot expected losses
  # Format to plot frontier
  strats <- 1:(length(unique(df.exp.loss.plot$Strategy)) - 1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot ELC
  elc.gg <- ggplot(data = df.exp.loss.plot, aes(x = WTP/1000, y = -value)) +
    geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
    geom_line(aes(color = Strategy)) +
    ggtitle("Expected Loss Curves") + 
    #scale_colour_hue(l = 50, values = colors) +
    scale_x_continuous(breaks = number_ticks(20))+
    scale_y_continuous(labels = comma, breaks = number_ticks(8)) +
    xlab("Willingness to Pay (Thousand $/QALY)") +
    ylab("Expected loss ($)") +
    scale_shape_manual(values = point.shapes) +
    scale_color_manual(values = colors) + 
    scale_size_manual(values = point.size) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  return(elc.gg)
}

#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#
samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}