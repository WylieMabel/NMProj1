# Network Modeling - HS 2024
# C. Stadtfeld, A. Uzaheta and I. Smokovic
# Social Networks Lab
# Department of Humanities, Social and Political Sciences
# ETH Zurich
# 14 October 2024
#
# Assignment 1 - Task 2

# computeStats ------------------------------------------------------------------
#' Calculates the three statistics: edge count, reciprocal edges, star2
#' and outputs them as a vector.
#'
#' @param net adiacency matrix (network) to compute
#'
#' @return a vector containing the statistics (edge count, reciprocal edges, star2)
computeStats = function(net)
{
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Edge count
  # Not yet correct must remove reciprocal edges as they got double counted
  stat1 = sum(net);
  
  # Sum of matrix representing reciprocal ties. Warning, Hadamard product not matrix prod.
  stat2 = sum(net * t(net)) / 2;
  # Now stat1 is correct
  stat1 = stat1 - stat2;
  
  # Vector where each row is the indegree of the correspoding node
  # trnaspose is necessary to ensure each row is a receiver instead of a sender
  indegree = as.vector(t(net) %*% rep(1,dim(net)[1]));
  # New vector that filters only indegrees of 2
  star_2 = ifelse(indegree[1:nvertices] == 2, 1, 0);
  stat3 = sum(star_2);
  
  return(c(stat1,stat2,stat3))
}


# MHstep ------------------------------------------------------------------
#' Simulate the next step of a network in Markov chain using Metropolis-Hasting
#' 
#' The function `MHstep` simulates the the Metropolis-Hastings step that defines
#' the Markov chain whose stationary distribution is the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param theta statistics vector, for task 2 this are the following values:
#' edge count, reciprocal tie count, and 2-istar count.
#'
#' @return next state of the Markov Chain
MHstep <- function(net, theta){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace = FALSE
  tie <- sample(1:nvertices, 2, replace = FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  new_net = net;
  new_net[i,j] = ifelse(new_net[i,j]== 0, 1, 0);
  
  # Compute the change statistics
  
  current_stat = computeStats(net);
  new_stat = computeStats(new_net);
  delta_stat = new_stat - current_stat;
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  p = min(1,exp(sum(theta * delta_stat)));
  
  # Select the next state: 
  # Return the next state of the chain
  
  r = runif(1);
  if(r <= p)
  {
     next_net = new_net;
  }
  else
  {
    next_net = net;
  }
  return(next_net)
}

# Markov Chain simulation -------------------------------------------------
#' The function MarkovChain simulates the networks from the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param theta statistics vector, for task 2 this are the following values:
#' edge count, reciprocal tie count, and 2-istar count.
#' @param burnin an integer value.
#'   Number of steps to reach the stationary distribution.
#' @param thinning an integer value. Number of steps between simulated networks.
#' @param nNet an integer value. Number of simulated networks to return as output.
#'
#' @return a named list:
#'   - netSim: an `array` with the adjancency matrices of the simulated networks.
#'   - statSim: a `matrix` with the value of the statistic defining the ERGM.
MarkovChain <- function(
    net,
    theta,
    burnin = 10000, thinning = 1000, nNet = 1000){
  
  # Burnin phase: repeating the steps of the chain "burnin" times  
  nvertices <- nrow(net)
  burninStep <- 1 # counter for the number of burnin steps
  
  net_t = net;
  # Perform the burnin steps
  for (i in 1:burnin)
  {
    net_t = MHstep(net_t,theta);
  }
  
  # After the burnin phase we draw the networks
  # The simulated networks and statistics are stored in the objects
  # netSim and statSim
  netSim <- array(0, dim = c(nvertices, nvertices, nNet))
  statSim <- matrix(0, nNet, 3)
  thinningSteps <- 0 # counter for the number of thinning steps
  netCounter <- 1 # counter for the number of simulated network
  
  for(i in 1:nNet)
  {
    net_t = MHstep(net_t,theta);
    netSim[,,i] = net_t;
    statSim[i,] = computeStats(net_t);
  }
  
  # Return the simulated networks and the statistics
  return(list(netSim = netSim, statSim = statSim))
}
