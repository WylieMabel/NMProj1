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
  # Unit vector of size equal to matix row count
  unit = rep(1,nvertices);
  # Not yet correct must remove reciprocal edges as they got double counted
  stat1 = sum((net %*% unit) * unit);
  
  # Matrix representing reciprocal ties. Warning, Hadamard product not matrix prod.
  recip = net * t(net)
  # Symmetric matrix, /2 to avoid double count
  stat2 = sum((recip %*% unit) * unit) / 2;
  # Now stat1 is correct
  stat1 = stat1 - stat2;
  
  # Vector where each row is the indegree of the correspoding node
  # trnaspose is necessary to ensure each row is a receiver instead of a sender
  indegree = as.vector(t(net) %*% unit);
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
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the istar(2) parameter of the ERGM.
#'
#' @return next state of the Markov Chain
#'
#' @examples
#' MHstep(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   -log(0.5), log(0.4), log(.8)
#' )
MHstep <- function(net, theta1, theta2, theta3){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace = FALSE
  tie <- sample(1:nvertices, 2, replace = FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  new_net = net;
  new_net[i,j] = ifelse(new_net[i,j]== 0, 1, 0);
  
  # Compute the change statistics
  
  theta = c(theta1,theta2,theta3);
  current_stat = computeStats(net);
  new_stat = computeStats(new_net);
  delta_stat = new_stat - current_stat;
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  p = min(1,exp(sum(theta * delta_stat)));
  
  # Select the next state: 
  # Return the next state of the chain
  
  r = runif(1);
  if (r < p) return(new_net)
  else return(net)
}

# Markov Chain simulation -------------------------------------------------
#' The function MarkovChain simulates the networks from the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the istar(2) parameter of the ERGM.
#' @param burnin an integer value.
#'   Number of steps to reach the stationary distribution.
#' @param thinning an integer value. Number of steps between simulated networks.
#' @param nNet an integer value. Number of simulated networks to return as output.
#'
#' @return a named list:
#'   - netSim: an `array` with the adjancency matrices of the simulated networks.
#'   - statSim: a `matrix` with the value of the statistic defining the ERGM.
#'
#' @examples
#' MarkovChain(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   -log(0.5), log(0.4), log(.8)
#' )
MarkovChain <- function(
    net,
    theta1, theta2, theta3,
    burnin = 10000, thinning = 1000, nNet = 1000){
  
  # Burnin phase: repeating the steps of the chain "burnin" times  
  nvertices <- nrow(net)
  burninStep <- 1 # counter for the number of burnin steps
  
  # Perform the burnin steps
  for (i in 1:burnin)
  {
    net = MHstep(net,theta1,theta2,theta3);
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
    net = MHstep(net,theta1,theta2,theta3);
    netSim[,,i] = net;
    statSim[i] = computeStats(net);
  }
  
  # Return the simulated networks and the statistics
  return(list(netSim = netSim, statSim = statSim))
}
