# task3Stat ------------------------------------------------------------------
#' Computes all of the necessary statistics for task 3 from the given network
#'
#' @param net network to evaluate
#'
#' @return nothing, it just prints out the evaluated statistics (see evalStat)
task3Stat = function(net)
{
  recip = sum(net * t(net)) / 2;
  edges = sum(net) - recip;
  
  return(c(edges,recip))
}

# covariance ------------------------------------------------------------------
#' Calculates the covariance matrix
#'
#' @param statsSim simulated statistics
#' @param mean mean of simulated statistics
#'
#' @return covariance matrix
covariance = function(statSim, mean)
{
  iter = dim(statSim)[1];
  s = dim(statSim)[2];
  
  m = matrix(0,s,s);
  
  for(i in 1:s)
  {
    for(j in 1:s)
    {
      m_i = mean[i];
      m_j = mean[j];
      cov = 0;
      
      for(n in 1:iter) 
        cov = cov + (statSim[n,i] - m_i) * (statSim[n,j] - m_j);
      
      cov = cov / iter;
      m[i,j] = cov;
    }
  }
  return(m)
}

# MCMC ------------------------------------------------------------------
#' Perform the Monte Carlo Markov Chain estimation, the statistics used are
#' defined in `task3Stat`.
#'
#' @param obs_net the observed network
#' @param iter number of iterations before stopping
#' @param chain_iter how many iterations for each chain (on top of burnin)
#' @param burnin burnin iterations
#' @param thinning still don't know what it is
#'
#' @return data frame with final iteration statistics and parameters
MCMC = function(obs_net, iter, chain_iter = 1000, burnin = 10000, thinning = 1000)
{
  observed = task3Stat(obs_net);
  s = length(observed);
  n = dim(obs_net)[1];
  
  # Start net
  net = matrix(0,n,n);
  #Parameters
  theta = vector(length = s);
  # Gain
  alpha_0 = 0.01;
  alpha = alpha_0;
  # Evaluated statistics
  eval = 0;
  # Sensitivity matrix
  D = diag(s);
  
  for(i in 1:iter)
  {
    print(paste("Iteration", toString(i), ". Parameters:", toString(theta), sep = " "));
    result = MarkovChain(net = net, theta = theta, statFunc = task3Stat, burnin = burnin,
                         thinning = thinning, nNet = chain_iter);
    eval = evalStat(result$statSim, observed = observed);
    
    # Avoid updating theta at last iterations otherwise the evaluated statistics
    # will not match the parameters
    if(i != iter)
    {
      # Robbins-Monro algorithm
      theta = theta - alpha * (D %*% (eval$mean - observed));
      alpha = alpha_0 / (i + 1);
      D = solve(covariance(result$statSim, eval$mean));
    }
  }
  eval$theta = theta;
  
  return(eval)
}
