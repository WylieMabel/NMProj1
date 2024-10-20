# task_2_2 ------------------------------------------------------------------
#' Executes task 2.2
#'
#' @param advice advice adiancency matrix
#'
#' @return nothing, it just prints out the evaluated statistics (see evalStat)
task_2_2 = function(advice)
{
  iterations = 1000
  
  obs_stats = computeStats(advice)
  n = dim(advice)[1]
  
  result = MarkovChain(matrix(0,n,n),c(-2.76,0.68,0.05), statFunc = computeStats, nNet = iterations)
  
  out = evalStat(result$statSim,obs_stats)
  print(out)
}

# evalStat ------------------------------------------------------------------
#' Calculates the mean, standard deviation, p_values, and t_ratios from a
#' given matrix of statistics vectors, each row should be a statistic from a
#' simulated network.
#'
#' @param stats the statistics matrix collected from the simulation
#' @param observed the statistics observed in the real network (should be a vector)
#'
#' @return a data frame containing the following columns in this order: observed,
#' mean, std_dev, p_value, t_ratio
evalStat = function(stats, observed)
{
  iter = dim(stats)[1];
  s = dim(stats)[2];
  mean = vector(length = s)
  dev = vector(length = s)
  p = vector(length = s)
  t = vector(length = s)
  
  #Calculate mean and p value
  for(i in 1:iter)
  {
    st = stats[i,];
    for(j in 1:s)
    {
      stat = st[j];
      mean[j] = mean[j] + stat;
      if(observed[j] <= stat) p[j] = p[j] + 1;
    }
  }
  mean = mean / iter;
  p = p / iter;
  
  #Calculate standard deviation
  for(i in 1:iter)
  {
    st = stats[i,];
    for(j in 1:s) dev[j] = dev[j] + abs(mean[j] - st[j])^2;
  }
  dev = sqrt(dev / iter);
  
  #Calculate t-ratios
  for(i in 1:s) t[i] = abs(observed[i] - mean[i]) / dev[i];
  
  return(data.frame(observed = observed, mean = mean, std_dev = dev, p_value = p, t_ratio = t))
}
