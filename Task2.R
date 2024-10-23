# task_2_2 ------------------------------------------------------------------
#' Executes task 2.2
#'
#' @param advice advice adiancency matrix
#'
#' @return nothing, it just prints out the evaluated statistics (see evalStat)
task2_2 = function(advice)
{
  outdegree = rowSums(advice);
  indegree = colSums(advice);
  par(mfrow=c(2,1),mar=c(4,3,1,3));
  hist(outdegree, xlab="Outdegree",col="grey", main="", breaks = 21);
  hist(indegree, xlab="Indegree",col="grey", main="", breaks = 21);
  
  iterations = 100
  
  obs_stats = computeStats(advice)
  n = dim(advice)[1]
  
  #' Task 2.2 the given parameters deviate too much from the observed values:
  #' observed  mean   std_dev p_value
  # 1      190 32.22  6.606669       0
  # 2       90  6.32  3.637293       0
  # 3      930 25.86 11.646320       0
  #' the main issue seems to lie in the first parameter, as way less edges
  #' appear than the observed quantity, as such for task 2.3 it makes sense
  #' to increase that parameter, this will also inevitably increase the amount
  #' of 2-istars so we keep the third parameter constant or only adjust it
  #' slightly as we see fit.
  result = MarkovChain(matrix(0,n,n),c(-2.76,0.68,0.05), nNet = iterations)
  
  out = evalStat(result$statSim,obs_stats)
  print(out)
}

#' Task 2.3 the estimates for the parameters are obtained by trying random
#' values until it improved the simulated results. At the moment I have
#' obtained better results with 0,0,0 compared than -2.76,0.68,0.05.
#' For sure the third parameter needs to be very small as 2 istars grow quickly,
task2_3 = function(advice, theta1, theta2, theta3)
{
  iterations = 20
  
  obs_stats = computeStats(advice)
  n = dim(advice)[1]
  
  result = MarkovChain(matrix(0,n,n),c(theta1,theta2,theta3), nNet = iterations)
  
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

  #Calculate mean and p value, though I am not sure it works
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
  for(i in 1:s) dev[i] = sd(stats[,i])
  
  return(data.frame(observed = observed, mean = mean, std_dev = dev, p_value = p))
}
