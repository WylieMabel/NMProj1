task_2_2 = function(advice, friendship, report_to)
{
  iterations = 1000
  results = MarkovChain(advice,-2.76,0.68,0.05, nNet = iterations)
  chain_edges = sum(results$statSim[,1]) / iterations
  chain_reciprocal = sum(results$statSim[,2]) / iterations
  chain_star_2 = sum(results$statSim[,3]) / iterations
  
  print(c(chain_edges,chain_reciprocal,chain_star_2))
}
