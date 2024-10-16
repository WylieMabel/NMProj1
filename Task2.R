# task_2_2 ------------------------------------------------------------------
#' Executes task 2.2
#'
#' @param n number of vertices
#' @param advice advice adiancency matrix
#' @param friendship friendship adiancency matrix
#' @param report_to report_to adiancency matrix
#'
#' @return nothing, it just prints out the results at the moment, it first prints
#' out the observed stats (edges,reciprocal,star2 in this order), then it outputs
#' the average of those stats across all the simulated networks, right after stat
#' the corrsponding p value is also printed out.
#' TODO better output formatting
task_2_2 = function(n, advice, friendship, report_to)
{
  iterations = 1000
  labels = c("ADVICE","FRIENDSHIP","REPORTS_TO")
  simulate = list(advice,friendship,report_to)
  
  for(i in 1:length(simulate))
  {
    obs_stats = computeStats(simulate[[i]])
    print(c("OBSERVED",labels[i],obs_stats[1],obs_stats[2],obs_stats[3]))
    result = MarkovChain(matrix(0,n,n),-2.76,0.68,0.05, nNet = iterations)
    sim_edges = 0
    p_edges = 0
    sim_reciprocal = 0
    p_reciprocal = 0
    sim_star_2 = 0
    p_star_2 = 0
    for(j in 1:iterations)
    {
      stat = result$statSim[j,]
      sim_edges = sim_edges + stat[1]
      if(obs_stats[1] <= stat[1]) p_edges = p_edges + 1
      sim_reciprocal = sim_reciprocal + stat[2]
      if(obs_stats[2] <= stat[2]) p_reciprocal = p_reciprocal + 1
      sim_star_2 = sim_star_2 + stat[3]
      if(obs_stats[3] <= stat[3]) p_star_2 = p_star_2 + 1
    }
    sim_edges = sim_edges / iterations
    p_edges = p_edges / iterations
    sim_reciprocal = sim_reciprocal / iterations
    p_reciprocal = p_reciprocal / iterations
    sim_star_2 = sim_star_2 / iterations
    p_star_2 = p_star_2 / iterations
    
    print(c("SIMULATED",labels[i],sim_edges,p_edges,sim_reciprocal,p_reciprocal,sim_star_2,p_star_2))
  }
}
