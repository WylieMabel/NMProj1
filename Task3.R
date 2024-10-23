
plotNetwork = function(net, nodes)
{
  net = network(friendship,directed=TRUE);
  net %v% "level" = nodes$nodeLevel;
  net %v% "department" = nodes$nodeDepartment;
  
  ggraph(net) +
    geom_edge_link0(
                    arrow = arrow(
                      angle = 10,
                      length = unit(4, "mm"),
                      type = "closed"
                    )
    ) +
    geom_node_point(
      size = 5,
      aes(
        shape = as.factor(level),
        fill = as.factor(department)
      )
    ) +
    scale_fill_discrete(
      name = "Department"
    ) +
    guides(fill = guide_legend(
      override.aes = list(shape = 21)
    )) +
    scale_shape_manual(
      name = "Level",
      values = c("1" = 21, "2" = 22, "3" = 23),
      labels = c("1" = "CEO", "2" = "Vice President", "3" = "Manager")
    );
}


task3_1 = function(friendship, nodes)
{
  net = network(friendship,directed=TRUE);
  net %v% "level" = nodes$nodeLevel;
  net %v% "department" = nodes$nodeDepartment;
  
  outdegree = rowSums(friendship);
  indegree = colSums(friendship);
  par(mfrow=c(2,1),mar=c(4,3,1,3));
  hist(outdegree, xlab="Outdegree",col="grey", main="", breaks = 21);
  hist(indegree, xlab="Indegree",col="grey", main="", breaks = 21);
  
  model = ergm(net ~ edges + nodematch('department'));
  print(summary(model))
  
  theta_1 = coef(model)[1]
  theta_2 = coef(model)[2]
  
  # Odds of tie existing given it is not homophile
  odds1 = exp(theta_1)
  # Probability
  prob1 = exp(odds1) / (1+exp(odds1))
  print(paste("Tie existing given it is not going homophile.","Odds:",odds1,"Probability:",prob1))
  
  # Odds of a tie being homophile given it exists
  odds2 = exp(theta_2)
  # Probability
  prob2 = exp(odds2) / (1+exp(odds2))
  print(paste("Tie being homophile given it already exists.","Odds:",odds2,"Probability:",prob2))
  
  # Odds of a tie existing given it is homophile
  odds12 = exp(theta_1 + theta_2)
  # Probability
  prob12 = exp(odds12) / (1+exp(odds12))
  print(paste("Tie existing given it is homophile","Odds:",odds12,"Probability:",prob12))
  
}

task3_234 = function(friendship, nodes)
{
  net = network(friendship,directed=TRUE)
  net %v% "level" = nodes$nodeLevel
  net %v% "department" = nodes$nodeDepartment
  
  # Task 3.2
  #' Motivation for statistics choice: as asked in task 3.1 we are carrying over the edges and
  #' nodematch('department') statistics, which measure edge count and department homophile
  #' edges respectively. The new additions are mutual, ttriple/gwesp, and istar(2). Mutual counts
  #' the number of mutual ties, this tests for hypothesis I. Ttriple (transitive triple) and 
  #' gwesp (geometrically weighted edgewise shared partner distribution) are both used to count
  #' the number of transitive triangles, ttriple will not be discarded due to poor convergence
  #' (see task 3.3). This statistics are used to test hypothesis II. For hypothesis III we are
  #' using the istar(2) statistic, which counts the amount istar of degree 2. As pointed out
  #' in the assignment description for task 2 and by [ERGMs (Koskinen & Daraganova, 2013) Datei]
  #' a skewed indegree distribution that favours high degree nodes will lead to a higher count of
  #' istar(2).
  #'
  # Does not converge (near degeneracy)
  # model1 = ergm(net ~ edges + nodematch('department') + mutual + ttriple + istar(2))
  # Does converge
  model = ergm(net ~ edges + nodematch('department') + mutual + gwesp(decay=0.3,fixed=TRUE) + istar(2))
  print(summary(model))
  
  #' Comments of convergence (3.3): as seen in the diagnostic image the trace of the statistics
  #' fluctuate around zero and the density approximates a normal distribution for all
  #' used statistics. This implies that the model converged but not necessarly that it
  #' is a good fit, for that GOF is going to be necessary (task 3.4). This results are
  #' only obtained when using the gwesp statistic, its counterpart ttriple does not mix
  #' well and prevents the model from converging. This as seen in class is most likely
  #' due to the phenomenon of near degeneracy. Gwesp counteracts this effect.
  mcmc.diagnostics(model)
  
  #' Comments on GOF (3.4): as seen in the GOF plots all the custom defined statistics as well
  #' as the four auxiliary statistics provided by the statnet library gof function 
  #' (odegree, idegree, edge-wise shared partners and minimum geodesic distance) fall within
  #' the lower and upper quantiles therefore showing a good fit. However the mean of the odegree
  #' does deviate in some cases from the observed one. While it does not leave the acceptable
  #' region it does point towards a discrepancy from the observed network. A possible origin for
  #' this discrepancy is the skewed distribution of the odegree in the observed network (see image)
  #' lower degree nodes seem to appear far more often than all others. This is even more evident
  #' when comparing it with the indegree distribution, it will be discussed in more detail in task 3.5.
  modelgof = gof(model)
  par(mfrow=c(2,3))
  plot(modelgof)
}

task3_5 = function(friendship, nodes)
{
  nvertices = dim(nodes)[1]
  superiority = matrix(0,nvertices,nvertices)
  for( i in 1:nvertices)
  {
    for( j in  1:nvertices)
    {
      superiority[i,j] = -nodes$nodeLevel[i] + nodes$nodeLevel[j];
    }
  }
  
  net = network(friendship,directed=TRUE)
  net %v% "level" = nodes$nodeLevel
  net %v% "department" = nodes$nodeDepartment
  
  model = ergm(net ~ edges + nodematch('department') + mutual +
                  gwesp(decay=0.3,fixed=TRUE) + istar(2) + edgecov(superiority))
  print(summary(model))
  mcmc.diagnostics(model)
  modelgof = gof(model)
  par(mfrow=c(2,3))
  plot(modelgof)
  
  #' Task 3.5
  #' Based on the results from the model of task 3.2 we can see a significant
  #' estimate for department homophily and transitive triangles hinting that
  #' both have an impact on tie formation, both of which improve the odds of
  #' a tie forming.
  #' Mutual ties while still show some significance but it is lower than almost
  #' all other statistics, possibly showing that reciprocity, while present, is
  #' comparatively less impactful. 
  #' This is particularly strange, as reciprocity has been observed
  #' to be particularly significant in most friendship
  #' networks and our results go against this pattern, possibly hinting at an
  #' additional mechanic that is affecting the formation of mutual ties, this 
  #' will be discussed in a second, first we want to observe the last estimate
  #' tied to istars.
  #' Our result show no significant effect and a value close to
  #' zero, this could imply no particular skew in the indegree distribution,
  #' which does indeed match our observations (see friendship_degree_distribution
  #' image). 
  #' The odds of a tie forming if not helped by any of effects mentioned before
  #' (homophily, reciprocity, transitivity) seem to be low as shown by the
  #' negative estimate of the edge statistic.
  #' Upon visual observation of the network (see friendship network plot)
  #' one can observe that ties that go from a subordinate to their superior
  #' (manager to vice president and vice president to CEO) are almost
  #' never reciprocated. One potential real-life explanation for this could
  #' be that while some subordinates might like their superior and consider
  #' them a friend, the superior might not want to/ have the chance to
  #' familiarize with a subordinate because they want to remain professional
  #' or might be/feel detached from the subordinates. To explore this option
  #' a new edge covariate is introduced: "superiority" which describes the
  #' difference in level between two nodes: a tie from a vice president to a manager
  #' has a superiority of one, or -1 in the opposite direction. Same applies
  #' in a tie from the CEO to a vice president and vice-versa. A tie from the CEO
  #' to a manager has a superiority of 2, or -2 in the opposite direction.
  #' We added this covariate to the model and performed the fit once again.
  #' The diagnostics show good mixing in all statistics, and goodness of fit
  #' also is within the acceptable range of all previous statistics as well
  #' as the newly introduced one.
  #' In the results we can see that superiority shows significance and that
  #' reciprocity has become more significant while the other values stayed
  #' nearly unchanged. Additionally the superiority parameter is negative, as
  #' expected, meaning that a positive level difference (superior to subordinate)
  #' reduces the odds of a tie existing. The mutual parameter has also increased
  #' which also makes sense as we are now controlling the level difference
  #' with the superiority statistic so we can now isolate the two phenomena.
  #' All of this does support our original theory that superiority affects
  #' reciprocity in this network.
}
