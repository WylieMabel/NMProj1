task4_1 = function(friendship, nodes)
{
  par(mfrow=c(1,1),mar=c(4,3,1,3));
  hist(nodes$nodeAge, xlab="Age",col="grey", main="", breaks = max(nodes$nodeAge));
  
  net = network(friendship,directed=TRUE)
  net %v% "level" = nodes$nodeLevel;
  net %v% "department" = nodes$nodeDepartment
  net %v% "age" = nodes$nodeAge
  
  #' ERGM parameter explanation
  #' nodematch tests for homophily
  #' absdiff tests for the absolute age difference
  #' nodeocov for the effect of age on the outdegree, not to be confused with nodeofactor
  #' which does a 'categorical' test, that means: it tests for example being specifically
  #' 27,28,29.. has an effect on tie formation, it creates one paramater for each age group
  #' nodeocov is 'continuous' and tests for the effect on tie formation when increasing
  #' the age no matter where you started from.
 
  model0 = ergm(net ~ edges + nodematch('department') + absdiff('age') + nodeocov('age'))
  print(summary(model0))
  
  #' From Task 1.3 MR-QAP
  #'                 Estimate    Exp(b)    Pr(<=b) Pr(>=b) Pr(>=|b|)
  #' intercept       -1.97460863 0.1388156 0.0056  0.9944  0.0056  *
  #' adviceMatrix     0.78779191 2.1985365 0.9822  0.0178  0.0292  * 
  #' sameDep          1.18388201 3.2670323 0.9996  0.0004  0.0004  **
  #' senioritySender  0.04876536 1.0499740 0.8688  0.1312  0.2336   
  #' ageDiff         -0.04722080 0.9538768 0.0340  0.9660  0.0954   
  #' 
  #' By comparing the results of the ERGM fit (see results.txt) we can see that hypothesis
  #' conincides with both methods, with homophily having a positive impact on tie formation
  #' likelyhood and being significant, even more so in the ERGM. Hypothesis 3 is also unchanged
  #' with the estimate being negative values close to zero and showing little significance
  #' therefore the null hypothesis is not rejected even with the ERGM. The differences arise
  #' with hypothesis 2: the ERGM estimates a small negative impact on tie formation with increasing
  #' age, while not extremely significant it still results in a p-value below 0.05 and therefore
  #' cannot be completely ignored. This is inconsistent with the previous results obtained with
  #' MR-QAP which showed an insignificant positive estimate. Another discrepancy appears in the
  #' baseline odds of forming a tie, in the ERGM this corresponds to the edges parameter, while
  #' in MR-QAP this corresponds to the intercept. We can see a large difference in the estimate
  #' that implies that ties are much less likely to occur unaided (no help from homophily or other
  #' phenomena) in the MR-QAP estimate compared to the ERGM.
   
  model1 = ergm(net ~ edges + nodematch('department') + absdiff('age') + nodeocov('age') +
            mutual + gwesp(decay=0.3,fixed=TRUE) + istar(2))
  print(summary(model1))
  
  #'After this second we that the results match the ones of the MR-QAP much more closely:
  #'all the results concerning the hypotheses in task 1.3 correspond with the estimates
  #'of ERGM, the age difference and age contribution to the outdegree both show small
  #'estimates close to zero with no significance therefore we cannot reject the null hypotheses.
  #'Department homophily remained a positive and very significant contribution much like in the
  #'MR-QAP. Additionally the baseline odds of forming a tie (edges for ERGM and intercept for MR-QAP)
  #'are now much more closer to each other. One possible explanation for the original discrepancy
  #'is that the ERGM does not take into account any effect that has not been specified in the model,
  #'as such, a simplistic ERGM will not succeed in modeling the original network. MR-QAP on the other
  #'hand work by switching the labels of nodes without changing the structure of the network, as such
  #'structures like reciprocal ties and transitive triangles will persist and affect the final result.
  #'This was however resolved as soon as we introduced the missing structures in the second ERGM.
}
