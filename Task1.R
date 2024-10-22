# ------------------------------- Loading data ------------------------------- #
# Setting the WORKING DIRECTORY
# This is done using the function setwd and the path of the folder
# in which data are and results will be saved

getwd()

list.files() # Listing the files in the working directory

# Installing and loading the necessary R packages 
# (installation to be done only the first time you use a library):
# install.packages("sna")
# install.packages("network")
library(sna)
library(network)


# Loading data: 
layersPath <- "Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_layers.txt"
nodesPath <- "Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_nodes.txt"
edgesPath <- "Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_multiplex.edges"
layers <- read.csv(layersPath,header=TRUE, sep = " ")
nodes <- read.csv(nodesPath,header=TRUE, sep = " ")
edges <- read.table(edgesPath,header=FALSE, sep = " ")
edgesCols <- c ("layerID", "sendingNodeID", "recievingNodeID", "weight")
colnames(edges) <- edgesCols

# Setting Adjacency Matrices for each layer
adviceTies <- edges[edges$layerID==1, ]
friendshipTies <- edges[edges$layerID==2, ]
reportsToTies <- edges[edges$layerID==3, ]

adviceMatrix <- matrix(0L, nrow=21, ncol=21)
friendshipTiesMatrix <- matrix(0L, nrow=21, ncol=21)
reportsToTiesMatrix <- matrix(0L, nrow=21, ncol=21)
rownames(adviceMatrix) <- 1:21
colnames(adviceMatrix) <- 1:21
rownames(friendshipTiesMatrix) <- 1:21
colnames(friendshipTiesMatrix) <- 1:21
rownames(reportsToTiesMatrix) <- 1:21
colnames(reportsToTiesMatrix) <- 1:21


for(i in 1:dim(edges)) {
  layer <- edges[i, 1]
  sendingNode <- edges[i, 2]
  recievingNode <- edges[i, 3]
  weightTie <- edges[i, 4]
  if (layer == "1") {
    adviceMatrix[sendingNode, recievingNode] <- weightTie
  } else if (layer == "2") {
    friendshipTiesMatrix[sendingNode, recievingNode] <- weightTie
  } else {
    reportsToTiesMatrix[sendingNode, recievingNode] <- weightTie
  }
}
# 1 #
# ------------------------------ QAP regression ------------------------------ #
# We use the QAP regression to test if lawyers seek out their personal friends 
# for work-related advice
set.seed(18) #To reproduce the results 
permutations <- 5000 # Number of permutations
nl0 <- netlogit(adviceMatrix, friendshipTiesMatrix, rep=permutations, nullhyp="qapy") #permute advice network labels

# Is qapy = qapspp in this case? Why?
nl0b <- netlogit(adviceMatrix, friendshipTiesMatrix, rep=permutations, nullhyp="qapspp")
table(nl0$coefficients == nl0b$coefficients)

# netlm for the linear regression model 
nl0$names <- c("intercept", "friendship")

summary(nl0)

# 2 #
# ---------------------------- MR-QAP regression ----------------------------- #
# We use MR-QAP regression to test whether the following hypotheses are supported
# by the data:
# - Hp. 1: A friendship nomination is more likely between a pair of managers within the same deparment.
# - Hp. 2: Senior managers are less likely to nominate friends.
# - Hp. 3: A friendship nomination is more likely between a pair of managers of a similar age.

# Step 1. Model specification:
# Creating the matrices Z with the information on the explanatory variables.
# Combining these matrices in a list or an array

# Hp 1: Same department
departments <- nodes[,5]
sameDepartment <- abs(outer(departments,departments,"=="))

# Hp 2: Senior Manager (older the manager, less likely)
seniority <- nodes[,3]
senioritySender <- matrix(seniority,21,21,byrow=FALSE)

# Hp 3: Similar Age (as diff increases then less likely)
age <- nodes[,2]
ageDiff <- abs(outer(age,age,"-"))

# The explanatory variables must be combined in a list
zm <- list(sameDepartment, senioritySender, ageDiff)

# 3 #
# Step 2: running the MR-QAP regression
set.seed(18) #To reproduce the results 
permutations <- 5000 # Number of permutations
nl <- netlogit(friendshipTiesMatrix, zm, rep=permutations, nullhyp="qapspp")
nl$names <- c("intercept","sameDep", "senioritySender", "ageDiff")
summary(nl)

# Step 3: Model check and interpretation
# Understanding the empirical p-values
z.values <- rbind(nl$dist,nl$tstat)
p.values <- function(x,permutations){
  sum(abs(x[1:permutations]) > abs(x[permutations+1]))/permutations}
empirical.p.values <- apply(z.values,2,p.values,permutations)
empirical.p.values

# Visualizing the empirical p-values
par(mfrow=c(2,2))
for (i in 1:4)
{
  hist(nl$dist[,i],
       breaks = 30,
       xlim = c(min(c(nl$tstat[i],nl$dist[,i]))-1,
                max(c(nl$tstat[i],nl$dist[,i]))+1),
       main = nl$names[i],xlab="z-values")
  abline(v = nl$tstat[i],col="red",lwd=3,lty=2)
}

# Interpreting the parameters (see the commented output in the lecture notes)

# Step 4: Format and export the results
res <- summary(nl)
expRes <- cbind(res$coefficients, exp(res$coefficients), res$se, res$pgreqabs)
colnames(expRes) <- c("ESt.", "exp(Est.)", "s.e.", "p-value")
rownames(expRes) <- res$names
expRes
write.csv(expRes,"resQAP.csv")

# Exporting results in tex
library(xtable)
xtable(expRes,digits=3)

# 4 #
# Another Hypothesis could be that Nodes on the same level are more likely to be friends #
# (e.g managers with managers, vice pres with vice pres) #

# 5 #

# Hp 4: Same level 
levels <- nodes[,4]
sameLevel <- abs(outer(levels,levels,"=="))

# The explanatory variables must be combined in a list
zm1 <- list(sameDepartment, senioritySender, ageDiff, sameLevel)

# 3 #
# Step 2: running the MR-QAP regression
set.seed(18) #To reproduce the results 
permutations <- 5000 # Number of permutations
nl <- netlogit(friendshipTiesMatrix, zm1, rep=permutations, nullhyp="qapspp")
nl$names <- c("intercept","sameDep", "senioritySender", "ageDiff", "sameLevel")
summary(nl)


# Step 3: Model check and interpretation
# Understanding the empirical p-values
z.values <- rbind(nl$dist,nl$tstat)
p.values <- function(x,permutations){
  sum(abs(x[1:permutations]) > abs(x[permutations+1]))/permutations}
empirical.p.values <- apply(z.values,2,p.values,permutations)
empirical.p.values

# Visualizing the empirical p-values
par(mfrow=c(2,3))
for (i in 1:5)
{
  hist(nl$dist[,i],
       breaks = 30,
       xlim = c(min(c(nl$tstat[i],nl$dist[,i]))-1,
                max(c(nl$tstat[i],nl$dist[,i]))+1),
       main = nl$names[i],xlab="z-values")
  abline(v = nl$tstat[i],col="red",lwd=3,lty=2)
}

# Interpreting the parameters (see the commented output in the lecture notes)

# Step 4: Format and export the results
res <- summary(nl)
expRes <- cbind(res$coefficients, exp(res$coefficients), res$se, res$pgreqabs)
colnames(expRes) <- c("ESt.", "exp(Est.)", "s.e.", "p-value")
rownames(expRes) <- res$names
expRes
write.csv(expRes,"res1_5QAP.csv")

# Exporting results in tex
library(xtable)
xtable(expRes,digits=3)

# Still to display the results and comment on them #

