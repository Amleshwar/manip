library(igraph)
library(optparse)
library(Matrix)
option_list = list(
  make_option(c("-n", "--vertices"), type="integer", default=NULL, 
              help="number of vertices", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.pdf", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$vertices)){
  print(opt$vertices)
  print_help(opt_parser)
  stop("At least one argument must be supplied integer", call.=FALSE)
}

n = opt$vertices
mat = matrix(0, nrow = n, ncol = n)
qualities = vector("integer", n)
f1 = ceiling(0.05*n)
f2 = ceiling(0.15*n)
f3 = ceiling(0.3*n)
f4 = ceiling(0.5*n)

qualities[1:f1] = runif(f1, 0.9,1)
qualities[(f1 + 1):(f1 + f2)] = runif(f2 , 0.7,0.9)
qualities[(f1+f2+1):(f1+f2+f3)] = runif(f3 , 0.5, 0.7)
qualities[(f1+f2+f3+1): n] = runif((n - (f1+f2+f3)) , 0, 0.5)                                        
for(i in 1:n)
  for(j in 1:n)
    if(i != j)
      mat[i,j] = rbinom(1,1, qualities[j])


rand_net = graph_from_adjacency_matrix(mat, mode = "directed") 

mat_coal = matrix(0, nrow = n, ncol = n)
for(i in 1:n)
  for(j in 1:n)
    if(i != j)
      mat_coal[i,j] = rbinom(1,1, qualities[j])

ncom = ceiling(runif(1, 1, 0.15*n))
total = c(1:n)
community = list()
for(i in 1:ncom)
{
  community[[i]] = sample(total, ceiling(runif(1, 2, (0.15*n))))
  total = total[-community[[i]]]
}

for(i in 1:length(community))
  for(j in 1:length(community[[i]]))
    for(k in 1:length(community[[i]]))
      if(community[[i]][j] != community[[i]][k])
        mat_coal[community[[i]][j], community[[i]][k]] = 1

coal_net = graph_from_adjacency_matrix(mat_coal, mode = "directed")
l <- layout_in_circle(coal_net)
plot(coal_net, layout = l)

####Symmetrization


##########################
rand_clust = cluster_walktrap(rand_net)
rand_clust$membership
coal_clust = cluster_walktrap(coal_net)
coal_clust$membership
modularity(rand_clust)
modularity(coal_clust)
rand_mod  = modularity_matrix(rand_net, rand_clust$membership)
coal_mod = modularity_matrix(coal_net, coal_clust$membership)


credit_rand = degree(rand_net, mode = "in")
credit_coal = degree(coal_net, mode = "in")
total = c(1:n)
####SUBSTITUTE 10 WITH BETA/2M
for(i in 1:length(table(coal_clust$membership)))
{
  com = total[coal_clust$membership == i]
  if(length(com) != 1)
  {
    for(k in 1:length(com))
      credit_coal[com[k]] = credit_coal[com[k]] - 0.0001*(sum(coal_mod[com[k], com[-k]]) + sum(coal_mod[com[-k], com[k]]))
    
  }
  
}


for(i in 1:length(table(rand_clust$membership)))
{
  com = total[rand_clust$membership == j]
  if(length(com) != 1)
    for(k in 1:length(com))
      credit_rand[com[k]] = credit_rand[com[k]] + 0.01*(sum(rand_mod[com[k], com[-k]]) + sum(rand_mod[com[-k], com[k]]))
}
write.csv(credit_rand, "score_random.csv" , row.names = F)
write.csv(credit_coal, "score_coalition.csv", row.names = F)

