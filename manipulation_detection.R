setwd("~/aamas/manip/data")

library(igraph)
e1 = read.csv("socialnetwork.txt", header  =  FALSE, sep = ' ')
e2 = read.csv("socialnetwork2.txt", header = FALSE, sep = ' ')
e1 = as.matrix(e1)
sn1 = graph_from_edgelist(e1, directed = TRUE)
e2 = as.matrix(e2)
sn2 = graph_from_edgelist(e2, directed = TRUE)
er <- sample_gnm(n=53, m=265, directed = TRUE)# Erdos Renyi grpah with same number of nodes and edges as of sn2
plot(sn2, edge.arrow.size=0.3, vertex.size = 3, edge.curved = 0.1, layout = layout_with_fr(sn2), vertex.label = NA)
plot(er, edge.arrow.size=0.3, vertex.size = 3, edge.curved = 0.1, layout = layout_with_fr(er), vertex.label = NA)
plot(er, edge.arrow.size=0.3, vertex.size = 3, edge.curved = 0.1, layout = layout_in_circle(er))
plot(sn2, edge.arrow.size=0.3, vertex.size = 3, edge.curved = 0.1, layout = layout_in_circle(sn2), vertex.label = NA)


layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(sn2)) 
  plot(sn2, edge.arrow.mode=0,edge.arrow.size=0.3, vertex.size = 3, layout=l, main=layout, vertex.label = NA)
   }
#Calculating the reciprocity of The proportion of reciprocated ties (for a directed network).
#Classify dyads in a directed graphs. The relationship between each pair of vertices is measured. It can be in three states: mutual, asymmetric or non-existent.
# The measure of reciprocity defines the proporsion of mutual connections, in a directed graph.
# It is most commonly defined as the probability that the opposite counterpart of a 
# directed edge is also included in the graph. Or in adjacency matrix notation:
#   sum(i, j, (A.*A')ij) / sum(i, j, Aij), where A.*A'  is the element-wise product of matrix A and its transpose. 
#   This measure is calculated if the mode argument is default.
#                                                                                                                                                                                                                                                                         
reciprocity(sn1)
dyad_census(sn1) # Mutual, asymmetric, and nyll node pairs
2*dyad_census(sn1)$mut/ecount(sn1) # Calculating reciprocity

reciprocity(sn2)
dyad_census(sn2) # Mutual, asymmetric, and nyll node pairs
2*dyad_census(sn2)$mut/ecount(sn2) # Calculating reciprocity

reciprocity(er)
dyad_census(er) # Mutual, asymmetric, and nyll node pairs
2*dyad_census(er)$mut/ecount(er) # Calculating reciprocity

#Hubs and authorities
hs1 <- hub_score(sn1, weights=NA)$vector
as1 <- authority_score(sn1, weights=NA)$vector
par(mfrow=c(1,2))
plot(sn1, vertex.size=hs1*25, main="Hubs", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(sn1))
plot(sn1, vertex.size=as1*15, main="Authorities", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(sn1))

hs2 <- hub_score(sn2, weights=NA)$vector
as2 <- authority_score(sn2, weights=NA)$vector
par(mfrow=c(1,2))
plot(sn2, vertex.size=hs2*25, main="Hubs", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(sn2))
plot(sn2, vertex.size=as2*15, main="Authorities", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(sn2))


hs_er <- hub_score(er, weights=NA)$vector
as_er <- authority_score(er, weights=NA)$vector
par(mfrow=c(1,2))
plot(er, vertex.size=hs_er*25, main="Hubs", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(er))
plot(er, vertex.size=as_er*15, main="Authorities", vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(er))

net1.sym <- as.undirected(sn1, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))
cliques(net1.sym) # list of cliques       
sapply(cliques(net1.sym), length) # clique sizes
largest_cliques(net1.sym) # cliques with max number of nodes

net2.sym <- as.undirected(sn2, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))
cliques(net2.sym) # list of cliques       
sapply(cliques(net2.sym), length) # clique sizes
largest_cliques(net2.sym) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net2.sym))
vcol[unlist(largest_cliques(net2.sym))] <- "gold"
plot(as.undirected(net2.sym), vertex.label=V(net2.sym)$name, vertex.color=vcol)
plot(as.undirected(net2.sym), vertex.label=V(net1.sym)$name, vertex.color=vcol, vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_with_fr(er))


net_er.sym <- as.undirected(er, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))
cliques(net_er.sym) # list of cliques       
sapply(cliques(net_er.sym), length) # clique sizes
largest_cliques(net_er.sym) # cliques with max number of nodes

######Community detection using edge betweenness(Newman Girvan)
ceb2 <- cluster_edge_betweenness(sn2) 
dendPlot(ceb2, mode="hclust") 
plot(ceb2, sn2, vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_in_circle(er)) 
length(ceb2)
modularity(ceb2)
#Modularity comes out to be 0.09505162
ceb_er <- cluster_edge_betweenness(er) 
dendPlot(ceb_er, mode="hclust") 
plot(ceb_er, er, vertex.label = NA, edge.arrow.size=0.3, vertex.size = 3, layout = layout_in_circle(er)) 
length(ceb_er)
modularity(ceb_er)


#Modularity comes out to be 0.02173015


#########Community detection based on propagating labels
clp2 <- cluster_label_prop(sn2)
plot(clp2, sn2, vertex.label = NA, edge.arrow.size=0.01, vertex.size = 3) 

str(V(sn2)$name)
V(sn2)$id = c(1:53)
str(V(sn1)$name)
V(sn1)$id = c(1:59)
str(V(sn1)$id)
write.graph(sn1, "anon_follow.net" , format = "pajek")
write.graph(sn2, "anon_upvote.net" , format = "pajek")
