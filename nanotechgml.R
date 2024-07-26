rm(list=ls())
library(igraph)
library(tidyverse)
library(readr)
library(sna)
library(influenceR)
library(RColorBrewer)

g1 = read_graph("Nanotech_UTF8.gml", format = "gml")
summary(g1)
vcount(g)
ecount(g)
V(g)
E(g)[sample(1:ecount(g), 10)]

coords = layout.kamada.kawai(g)
plot(g, layout=coords, vertex.label=NA, vertex.size=10)

V(g)$OrgType

#allocating vectors
industry=V(g)$OrgType=="IND"
industry
sum(industry)

education=V(g)$OrgType=="EDU"
education
sum(education)

research_org=V(g)$OrgType=="ROR"
research_org
sum(research_org)



coords = layout.kamada.kawai(g)
plot(g, layout=coords,vertex.color=colrs ,vertex.label=NA, vertex.size=5)
legend(x= 1.5,y= 1.1, c("EDU","IND","ROR","GOV","CON","OTH","RH","PNP"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
##############
#education subgraph
education=V(g)$OrgType=="EDU"
education
sum(education)
h1 = induced.subgraph(g, V(g)[education])
plot(h1, vertex.label=NA,vertex.size=5)
transitivity_g1 <- transitivity(h1, type = "globalundirected")
#write_graph(h1, "NanotechEDU.gml", format = "gml")

#industry subgraph
industry=V(g)$OrgType=="IND"
industry
sum(industry)
h2 = induced.subgraph(g, V(g)[industry])
plot(h2, vertex.label=NA,vertex.size=5)
transitivity_g2<- transitivity(h2, type = "globalundirected")
write_graph(h2, "NanotechIND1.gml", format = "gml")

h2.new <- delete_vertices(h2, V(h2)$degree_w==0)
plot(h2.new, vertex.label=NA,vertex.size=5)
write_graph(h2.new, "NanotechIND2.gml", format = "gml")
#research organisation subgraph

research_org=V(g)$OrgType=="ROR"
research_org
sum(research_org)

h3 = induced.subgraph(g, V(g)[research_org])
plot(h3, vertex.label=NA,vertex.size=5)
transitivity_g3 <- transitivity(h3, type = "globalundirected")
#write_graph(h3, "NanotechROR2.gml", format = "gml")


#EDU IND RES Subgraph
h4 = induced.subgraph(g, V(g)[industry | education | research_org])
plot(h4, vertex.label=NA,vertex.size=5)
write_graph(h4, "NanotechERI1.gml", format = "gml")
summary(h4)

legend(x= 1.5,y= 1.1, c("EDU","IND","ROR","GOV","CON","OTH","RH","PNP"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
write_graph(h, "nanotechF1.gml", format = "gml")

#EDu IND
h5 = induced.subgraph(g, V(g)[industry | education ])
plot(h5, vertex.label=NA,vertex.size=5,vertex.color=c("tomato","lightblue"))
write_graph(h5, "NanotechEI.gml", format = "gml")
transitivity_g5 <- transitivity(h5, type = "globalundirected")


#EDU RES
h6 = induced.subgraph(g, V(g)[education | research_org ])
summary(h6)
plot(h6, vertex.label=NA,vertex.size=5,vertex.color=c("lightblue","darksalmon"))
write_graph(h6, "NanotechER.gml", format = "gml")
transitivity_g6 <- transitivity(h6, type = "globalundirected")


#IND RES
h7 = induced.subgraph(g, V(g)[industry | research_org ])
plot(h7, vertex.label=NA,vertex.size=5,vertex.color=c("tomato","darksalmon"))
write_graph(h7, "NanotechIR.gml", format = "gml")


V(h2)$degree_w <- strength(h2, mode = "all", loops = F)
V(h4)$degree_w

V(h4)$closeness_w <- igraph::closeness(h4)
V(h4)$closeness_w

V(h4)$betweenness_w <- betweenness(h4)
V(h4)$betweenness_w

hist(igraph::degree(h4), col= "pink",
     xlab = "degree",
     ylab = "frequency",
     main = "")

hist(betweenness(h4), col= "light blue",
     xlab = "Betweenness",
     ylab = "Frequency",
     main = "")




G <- get.adjacency(h4, sparse = F)                                   #Get the adjacency matrix

br <- sna::brokerage(G, V(h4)$betweenness_w)                                  #Calculate brokerage measures
summary(br)

brraw_nli_csv <- write.csv(br$raw.nli,"br_raw1nanotech.csv")
centr_degree(h4, mode = "total", loops = F)                         #Degree centralisation  
centr_clo(h4, mode = "total")                                       #Closeness centralisation
centr_betw(h4, directed = F)                                        #Betweenness centralisation


d_g <- diameter(h5, directed = FALSE, unconnected = FALSE)          #diameter
apl_g  <- mean_distance(h4, directed = FALSE, unconnected = FALSE)  #APL
ed_g <- edge_density(h4)                                            #Calculate density
cp_g <- articulation_points(h4)                                     #Cutpoints
cliques_g <- cliques(h4, min = 3)                                   #List of cliques
numcliques_g <- count_max_cliques(h4, min = 3)                      #Number of cliques
transitivity_g <- transitivity(h4, type = "globalundirected")       #Calculate transitivity

numisolates_g <- sum(igraph::degree(h4)==0)                                 #Number of isolates
isolates_g <- V(h4)[igraph::degree(h4)==0]                                   #List of isolates
inclusiveness_g <- (vcount(h4)-numisolates_g)/vcount(h4)

effective_net <- influenceR::ens(h4) 


const <- constraint(h4) 
invConstraint <-1.125 - const
names = V(h4)$Label
d <-data.frame(node.name=names, constraint=const,Inverse_constraint= invConstraint,effectiveNetworksize=effective_net)
write.csv(d,"constraint1.csv")
