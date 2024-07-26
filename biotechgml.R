rm(list=ls())
library(igraph)
library(tidyverse)
library(readr)
library(sna)
library(influenceR)
library(RColorBrewer)

g = read_graph("Biotech_UTF8.gml", format = "gml")

summary(g)
vcount(g)
ecount(g)
V(g)
E(g)[sample(1:ecount(g), 10)]

coords = layout.kamada.kawai(g)
plot(g, layout=coords, vertex.label=NA, vertex.size=10,vertex.color=colrs)

title(main = "Collaborations in different sectors in Biotechnology",cex.main=1)
legend(x=-1.5,y=-1.1, c("EDU","IND","ROR","GOV","CON","OTH","RH","PNP"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)



V(g)$OrgType





##subgraphs
#education subgraph for biotech
education=V(g)$OrgType=="EDU"
education
sum(education)
h1 = induced.subgraph(g, V(g)[education])
plot(h1, vertex.label=NA,vertex.size=5,vertex.color="lightblue")
transitivity_g1 <- transitivity(h1, type = "globalundirected")
#write_graph(h1, "BiotechEDU.gml", format = "gml")

#industry subgraph
industry=V(g)$OrgType=="IND"
industry
sum(industry)

h2 = induced.subgraph(g, V(g)[industry])
plot(h2, vertex.label=NA,vertex.size=5)
write_graph(h2, "BiotechIND.gml", format = "gml")
transitivity_g2 <- transitivity(h2, type = "globalundirected")

#res
research_org=V(g)$OrgType=="ROR"
research_org
sum(research_org)
h3 = induced.subgraph(g, V(g)[research_org])
plot(h3, vertex.label=NA,vertex.size=5)
write_graph(h3, "BiotechROR2.gml", format = "gml")
transitivity_g3 <- transitivity(h3, type = "globalundirected")




#EDu IND
h5 = induced.subgraph(g, V(g)[industry | education ])
plot(h5, vertex.label=NA,vertex.size=5,vertex.color=c("tomato","lightblue"))
write_graph(h5, "BiotechEI.gml", format = "gml")

transitivity_g5 <- transitivity(h5, type = "globalundirected")
#EDU RES
h6 = induced.subgraph(g, V(g)[education | research_org ])
summary(h6)
plot(h6, vertex.label=NA,vertex.size=5,vertex.color=c("lightblue","darksalmon"))
write_graph(h6, "BiotechER.gml", format = "gml")
transitivity_g6 <- transitivity(h6, type = "globalundirected")

#IND RES
h7 = induced.subgraph(g, V(g)[industry | research_org ])
plot(h7, vertex.label=NA,vertex.size=5,vertex.color=c("tomato","darksalmon"))
write_graph(h7, "BiotechIR.gml", format = "gml")

#EDU IND ROR subgraph
h4 = induced.subgraph(g, V(g)[industry | education | research_org])
plot(h4, vertex.label=NA,vertex.size=5,vertex.color=c("tomato","lightblue","darksalmon"))
#write_graph(h4, "BiotechERI.gml", format = "gml")
summary(h4)

h4.new <- delete_vertices(h4, V(h4)$degree_w==0)


par(mfrow=c(1, 1), mar=c(0,0,0,0))                      
plot(h4.new, vertex.size=igraph::degree(h4.new)*0.5, vertex.label=NA, layout=layout_with_kk,vertex.color=c("tomato","lightblue","darksalmon"))
title(main = "Collaborations between Industry, Educationa and Research Organisations in Biotechnology",cex.main=1)
legend(x=-1.5,y=-1.1, c("IND","EDU","ROR"), pch=21,
       col="#777777", pt.bg=c("tomato","lightblue","darksalmon"), pt.cex=2, cex=.8, bty="n", ncol=1)
#pdf('gmlrbiotech1.pdf', width = 75, height = 75)
dev.off()

#Centrality measures for weighted graphs node level measures
V(h4)$degree_w <- strength(h4, mode = "all", loops = F)
V(h4)$degree_w

V(h4)$closeness_w <- igraph::closeness(h4.new)
V(h4)$closeness_w

V(h4)$betweenness_w <- betweenness(h4)
V(h4)$betweenness_w

hist(igraph::degree(h4.new), col= "pink",
     xlab = "degree",
     ylab = "frequency",
     main = "Number of Collaborations")
dev.off()
hist(betweenness(h4.new), col= "light blue",
     xlab = "Betweenness",
     ylab = "Frequency",
     main = "")


#brokerage
G <- get.adjacency(h4, sparse = F)                                   #Get the adjacency matrix

br <- sna::brokerage(G, V(h4.new)$betweenness_w)                                  #Calculate brokerage measures
summary(br)

brraw_nli_csv <- write.csv(br$raw.nli,"br_raw1biotech.csv")
centr_degree(h4.new, mode = "total", loops = F)                         #Degree centralisation  
centr_clo(h4.new, mode = "total")                                       #Closeness centralisation
centr_betw(h4.new, directed = F)                                        #Betweenness centralisation


d_g <- diameter(h4.new, directed = FALSE, unconnected = FALSE)          #diameter
apl_g  <- mean_distance(h4, directed = FALSE, unconnected = FALSE)  #APL
ed_g <- edge_density(h4)                                            #Calculate density
cp_g <- articulation_points(h4.new)                                     #Cutpoints
cliques_g <- cliques(h4.new, min = 3)                                   #List of cliques
numcliques_g <- count_max_cliques(h4, min = 3)                      #Number of cliques
transitivity_g <- transitivity(h4.new, type = "globalundirected")       #Calculate transitivity

numisolates_g <- sum(igraph::degree(h4)==0)                                 #Number of isolates
isolates_g <- V(h4)[igraph::degree(h4)==0]                                   #List of isolates
inclusiveness_g <- (vcount(h4)-numisolates_g)/vcount(h4) 


V(h4.new)$constraint <- constraint(h4.new)    

V(h4.new)$effective_net <- influenceR::ens(h4.new) 
#export data in csv
V(h4.new)
summary(h4.new)

const <- constraint(h4) 
invConstraint <-1.125 - const
names = V(h4)$Label
d <-data.frame(node.name=names, constraint=const,Inverse_constraint= invConstraint)
write.csv(d,"constraint.csv")

