# Looking at how many are dropped based on distance to edge

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste0(root,"/driversData/data/idaho/speciesData")); 

setwd("./ARTR");
X = read.csv("survD.csv");
1- mean(X$allEdge);
Y = read.csv("ARTR_nbhood_rings.csv");
mean(Y$propQ.8.10>0.85)

setwd("../HECO");
X = read.csv("survD.csv");
1- mean(X$allEdge);
Y = read.csv("HECO_nbhood_rings.csv");
mean(Y$propQ.8.10>0.85)

setwd("../POSE");
X = read.csv("survD.csv");
1- mean(X$allEdge);
Y = read.csv("POSE_nbhood_rings.csv");
mean(Y$propQ.8.10>0.85)

setwd("../PSSP");
X = read.csv("survD.csv");
1- mean(X$allEdge);
Y = read.csv("PSSP_nbhood_rings.csv");
mean(Y$propQ.8.10>0.85)