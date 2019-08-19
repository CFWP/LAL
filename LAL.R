################################################################################
################################################################################
################################################################################
##
## Name:   LAL.R
## Author: Carel F.W. Peeters
##         Statistics for Omics Research Unit
##         Dept. of Epidemiology & Biostatistics
##         Amsterdam Public Health research institute
##         VU University medical center
##         Amsterdam, the Netherlands
## Email:	 cf.peeters@vumc.nl
##
## Last Update:	19/08/2019
## Description:	R script for the analyses contained in the manuscript:
##              - Peeters, C.F.W. (2019). "Social network analysis of corruption
##                structures: Adjacency matrices supporting the visualization 
##                and quantification of layeredness".
##
################################################################################
################################################################################
################################################################################


#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Preliminaries**
#' **------------------------------------------------------------------------**

## Set working directory to convenience
setwd("")

## Load packages
require("rags2ridges")
require("igraph")





#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Convenience functions**
#' **------------------------------------------------------------------------**

SharedAM <- function(AM1, AM2){
  ##############################################################################
  # Function that determines the shared elements in two adjacency matrices
  # - AM1 > Adjacency matrix 1
  # - AM2 > Adjacency matrix 2
  #
  # Notes:
  # - AM1 and AM2 should be square and of the same dimension
  # - Assumes both AM1 and AM2 are unweighted adjacency matrices
  ##############################################################################
  
  # Dependencies:
  # require("base")
  
  # Determine shared elements
  AMshare <- mat.or.vec(ncol(AM1),nrow(AM1))
  colnames(AMshare) = rownames(AMshare) <- colnames(AM1)
  for (i in 1:nrow(AM1)){
    for(j in 1:ncol(AM1)){
      if (AM1[i,j] == 1 & AM2[i,j] == 1){
        AMshare[i,j] <- 1
      }
    }
  }
  diag(AMshare) <- 0
  
  # Return
  return(AMshare)
}





#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Producing example 1**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Toy Data*\
#'**---------------------------------**\

## Relation-type 1
A1 <- cbind(
  c(0,1,1,1,1,1),
  c(1,0,1,0,0,1),
  c(1,1,0,1,0,0),
  c(1,0,1,0,1,0),
  c(1,0,0,1,0,1),
  c(1,1,0,0,1,0)
)

## Relation-type 2
A2 <- cbind(
  c(0,0,0,0,0,0),
  c(0,0,1,0,0,1),
  c(0,1,0,1,0,0),
  c(0,0,1,0,1,0),
  c(0,0,0,1,0,1),
  c(0,1,0,0,1,0)
)

## Labeling
ToyLabels <- c("A","B","C","D","E","F")
rownames(A1) = colnames(A1) <- ToyLabels
rownames(A2) = colnames(A2) <- ToyLabels


#'**---------------------------------**\
#'**Visualization*\
#'**---------------------------------**\

pdf("ToyExampleNetworks.pdf", width = 9, height = 5)
  par(mfrow=c(1,2))
  ## Relation-type 1
  CoordsToy <- Ugraph(A1, 
                      type = "weighted", 
                      lay = "layout_with_fr",
                      Vcolor = "white",
                      VBcolor = "black",
                      main = "Relation-type 1",
                      scale = 1)
  
  ## Relation-type 2
  Ugraph(A2, 
         type = "weighted", 
         lay = NULL,
         coords = CoordsToy,
         Vcolor = "white",
         VBcolor = "black",
         main = "Relation-type 2",
         scale = 1)
  par(mfrow=c(1,1))
dev.off()


#'**---------------------------------**\
#'**Calculating degrees*\
#'**---------------------------------**\
  
## Relation-Type 1 network
degreesA1 <- degree(graph.adjacency(A1, mode = "undirected"))
degreesA1

## Relation-Type 2 network
degreesA2 <- degree(graph.adjacency(A2, mode = "undirected"))
degreesA2





#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Producing example 2**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Visualization*\
#'**---------------------------------**\

## Determine shared adjacency matrices for relation-type combinations
shareA12 <- SharedAM(A1, A2)

## Visualization of differential and shared graphs
pdf("Example_DiffShared_Networks.pdf", width = 9, height = 5)
  ## Plot combination
  par(mfrow=c(1,2))
  
  ## Differential graph
  DiffGraph(A1, A2, 
            lay = NULL, 
            Vcolor = "white",
            VBcolor = "black",
            coords = CoordsToy,
            P1color = "red",
            P2color = "green",
            main = "Differential")
  
  ## Shared graph
  Ugraph(shareA12, 
         type = "weighted", 
         lay = NULL,
         Vcolor = "white",
         VBcolor = "black",
         coords = CoordsToy,
         scale = 1,
         pEcolor = "khaki4",
         main = "Shared") 
  par(mfrow=c(1,1))
dev.off()


#'**---------------------------------**\
#'**Calculating degrees*\
#'**---------------------------------**\

## Degrees from differential network
DiffA12 <- A1 - A2
degreeDiffA12 <- rowSums(abs(DiffA12))

## Degrees from shared network
degreeShareA12 <- rowSums(shareA12)

## (Ordered) Table
DegreeTableEx <- cbind(degreesA1,
                       degreesA2,
                       degreeDiffA12,
                       degreeShareA12)
DegreeTableEx <- DegreeTableEx[order(DegreeTableEx[,1], decreasing = TRUE),]
DegreeTableEx





#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Analysis Rath Affair**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Data*\
#'**---------------------------------**\

## Full network adjacency matrix
FullAM <- cbind(
  c(0,1,1,1,1,1,0,1,0,0,0),
  c(1,0,1,0,1,1,0,1,0,0,0),
  c(1,1,0,1,1,1,0,1,0,0,1),
  c(1,0,1,0,1,1,0,1,0,1,0),
  c(1,1,1,1,0,1,1,1,0,0,1),
  c(1,1,1,1,1,0,1,1,1,1,1),
  c(0,0,0,0,1,1,0,1,0,0,0),
  c(1,1,1,1,1,1,1,0,1,1,1),
  c(0,0,0,0,0,1,0,1,0,0,0),
  c(0,0,0,1,0,1,0,1,0,0,0),
  c(0,0,1,0,1,1,0,1,0,0,0)
)

## Pre-existing ties network adjacency matrix
PreAM <- cbind(
  c(0,1,0,0,0,0,0,0,0,0,0),
  c(1,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,0,1,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,1,0,0,0,0,1,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,1,0,0,1,0,0,0,0,1),
  c(0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,0,0,1,0,0,0)
)

## Resource transfer network adjacency matrix
RTAM <- cbind(
  c(0,0,1,0,0,1,0,0,0,0,0),
  c(0,0,1,0,0,0,0,0,0,0,0),
  c(1,1,0,0,1,0,0,0,0,0,1),
  c(0,0,0,0,1,0,0,0,0,0,0),
  c(0,0,1,1,0,1,0,1,0,0,0),
  c(1,0,0,0,1,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,1,0,0,0),
  c(0,0,0,0,1,0,1,0,1,1,1),
  c(0,0,0,0,0,1,0,1,0,0,0),
  c(0,0,0,0,0,1,0,1,0,0,0),
  c(0,0,1,0,0,1,0,1,0,0,0)
)

## Resource transfer network adjacency matrix
CoAM <- cbind(
  c(0,1,0,1,1,1,0,1,0,0,0),
  c(1,0,1,0,1,1,0,1,0,0,0),
  c(0,1,0,1,0,1,0,1,0,0,0),
  c(1,0,1,0,0,1,0,1,0,1,0),
  c(1,1,0,0,0,1,1,1,0,0,1),
  c(1,1,1,1,1,0,1,1,0,1,0),
  c(0,0,0,0,1,1,0,1,0,0,0),
  c(1,1,1,1,1,1,1,0,1,0,0),
  c(0,0,0,0,0,0,0,1,0,0,0),
  c(0,0,0,1,0,1,0,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0)
)

## Assign labels to rows and columns
LabelsNumeric <- c("1","2","3","4","5","6","7","8","9","10","11")
LabelsName    <- c("Drazdansky",
                   "Mlady",
                   "Rath",
                   "Salacova",
                   "Novanska",
                   "Pancova",
                   "Hajek",
                   "Kott",
                   "Jires",
                   "Kovanda",
                   "Rehak")
rownames(FullAM) = colnames(FullAM) <- LabelsName
rownames(PreAM)  = colnames(PreAM)  <- LabelsName
rownames(RTAM)   = colnames(RTAM)   <- LabelsName
rownames(CoAM)   = colnames(CoAM)   <- LabelsName


#'**---------------------------------**\
#'**Visualizations*\
#'**---------------------------------**\

## Node colors
## Coloring according to politician (blue) and business people (red)
COL <- c(rep("orangered",2), 
         "lightblue", 
         rep("orangered",2), 
         "lightblue", 
         "orangered", 
         "lightblue", 
         rep("orangered",3))


## First run with layout according to FR algorithm
## Retaining node-coordinates
CoordsFR <- Ugraph(FullAM, 
                   type = "weighted", 
                   Vcolor = COL,
                   lay = "layout_with_fr",
                   scale = 1,
                   pEcolor = "gray21",
                   main = "All relations")
save(CoordsFR, file = "CoordsFR.Rdata")
#load("CoordsFR.Rdata")


## Visualization of full and relation-specific networks
pdf("Rath_Relational_Networks.pdf", width = 15, height = 15)
  ## Plot combination
  par(mfrow=c(2,2))
  
  ## Full network
  Ugraph(FullAM, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "gray21",
         main = "Full network") 
  ## Relation-specific networks 
  ## All in same node-coordinates
  Ugraph(PreAM, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "gray21",
         main = "Pre-existing relations") 
  Ugraph(RTAM, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "gray21",
         main = "Resource transfer relations") 
  Ugraph(CoAM, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "gray21",
         main = "Collaboration relations") 
dev.off()


## Determine shared adjacency matrices for relation-type combinations
sharePreRT <- SharedAM(PreAM, RTAM)
sharePreCo <- SharedAM(PreAM, CoAM)
shareRTCo  <- SharedAM(RTAM, CoAM)

## Visualization of differential and shared graphs
pdf("Rath_DiffShared_Networks.pdf", width = 15, height = 20)
  ## Plot combination
  par(mfrow=c(3,2))
  
  ## Differential and shared: Pre-existing and resource transfer
  DiffGraph(PreAM, RTAM, 
            lay = NULL, 
            Vcolor = COL,
            coords = CoordsFR,
            P1color = "red",
            P2color = "green",
            main = "Differential: Pre-existing and resource transfer")
  Ugraph(sharePreRT, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "khaki4",
         main = "Shared: Pre-existing and resource transfer") 
  
  ## Differential and shared: Pre-existing and collaboration
  DiffGraph(PreAM, CoAM, 
            lay = NULL, 
            Vcolor = COL,
            coords = CoordsFR,
            P1color = "red",
            P2color = "blue",
            main = "Differential: Pre-existing and collaboration")
  Ugraph(sharePreCo, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "khaki4",
         main = "Shared: Pre-existing and collaboration") 
  
  ## Differential and shared: Resource transfer and collaboration
  DiffGraph(RTAM, CoAM, 
            lay = NULL, 
            Vcolor = COL,
            coords = CoordsFR,
            P1color = "green",
            P2color = "blue",
            main = "Differential: Resource transfer and collaboration")
  Ugraph(shareRTCo, 
         type = "weighted", 
         lay = NULL,
         Vcolor = COL,
         coords = CoordsFR,
         scale = 1,
         pEcolor = "khaki4",
         main = "Shared: Resource transfer and collaboration") 
dev.off()


#'**---------------------------------**\
#'**Analysis*\
#'**---------------------------------**\

## Raw degrees
degreeFull <- rowSums(FullAM)
degreePre  <- rowSums(PreAM)
degreeRT   <- rowSums(RTAM)
degreeCo   <- rowSums(CoAM)

## Degrees from differential network
diffPreRT <- PreAM - RTAM
diffPreCo <- PreAM - CoAM
diffRTCo  <- RTAM - CoAM
degreeDiffPreRT <- rowSums(abs(diffPreRT))
degreeDiffPreCo <- rowSums(abs(diffPreCo))
degreeDiffRTCo  <- rowSums(abs(diffRTCo))

## Degrees from shared network
degreeSharePreRT <- rowSums(sharePreRT)
degreeSharePreCo <- rowSums(sharePreCo)
degreeShareRTCo  <- rowSums(shareRTCo)

## (Ordered) Table
DegreeTable <- cbind(degreeFull,
                     degreePre,
                     degreeRT,
                     degreeCo,
                     degreeDiffPreRT,
                     degreeDiffPreCo,
                     degreeDiffRTCo,
                     degreeSharePreRT,
                     degreeSharePreCo,
                     degreeShareRTCo)
DegreeTable <- DegreeTable[order(DegreeTable[,1], decreasing = TRUE),]
DegreeTable


#'**---------------------------------**\
#'**Further Analysis*\
#'**---------------------------------**\

## Communities
Communities(FullAM)
Communities(PreAM)
Communities(RTAM)
Communities(CoAM)

## From a (multi-dimensional) scaling perspective:
## Pancova and Kott pair: Single operational unit in full graph
par(mfrow=c(1,1))
Ugraph(FullAM, 
       type = "weighted", 
       lay = "layout_with_mds",
       Vcolor = COL,
       scale = 1,
       pEcolor = "gray21",
       main = "Full network") 
Ugraph(PreAM, 
       type = "weighted", 
       lay = "layout_with_mds",
       Vcolor = COL,
       scale = 1,
       pEcolor = "gray21",
       main = "Pre-existing relations") 
Ugraph(RTAM, 
       type = "weighted", 
       lay = "layout_with_mds",
       Vcolor = COL,
       scale = 1,
       pEcolor = "gray21",
       main = "Resource transfer relations") 
Ugraph(CoAM, 
       type = "weighted", 
       lay = "layout_with_mds",
       Vcolor = COL,
       scale = 1,
       pEcolor = "gray21",
       main = "Collaboration relations") 

