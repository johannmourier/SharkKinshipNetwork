#### LOAD PACKAGES ####

if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(asnipe)){install.packages('asnipe'); library(asnipe)}
if(!require(assortnet)){install.packages('assortnet'); library(assortnet)}
if(!require(sna)){install.packages('sna'); library(sna)}
if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)}
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)}
if(!require(vegan)){install.packages('vegan'); library(vegan)}
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)}
if(!require(magrittr)){install.packages('magrittr'); library(magrittr)}
if(!require(stringr)){install.packages('stringr'); library(stringr)}
if(!require(corrplot)){install.packages('corrplot'); library(corrplot)}
if(!require(devtools)){install.packages("devtools"); library(devtools)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}
if(!require(ade4)){install.packages('ade4'); library(ade4)}
if(!require(pander)){install.packages('pander'); library(pander)}
if(!require(related)){install.packages('related'); library(related)}
if(!require(lattice)){install.packages('lattice'); library(lattice)}
if(!require(plyr)){install.packages('plyr'); library(plyr)}
if(!require(ggridges)){install.packages('ggridges'); library(ggridges)}
if(!require(reshape)){install.packages('reshape'); library(reshape)}
if(!require(boot)){install.packages('boot'); library(boot)}
if(!require(boot)){install.packages('reshape2'); library(reshape2)}
#Color Palettes
if(!require(fishualize)){install.packages('fishualize'); library(fishualize)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}

#setup your repository
setwd("")
# load functions ====
source("02_functions.R")
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

#### IMPORT ALL DATA ####
individuals<-read.table("DataNetwork.csv",sep=",",header=T)
gbitot<-get_group_by_individual(individuals,location=individuals[,7], data_format="individuals")
gbitot
## mean group size
mean(rowSums(gbitot))
sd(rowSums(gbitot))

### Import data from individual sighted >14 times
data.raw<-read.csv("DataBTRS_15sightings.csv",
                   header=TRUE, 
                   sep = ",",
                   as.is = c("IDs"),
                   col.names = c("date", "Year", "Month", "Period","Site","Lat","Long","SP", "IDs")
)

data.raw$ID.list <- strsplit(data.raw$IDs, split = " ")
data.raw$identities <- lapply(data.raw$ID.list, unlist)     

# Behavior and Coop status ----
rawINFO <- data.raw[,c("date", "Year", "Month", "Period", "Site", "Lat","Long","SP")]


#### ID INFO ####

idcov <- read.csv("Attribute.csv",
                  header = TRUE,
                  sep = ",",
                  na.strings = "NA")

rownames(idcov) <- idcov$ID

# distribution of number of sightings per individuals
sightings<- read.csv("No.sightings.csv", header = TRUE, sep = ",", na.strings = "NA")
sightings<-sightings[which(sightings$Context=="All"),]

cdat <- sightings %>%
  summarise(mean = mean(Sightings),median(Sightings))
cdat.dens <- ggplot_build(ggplot(sightings, aes(x=Sightings)) + geom_density())$data[[1]] 

cdat.dens <- ggplot_build(ggplot(sightings, aes(x=Sightings, fill=Context, colour = Context)) + geom_density())$data[[1]] 
# plot density of sightings
ggplot(sightings, aes(x=Sightings))+
  geom_density(color="darkblue", fill="lightblue",alpha=0.4)+
  geom_segment(data = cdat.dens, aes(x = 14, xend = 14, y = 0, yend = mean(cdat.dens[which(round(cdat.dens$x)==14),3]), colour = "blue"),
               linetype = "dashed", size = 1) +
  geom_segment(data = cdat.dens, aes(x = mean(sightings$Sightings), xend = mean(sightings$Sightings), y = 0, yend = mean(cdat.dens[which(round(cdat.dens$x,1)==round(mean(sightings$Sightings),1)),3]), colour = "black"),
               linetype = "dotted", size = 1)+
  scale_color_manual(labels = c("Mean", "Median"), values = c("blue", "black"))

#mean signthings
mean(sightings$Sightings) # mean = 14.92381
median(sightings$Sightings) # median = 14
sd(sightings$Sightings) # sd = 8.048764



#### DATA ####

##### Genetic relatedness ####

related.matrix=read.csv("MatrixGeneticRelatedness_Trioml.csv",sep=",",header=T)
related.matrix<-related.matrix[,-1]
rownames(related.matrix)=colnames(related.matrix)
related.matrix<-as.matrix(related.matrix)
related.matrix <- related.matrix + t(related.matrix) - diag(diag(related.matrix))


# All data ----
rawGBI <- get_group_by_individual(data.raw$identities, 
                                  data_format = "groups")
rownames(rawGBI) <- paste("G",1:nrow(rawGBI), sep = "")
rownames(rawINFO) <- paste("G",1:nrow(rawGBI), sep = "")


#### Structural factors ####

####  Home Range and Spatial overlap matrix #### 
###

SpatialProfile=read.table("SpatialProfile.csv",sep=",",header=T)
SpatialOverlap<-vegdist(SpatialProfile[,2:8], method="bray",diag=TRUE, upper=TRUE)
SpatialOverlap<-as.matrix(1-SpatialOverlap)
rownames(SpatialOverlap)<-SpatialProfile[,1]
colnames(SpatialOverlap)<-rownames(SpatialOverlap)
#write.csv(SpatialOverlap,file = "MatrixSpatialOverlap.csv")

### Spatial profile graph for supplementary materials
SP<- SpatialProfile[order(SpatialProfile[,1]),]    
colnames(SP)<-c("ID", "1", "2", "3", "4","5","6","7")
rownames(SP)<-SP[,1]

SP<-SP[which(rownames(SP) %in% rownames(traits_by_related$related)),]

# Plot heatmap of spatial profiles (site encounter rate)

cols2 <- colorRampPalette(c("white","#F3D617","#FF7F00","#FF0000","#ED0000"))(256)
#cols2 <- colorRampPalette(c("white","#FDF2D6","#FBD04F","#FF7F41","#D74850"))(256)
P1<-levelplot(t(as.matrix(SP[,-1])), at=seq(0, 1, 0.01), col.regions=cols2,grid=TRUE)


### Temporal overlap ###
### SRI with SP = month

individuals=read.csv("TempData.csv",sep=",",header=T)
SPs<-get_sampling_periods(individuals[,c(1,2)],individuals[,5],1,data_format="individuals")
SPs[1,,]
dim(SPs)
TempOverlap<-get_network(SPs, association_index = "SRI", data_format="SP")
#write.csv(TempOverlap,file = "MatrixTempOverlap.csv")

#### Sex similarity matrix ###
#Create an empty N x N matrix to store species similarities
N=nrow(idcov)
sex_sim <- matrix(0, nrow=N, ncol=N)
# Loop through each row and each column in the data
for (row in c(1:N)) {
  for (col in c(1:N)) {
    # Test if the sexes are the same
    if (idcov$Sex[row] == idcov$Sex[col]) {
      sex_sim[row,col] <- 1
    } else {
      sex_sim[row,col] <- 0
    }
  }
}
rownames(sex_sim)<-rownames(idcov)
colnames(sex_sim)<-rownames(idcov)

######## Size similarity #########
idcov$SizeClass<-rep(NA, nrow(idcov))

for (i in 1:nrow(idcov)){
  if(idcov$Size[i]<110){
    idcov$SizeClass[i]<-"A"
  }
  else {
    if(idcov$Size[i]>=110 & (idcov$Size[i]<120)){
      idcov$SizeClass[i]<-"B"
    }
    else {
      if(idcov$Size[i]>=120 & (idcov$Size[i]<130)){
        idcov$SizeClass[i]<-"C"
      }
      else {
        if(idcov$Size[i]>=130 & (idcov$Size[i]<140)){
          idcov$SizeClass[i]<-"D"
        }
        else {
          if(idcov$Size[i]>=140 & (idcov$Size[i]<150)){
            idcov$SizeClass[i]<-"E"
          }
          else {
            if(idcov$Size[i]>=150){
              idcov$SizeClass[i]<-"F"
            }
          }
        }
      }
    }
  }
}
g=nrow(idcov)
size_sim <- matrix(0, nrow=g, ncol=g)
# Loop through each row and each column in the data
for (row in c(1:g)) {
  for (col in c(1:g)) {
    # Test if the sexes are the same
    if (idcov$SizeClass[row] == idcov$SizeClass[col]) {
      size_sim[row,col] <- 1
    } else {
      size_sim[row,col] <- 0
    }
  }
}
rownames(size_sim)<-rownames(idcov)
colnames(size_sim)<-rownames(idcov)


####### Gregariousness ######
rawGBI <- get_group_by_individual(data.raw$identities, 
                                  data_format = "groups")
rownames(rawGBI) <- paste("G",1:nrow(rawGBI), sep = "")
rownames(rawINFO) <- paste("G",1:nrow(rawGBI), sep = "")

network <- get_network(rawGBI, association_index = "SRI")

Gregariousness<-network
for (i in 1:nrow(network)){
  for (j in 1:ncol(network)){
    Gregariousness[i,j]<-log((rowSums(network)[i]-network[i,j])*(colSums(network)[j]-network[i,j]))
  }
}
diag(Gregariousness)<-0


#########################
# All data ----
rawGBI

# Subset by home range data ---- 
idsSUBSET <- colnames(rawGBI)[which(colnames(rawGBI) %in% 
                                      colnames(related.matrix))]

# filter Numrec 5% ----
GBI_fullSET <- rawGBI[,which(colnames(rawGBI) %in% idsSUBSET)]

GBI_fullSET <- cbind(GBI_fullSET, 
                     rawINFO[which(rownames(rawINFO) %in% 
                                     rownames(GBI_fullSET)),])

gbi_ALL <- GBI_fullSET %>% 
  dplyr::select(-date, -Year, -Month, -Period, -Site, -Lat, -Long, -SP)


#### COVARIABLES ####

# set ID labels
IDlabels <- sort(colnames(gbi_ALL))
info_IDrestrict <- subset(idcov, ID %in% IDlabels) 

# Home Range Overlap ----
HROmatrix <- matrix_subset(SpatialOverlap, IDlabels) %>% 
  as.matrix() %>% 
  order_matrix()

dim(HROmatrix)

# Temporal Overlap ----
TOmatrix <- matrix_subset(TempOverlap, IDlabels) %>% 
  as.matrix() %>% 
  order_matrix()

dim(TOmatrix)

# Sex ----
SEXmatrix <- matrix_subset(sex_sim, IDlabels) %>% 
  as.matrix() %>% 
  order_matrix()

dim(SEXmatrix)

# Size ----
AGEmatrix <- matrix_subset(size_sim, IDlabels) %>% 
  as.matrix() %>% 
  order_matrix()

dim(AGEmatrix)

# Gregariousmess ----
Gregmatrix <- matrix_subset(Gregariousness, IDlabels) %>% 
  as.matrix() %>% 
  order_matrix()

dim(Gregmatrix)


#### Subset by individual traits ####

# ID by relatedness ----
subsetID_related <- Reduce(intersect, list(rownames(TOmatrix), rownames(HROmatrix),rownames(Gregmatrix), rownames(AGEmatrix), rownames(SEXmatrix), rownames(related.matrix)))

GBI_related <- rawGBI[,which(colnames(rawGBI) %in% subsetID_related)]
GBI_related_filter <-  GBI_related[which(rowSums(GBI_related) > 0),
                                   which(colSums(GBI_related) > nrow(GBI_related)*0.05)]
GBI_related_filter <- cbind(GBI_related_filter, rawINFO[which(rownames(rawINFO) %in% rownames(GBI_related_filter)),])


# List GBI by Traits ----
listTRAIT <- list(related = GBI_related_filter)

gbi_ALL <- listTRAIT[[1]] %>% 
  dplyr::select(-date, -Year, -Month, -Period, -Site, -Lat, -Long, -SP)


# Covariables subset ----
traits <- list(to = TOmatrix,
               hro = HROmatrix,
               greg = Gregmatrix,
               size = AGEmatrix,
               sex = SEXmatrix,
               related = related.matrix)

# traits by related ----
traits_by_related <- list()
for(i in 1:length(traits)){
  traits_by_related[[i]] <- matrix_subset(traits[[i]], subsetID_related)
}
names(traits_by_related) <- names(traits)


# Plot heatmap of spatial profiles (site encounter rate)

SP<-SP[which(rownames(SP) %in% rownames(traits_by_related$related)),]
cols2 <- colorRampPalette(c("white","#F3D617","#FF7F00","#FF0000","#ED0000"))(256)
#cols2 <- colorRampPalette(c("white","#FDF2D6","#FBD04F","#FF7F41","#D74850"))(256)
P1<-levelplot(t(as.matrix(SP[,-1])), at=seq(0, 1, 0.01), col.regions=cols2,grid=TRUE)



############# SRI + DENMAT ###

SRI_real <- SRI(gbi_ALL)$SRI %>% 
  order_matrix() %>% 
  as.matrix()

INMAT_real <- SRI(gbi_ALL)$SRI.numerator %>% 
  order_matrix() %>% 
  as.matrix()

DENMAT_real <- SRI(gbi_ALL)$SRI.denominator %>% 
  order_matrix() %>% 
  as.matrix()

SD_SRI_real <-sd(matrix_unfold(SRI_real))


# RANDOM GBI 
Npermute_GBI <- 25000
list_RANDOM_GBI <- list()
list_RANDOM_GBI <- null_checkerboard(gbi_ALL,
                                     iter = Npermute_GBI)

#### RANDOM SRI + INMAT + DENMAT ####
list_RANDOM_SRI <- list()
list_RANDOM_INMAT <- list()
list_RANDOM_DENMAT <- list()

list_SRI <- list()
list_inmat <- list()
list_denmat <- list()

for(j in 5001:Npermute_GBI){ # Remove the first 5000 random GBI
  tmp.list <- SRI(list_RANDOM_GBI$perm[j][[1]])
  
  list_SRI[[j]] <- order_matrix(tmp.list$SRI)
  list_inmat[[j]] <- order_matrix(tmp.list$SRI.numerator)
  list_denmat[[j]] <- order_matrix(tmp.list$SRI.denominator)
}

list_RANDOM_SRI <- list_SRI
list_RANDOM_INMAT <- list_inmat
list_RANDOM_DENMAT <- list_denmat


# Remove empty levels of the list
list_RANDOM_SRI <- list_RANDOM_SRI[lapply(list_RANDOM_SRI,length)>0]
list_RANDOM_INMAT <- list_RANDOM_INMAT[lapply(list_RANDOM_INMAT,length)>0]
list_RANDOM_DENMAT <- list_RANDOM_DENMAT[lapply(list_RANDOM_DENMAT,length)>0]


# set number of permutations
Npermute <- 20000

#### Modularity SRI Real ----
modularity_SRI_real <- vector()
modularity_size_SRI_real <- list()
modularity_memb_SRI_real <- list()

tmp1 <- graph.adjacency(SRI_real, 
                        mode="undirected",
                        weighted=TRUE,
                        diag=FALSE)

tmpmod <- cluster_leading_eigen(graph = tmp1,
                                weights = E(tmp1)$weight)

modularity_SRI_real <- modularity(tmpmod)
modularity_size_SRI_real <- table(tmpmod$membership)
modularity_memb_SRI_real <- data.frame(Id=row.names(SRI_real), module=tmpmod$membership)


#### Modularity SRI Random ----
modularity_SRI_random <- NULL
listMOD_SRI_random <- list()

for(j in 1:Npermute){
  
  tmp.graph <- graph.adjacency(as.matrix(list_RANDOM_SRI[[j]]),
                               mode = "undirected",
                               diag = FALSE,
                               weighted = TRUE)
  #tmp.graph <- giant.component(tmp.graph)
  
  tmpmod2 <- cluster_leading_eigen(graph = tmp.graph,
                                   weights = E(tmp.graph)$weight,
                                   options=list(maxiter=1000000))
  
  modularity_SRI_random[j] <- modularity(tmpmod2)
  
}

listMOD_SRI_random <- modularity_SRI_random

#### SD random SRI ----
SDSRI_random <- list()

sdtmp <- vector()

for(k in 1:Npermute){
  sdtmp[k] <- sd(matrix_unfold(list_RANDOM_SRI[[k]]))
}

SDSRI_random <- sdtmp


#### Descriptive results ----

# Sex ----
info_IDrestrict[which(info_IDrestrict$DNA=="Y"),]$Sex %>% 
  na.omit() %>% 
  table()

# Genetic Relatedness ----
# subset data by individuals in this study
#function standar error
se <- function(x) {sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}

mean(matrix_unfold(related.matrix))
sd(matrix_unfold(related.matrix))
se(matrix_unfold(related.matrix))

#### CORRELOGRAM
rgb.palette3 <- colorRampPalette(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "white", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"), space = "rgb")
rgb.palette4 <- colorRampPalette(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "white", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"), space = "rgb")

corrplot(traits_by_related$related, method = "color",type="lower", order="original",tl.cex=0.7,tl.pos = "lt",
         col=rev(rgb.palette3(100)), diag = TRUE, cl.lim = c(0,1),tl.col="black",is.corr = FALSE)
corrplot(traits_by_related$hro, type="upper",method = "color", order="original",tl.cex=0.7,tl.pos = "n",
         col=rgb.palette4(100), diag = TRUE,cl.lim = c(0,1),add=TRUE,tl.col="black",is.corr = FALSE)


### Matel test between relatedness and spatial overlap matrices
mantel(traits_by_related$related, traits_by_related$hro, permutations = 20000)# not different


#### MRQAP ####

# MRQAP subset 43 individuals with relatedness info ----
mrqap.related <- list()

### MRQAP function with custom permutation networks

# convert list list_RANDOM_SRI in 3x20000x43x43 array
randomSRI<-array(data = NA, dim = c(Npermute,dim(traits_by_related$related)))

for (j in 1:Npermute){
  randomSRI[j,,]<-as.matrix(list_RANDOM_SRI[[j]])
  colnames(randomSRI[j,,])<-colnames(traits_by_related$related)
  rownames(randomSRI[j,,])<-colnames(traits_by_related$related)
}


mrqap.related <- mrqap.custom.null(SRI_real ~
                                     traits_by_related$to +
                                     traits_by_related$hro + 
                                     traits_by_related$greg + 
                                     traits_by_related$size +
                                     traits_by_related$sex,
                                   random.y=randomSRI, 
                                   intercept = FALSE,
                                   test.statistic = "beta")

#Results of MRQAP
mrqap.related



# Mantel SRI ~ R ----
rl1 <- ade4::mantel.rtest(as.dist(SRI_real), as.dist(traits_by_related$related), 9999)

pdf(file = "mantel_SRI_relatedness.pdf", 
    height = 6, width = 7.25, units = "in", res = 400)

plot(rl1, main = "")
text(rl1$plot$hist$breaks[1] * 0.7,
     max(rl1$plot$hist$counts) * 0.9, 
     paste("Observation =", format(rl1$obs, digits = 3)))
text(rl1$plot$hist$breaks[1] * 0.7,
     max(rl1$plot$hist$counts) * 0.8, 
     paste("Simulated p-value =", format(rl1$pvalue, digits = 3)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.9, 
     paste("Std. Obs. =", format(rl1$expvar[1], digits = 4)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.8, 
     paste("Expectation =", format(rl1$expvar[2], digits = 4)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.7, 
     paste("Variance =", format(rl1$expvar[3], digits = 4)))

dev.off()

#### GLM data for relatedness n=43 ####

#GLMdata <- data.frame(HROunfold = matrix_unfold(traits_by_related$hro),
#                      AGEunfold = matrix_unfold(traits_by_related$size),
#                      TOunfold = matrix_unfold(traits_by_related$to))

GLMdata <- data.frame(HROunfold = matrix_unfold(traits_by_related$hro),
                      AGEunfold = matrix_unfold(traits_by_related$size),
                      TOunfold = matrix_unfold(traits_by_related$to),
                      SEXunfold = matrix_unfold(traits_by_related$sex),
                      GREGunfold = matrix_unfold(traits_by_related$greg))


### GAI #####
#### GAI Real ####
listGAI_real <- list()
SDgai_real <- list()


resids <- vector()
sd.gai <- vector()
sd.list <- list()
residsLIST <- list()

tmp.SRI <- matrix_unfold(SRI_real)
#tmp.inmat <- matrix_unfold(list_INMAT_real)
#tmp.denmat <- matrix_unfold(list_DENMAT_real)
#tmp.glm <- glm(cbind(tmp.inmat, (tmp.denmat - tmp.inmat)) ~ TOunfold + HROunfold, 
#               data = GLMdata,
#               family = binomial(link = "logit"))

#tmp.glm <- glm(tmp.SRI ~ TOunfold + HROunfold, 
#               data = GLMdata,
#               family = binomial(link = "logit"))

tmp.glm <- glm(tmp.SRI ~ TOunfold + HROunfold + AGEunfold + SEXunfold + GREGunfold, 
               data = GLMdata,
               family = binomial(link = "logit"))

resids <- unname(residuals.glm(tmp.glm,
                               type = "deviance"))

listGAI_real <- resids
SDgai_real <- sd(resids)


# Modularity GAI real ----
modularity_GAI_real <- vector()
modularity_size_GAI_real <- list()
modularity_memb_GAI_real <- list()
listMATRIX_GAI_real <- list()

tmpmat <- matrix_folding(matrix.unfold = listGAI_real,
                         matrix.labels = rownames(SRI_real),
                         matrix.reference = SRI_real)

tmp1 <- graph.adjacency(as.matrix(tmpmat), 
                        mode="undirected",
                        weighted=TRUE,
                        diag=FALSE)

tmpmod3 <- cluster_leading_eigen(graph = tmp1,
                                 weights = E(tmp1)$weight,
                                 options=list(maxiter=1000000))

listMATRIX_GAI_real <- tmpmat
modularity_GAI_real <- modularity(tmpmod3)
modularity_size_GAI_real <- table(tmpmod3$membership)
modularity_memb_GAI_real <- data.frame(Id=row.names(tmpmat), module=tmpmod3$membership)

#### GAI RANDOM ####

listGAI_random <- list()
SDgai_random <- list()
listMATRIX_GAI_random <- list()

resids <- vector()
sd.gai <- vector()
sd.list <- list()
residsLIST <- list()
tmp.GAImat <- list()

for(j in 1:Npermute){
  
  tmp.rand.SRI <- matrix_unfold(list_RANDOM_SRI[[j]])
  #tmp.rand.inmat <- matrix_unfold(list_RANDOM_INMAT[[i]][[j]])
  # tmp.rand.denmat <- matrix_unfold(list_RANDOM_DENMAT[[i]][[j]])
  
  # tmp.glm <- glm(cbind(tmp.rand.inmat, (tmp.rand.denmat - tmp.rand.inmat)) ~ TOunfold + HROunfold, 
  #                data = GLMdata,
  #                family = binomial(link = "logit"))
  
  
  #tmp.glm <- glm(tmp.rand.SRI ~ TOunfold + HROunfold, 
  #              data = GLMdata,
  #              family = binomial(link = "logit"))
  
  tmp.glm <- glm(tmp.rand.SRI ~ TOunfold + HROunfold + AGEunfold + SEXunfold + GREGunfold, 
                 data = GLMdata,
                 family = binomial(link = "logit"))
  
  resids <- unname(residuals.glm(tmp.glm, type = "deviance"))
  
  residsLIST[[j]] <- resids
  sd.gai[j] <- sd(resids)
  
  tmp.GAImat[[j]] <- matrix_folding(matrix.unfold = residsLIST[[j]],
                                    matrix.labels = rownames(SRI_real),   
                                    matrix.reference = SRI_real)        
  
}

listGAI_random <- residsLIST
SDgai_random <- sd.gai
listMATRIX_GAI_random <- tmp.GAImat


#### Modularity GAI Random ####

listMOD_GAI_random <- list()
modularity_GAI_random <- vector()
sizemod_gai_rand_vec <- vector()
sizemod_GAI_random <- list()

modularity_GAI_random <- vector()
assort_random <- vector()
for(j in 1:Npermute){
  
  tmp1 <- graph.adjacency(listMATRIX_GAI_random[[j]], 
                          mode="undirected",
                          weighted=TRUE,
                          diag=FALSE)
  tmpmod4 <- cluster_leading_eigen(graph = tmp1,
                                   weights = E(tmp1)$weight,
                                   options = list(maxiter = 1000000))
  modularity_GAI_random[j] <- modularity(tmpmod4)
  sizemod_gai_rand_vec[j] <- length(communities(tmpmod4))
}
listMOD_GAI_random <- modularity_GAI_random

### Make GAI matrices symmetrics ####

list_SYM_MATRIX_GAI_real <- listMATRIX_GAI_real

list_SYM_MATRIX_GAI_real[upper.tri(list_SYM_MATRIX_GAI_real)] = t(list_SYM_MATRIX_GAI_real)[upper.tri(t(list_SYM_MATRIX_GAI_real))]


list_SYM_MATRIX_GAI_random <- listMATRIX_GAI_random

for(j in 1:Npermute){
  list_SYM_MATRIX_GAI_random[[j]][upper.tri(list_SYM_MATRIX_GAI_random[[j]])] = t(list_SYM_MATRIX_GAI_random[[j]])[upper.tri(t(list_SYM_MATRIX_GAI_random[[j]]))]
}


# Mantel SRI ~ R ----
rl1 <- ade4::mantel.rtest(as.dist(list_SRI_related), as.dist(traits_by_related$related), 9999)

png(file = "mantel_SRI_relatedness.png", 
    height = 6, width = 7.25, units = "in", res = 400)

plot(rl1, main = "")
text(rl1$plot$hist$breaks[1] * 0.7,
     max(rl1$plot$hist$counts) * 0.9, 
     paste("Observation =", format(rl1$obs, digits = 3)))
text(rl1$plot$hist$breaks[1] * 0.7,
     max(rl1$plot$hist$counts) * 0.8, 
     paste("Simulated p-value =", format(rl1$pvalue, digits = 3)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.9, 
     paste("Std. Obs. =", format(rl1$expvar[1], digits = 4)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.8, 
     paste("Expectation =", format(rl1$expvar[2], digits = 4)))
text(max(rl1$plot$hist$breaks) * 0.7,
     max(rl1$plot$hist$counts) * 0.7, 
     paste("Variance =", format(rl1$expvar[3], digits = 4)))

dev.off()

# Mantel GAI ~ R ----
rl2 <- ade4::mantel.rtest(as.dist(list_SYM_MATRIX_GAI_real), as.dist(traits_by_related$related), 9999)
plot(rl2, main = "")
text(rl2$plot$hist$breaks[1] * 0.7,
     max(rl2$plot$hist$counts) * 0.9, 
     paste("Observation =", format(rl2$obs, digits = 3)))
text(rl2$plot$hist$breaks[1] * 0.7,
     max(rl2$plot$hist$counts) * 0.8, 
     paste("Simulated p-value =", format(rl2$pvalue, digits = 3)))
text(max(rl2$plot$hist$breaks) * 0.7,
     max(rl2$plot$hist$counts) * 0.9, 
     paste("Std. Obs. =", format(rl2$expvar[1], digits = 4)))
text(max(rl2$plot$hist$breaks) * 0.7,
     max(rl2$plot$hist$counts) * 0.8, 
     paste("Expectation =", format(rl2$expvar[2], digits = 4)))
text(max(rl2$plot$hist$breaks) * 0.7,
     max(rl2$plot$hist$counts) * 0.7, 
     paste("Variance =", format(rl2$expvar[3], digits = 4)))



#### FIGURES ####

# Social preferences ====

dfsd.SRI <- data.frame()
dfsd.gai <- data.frame()

tmp.ci.sd.SRI <- quantile(SDSRI_random,
                          probs = c(0.025, 0.975),
                          type = 2)  

tmp.ci.sd.gai <- quantile(SDgai_random,
                          probs = c(0.025, 0.975),
                          type = 2)

dfsd.SRI[1, "SD"] <- SD_SRI_real
dfsd.SRI[1, "context"] <- "Global"
dfsd.SRI[1, "lowCI"] <- tmp.ci.sd.SRI[1]
dfsd.SRI[1, "uprCI"] <- tmp.ci.sd.SRI[2]
dfsd.SRI[1, "index"] <- "SRI"

dfsd.gai[1, "SD"] <- SDgai_real
dfsd.gai[1, "context"] <- "Global"
dfsd.gai[1, "lowCI"] <- tmp.ci.sd.gai[1]
dfsd.gai[1, "uprCI"] <- tmp.ci.sd.gai[2]
dfsd.gai[1, "index"] <- "GAI"


dfSD_real <- rbind(dfsd.SRI, dfsd.gai) %>% 
  mutate(index = factor(index, levels = c("SRI", "GAI")))

dfSD_random <- cbind(
  rbind(reshape2::melt(SDgai_random),
        reshape2::melt(SDSRI_random)),
  rep(c("GAI", "SRI"), 
      each = (dim(reshape2::melt(SDSRI_random))[1])),
  rep("Global",2*dim(reshape2::melt(SDSRI_random))[1])) 
names(dfSD_random) <- c("SD", "index","context")
dfSD_random$index <- factor(dfSD_random$index, 
                            levels = c("SRI", "GAI"))


# Plot standard deviation ====

gplot.sd <- ggplot(dfSD_random,
                   aes(x = context, y = SD,
                       fill = index)) +
  geom_flat_violin(trim = FALSE,
                   scale = "width",
                   adjust = 2,
                   colour = "NA") +
  #facet_grid(index ~ .)+
  facet_wrap(.~index, ncol = 2, scales = "free", dir = "v") +
  scale_fill_manual(values=c(alpha("grey80", 0.5), 
                             alpha("#D55E00", 0.25))) +
  #scale_x_discrete(labels = c("SRI", "GAI")) +
  scale_y_continuous(breaks = equal_breaks(4, 0.05),
                     minor_breaks = FALSE,
                     expand = c(0.05, 0),
                     #label = scientific_10,
                     label = function(x){round(x, digits = 3)}) +
  geom_errorbar(data = dfSD_real, 
                aes(x = context,
                    ymin = lowCI, 
                    ymax = uprCI),
                size = 0.5,
                width = 0.075, 
                color = "#0000ff96") +
  geom_point(data = dfSD_real, 
             aes(x = context, 
                 y = SD),
             size = 2,
             color = "#ff000096") +
  theme_violin + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  coord_flip() 


# Social division ====

dfmod.SRI <- data.frame()
dfmod.gai <- data.frame()

tmp.ci.Q.SRI <- quantile(listMOD_SRI_random,
                         probs = c(0.025, 0.975),
                         type = 2)  

tmp.ci.Q.gai <- quantile(listMOD_GAI_random,
                         probs = c(0.025, 0.975),
                         type = 2)

dfmod.SRI[1, "Q"] <- modularity_SRI_real
dfmod.SRI[1, "context"] <- "Global"
dfmod.SRI[1, "lowCI"] <- tmp.ci.Q.SRI[1]
dfmod.SRI[1, "uprCI"] <- tmp.ci.Q.SRI[2]
dfmod.SRI[1, "index"] <- "SRI"

dfmod.gai[1, "Q"] <- modularity_GAI_real
dfmod.gai[1, "context"] <- "Global"
dfmod.gai[1, "lowCI"] <- tmp.ci.Q.gai[1]
dfmod.gai[1, "uprCI"] <- tmp.ci.Q.gai[2]
dfmod.gai[1, "index"] <- "GAI"

dfmod_real <- rbind(dfmod.SRI, dfmod.gai) %>% 
  mutate(index = factor(index, levels = c("SRI", "GAI")))

dfmod_random <- cbind(
  rbind(reshape2::melt(listMOD_SRI_random),
        reshape2::melt(listMOD_GAI_random)),
  rep(c("SRI", "GAI"), 
      each = (dim(reshape2::melt(listMOD_SRI_random))[1])),
  rep("Global",2*dim(reshape2::melt(listMOD_SRI_random))[1])) 
names(dfmod_random) <- c("Q", "index", "context")
dfmod_random$index <- factor(dfmod_random$index,
                             levels = c("SRI", "GAI"))


# Plot Modularity ====

gplot.mod <- ggplot(dfmod_random,
                    aes(x = context, y = Q,
                        fill = index)) +
  geom_flat_violin(trim = FALSE,
                   scale = "width",
                   adjust = 2,
                   colour = "NA") +
  #facet_grid(index ~ .)+
  facet_wrap(.~index, ncol = 2, scales = "free", dir = "v") +
  scale_fill_manual(values=c(alpha("grey80", 0.5), 
                             alpha("#D55E00", 0.25))) +
  #scale_x_discrete(labels = c("SRI", "GAI")) +
  scale_y_continuous(breaks = equal_breaks(4, 0.05),
                     minor_breaks = FALSE,
                     expand = c(0.05, 0),
                     #label = scientific_10,
                     label = function(x){round(x, digits = 3)}) +
  geom_errorbar(data = dfmod_real, 
                aes(x = context,
                    ymin = lowCI, 
                    ymax = uprCI),
                size = 0.5,
                width = 0.075, 
                color = "#0000ff96") +
  geom_point(data = dfmod_real, 
             aes(x = context, 
                 y = Q),
             size = 2,
             color = "#ff000096") +
  theme_violin + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  coord_flip() 


# TABLE RESULTS ====
pander(dfSD_real)
pander(dfmod_real)

# PLOT FIGURE 01 ====


pdf(file = "Figure 2.pdf",
    height = 7, width = 8)

(gplot.sd + gplot.mod) + 
  plot_layout(nrow = 2, heights = c(0.45, 0.45))
#plot_layout(ncol = 2, widths = c(0.6, 0.6))

dev.off()


#### PLOT NETWORKS ====

#color palettes
pal<-fish(5, option = "Sparisoma_viride")
pal2<-rev(fish(100, option = "Prionace_glauca"))


# Colorblind palette
color.mod <- c("#E69F00", # orange
               "#56B4E9", # sky blue
               "#009E73", # bluish green
               "#F0E442", # yellow
               "#D55E00", # vermilion
               "#CC79A7", # reddish purple
               "#0072B2", # blue
               "#999999" # grey
)
svg(file = "Networks_Global1.svg",
    height = 9.2, width = 7.2)

# SRI Networks 
par(mfrow = c(1,2), mar = c(0,1,0,1))

layout.net <- list()
community.SRI<-list()

SRInet <- graph.adjacency(as.matrix(SRI_real),
                          mode = "undirected",
                          weighted = TRUE,
                          diag = FALSE)

mod.SRI <- cluster_leading_eigen(SRInet)
community.SRI<-mod.SRI

# Set node size
node.size.fp <- idcov$Size

layout.net <- igraph::layout.fruchterman.reingold(SRInet)

#layout.net[[i]] <- layout.forceatlas2(SRInet,directed=FALSE, iterations = 500, plotstep = 10)

plot.igraph(SRInet,
            layout = layout.net,
            main = "SRI",
            edge.color = alpha('grey80', 0.9),
            vertex.label = NA,
            vertex.shape = "sphere",
            vertex.size = rescale(node.size.fp, 
                                  r.out = range(5,15)),
            edge.width = E(SRInet)$weight * 10,
            vertex.color = pal[mod.SRI$membership]) 

# GAI NETWORKS ----
community.GAI<-list()

gainet <- graph.adjacency(as.matrix(listMATRIX_GAI_real),
                          mode = "undirected",
                          weighted = TRUE,
                          diag = FALSE)

gainet <- delete.edges(gainet, which(E(gainet)$weight < 0))
mod.gai <- cluster_leading_eigen(graph = gainet,
                                 weights = E(gainet)$weight,
                                 options = list(maxiter=1000000))
community.GAI<-mod.gai

plot.igraph(gainet,
            layout = layout.net,
            main = "GAI",
            edge.color = scales::alpha('#D55E00', 0.25),
            vertex.label = NA,
            vertex.shape = "sphere",
            vertex.size = rescale(node.size.fp, 
                                  r.out = range(5,15)),
            edge.width = E(gainet)$weight*5,
            #edge.width = 5 * rescale(E(gainet)$weight, 
            #                         r.out = range(E(SRInet)$weight)),
            vertex.color = pal[mod.gai$membership]) 


# Relatedness NETWORKS ----
idcov_genet<-idcov[which(idcov$DNA=="Y"),]

idcov_genet <- idcov_genet[which(rownames(idcov_genet) %in% 
                                   rownames(listMATRIX_GAI_real)),]
idcov_genet <- idcov_genet[order(idcov_genet$ID),]

Rnet <- graph.adjacency(as.matrix(traits_by_related$related),
                        mode = "undirected",
                        weighted = TRUE,
                        diag = FALSE)

#sub.Rnet <- delete.edges(as.network(Rnet), which(abs(E(Rnet)$weight) < 0.25))
sub.Rnet <- subgraph.edges(Rnet, E(Rnet)[E(Rnet)$weight>0.25], del=F)

c_scale <- colorRamp(pal2)
V(sub.Rnet)[idcov_genet$Sex=="M"]$color<-"grey50"
V(sub.Rnet)[idcov_genet$Sex=="F"]$color<-"seashell2"
plot.igraph(sub.Rnet,
            layout = layout.net,
            main = "Genetic relatedness",
            edge.color = apply(c_scale(E(sub.Rnet)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) ),
            vertex.label = NA,
            vertex.shape = "sphere",
            vertex.size = rescale(node.size.fp, 
                                  r.out = range(5,15)),
            #edge.width = E(gainet)$weight,
            edge.width = 5 * rescale(E(sub.Rnet)$weight, 
                                     r.out = range(E(SRInet)$weight))
) 

dev.off()


# Calculate Social differentiation S (Whitehead 2008) using boostraps

library(boot)
# function to get coefficient of variation 
get_cv <- function(x){ sd(x,na.rm=T)/mean(x,na.rm=T)}
fc <- function(d, i){
  d2 <- d[i,]
  network<-get_network(d2)
  cv<-get_cv(network)
  return(cv)
}
results=boot(data=GBI_related, statistic = fc, R = 1000)
results



#### Sex matrices ###
#Create an empty N x N matrix to store species similarities

#males
n=nrow(idcov_genet)
males <- matrix(0, nrow=n, ncol=n)
# Loop through each row and each column in the data
for (row in c(1:n)) {
  for (col in c(1:n)) {
    # Test if the sexes are the same
    if (idcov_genet$Sex[row] == idcov_genet$Sex[col] & idcov_genet$Sex[row] == "M") {
      males[row,col] <- 1
    } else {
      males[row,col] <- 0
    }
  }
}
rownames(males)<-rownames(idcov_genet)
colnames(males)<-rownames(idcov_genet)

#females
females <- matrix(0, nrow=n, ncol=n)
# Loop through each row and each column in the data
for (row in c(1:n)) {
  for (col in c(1:n)) {
    # Test if the sexes are the same
    if (idcov_genet$Sex[row] == idcov_genet$Sex[col] & idcov_genet$Sex[row] == "F") {
      females[row,col] <- 1
    } else {
      females[row,col] <- 0
    }
  }
}
rownames(females)<-rownames(idcov_genet)
colnames(females)<-rownames(idcov_genet)

#male-female
mf<- matrix(0, nrow=n, ncol=n)
# Loop through each row and each column in the data
for (row in c(1:n)) {
  for (col in c(1:n)) {
    # Test if the sexes are the same
    if (idcov_genet$Sex[row] == idcov_genet$Sex[col]) {
      mf[row,col] <- 0
    } else {
      mf[row,col] <- 1
    }
  }
}
rownames(mf)<-rownames(idcov_genet)
colnames(mf)<-rownames(idcov_genet)



t.relatedness.m<-traits_by_related$related
# are male-male dyads more or less related than other dyads? (p is one-tailed !!)
(temp <- mantel(t.relatedness.m, males, permutations = 20000)) # 
MM.relatedness.stat <- temp$statistic
MM.relatedness.p <- temp$signif
MM.relatedness.number.sharks <- max(colSums(males))

# are female-female dyads more or less related than other dyads? (p is one-tailed !!)
(temp <- mantel(t.relatedness.m, females, permutations = 20000)) # 
FF.relatedness.stat <- temp$statistic
FF.relatedness.p <- temp$signif
FF.relatedness.number.sharks <- max(colSums(females))

# are mixed-sex dyads more or less related than other dyads? (p is one-tailed !!)
(temp <- mantel(t.relatedness.m, mf, permutations = 20000))# not different
FM.relatedness.stat <- temp$statistic
FM.relatedness.p <- temp$signif
FM.relatedness.number.sharks <- sum(unique(colSums(mf)))



### Subset matrix for MM, FF, MF
#males
onlymales<-males[rowSums(males[,-1]) != 0,]
onlymales<-onlymales[,colSums(onlymales[-1,]) != 0]
onlymales_R <- matrix_subset(t.relatedness.m, rownames(onlymales)) %>% 
  as.matrix() %>% 
  order_matrix()
onlymales_SRI <- matrix_subset(SRI_real, rownames(onlymales)) %>% 
  as.matrix() %>% 
  order_matrix()
onlymales_GAI <- matrix_subset(list_SYM_MATRIX_GAI_real, rownames(onlymales)) %>% 
  as.matrix() %>% 
  order_matrix()

#females
onlyfemales<-females[rowSums(females[,-1]) != 0,]
onlyfemales<-onlyfemales[,colSums(onlyfemales[-1,]) != 0]
onlyfemales_R <- matrix_subset(t.relatedness.m, rownames(onlyfemales)) %>% 
  as.matrix() %>% 
  order_matrix()
onlyfemales_SRI <- matrix_subset(SRI_real, rownames(onlyfemales)) %>% 
  as.matrix() %>% 
  order_matrix()
onlyfemales_GAI <- matrix_subset(list_SYM_MATRIX_GAI_real, rownames(onlyfemales)) %>% 
  as.matrix() %>% 
  order_matrix()

#male-female
onlymf<-males[rowSums(males[,-1]) != 0,]
onlymf<-onlymf[,colSums(onlymf) == 0]
onlymf_R <- t.relatedness.m[rownames(t.relatedness.m) %in% rownames(onlymf),colnames(t.relatedness.m) %in% colnames(onlymf)] %>% 
  as.matrix() %>% 
  order_matrix()
onlymf_SRI <- SRI_real[rownames(SRI_real) %in% rownames(onlymf),colnames(SRI_real) %in% colnames(onlymf)] %>% 
  as.matrix() %>% 
  order_matrix()
onlymf_GAI <- list_SYM_MATRIX_GAI_real[rownames(list_SYM_MATRIX_GAI_real) %in% rownames(onlymf),colnames(list_SYM_MATRIX_GAI_real) %in% colnames(onlymf)] %>% 
  as.matrix() %>% 
  order_matrix()

t.relatedness.m <- t.relatedness.m %>% 
  as.matrix() %>% 
  order_matrix()


mm_SRI_R_test_obs<-cor.test(matrix_unfold(onlymales_SRI),matrix_unfold(onlymales_R), method = "spearm")
ff_SRI_R_test_obs<-cor.test(matrix_unfold(onlyfemales_SRI),matrix_unfold(onlyfemales_R), method = "spearm")
mf_SRI_R_test_obs<-cor.test(matrix_unfold(onlymf_SRI),matrix_unfold(onlymf_R), method = "spearm")

mm_estimate_obs<-mm_SRI_R_test_obs$estimate
ff_estimate_obs<-ff_SRI_R_test_obs$estimate
mf_estimate_obs<-mf_SRI_R_test_obs$estimate

#subset sex-specific random SRI matrices
onlymales_SRI_rand<-array(0, dim=c(nrow(onlymales),ncol(onlymales),Npermute))
onlyfemales_SRI_rand<-array(0, dim=c(nrow(onlyfemales),ncol(onlyfemales),Npermute))
onlymf_SRI_rand<-array(0, dim=c(nrow(onlymf),ncol(onlymf),Npermute))
for (i in c(1:Npermute)) {
  onlymales_SRI_rand[,,i]<- list_RANDOM_SRI[[i]][rownames(list_RANDOM_SRI[[i]]) %in% rownames(onlymales),colnames(list_RANDOM_SRI[[i]]) %in% colnames(onlymales)] %>% 
    as.matrix() %>% 
    order_matrix()
  onlyfemales_SRI_rand[,,i]<- list_RANDOM_SRI[[i]][rownames(list_RANDOM_SRI[[i]]) %in% rownames(onlyfemales),colnames(list_RANDOM_SRI[[i]]) %in% colnames(onlyfemales)] %>% 
    as.matrix() %>% 
    order_matrix()
  onlymf_SRI_rand[,,i]<- list_RANDOM_SRI[[i]][rownames(list_RANDOM_SRI[[i]]) %in% rownames(onlymf),colnames(list_RANDOM_SRI[[i]]) %in% colnames(onlymf)] %>% 
    as.matrix() %>% 
    order_matrix()
}
rownames(onlymales_SRI_rand)<-rownames(onlymales)
colnames(onlymales_SRI_rand)<-colnames(onlymales)
rownames(onlyfemales_SRI_rand)<-rownames(onlyfemales)
colnames(onlyfemales_SRI_rand)<-colnames(onlyfemales)
rownames(onlymf_SRI_rand)<-rownames(onlymf)
colnames(onlymf_SRI_rand)<-colnames(onlymf)

#random correlation tests
mm_estimate_rand <- rep(0,Npermute)
ff_estimate_rand <- rep(0,Npermute)
mf_estimate_rand <- rep(0,Npermute)
for (i in c(1:Npermute)) {
  mm_estimate_rand[i]<-cor.test(matrix_unfold(onlymales_SRI_rand[,,i]),matrix_unfold(onlymales_R), method = "spearm")$estimate
  ff_estimate_rand[i]<-cor.test(matrix_unfold(onlyfemales_SRI_rand[,,i]),matrix_unfold(onlyfemales_R), method = "spearm")$estimate
  mf_estimate_rand[i]<-cor.test(matrix_unfold(onlymf_SRI_rand[,,i]),matrix_unfold(onlymf_R), method = "spearm")$estimate
}

## estimate and p value mantel SRI ~ R for MM
mm_estimate_obs
sum(abs(mm_estimate_obs) < abs(mm_estimate_rand))/Npermute
## estimate andp value mantel SRI ~ R for MM
ff_estimate_obs
sum(abs(ff_estimate_obs) < abs(ff_estimate_rand))/Npermute
## estimate andp value mantel SRI ~ R for MM
mf_estimate_obs
sum(abs(mf_estimate_obs) < abs(mf_estimate_rand))/Npermute



#### GAI
mm_GAI_R_test_obs<-cor.test(matrix_unfold(onlymales_GAI),matrix_unfold(onlymales_R), method = "spearm")
ff_GAI_R_test_obs<-cor.test(matrix_unfold(onlyfemales_GAI),matrix_unfold(onlyfemales_R), method = "spearm")
mf_GAI_R_test_obs<-cor.test(matrix_unfold(onlymf_GAI),matrix_unfold(onlymf_R), method = "spearm")

mm_estimate_obs_GAI<-mm_GAI_R_test_obs$estimate
ff_estimate_obs_GAI<-ff_GAI_R_test_obs$estimate
mf_estimate_obs_GAI<-mf_GAI_R_test_obs$estimate

#subset sex-specific random GAI matrices
onlymales_GAI_rand<-array(0, dim=c(nrow(onlymales),ncol(onlymales),Npermute))
onlyfemales_GAI_rand<-array(0, dim=c(nrow(onlyfemales),ncol(onlyfemales),Npermute))
onlymf_GAI_rand<-array(0, dim=c(nrow(onlymf),ncol(onlymf),Npermute))
for (i in c(1:Npermute)) {
  onlymales_GAI_rand[,,i]<- list_SYM_MATRIX_GAI_random[[i]][rownames(list_SYM_MATRIX_GAI_random[[i]]) %in% rownames(onlymales),colnames(list_SYM_MATRIX_GAI_random[[i]]) %in% colnames(onlymales)] %>% 
    as.matrix() %>% 
    order_matrix()
  onlyfemales_GAI_rand[,,i]<- list_SYM_MATRIX_GAI_random[[i]][rownames(list_SYM_MATRIX_GAI_random[[i]]) %in% rownames(onlyfemales),colnames(list_SYM_MATRIX_GAI_random[[i]]) %in% colnames(onlyfemales)] %>% 
    as.matrix() %>% 
    order_matrix()
  onlymf_GAI_rand[,,i]<- list_SYM_MATRIX_GAI_random[[i]][rownames(list_SYM_MATRIX_GAI_random[[i]]) %in% rownames(onlymf),colnames(list_SYM_MATRIX_GAI_random[[i]]) %in% colnames(onlymf)] %>% 
    as.matrix() %>% 
    order_matrix()
}
rownames(onlymales_GAI_rand)<-rownames(onlymales)
colnames(onlymales_GAI_rand)<-colnames(onlymales)
rownames(onlyfemales_GAI_rand)<-rownames(onlyfemales)
colnames(onlyfemales_GAI_rand)<-colnames(onlyfemales)
rownames(onlymf_GAI_rand)<-rownames(onlymf)
colnames(onlymf_GAI_rand)<-colnames(onlymf)

#random correlation tests
mm_estimate_rand_GAI <- rep(0,Npermute)
ff_estimate_rand_GAI <- rep(0,Npermute)
mf_estimate_rand_GAI <- rep(0,Npermute)
for (i in c(1:Npermute)) {
  mm_estimate_rand_GAI[i]<-cor.test(matrix_unfold(onlymales_GAI_rand[,,i]),matrix_unfold(onlymales_R), method = "spearm")$estimate
  ff_estimate_rand_GAI[i]<-cor.test(matrix_unfold(onlyfemales_GAI_rand[,,i]),matrix_unfold(onlyfemales_R), method = "spearm")$estimate
  mf_estimate_rand_GAI[i]<-cor.test(matrix_unfold(onlymf_GAI_rand[,,i]),matrix_unfold(onlymf_R), method = "spearm")$estimate
}

## estimate and p value mantel GAI ~ R for MM
mm_estimate_obs_GAI
sum(abs(mm_estimate_obs_GAI) < abs(mm_estimate_rand_GAI))/Npermute
## estimate andp value mantel GAI ~ R for MM
ff_estimate_obs_GAI
sum(abs(ff_estimate_obs_GAI) < abs(ff_estimate_rand_GAI))/Npermute
## estimate andp value mantel GAI ~ R for MM
mf_estimate_obs_GAI
sum(abs(mf_estimate_obs_GAI) < abs(mf_estimate_rand_GAI))/Npermute


##### Correlations R~SRI or GAI by sex dyads
GAIvsRmm<-data.frame(matrix_unfold(onlymales_R),matrix_unfold(onlymales_GAI),rep("MM",length(matrix_unfold(onlymales_R))),rep("GAI",length(matrix_unfold(onlymales_R))))
colnames(GAIvsRmm)<-c("Relatedness","Sociality","Sex.dyad","Index")
GAIvsRff<-data.frame(matrix_unfold(onlyfemales_R),matrix_unfold(onlyfemales_GAI),rep("FF",length(matrix_unfold(onlyfemales_R))),rep("GAI",length(matrix_unfold(onlyfemales_R))))
colnames(GAIvsRff)<-c("Relatedness","Sociality","Sex.dyad","Index")
GAIvsRmf<-data.frame(matrix_unfold(onlymf_R),matrix_unfold(onlymf_GAI),rep("MF",length(matrix_unfold(onlymf_R))),rep("GAI",length(matrix_unfold(onlymf_R))))
colnames(GAIvsRmf)<-c("Relatedness","Sociality","Sex.dyad","Index")

GAIvsR<-rbind(GAIvsRmm,GAIvsRff,GAIvsRmf)

SRIvsRmm<-data.frame(matrix_unfold(onlymales_R),matrix_unfold(onlymales_SRI),rep("MM",length(matrix_unfold(onlymales_R))),rep("SRI",length(matrix_unfold(onlymales_R))))
colnames(SRIvsRmm)<-c("Relatedness","Sociality","Sex.dyad","Index")
SRIvsRff<-data.frame(matrix_unfold(onlyfemales_R),matrix_unfold(onlyfemales_SRI),rep("FF",length(matrix_unfold(onlyfemales_R))),rep("SRI",length(matrix_unfold(onlyfemales_R))))
colnames(SRIvsRff)<-c("Relatedness","Sociality","Sex.dyad","Index")
SRIvsRmf<-data.frame(matrix_unfold(onlymf_R),matrix_unfold(onlymf_SRI),rep("MF",length(matrix_unfold(onlymf_R))),rep("SRI",length(matrix_unfold(onlymf_R))))
colnames(SRIvsRmf)<-c("Relatedness","Sociality","Sex.dyad","Index")

SRIvsR<-rbind(SRIvsRmm,SRIvsRff,SRIvsRmf)

SocialityR<-rbind(GAIvsR,SRIvsR)

#############
#Plot correlations
pdf(file = "Figure 3.pdf",
    height = 3.5, width = 8)
SocialityR$index = factor(SocialityR$Index, levels=c('SRI','GAI'))
ggplot(SocialityR, aes(x=Relatedness, y=Sociality, color=Sex.dyad)) +
  geom_point(alpha = 0.5, size = 3) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07", "#E7B800"))+
  geom_abline(slope = c(0.1028,-0.0653,0.0136,0.0404,0.0361,-0.0240), intercept = c(0,0,0,0,0,0))+
  facet_wrap(.~index, ncol = 2, scales = "free", dir = "h")+
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

dev.off()



##########                                  ############
############ Genetic relatedness analyses ##############
##########                                  ############
#  install.packages("related_1.0.tgz", repos=NULL, type="mac.binary")


#Comparing estimators
input<-readgenotypedata("Genotype(16markers).txt")

simdata <- familysim ( input$freqs , 100)
output <- coancestry ( simdata , dyadml =1 , trioml =1 , quellergt =1 , wang =1)
simrel <- cleanuprvals ( output$relatedness , 100)

triomlpo <- simrel [1:100 , 5]
triomlfs <- simrel [(100 + 1) : (2 * 100) , 5]
triomlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 5]
triomlur <- simrel [((3 * 100) + 1) : (4 * 100) , 5]
wangpo <- simrel [1:100 , 6]
wangfs <- simrel [(100 + 1) : (2 * 100) , 6]
wanghs <- simrel [((2 * 100) + 1) : (3 * 100) , 6]
wangur <- simrel [((3 * 100) + 1) : (4 * 100) , 6]
quellergtpo <- simrel [1:100 , 10]
quellergtfs <- simrel [(100 + 1) : (2 * 100) , 10]
quellergths <- simrel [((2 * 100) + 1) : (3 * 100) , 10]
quellergtur <- simrel [((3 * 100) + 1) : (4 * 100) , 10]
dyadmlpo <- simrel [1:100 , 11]
dyadmlfs <- simrel [(100 + 1) : (2 * 100) , 11]
dyadmlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 11]
dyadmlur <- simrel [((3 * 100) + 1) : (4 * 100) , 11]


trioml <- rep ("tri", 100)
wang <- rep ("W", 100)
quellergt <- rep ("QG", 100)
dyadml <- rep ("di", 100)
estimator2 <- c( trioml , wang , quellergt , dyadml )
Estimator <- rep ( estimator2 , 4)
po <- rep (" Parent - Offspring ", (4 * 100) )
fs <- rep ("Full - Sibs ", (4 * 100) )
hs <- rep ("Half - Sibs ", (4 * 100) )
ur <- rep (" Unrelated ", (4 * 100) )
relationship <- c(po , fs , hs , ur )

relatednesspo <- c( triomlpo , wangpo , quellergtpo , dyadmlpo )
relatednessfs <- c( triomlfs , wangfs , quellergtfs , dyadmlfs )
relatednesshs <- c( triomlhs , wanghs , quellergths , dyadmlhs )
relatednessur <- c( triomlur , wangur , quellergtur , dyadmlur )
Relatedness_Value <- c( relatednesspo , relatednessfs , relatednesshs , relatednessur )

combineddata <- as.data.frame ( cbind ( Estimator , relationship , Relatedness_Value ))
combineddata$Relatedness_Value <- as.numeric ( as.character (
  combineddata$Relatedness_Value ))

ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 ,
                                                                                   1.0) ) +
  geom_boxplot () +
  facet_wrap (~ relationship )

urval <- rep (0 , 100)
hsval <- rep (0.25 , 100)
fsval <- rep (0.5 , 100)
poval <- rep (0.5 , 100)
relvals <- c( poval , fsval , hsval , urval )

cor ( relvals , simrel [ , 5]) #Trioml
cor ( relvals , simrel [ , 6]) #
cor ( relvals , simrel [ , 10]) #
cor ( relvals , simrel [ , 11]) #dyadml

sd(simrel [ , 5]) #Trioml
sd(simrel [ , 6]) #Wang
sd(simrel [ , 10]) #QG
sd(simrel [ , 11]) #dyadml

# SD of relatedness values for each index


ddply(combineddata, .(Estimator, relationship), summarize,
      mean = round(mean(Relatedness_Value), 2),
      sd = round(sd(Relatedness_Value), 2))
###



relAll <- coancestry(input$gdata, trioml = 1)

#### Build Relatedness Matrix total
relmatrixAll <- matrix(0, nrow = 83, ncol = 83)
relmatrixAll[lower.tri(relmatrixAll, diag = FALSE)] <- relAll$relatedness[, 5]
individuals <- as.vector(input$gdata[, 1]) 
colnames(relmatrixAll) <- individuals
rownames(relmatrixAll) <- individuals
write.csv(relmatrixAll,"MatrixGeneticRelatedness_Trioml.csv")

#### Analazing relatedness within size classes
## create labels for groups
idcov_genet$label_size<-paste0(idcov_genet$SizeClass,idcov_genet$SizeClass,"_",rownames(idcov_genet))

datasize<-input$gdata

for (i in 1:nrow(datasize)){
  if (datasize[i,1] %in% rownames(idcov_genet)){
    datasize[i,1]<-as.character(idcov_genet[which(rownames(idcov_genet)==datasize[i,1]),7])
  }
  else {
    datasize[i,1]<-paste0("XX","_",datasize[i,1])
  }
}

grouprel(genotype = datasize, estimatorname = "trioml", usedgroups= c("BB","CC","DD","EE","FF"), iterations = 100)


Size_relatedness <- coancestry(datasize, trioml = 1)
Size_rel<-Size_relatedness$relatedness
Size_rel<-Size_rel[,c(2,3,4,5)]
substring="XX"
Size_rel<-Size_rel[-which(grepl(substring,Size_rel[,3])==TRUE),] #remove all "XX"

#mean within group relatedness 
WithinGroup<-Size_rel[which(Size_rel$group %in% c("AAAA","BBBB","CCCC","DDDD","EEEE","FFFF")),]
BetweenGroup<-Size_rel[-which(Size_rel$group %in% c("AAAA","BBBB","CCCC","DDDD","EEEE","FFFF")),]
#mean within group relatedness 
mean(WithinGroup$trioml)
sd(WithinGroup$trioml)
se(WithinGroup$trioml)
nrow(WithinGroup)
#mean between group relatedness 
mean(BetweenGroup$trioml)
sd(BetweenGroup$trioml)
se(BetweenGroup$trioml)
nrow(BetweenGroup)


#mean withing group BB 110-119
WithinSize_BB<-Size_rel[which(Size_rel$group == "BBBB"),]
mean(WithinSize_BB$trioml)
sd(WithinSize_BB$trioml)
length(unique(c(WithinSize_BB$ind1.id,WithinSize_BB$ind2.id)))
#mean withing group CC 120-129
WithinSize_CC<-Size_rel[which(Size_rel$group == "CCCC"),]
mean(WithinSize_CC$trioml)
sd(WithinSize_CC$trioml)
length(unique(c(WithinSize_CC$ind1.id,WithinSize_CC$ind2.id)))
#mean withing group DD 130-139
WithinSize_DD<-Size_rel[which(Size_rel$group == "DDDD"),]
mean(WithinSize_DD$trioml)
sd(WithinSize_DD$trioml)
length(unique(c(WithinSize_DD$ind1.id,WithinSize_DD$ind2.id)))
#mean withing group EE 140-149
WithinSize_EE<-Size_rel[which(Size_rel$group == "EEEE"),]
mean(WithinSize_EE$trioml)
sd(WithinSize_EE$trioml)
length(unique(c(WithinSize_EE$ind1.id,WithinSize_EE$ind2.id)))
#mean withing group FF >150
WithinSize_FF<-Size_rel[which(Size_rel$group == "FFFF"),]
mean(WithinSize_FF$trioml)
sd(WithinSize_FF$trioml)
length(unique(c(WithinSize_FF$ind1.id,WithinSize_FF$ind2.id)))

## Analazing relatedness within communities global

# with SRI
length(community.SRI[[1]])
#Create labels for groups (AA = Community 1, BB = Community 2 and CC = Community 3)
Com1<-as.vector(unlist(community.SRI[1]))
Com2<-as.vector(unlist(community.SRI[2]))
Com3<-as.vector(unlist(community.SRI[3]))
datacomm<-input$gdata
for (i in 1:nrow(datacomm)){
  if (datacomm[i,1] %in% Com1){
    datacomm[i,1]<-paste0("AA","_",datacomm[i,1])
  }
  else {
    if (datacomm[i,1] %in% Com2){
      datacomm[i,1]<-paste0("BB","_",datacomm[i,1])
    }
    else{
      if (datacomm[i,1] %in% Com3){
        datacomm[i,1]<-paste0("CC","_",datacomm[i,1])
      }
      else{
        datacomm[i,1]<-paste0("XX","_",datacomm[i,1])
      }
    }
  }
}

substring="XX"
datacomm<-datacomm[-which(grepl(substring,datacomm[,1])==TRUE),] #remove all "XX"


grouprel(genotype = datacomm, estimatorname = "trioml", usedgroups= c("AA","BB","CC"), iterations = 1000)

SRIwithinCom1<-(SRI_real[which(rownames(SRI_real) %in% Com1),which(colnames(SRI_real) %in% Com1)])
SRIwithinCom2<-(SRI_real[which(rownames(SRI_real) %in% Com2),which(colnames(SRI_real) %in% Com2)])
SRIwithinCom3<-(SRI_real[which(rownames(SRI_real) %in% Com3),which(colnames(SRI_real) %in% Com3)])

se <- function(x) {sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}

#mean SRI 
mean(SRI_real)
sd(SRI_real)
nrow(SRI_real)

#mean SRI within com1 
mean(SRIwithinCom1)
sd(SRIwithinCom1)
nrow(SRIwithinCom1)
#mean SRI within com2 
mean(SRIwithinCom2)
sd(SRIwithinCom2)
nrow(SRIwithinCom2)
#mean SRI within com3 
mean(SRIwithinCom3)
sd(SRIwithinCom3)
nrow(SRIwithinCom3)

SRI_relatedness <- coancestry(datacomm, trioml = 1)
SRI_rel<-SRI_relatedness$relatedness
SRI_rel<-SRI_rel[,c(2,3,4,5)]
substring="XX"
SRI_rel<-SRI_rel[-which(grepl(substring,SRI_rel[,3])==TRUE),] #remove all "XX"

#mean within group relatedness 
WithinGroup<-SRI_rel[which(SRI_rel$group %in% c("AAAA","BBBB","CCCC")),]
BetweenGroup<-SRI_rel[-which(SRI_rel$group %in% c("AAAA","BBBB","CCCC")),]
#mean within group relatedness 
mean(WithinGroup$trioml)
sd(WithinGroup$trioml)
se(WithinGroup$trioml)
nrow(WithinGroup)
#mean between group relatedness 
mean(BetweenGroup$trioml)
sd(BetweenGroup$trioml)
se(BetweenGroup$trioml)
nrow(BetweenGroup)

#mean withing group AA
WithinSRI_AA<-SRI_rel[which(SRI_rel$group == "AAAA"),]
mean(WithinSRI_AA$trioml)
se(WithinSRI_AA$trioml)
#mean withing group BB
WithinSRI_BB<-SRI_rel[which(SRI_rel$group == "BBBB"),]
mean(WithinSRI_BB$trioml)
se(WithinSRI_BB$trioml)
#mean withing group CC
WithinSRI_CC<-SRI_rel[which(SRI_rel$group == "CCCC"),]
mean(WithinSRI_CC$trioml)
se(WithinSRI_CC$trioml)

# with GAI
length(community.GAI[[1]])
#Create labels for groups (AA = Community 1, BB = Community 2 and CC = Community 3)
Com1<-as.vector(unlist(community.GAI[1]))
Com2<-as.vector(unlist(community.GAI[2]))
Com3<-as.vector(unlist(community.GAI[3]))
Com4<-as.vector(unlist(community.GAI[4]))
#Com5<-as.vector(unlist(community.GAI[5]))

datacomm<-input$gdata
for (i in 1:nrow(datacomm)){
  if (datacomm[i,1] %in% Com1){
    datacomm[i,1]<-paste0("AA","_",datacomm[i,1])
  }
  else {
    if (datacomm[i,1] %in% Com2){
      datacomm[i,1]<-paste0("BB","_",datacomm[i,1])
    }
    else{
      if (datacomm[i,1] %in% Com3){
        datacomm[i,1]<-paste0("CC","_",datacomm[i,1])
      }
      else{
        if (datacomm[i,1] %in% Com4){
          datacomm[i,1]<-paste0("DD","_",datacomm[i,1])
        }
        #else{
        #  if (datacomm[i,1] %in% Com5){
        #    datacomm[i,1]<-paste0("EE","_",datacomm[i,1])
        #  }
        else{
          datacomm[i,1]<-paste0("XX","_",datacomm[i,1])
        }
        # }
      }
    }
  }
}

substring="XX"
datacomm1<-datacomm[-which(grepl(substring,datacomm[,1])==TRUE),] #remove all "XX"

grouprel(genotype = datacomm1, estimatorname = "trioml", usedgroups= c("AA","BB","CC","DD"), iterations = 1000)

GAIwithinCom1<-(listMATRIX_GAI_real[which(rownames(listMATRIX_GAI_real) %in% Com1),which(colnames(listMATRIX_GAI_real) %in% Com1)])
GAIwithinCom2<-(listMATRIX_GAI_real[which(rownames(listMATRIX_GAI_real) %in% Com2),which(colnames(listMATRIX_GAI_real) %in% Com2)])
GAIwithinCom3<-(listMATRIX_GAI_real[which(rownames(listMATRIX_GAI_real) %in% Com3),which(colnames(listMATRIX_GAI_real) %in% Com3)])
GAIwithinCom4<-(listMATRIX_GAI_real[which(rownames(listMATRIX_GAI_real) %in% Com4),which(colnames(listMATRIX_GAI_real) %in% Com4)])
#GAIwithinCom5<-(listMATRIX_GAI_real[which(rownames(listMATRIX_GAI_real) %in% Com5),which(colnames(listMATRIX_GAI_real) %in% Com5)])

#mean GAI 
mean(listMATRIX_GAI_real)
sd(listMATRIX_GAI_real)
nrow(listMATRIX_GAI_real)

#mean GAI within com1 
mean(GAIwithinCom1)
sd(GAIwithinCom1)
nrow(GAIwithinCom1)
#mean GAI within com2 
mean(GAIwithinCom2)
sd(GAIwithinCom2)
nrow(GAIwithinCom2)
#mean GAI within com3 
mean(GAIwithinCom3)
sd(GAIwithinCom3)
nrow(GAIwithinCom3)
#mean GAI within com4 
mean(GAIwithinCom4)
sd(GAIwithinCom4)
nrow(GAIwithinCom4)
#mean GAI within com5 
#mean(GAIwithinCom5)
#sd(GAIwithinCom5)
#nrow(GAIwithinCom5)

GAI_relatedness <- coancestry(datacomm, trioml = 1)
GAI_rel<-GAI_relatedness$relatedness
GAI_rel<-GAI_rel[,c(2,3,4,5)]
substring="XX"
GAI_rel<-GAI_rel[-which(grepl(substring,GAI_rel[,3])==TRUE),] #remove all "XX"

#mean within group relatedness 
WithinGroupGAI<-GAI_rel[which(GAI_rel$group %in% c("AAAA","BBBB","CCCC","DDDD")),]
BetweenGroupGAI<-GAI_rel[-which(GAI_rel$group %in% c("AAAA","BBBB","CCCC","DDDD")),]
#mean within group relatedness 
mean(WithinGroupGAI$trioml)
sd(WithinGroupGAI$trioml)
se(WithinGroupGAI$trioml)
nrow(WithinGroupGAI)
#mean between group relatedness 
mean(BetweenGroupGAI$trioml)
sd(BetweenGroupGAI$trioml)
se(WithinGroupGAI$trioml)
nrow(BetweenGroupGAI)

#mean withing group AA
WithinGAI_BB<-GAI_rel[which(GAI_rel$group == "AAAA"),]
mean(WithinGAI_BB$trioml)
se(WithinGAI_BB$trioml)
#mean withing group BB
WithinGAI_BB<-GAI_rel[which(GAI_rel$group == "BBBB"),]
mean(WithinGAI_BB$trioml)
se(WithinGAI_BB$trioml)
#mean withing group CC
WithinGAI_CC<-GAI_rel[which(GAI_rel$group == "CCCC"),]
mean(WithinGAI_CC$trioml)
se(WithinGAI_CC$trioml)
#mean withing group DD
WithinGAI_DD<-GAI_rel[which(GAI_rel$group == "DDDD"),]
mean(WithinGAI_DD$trioml)
se(WithinGAI_DD$trioml)
#mean withing group EE
#WithinGAI_EE<-GAI_rel[which(GAI_rel$group == "EEEE"),]
#mean(WithinGAI_EE$trioml)
#se(WithinGAI_EE$trioml)



KinFreq<-function (genotypes, estimatorname, usedgroups, iterations, threshold) 
{
  if (estimatorname == "trioml") {
    estimator = 5
  }
  if (estimatorname == "wang") {
    estimator = 6
  }
  if (estimatorname == "lynchli") {
    estimator = 7
  }
  if (estimatorname == "lynchrd") {
    estimator = 8
  }
  if (estimatorname == "ritland") {
    estimator = 9
  }
  if (estimatorname == "quellergt") {
    estimator = 10
  }
  if (estimatorname == "dyadml") {
    estimator = 11
  }
  if (estimatorname == "trioml") {
    relatives <- coancestry(genotypes, trioml = 1)
  }
  if (estimatorname == "wang") {
    relatives <- coancestry(genotypes, wang = 1)
  }
  if (estimatorname == "lynchli") {
    relatives <- coancestry(genotypes, lynchli = 1)
  }
  if (estimatorname == "lynchrd") {
    relatives <- coancestry(genotypes, lynchrd = 1)
  }
  if (estimatorname == "ritland") {
    relatives <- coancestry(genotypes, ritland = 1)
  }
  if (estimatorname == "quellergt") {
    relatives <- coancestry(genotypes, quellergt = 1)
  }
  if (estimatorname == "dyadml") {
    relatives <- coancestry(genotypes, dyadml = 1)
  }
  rels <- relatives$relatedness[, c(2,3,4,estimator)]
  
  groups <- usedgroups
  within <- paste(groups, groups, sep = "")
  relvalues <- 1:length(within)
  relvaluesinf<-1:length(within)
  freqvalues<- 1:length(within)
  cat("\n Calculating within-group r-values above thresold...\n")
  for (i in 1:length(within)) {
    relvalues[i]<-dim(subset(rels, group == within[i] & rels[,estimator-1] >= threshold))[1]
  }
  for (i in 1:length(within)) {
    relvaluesinf[i]<-dim(subset(rels, group == within[i] & rels[,estimator-1] < threshold))[1]
  }
  
  
  withingroup_sup<-sum(relvalues)
  withingroup_inf<-sum(relvaluesinf)
  betweengroup_sup<-dim(subset(rels, rels[,estimator-1] >= threshold))[1]-sum(relvalues)
  betweengroup_inf<-dim(subset(rels, rels[,estimator-1] < threshold))[1]-sum(relvaluesinf)
  tab = matrix(c(c(withingroup_sup,withingroup_inf),c(betweengroup_sup,betweengroup_inf)),2,2,byrow=T)
  test_obs<-chisq.test(tab)$statistic
  df<-chisq.test(tab)$parameter
  p<- chisq.test(tab)$p.value
  cat("\n Number of between groups...\n", test_obs)
  cat("\n Number of within groups above thresold...\n",withingroup_sup)
  cat("\n Number of between groups above thresold...\n",betweengroup_sup)
  cat("\n Number of within groups under thresold...\n",withingroup_inf)
  cat("\n Number of between groups under thresold...\n",betweengroup_inf)
  cat("\n Total Number of relationships...\n", dim(rels)[1])
  cat("\n Statistic of Chi test ...\n",  test_obs)
  cat("\n Degree of freedom of Chi test ...\n",  df)
  cat("\n p value of Chi test ...\n",  p)
  
  
  simresults_sup <- data.frame(matrix(nrow = iterations, ncol = (length(within))))
  simresults_inf <- data.frame(matrix(nrow = iterations, ncol = (length(within))))
  withingroup_sim_sup<-data.frame(matrix(nrow = iterations, ncol = 1))
  withingroup_sim_inf <- data.frame(matrix(nrow = iterations, ncol = 1))
  betweengroup_sim_sup<-data.frame(matrix(nrow = iterations, ncol = 1))
  betweengroup_sim_inf <- data.frame(matrix(nrow = iterations, ncol = 1))
  tab2 <- array(rep(1, 2*2*iterations), dim=c(2, 2, iterations))
  test_sim<-1:iterations
  for (j in 1:iterations) {
    cat(sprintf("Iteration %d\n", j))
    randlist <- 1:length(genotypes[, 1])
    randlist2 <- sample(randlist, length(randlist), replace = FALSE)
    randgenos <- data.frame(matrix(nrow = length(randlist), 
                                   ncol = length(genotypes[1, ])))
    for (i in 1:length(randlist)) {
      randgenos[i, ] <- genotypes[randlist2[i], ]
      randgenos[,1]<-genotypes[,1]
    }
    if (estimatorname == "trioml") {
      simrels <- coancestry(randgenos, trioml = 1)
    }
    if (estimatorname == "wang") {
      simrels <- coancestry(randgenos, wang = 1)
    }
    if (estimatorname == "lynchli") {
      simrels <- coancestry(randgenos, lynchli = 1)
    }
    if (estimatorname == "lynchrd") {
      simrels <- coancestry(randgenos, lynchrd = 1)
    }
    if (estimatorname == "ritland") {
      simrels <- coancestry(randgenos, ritland = 1)
    }
    if (estimatorname == "quellergt") {
      simrels <- coancestry(randgenos, quellergt = 1)
    }
    if (estimatorname == "dyadml") {
      simrels <- coancestry(randgenos, dyadml = 1)
    }
    
    relsim<-simrels$relatedness[,c(2,3,4,estimator)]
    
    for (k in 1:length(within)) {
      simresults_sup[j,k]<-dim(subset(relsim,group==within[k] & relsim[,estimator-1] >=threshold))[1]
    }
    for (k in 1:length(within)) {
      simresults_inf[j,k]<-dim(subset(relsim,group==within[k] & relsim[,estimator-1] < threshold))[1]
    }
    withingroup_sim_sup[j,]<-rowSums(simresults_sup[j,])
    withingroup_sim_inf[j,]<-rowSums(simresults_inf[j,])
    betweengroup_sim_sup[j,]<-dim(subset(relsim,relsim[,estimator-1] >=threshold))[1]
    betweengroup_sim_inf[j,]<-dim(subset(relsim,relsim[,estimator-1] < threshold))[1]
    
    tab2[,,j] = matrix(c(c(withingroup_sim_sup[j,],withingroup_sim_inf[j,]),c(betweengroup_sim_sup[j,],betweengroup_sim_inf[j,])),2,2,byrow=T)
    test_sim[j]<-chisq.test(tab2[,,j])$statistic
    
  }
  
  #write.csv(test_sim, file="/Users/johannmourier/Documents/Scripts/BTRS SocialNetwork Genetic/Chi2_sim.csv")
  
  
  cat("\n Number of within groups above thresold...\n",withingroup_sup)
  cat("\n Number of between groups above thresold...\n",betweengroup_sup)
  cat("\n Number of within groups under thresold...\n",withingroup_inf)
  cat("\n Number of between groups under thresold...\n",betweengroup_inf)
  cat("\n Total Number of relationships...\n", dim(rels)[1])
  cat("\n Statistic of Chi test ...\n",  test_obs)
  cat("\n Degree of freedom of Chi test ...\n",  df)
  cat("\n p value of Chi test ...\n",  p)
  
  
  nb<-length(test_sim[test_sim>=test_obs])
  cat("\n number of simulated statistics superior to observed statistics ...\n",  nb)
  #P-value
  p_val<-nb/iterations
  cat("\n p value of Chi test from simulations...\n",  p_val)
  
  #tab2 = matrix(c(c(withingroup_sup,withingroup_inf),c(betweengroup_sup,betweengroup_inf)),2,2,byrow=T)
  #test_obs<-chisq.test(tab2)$statistic
  #df<-chisq.test(tab)$parameter
  #p<- chisq.test(tab)$p.value
  
}



##############
## Analazing frequency of r>0.25 within communities global

# with SRI
length(community.SRI[[1]])
#Create labels for groups (AA = Community 1, BB = Community 2 and CC = Community 3)
Com1<-as.vector(unlist(community.SRI[1]))
Com2<-as.vector(unlist(community.SRI[2]))
Com3<-as.vector(unlist(community.SRI[3]))
datacomm<-input$gdata
for (i in 1:nrow(datacomm)){
  if (datacomm[i,1] %in% Com1){
    datacomm[i,1]<-paste0("AA","_",datacomm[i,1])
  }
  else {
    if (datacomm[i,1] %in% Com2){
      datacomm[i,1]<-paste0("BB","_",datacomm[i,1])
    }
    else{
      if (datacomm[i,1] %in% Com3){
        datacomm[i,1]<-paste0("CC","_",datacomm[i,1])
      }
      else{
        datacomm[i,1]<-paste0("XX","_",datacomm[i,1])
      }
    }
  }
}


substring="XX"
datacomm1<-datacomm[-which(grepl(substring,datacomm[,1])==TRUE),] #remove all "XX"


#grouprel(genotype = datacomm, estimatorname = "trioml", usedgroups= c("AA","BB","CC"), iterations = 500)
KinFreq(genotypes = datacomm1, estimatorname = "trioml", usedgroups= c("AA","BB","CC"), iterations = 1000, threshold = 0.25)


# with GAI
length(community.GAI[[1]])
#Create labels for groups (AA = Community 1, BB = Community 2 and CC = Community 3)
Com1<-as.vector(unlist(community.GAI[1]))
Com2<-as.vector(unlist(community.GAI[2]))
Com3<-as.vector(unlist(community.GAI[3]))
Com4<-as.vector(unlist(community.GAI[4]))
#Com5<-as.vector(unlist(community.GAI[5]))

datacomm<-input$gdata
for (i in 1:nrow(datacomm)){
  if (datacomm[i,1] %in% Com1){
    datacomm[i,1]<-paste0("AA","_",datacomm[i,1])
  }
  else {
    if (datacomm[i,1] %in% Com2){
      datacomm[i,1]<-paste0("BB","_",datacomm[i,1])
    }
    else{
      if (datacomm[i,1] %in% Com3){
        datacomm[i,1]<-paste0("CC","_",datacomm[i,1])
      }
      else{
        if (datacomm[i,1] %in% Com4){
          datacomm[i,1]<-paste0("DD","_",datacomm[i,1])
        }
        #else{
        #if (datacomm[i,1] %in% Com5){
        #  datacomm[i,1]<-paste0("EE","_",datacomm[i,1])
        #}
        else{
          datacomm[i,1]<-paste0("XX","_",datacomm[i,1])
        }
        # }
      }
    }
  }
}


substring="XX"
datacomm1<-datacomm[-which(grepl(substring,datacomm[,1])==TRUE),] #remove all "XX"

#grouprel(genotype = datacomm, estimatorname = "trioml", usedgroups= c("BB","CC","DD","EE"), iterations = 500)
KinFreq(genotypes = datacomm1, estimatorname = "trioml", usedgroups= c("AA","BB","CC","DD"), iterations = 1000, threshold = 0.25)


####  Figure 4

table<-data.frame(rbind(c("SRI","Within community","R>0.25", NA),c("SRI","Between community","R>0.25", NA),
                        c("SRI","Within community","R<0.25", NA),c("SRI","Between community","R<0.25", NA),
                        c("GAI","Within community","R>0.25", NA),c("GAI","Between community","R>0.25", NA),
                        c("GAI","Within community","R<0.25", NA),c("GAI","Between community","R<0.25", NA)))
colnames(table)<-c("Index","Relationship","Relatedness","NumberPairs")
table$NumberPairs<-c(17,22,300,564,8,35,206,654)
table$Relationship<-factor(table$Relationship, levels = rev(levels(table$Relationship)))

table2<-table[which(table$Relatedness=="R>0.25"),]

#Plot
ggplot(data=table2, aes(x=Index, y=NumberPairs, fill=Relationship)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired",direction=-1)+
  scale_x_discrete(limits=c("SRI","GAI"))+
  labs(y = "Number of pairs with r > 0.25")




######### density plots of space use by communities
SpaceProfile=read.table("CountSightings.csv",sep=",",header=T)
colnames(SpaceProfile)<-c("ID","1","2","3","4","5","6","7")
library(reshape)
SpPro<-melt(SpaceProfile, id="ID")
rownames(SpPro)<-SpPro$ID
SpPro<-SpPro[SpPro$ID %in% rownames(traits_by_related$related),]
colnames(SpPro)<-c("ID","Site","Sightings")

## SRI community membership
Com1SRI<-as.vector(unlist(community.SRI[1]))
Com2SRI<-as.vector(unlist(community.SRI[2]))
Com3SRI<-as.vector(unlist(community.SRI[3]))

SRI_com_space<-SpPro
SRI_com_space$Community<-rep(0,nrow(SRI_com_space))
for (i in 1:nrow(SRI_com_space)){
  if (SRI_com_space[i,1] %in% Com1SRI){
    SRI_com_space[i,4]<-"Community 1"
  }
  else{
    if (SRI_com_space[i,1] %in% Com2SRI){
      SRI_com_space[i,4]<-"Community 2"
    }
    else{
      SRI_com_space[i,4]<-"Community 3"
    }
  }
}

## GAI community membership
Com1GAI<-as.vector(unlist(community.GAI[1]))
Com2GAI<-as.vector(unlist(community.GAI[2]))
Com3GAI<-as.vector(unlist(community.GAI[3]))
Com4GAI<-as.vector(unlist(community.GAI[4]))
#Com5GAI<-as.vector(unlist(community.GAI[5]))

GAI_com_space<-SpPro
GAI_com_space$Community<-rep(0,nrow(GAI_com_space))
for (i in 1:nrow(GAI_com_space)){
  if (GAI_com_space[i,1] %in% Com1GAI){
    GAI_com_space[i,4]<-"Community 1"
  }
  else{
    if (GAI_com_space[i,1] %in% Com2GAI){
      GAI_com_space[i,4]<-"Community 2"
    }
    else{
      if (GAI_com_space[i,1] %in% Com3GAI){
        GAI_com_space[i,4]<-"Community 3"
      }
      else{
        if (GAI_com_space[i,1] %in% Com4GAI){
          GAI_com_space[i,4]<-"Community 4"
        }
        #else{
        #  GAI_com_space[i,4]<-"Community 5"
        #}
      }
    }
  }
}

SRI_com_space2<-na.omit(SRI_com_space)
GAI_com_space2<-na.omit(GAI_com_space)
SRI_com_space2 <- SRI_com_space2[rep(row.names(GAI_com_space2), SRI_com_space2$Sightings),]
GAI_com_space2 <- GAI_com_space2[rep(row.names(GAI_com_space2), GAI_com_space2$Sightings),]


pal<-fish(5, option = "Sparisoma_viride")
# SRI
ggplot(SRI_com_space2, aes(x = as.numeric(Site), y = Community, fill = Community)) +
  geom_density_ridges() +
  scale_x_continuous(breaks=c(1:7))+
  scale_fill_manual(values= pal, guide = "none")+
  theme_ridges() + 
  labs(x = "Sites")+
  theme(legend.position = "none")

#GAI
ggplot(GAI_com_space2, aes(x = as.numeric(as.character(Site)), y = Community, fill = Community)) +
  geom_density_ridges() +
  scale_x_continuous(breaks=c(1:7))+
  scale_fill_manual(values= pal, guide = "none")+
  theme_ridges() + 
  labs(x = "Sites")+
  theme(legend.position = "none")

tab2 = matrix(c(c(8,206),c(35,654)),2,2,byrow=T)
test_obs<-chisq.test(tab2)$statistic
df<-chisq.test(tab2)$parameter
p<- chisq.test(tab2)$p.value

tab2 = matrix(c(c(11,218),c(29,645)),2,2,byrow=T)
test_obs<-chisq.test(tab2)$statistic
df<-chisq.test(tab2)$parameter
p<- chisq.test(tab2)$p.value

##########################################################################

####### Differences males vs females


#CALC MEAN MALES
network_ALL_M<-SRI_real[which(idcov_genet$Sex=="M"),which(idcov_genet$Sex=="M")]
dim(network_ALL_M)
mean(network_ALL_M[network_ALL_M>0])

#CALC MEAN FEMALES
network_ALL_F<-SRI_real[which(idcov_genet$Sex=="F"),which(idcov_genet$Sex=="F")]
dim(network_ALL_F)
mean(network_ALL_F[network_ALL_F>0])

t.test(network_ALL_M[network_ALL_M>0], network_ALL_F[network_ALL_F>0])
#Stored test statsitcs (t-value) real data
t.obs<-t.test(network_ALL_M[network_ALL_M>0], network_ALL_F[network_ALL_F>0])$statistic
t.obs

#create t-random
t.random<-rep(0,20000)

#doing loop

for (i in c(1:20000)) {
  network_male_rand <- randomSRI[i,which(idcov_genet$Sex=="M"), which(idcov_genet$Sex=="M")]
  network_female_rand <- randomSRI[i,which(idcov_genet$Sex=="F"), which(idcov_genet$Sex=="F")]
  t.random[i] <- t.test(network_female_rand[network_female_rand>0],
                        network_male_rand[network_male_rand>0])$statistic
}
hist(t.random,breaks=100)
abline(v=t.obs,col="red")

#get p-values
sum(t.obs<t.random)/20000

#Manually measure weighted degree
strengthALL<-rowSums(SRI_real)

#Measure binary degree
degreeALL<-SRI_real
degreeALL[degreeALL>0]<-1
deg_binaryALL<-rowSums(degreeALL)

#add degrees to attributes data
idcov_genet$strength<-as.numeric(strengthALL)
idcov_genet$degree<-deg_binaryALL

# Finally set up a 3 panel figure and plot each relationship ignoring NAs (from Farine & Whitehead 2015 sup materials p15)
#multi graph function
par(mfrow=c(1,3))
boxplot(idcov_genet$strength~idcov_genet$Sex, col=c("red","blue"))
boxplot(idcov_genet$degree~idcov_genet$Sex, col=c("red","blue"))
plot(idcov_genet$strength~idcov_genet$degree, col=c("red","blue")[idcov_genet$Sex])
legend("bottomright",c("Female","Male"),col=c("red","blue"),pch=1)

# GLM...(3rd line degree is strength, but refered to as degree within program function)

########### SRI ##############
#STRENGTH

netsex<-SRI_real
rownames(netsex)<-idcov_genet$Sex
colnames(netsex)<-idcov_genet$Sex
idcov_genet$strength<-degree(netsex,gmode="graph")
coef<-coef(glm(strength~Sex, data = idcov_genet))[2]
strength_rand<-apply(randomSRI,1,function(x) {degree(x,gmode = "graph")})
coefs<-apply(strength_rand,2,function(x) {coef(glm(x~Sex, data = idcov_genet))[2]})

#plot results and save as 4x11 inches (landscape) (same as box 4 in Farine & Whitehead 2015)
par(mfrow=c(1,2),cex.lab=1.5)
#plot network
nodecolor=as.character(idcov_genet$Sex)
nodecolor=gsub("M","red",nodecolor)
nodecolor=gsub("F","blue",nodecolor)
#plot(graph.adjacency(netsex,mode="undirected", diag=FALSE, weighted = TRUE, add.rownames = "code"),vertex.color=nodecolor,vertex.size=attribute$strength*4)
#plot observed difference
plot(idcov_genet$strength~factor(idcov_genet$Sex),col=c("blue","red"),xlab="Sex",ylab="Strength (weighted degree)",ylim=c(0,max(idcov_genet$strength)))
#Plot resulting distribution
a<-hist(coefs,xlim=c(min(coefs),max(coefs)), col="black",main="",xlab="Coefficient value",ylab="Frequency",breaks=100)
abline(v=coef, col="red")
#Model results
summary(glm(strength~Sex, data = idcov_genet))
#P value
sum(abs(coef)<abs(coefs))/20000


##########################################################
# double permutation procedure by Farine and Carter 2020 #
##########################################################

# Calculate effects
model <- lm(strength~Sex,data=idcov_genet)
coef.trait <- coefficients(model)[2]
p.model <- summary(model)$coefficients[2,4]
t.trait <- summary(model)$coefficients[2,3]

# Calculate effects
model <- lm(strength~Sex+OBS,data=idcov_genet)
coef.trait.obs.control <- coefficients(model)[2]
p.model.obs.control <- summary(model)$coefficients[2,4]



strength_rand <- apply(randomSRI,1,function(x) { rowSums(x)})

# function to get coefficients for each randomisation
get_res <- function(x) {
  model <- lm(x~Sex,data=idcov_genet)
  coef <- coefficients(model)[2]
  return(coef)
}
# function to get coefficients for each randomisation
get_res2 <- function(x) {
  model <- lm(x~Sex,data=idcov_genet)
  coef <- summary(model)$coefficients[2,3]
  return(coef)
}


coefs.data<-rep(NA, Npermute)
t.data<-rep(NA, Npermute)
p.data<-rep(NA, Npermute)
p.t.data<-rep(NA, Npermute)
# coefficients based on original degree
coefs.data1 <- apply(strength_rand,2,function(x) { get_res(x) })
t.data <- apply(strength_rand,2,function(x) { get_res2(x) })

# get controlled effect size
effect.size <- coef.trait - median(coefs.data)


# controlled degree (after pre-network permutation)
idcov_genet$DEGREE.cont <- idcov_genet$strength - apply(strength_rand,1,median)
# Calculate effect controlling for OBS
model <- lm(DEGREE.cont~Sex,data=idcov_genet)
coef.controlled <- coefficients(model)[2]
p.model.controlled <- summary(model)$coefficients[2,4]

idcov_genet$OBS<-c(18,17,22,22,17,17,23,19,20,19,20,18,23,26,32,24,22,18,17,35,16,20,35,16,20,25,28,18,32,19,18,16,22,16,39,30,27,30,15,23,24,16,23)

## Node permutations
coefs.tmp <- rep(NA,Npermute)
coefs.tmp.control <- rep(NA,Npermute)
coefs.tmp.obs.control <- rep(NA,Npermute)
for (i in 1:Npermute) {
  idcov_genet$Sex.tmp <- sample(idcov_genet$Sex)
  model.tmp <- lm(strength~Sex,data=idcov_genet)
  coefs.tmp[i] <- coefficients(model.tmp)[2]
  model.tmp <- lm(DEGREE.cont~Sex.tmp,data=idcov_genet)
  coefs.tmp.control[i] <- coefficients(model.tmp)[2]
  model.tmp <- lm(strength~Sex.tmp+OBS,data=idcov_genet)
  coefs.tmp.obs.control[i] <- coefficients(model.tmp)[2]
}

# P values
p.data <- sum(coef.trait<=coefs.data)/Npermute
p.node <- sum(coef.trait<=coefs.tmp)/Npermute
p.node.control <- sum(coef.controlled <= coefs.tmp.control)/Npermute
p.node.obs.control <- sum(coef.trait.obs.control <= coefs.tmp.obs.control)/Npermute
p.t <- sum(t.trait <= t.data)/Npermute




########## Centrality analyses ############
#DEGREE (binary)
#Measure binary degree

idcov_genet$degree_binary <- degree(netsex, gmode="graph",ignore.eval=TRUE)
coef2<-coef(glm(degree_binary~Sex, data = idcov_genet))[2]
degree_binary_rand <- apply(randomSRI,1,function(x) {degree(x,gmode = "graph",ignore.eval=TRUE)})
detach("package:sna")
coefs2<-apply(degree_binary_rand,2,function(x) {coef(glm(x~Sex, data = idcov_genet))[2]})

#plot results and save as 4x11 inches (landscape) (same as box 4 in Farine & Whitehead 2015)
par(mfrow=c(1,2),cex.lab=1.5)
#plot network
nodecolor=as.character(idcov_genet$Sex)
nodecolor=gsub("M","red",nodecolor)
nodecolor=gsub("F","blue",nodecolor)
#plot(graph.adjacency(netsex,mode="undirected", diag=FALSE, weighted = TRUE, add.rownames = "code"),vertex.color=nodecolor,vertex.size=attribute$degree_binary*0.5)
#plot observed difference
plot(idcov_genet$degree_binary~factor(idcov_genet$Sex),col=c("blue","red"),xlab="Sex",ylab="Degree (binary degree)",ylim=c(0,max(idcov_genet$degree_binary)))
#Plot resulting distribution
a<-hist(coefs2,xlim=c(min(coefs2),max(coef2)), col="black",main="",xlab="Coefficient value",ylab="Frequency",breaks=100)
abline(v=coef2, col="red")
#Model results
summary(glm(degree_binary~Sex, data = idcov_genet))
#P value
sum(abs(coef2)<abs(coefs2))/20000


#double permutation procedure

# Calculate effects
model <- lm(degree_binary~Sex,data=idcov_genet)
coef.trait <- coefficients(model)[2]
p.model <- summary(model)$coefficients[2,4]
t.trait <- summary(model)$coefficients[2,3]

# Calculate effects
model <- lm(degree_binary~Sex+OBS,data=idcov_genet)
coef.trait.obs.control <- coefficients(model)[2]
p.model.obs.control <- summary(model)$coefficients[2,4]



coefs.data<-rep(NA, Npermute)
t.data<-rep(NA, Npermute)
p.data<-rep(NA, Npermute)
p.t.data<-rep(NA, Npermute)
# coefficients based on original degree
coefs.data <- apply(degree_binary_rand,2,function(x) { get_res(x) })
t.data <- apply(degree_binary_rand,2,function(x) { get_res2(x) })

# get controlled effect size
effect.size <- coef.trait - median(coefs.data)


# controlled degree (after pre-network permutation)
idcov_genet$biDEGREE.cont <- idcov_genet$degree_binary - apply(degree_binary_rand,1,median)
# Calculate effect controlling for OBS
model <- lm(biDEGREE.cont~Sex,data=idcov_genet)
coef.controlled <- coefficients(model)[2]
p.model.controlled <- summary(model)$coefficients[2,4]

#idcov_genet$OBS<-c(18,17,22,22,17,17,23,19,20,19,20,18,23,26,32,24,22,18,17,35,16,20,35,16,20,25,28,18,32,19,18,16,22,16,39,30,27,30,15,23,24,16,23)

## Node permutations
coefs.tmp <- rep(NA,Npermute)
coefs.tmp.control <- rep(NA,Npermute)
coefs.tmp.obs.control <- rep(NA,Npermute)
for (i in 1:Npermute) {
  idcov_genet$Sex.tmp <- sample(idcov_genet$Sex)
  model.tmp <- lm(degree_binary~Sex,data=idcov_genet)
  coefs.tmp[i] <- coefficients(model.tmp)[2]
  model.tmp <- lm(biDEGREE.cont~Sex.tmp,data=idcov_genet)
  coefs.tmp.control[i] <- coefficients(model.tmp)[2]
  model.tmp <- lm(degree_binary~Sex.tmp+OBS,data=idcov_genet)
  coefs.tmp.obs.control[i] <- coefficients(model.tmp)[2]
}

# P values
p.data <- sum(coef.trait<=coefs.data)/Npermute
p.node <- sum(coef.trait<=coefs.tmp)/Npermute
p.node.control <- sum(coef.controlled <= coefs.tmp.control)/Npermute
p.node.obs.control <- sum(coef.trait.obs.control <= coefs.tmp.obs.control)/Npermute
p.t <- sum(t.trait <= t.data)/Npermute



########### GAI ##############
#STRENGTH
library(sna)
netsexGAI<-list_SYM_MATRIX_GAI_real
rownames(netsexGAI)<-idcov_genet$Sex
colnames(netsexGAI)<-idcov_genet$Sex
idcov_genet$strength<-degree(netsexGAI,gmode="graph")
coefGAI<-coef(glm(strength~Sex, data = idcov_genet))[2]

# convert list list_RANDOM_SRI in 3x20000x43x43 array
randomGAI<-array(data = NA, dim = c(Npermute,dim(traits_by_related$related)))

for (j in 1:Npermute){
  randomGAI[j,,]<-as.matrix(list_SYM_MATRIX_GAI_random[[j]])
  colnames(randomGAI[j,,])<-colnames(traits_by_related$related)
  rownames(randomGAI[j,,])<-colnames(traits_by_related$related)
}

strength_rand<-apply(randomGAI,1,function(x) {degree(x,gmode = "graph")})
coefsGAI<-apply(strength_rand,2,function(x) {coef(glm(x~Sex, data = idcov_genet))[2]})

#plot results and save as 4x11 inches (landscape) (same as box 4 in Farine & Whitehead 2015)
par(mfrow=c(1,2),cex.lab=1.5)
#plot network
nodecolor=as.character(idcov_genet$Sex)
nodecolor=gsub("M","red",nodecolor)
nodecolor=gsub("F","blue",nodecolor)
#plot(graph.adjacency(netsex,mode="undirected", diag=FALSE, weighted = TRUE, add.rownames = "code"),vertex.color=nodecolor,vertex.size=attribute$strength*4)
#plot observed difference
plot(idcov_genet$strength~factor(idcov_genet$Sex),col=c("blue","red"),xlab="Sex",ylab="Strength (weighted degree)",ylim=c(0,max(idcov_genet$strength)))
#Plot resulting distribution
a<-hist(coefsGAI,xlim=c(min(coefsGAI),max(coefsGAI)), col="black",main="",xlab="Coefficient value",ylab="Frequency",breaks=100)
abline(v=coefGAI, col="red")
#Model results
summary(glm(strength~Sex, data = idcov_genet))
#P value
sum(abs(coefGAI)<abs(coefsGAI))/20000


#DEGREE (binary)
#Measure binary degree


library(sna)
idcov_genet$degree_binary <- degree(netsexGAI, gmode="graph",ignore.eval=TRUE)
coef2GAI<-coef(glm(degree_binary~Sex, data = idcov_genet))[2]
degree_binary_randGAI <- apply(randomGAI,1,function(x) {degree(x,gmode = "graph",ignore.eval=TRUE)})
detach("package:sna")
coefs2GAI<-apply(degree_binary_randGAI,2,function(x) {coef(glm(x~Sex, data = idcov_genet))[2]})

#plot results and save as 4x11 inches (landscape) (same as box 4 in Farine & Whitehead 2015)
par(mfrow=c(1,2),cex.lab=1.5)
#plot network
nodecolor=as.character(idcov_genet$Sex)
nodecolor=gsub("M","red",nodecolor)
nodecolor=gsub("F","blue",nodecolor)
#plot(graph.adjacency(netsex,mode="undirected", diag=FALSE, weighted = TRUE, add.rownames = "code"),vertex.color=nodecolor,vertex.size=attribute$degree_binary*0.5)
#plot observed difference
plot(idcov_genet$degree_binary~factor(idcov_genet$Sex),col=c("blue","red"),xlab="Sex",ylab="Degree (binary degree)",ylim=c(0,max(idcov_genet$degree_binary)))
#Plot resulting distribution
a<-hist(coefs2GAI,xlim=c(min(coefs2GAI),max(coef2GAI)), col="black",main="",xlab="Coefficient value",ylab="Frequency",breaks=100)
abline(v=coef2GAI, col="red")
#Model results
summary(glm(degree_binary~Sex, data = idcov_genet))
#P value
sum(abs(coef2GAI)<abs(coefs2GAI))/20000




###########################################################
########## Farine's double permutation procedure ##########

rawGBI
head(rawINFO)

head(GBI_related_filter)

kinship.matrix <- related.matrix[rownames(networkDP),colnames(networkDP)]

networks_rand <- network_permutation(GBI_related, data_format = "GBI", permutations = 1000, association_matrix = NULL, within_day = time, within_location = loc)

out <- apply(networks_rand,1,function(x) { model <- mrqap.dsp(x~kinship.matrix,randomisations=1); return(c(model$coefficients[2],model$test.statistic[2])) })


# kinship matrix
kinship.matrix <- related.matrix[rownames(networkDP),colnames(networkDP)]

# Calculate network
#network <- get_network(rawGBI, data_format="GBI", association_index = "SRI", occurrences=occurrences, locations = rawINFO$Site)
networkDP <- get_network(GBI_related_filter[,1:43], data_format="GBI", association_index = "SRI")

# Calculate effects
model <- mrqap.dsp(networkDP~kinship.matrix,randomisations=Npermute)
coef.trait <- model$coefficients[2]
p.node <- model$P.greater[2]

# get t.statistic
model <- mrqap.dsp(networkDP~kinship.matrix,randomisations=1)
t.trait <- model$test.statistic[2]


coefs.data<-rep(NA, Npermute)
t.data<-rep(NA, Npermute)
p.data<-rep(NA, Npermute)
p.t.data<-rep(NA, Npermute)
## Data permutations
# Create random networks
networks_rand <- network_permutation(GBI_related_filter[,1:43], association_matrix=networkDP, data_format = "GBI", permutations = Npermute, locations = GBI_related_filter$Site, within_location = TRUE,days = GBI_related_filter$SP, within_day = TRUE)

# coefficients based on data permutations
out <- apply(networks_rand,1,function(x) { model <- mrqap.dsp(x~kinship.matrix,randomisations=1); return(c(model$coefficients[2],model$test.statistic[2])) })
coefs.data <- out[1,]
t.data <- out[2,]
p.data <- sum(coef.trait<=coefs.data)/Npermute
p.t <- sum(t.trait<=t.data)/Npermute

# get controlled effect size
effect.size <- coef.trait - median(coefs.data)

# control edge values by median of permuations
network_median <- apply(networks_rand,c(2,3),median)
network_controlled <- networkDP - network_median

# Calculate effect from controlled data
model <- mrqap.dsp(network_controlled~kinship.matrix,randomisations=Npermute)
coef.controlled <- model$coefficients[2]
p.node.control <- model$P.greater[2]


# sd
sds <- sd(networkDP)

# mantel
cors <- mantel(networkDP,kinship.matrix)$statistic


