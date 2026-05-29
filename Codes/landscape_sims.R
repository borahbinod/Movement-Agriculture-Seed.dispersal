l<-NLMR::nlm_randomcluster(ncol = 100,nrow = 100,resolution =50 ,p=0.5,ai=c(0.5,0.25),rescale = F)
raster::plot(l)

d<- as.data.frame(l,xy=T)
hist(d$clumps)
d[d$clumps==1,3]<-0
d[d$clumps==2,3]<-1
nrow(d)
length(which(d$clumps==1))/nrow(d)
hist(d$clumps)
lscp<- rasterFromXYZ(d)
plot(lscp)
writeRaster(x = lscp,filename = "Rasters/landscapes/frag_landscape_35.tif",overwrite=T)

### select start points
set.seed(100)
start_pops<- data.frame(matrix(NA,nrow=100,ncol=4))
start_pops$X1<-1000; start_pops$X2<- 1000; start_pops$X3<- 1; start_pops$X4<-1

lscp<- raster("Rasters/landscapes/frag_landscape_0.1_25.tif")
plot(lscp)
df<- as.data.frame(lscp,xy=T)
start_pops<- read.csv("Data_derived/arrival/start_pops.csv",header=T)


ds<- distanceFromPoints(lscp, c(6500,6500))
ds_s<- as.data.frame(ds,xy=T)
ds_s<- ds_s[ds_s$layer>=0 & ds_s$layer<= 4500,]
start_pops<- ds_s[sample(1:nrow(ds_s),size = 25,replace = F),]
points(start_pops[,c("x","y")],pch=19)
cell<- which(df$frag_landscape_0.1_50==1)
df$presence<- 0
#df[x==167.5 && y== 497.5,"presence"]<- 1
#start_pops[,]<- df[df$presence==1,]


disp_plant<- read.csv("Data_derived/arrival/Peru1_10percent.csv",header = T)
plant_sp<- unique(disp_plant$plant)
plant_list<- list()

## filter by plant species
for (i in 1:length(plant_sp)){
plant_subset<- disp_plant %>% filter(plant %in% c(plant_sp[i]))
## each bird has 100 interactions, one for each tree
bird_sp<- unique(plant_subset$bird)
bird_list<- list()
if (is.na(bird_sp)== FALSE){
for (j in 1:length(bird_sp)){
bird_subset<- plant_subset %>% filter(bird %in% c(bird_sp[j]))
bird_name<- unique(bird_subset$bird)
bird.filtered<- data.trait %>% filter(Scientific %in% bird_name)%>% select("Scientific","English","BodyMass.Value")
seeds_rm<- 0.808+0.311*bird.filtered$BodyMass.Value--0.246*plant_subset$size
cells_colonize<- data.frame(matrix(NA,nrow=100,ncol=1))
for (b in 1:nrow(bird_subset)){
  d1<- distanceFromPoints(lscp, start_pops[b,1:2])
  d<- as.data.frame(d1)
  df_new<- cbind(df,d)
  disp_cells<- df_new[df_new$layer>= bird_subset$X1[b] & df_new$layer<=bird_subset$X2[b],]
  if (nrow(disp_cells)>0){
    cells_arrived<- (rpois(1,lambda = seeds_rm))*0.01
    prop<-ifelse(cells_arrived/length(which(disp_cells$frag_landscape==1))!=Inf,cells_arrived/length(which(disp_cells$frag_landscape==1)),0)
  } else {
    prop<- 0
  }
  cells_colonize[b,1]<- prop
  cells_colonize$sim[b]<- b
  cells_colonize$bird[b]<- bird_name
  }
colnames(cells_colonize)<- c("proportion","sim","bird")
bird_list[[j]]<- cells_colonize
      }
    }
else {
 bird_list[[1]]<- data.frame(proportion=rep(0,100), sim=c(1:100),bird="extinct") 
}
list_birds<- do.call("rbind",bird_list)
list_birds$plant<- as.character(plant_sp[i])
plant_list[[i]]<- list_birds
}

plant_list[[4]]
p<- do.call("rbind",plant_list)
p$habitat<- 98
p$defaunation<- "10 percent"
p<- subset(p,sim <=100,)
e<- read.csv("Data_derived/arrival/Peru_arrival_10percent.csv",header = T)
p<- rbind(p,e)
write.csv(p, "Data_derived/arrival/Peru_arrival_10percent.csv")














































### other code
tree_sp<- read.csv("Data_derived/Tree_species.csv",header = T)
tree_species<- list()
for (i in 1: nrow(tree_sp)){
  sp<- as.character(tree_sp$Tree_species[i])
  type<- as.character(tree_sp$Fruit_type[i])
for (k in 1:100){
#init<- rasterFromXYZ(df_init)
d1<- distanceFromPoints(lscp, start_pops[k,1:2])
d<- as.data.frame(d1)
df_new<- cbind(df,d)
disp_cells<- df_new[df_new$layer>0 & df_new$layer<tree_sp$LDD_25[i],]
if (nrow(disp_cells)>0){
  prop<-length(which(disp_cells$frag_landscape==1))/length(cell)
} else {
  prop<- 0
  }
cells_colonize[k,1]<- prop
cells_colonize$sim[k]<- k
cells_colonize$name<- sp
cells_colonize$type<- type
  }
tree_species[[i]]<- cells_colonize
}

seed_arrival<- do.call("rbind",tree_species); 
write.csv(seed_arrival,"Data_derived/arrival_if_25percent.csv")


















### Dispersal kernels
f<- function(a,b,x){
  ((b/((2*pi*a^2)*gamma(2/b)))*exp(-x^b/a^b))
}
ncell<- 1:round(585/50)
ndist<- ncell*10
p<-f(10^-8,0.12,ndist)

### Generating dispersal events
disp_events<- MigClim::MigClim.migrate(iniDist = paste0("Rasters/Colombia1/init_",sim,".tif"),hsMap = "Rasters/frag_landscape",
                                  envChgSteps = 1, dispSteps = 99,dispKernel = p,
                                  testMode = F,overWrite = T)




x<- raster(xmn=0,xmx=1000,ymn=0,ymx=1000,res=100)
values(x)<-0
x[sample(1:ncell(x),size=50)]<- 1
plot(x)
df<- as.data.frame(x,xy=T)
df$id<- 1:nrow(df)

r<-distanceFromPoints(ra,c(2570,2570))
rr<-as.data.frame(r,xy=T)
rr$id<- 1:nrow(rr)
hist(rr$layer)
rrr<- subset(rr,rr$layer>400 & rr$layer<= 800,)
rrrr<- rrr[sample(1:nrow(rrr),size = 25,),]
rrrrr<-rrr[sample(1:nrow(rrr),size = 25,),]
rrrr<- rbind(rrrr,rrrrr)
points(rrrr[,c("x","y")])
mean(rrrr$layer)
rrrr$fragment<-1
hist(rrrr$layer)
median(rrrr$layer)


df$fragment<-rrrr$fragment[match(df$id,rrrr$id)]
df$fragment[which(is.na(df$fragment)==T)]<-0
ras<- rasterFromXYZ(df)
plot(ras$fragment)

#' Create neutral landscape maps
#' 
#' Use standard methods to generate fractal maps. Binary and continuous surfaces may be produced.
#' 
#' @param k integer. The extent of the map (2^k+1)^2 pixels
#' @param h numeric. Level of aggregation in the map.
#' @param p numeric (0,1). The proportion of map in habitat=1
#' @param binary logical. If TRUE, a 0/1 categorical landscape is produced.
#' @author Shannon Pittman, James Forester, modified by Lauren White
#' @export
#' @example examples/neutral.landscape_example.R
fracland_mod <- function(k, h, p, binary = TRUE) {
  ## Function for creating neutral landscapes Shannon Pittman University of Minnesota May, 2013 k = the extent of the map (2^k+1)^2 pixels h =
  ## how clumped the map should be (ranging from ?? to ??) -- weird behavior at higher values p = proportion of map in habitat 1 binary =
  ## plotflag == if TRUE will plot a filled contour version of the matrix
  
  ## function call: testmap=land(6,1,.5,FALSE,TRUE)
  A <- 2^k + 1  # Scalar-determines length of landscape matrix
  
  #Right now, as written (1-p) represents the amount of habitat listed as "1"
  B <- matrix(0, A, A)  # Creates landscape matrix
  
  B[1, 1] <- 0
  B[1, A] <- 0
  B[A, 1] <- 0
  B[A, A] <- 0
  
  
  iter <- 1
  for (iter in 1:k) {
    scalef <- (0.5 + (1 - h)/2)^(iter)
    
    d <- 2^(k - iter)
    
    # ALL SQUARE STEPS#
    for (i in seq(d + 1, A - d, 2 * d)) {
      for (j in seq(d + 1, A - d, 2 * d)) {
        B[i, j] <- mean(c(B[i - d, j - d], B[i - d, j + d], B[i + d, j - d], B[i + d, j + d])) + scalef * rnorm(n = 1)
      }
    }
    
    # OUTSIDE DIAMOND STEP#
    for (j in seq(d + 1, A - d, 2 * d)) {
      B[1, j] <- mean(c(B[1, j - d], B[1, j + d], B[1 + d, j])) + scalef * rnorm(n = 1)
      B[A, j] <- mean(c(B[A, j - d], B[A, j + d], B[A - d, j])) + scalef * rnorm(n = 1)
    }
    
    for (i in seq(d + 1, A - d, 2 * d)) {
      B[i, 1] <- mean(c(B[i - d, 1], B[i + d, 1], B[i, 1 + d])) + scalef * rnorm(n = 1)
      B[i, A] <- mean(c(B[i - d, A], B[i + d, A], B[i, A - d])) + scalef * rnorm(n = 1)
    }
    
    # INSIDE DIAMOND STEP#
    if (2 * d + 1 <= A - 2 * d) {
      for (i in seq(d + 1, A - d, 2 * d)) {
        for (j in seq(2 * d + 1, A - 2 * d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
      
      for (i in seq(2 * d + 1, A - 2 * d, 2 * d)) {
        for (j in seq(d + 1, A - d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
        }
      }
    }
    
    iter <- iter + 1
  }
  
  if (binary == T) {
    R <- sort(B)
    PosR <- (1 - p) * length(R)  #larger values become habitat, designated as 1
    pval <- R[PosR]
    T1 <- which(B > pval)
    T2 <- which(B <= pval)
    B[T1] <- 1  #habitat is 1
    B[T2] <- 0
  } 
  return(B)
}
set.seed(100)
r<- fracland_mod(k = 7,h = 0.001,p = 0.95,binary = T)
ra<- raster(ncol=129,nrow=129,xmn=0,ymn=0,xmx=129*120,ymx=129*120,vals=r)
plot(ra,col=c("light yellow","limegreen"))
ra
writeRaster(ra,"Rasters/landscapes/frag_landscape_0.9_75.tif")
ra_df<- as.data.frame(ra,xy=T)
length(which(ra_df$layer==1))
