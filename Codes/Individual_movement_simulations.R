euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
boots<- function(x){
  for (i in 1:1000){
  s1<- sample(x,size = 3,replace = T)
  ms<- mean(s1)
  s[i]<- ms
  }
  #m1<- mean(s)
  #k1<- sd(s)
  #r1<- c(m1,k1)
  return(s)
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  beta.params<- c(alpha,beta)
  return(beta.params)
}

estGammaParams<- function(mu,var){
  shp<- mu^2/var^2
  scl<- var^2/mu
  gamma.params<- c(shp,scl)
  return(gamma.params)
}

lscp<- raster("Rasters/landscapes/frag_landscape_0.1_50.tif")
plot(lscp)
df<- as.data.frame(lscp,xy=T)
colnames(df)<- c("x","y","habitat")
df$id<- 1:nrow(df)
start_pops<- read.csv("Rasters/landscapes/start_pops.csv",header = T)
xid<- start_pops$id
df[xid,"habitat"]<- 1


spg<- df
coordinates(spg)<- ~x+y
gridded(spg)<- TRUE
lscp<- raster(spg)
plot(lscp)
points(start_pops[,c("x","y")],pch=19)
est.data<- read.csv("Output/Estimates.csv",header=T)
est<- est.data %>% filter(Habitat %in% c("Specialist"))

### GRT in minutes, sd: 10 min
grt_data<- read.csv("Output/grt.csv",header=T)

## adjust the scale; h to reflect real movement;
## adjust the conc. parameter, r to reflect real movement
tree<- start_pops
sim_list<- list()
mp01<- list();mp02<- list(); mp03<- list()


### bootstrap means and sd
sl_en<- boots(est$sl.encamp)
sl_tvl<- boots(est$sl.travel)
ta_en<- boots(est$ta.encamp)
ta_tvl<- boots(est$ta.travel)
pr_en_hab<- boots(est$habitat.encamp); d3<- data.frame(pr.en=pr_en_hab,hab="frag",grp="spec")
pr_en_mat<- boots(est$habitat.encamp); d4<- data.frame(pr.en=pr_en_mat,hab="matrix",grp="spec")
pr_tvl_hab<- boots(est$habitat.travel);d2<- data.frame(pr.en=pr_tvl_hab,hab="frag",grp="gen")
pr_tvl_mat<- boots(est$matrix.travel);d4<- data.frame(pr.en=pr_tvl_mat,hab="matrix",grp="gen")
sel_hab<- boots(log(est$exp.coef.)); d2<- data.frame(rss=sel_hab,grp="spec")

dx1<- rbind(d1,d2)
obj<- glmer(rss~grp+(1|grp),data=dx1,family="gaussian")
summary(obj)
confint(obj) ## gen matrix travel: -0.05 [-0.06,-0.05]
## step lengths from gamma, conc. from beta, prob from beta and habitat selection from log normal
# gamma distribution parameters for step lengths in encamping/traveling mode
encamp.params<- estGammaParams(sl_en[1],sl_en[2])
travel.params<- estGammaParams(sl_tvl[1],sl_tvl[2])
# concentration parameters in encamping/traveling mode
encamp.conc<- estBetaParams(ta_en[1],ta_en[2])
travel.conc<- estBetaParams(ta_tvl[1],ta_tvl[2])
# state probabilities in habitats
hab.encamp<- estBetaParams(pr_en_hab[1],pr_en_hab[2])
hab.travel<- estBetaParams(pr_tvl_hab[1],pr_tvl_hab[2])
# state probabilities in matrix
mat.encamp<- estBetaParams(pr_en_mat[1],pr_en_mat[2])
mat.travel<- estBetaParams(pr_tvl_mat[1],pr_tvl_mat[2])



for (m in 1:100){
  trajectory_list<- list()
  for(t in 1:nrow(tree)){
    df.h<- data.frame(scl= c(rgamma(1,shape=encamp.params[1],scale = encamp.params[2]),rgamma(1,shape = travel.params[1],scale=travel.params[2]))/10,conc=c(rbeta(1,encamp.conc[1],encamp.conc[2]),rbeta(1,travel.conc[1],travel.conc[2])),pr=c(rbeta(1,hab.encamp[1],hab.encamp[2]),rbeta(1,hab.travel[1],hab.travel[2])))
    df.m<- data.frame(scl= c(rgamma(1,shape=encamp.params[1],scale = encamp.params[2]),rgamma(1,shape = travel.params[1],scale=travel.params[2]))/10,conc=c(rbeta(1,encamp.conc[1],encamp.conc[2]),rbeta(1,travel.conc[1],travel.conc[2])),pr=c(rbeta(1,mat.encamp[1],mat.encamp[2]),rbeta(1,mat.travel[1],mat.travel[2])))
    BodyMass.Value<- runif(1,min=min(est$Mass),max=max(est$Mass))   ## mass: 100 g
    grt<- 4.5* (BodyMass.Value/1000)^0.5 *60
    value.grt<-(rgamma(n = 1,shape = grt^2/(10)^2,scale = (10)^2/grt))
    start_point<- as.data.frame(matrix(NA,nrow=round(value.grt/10),ncol=2))
    colnames(start_point)<- c("x","y")
    start_point[1,]<- tree[t,c("x","y")]
    iter<- round(value.grt/10)
    for (i in 1:(iter-1)){
      source.attr<- cbind(raster::extract(lscp, start_point[i,1:2], start_point = T),start_point[i,1:2])
      colnames(source.attr)<- c("habitat","x","y")
      if (source.attr$habitat==1){ ## if bird in habitat
        c<- df.h[sample(c(1,2),size=1,prob=df.h$pr),] 
        traj1<-simm.crw(date=1:10, h = c$scl, r = c$conc,
                        x0=as.numeric(start_point[i,c("x","y")]), id="A1",
                        typeII=TRUE, proj4string=CRS())
        df<-as.data.frame(traj1[[1]])
        df<- df %>% dplyr::select(x=x,y=y)
        dest.attr<- cbind(raster::extract(lscp, df[10,1:2], df = T),df[10,1:2])
        colnames(dest.attr)<- c("ID","habitat","x","y")
        if (dest.attr$habitat==0){ ## if destination is matrix
          ## probability of residing in forest and not moving to matrix
          s<- sample(c(1,10),size = 1,prob = c(exp(rnorm(1,mean=sel_hab[1],sd=sel_hab[2])),1))# bird picks forest with prob p
          start_point[i+1,]<- df[s,]
        }
        else { ## if destination is habitat
          s<- sample(c(1,10),size = 1,prob = c(0.5,0.5)) # equal chance
          start_point[i+1,]<- df[s,] # moves
        }} else { ## if bird is in matrix
          e<- df.m[sample(c(1,2),size=1,prob=df.m$pr),] 
          traj2<-simm.crw(date=1:10, h = e$scl, r = e$conc,
                          x0=as.numeric(start_point[i,c("x","y")]), id="A1",
                          typeII=TRUE, proj4string=CRS())
          df<-as.data.frame(traj2[[1]])
          df<- df %>% dplyr::select(x=x,y=y)
          dest.attr<- cbind(raster::extract(lscp, df[10,1:2], df = T),df[10,1:2])
          colnames(dest.attr)<- c("ID","habitat","x","y") 
          if (dest.attr$habitat==1){ ## if destination is habitat
            s<- sample(c(1,10),size = 1,prob = c(1,rnorm(1,mean=sel_hab[1],sd=sel_hab[2])))# bird picks forest with prob p
            start_point[i+1,]<- df[s,]
          } 
          else { ## if destination is matrix
            s<- sample(c(1,10),size = 1,prob = c(0.5,0.5))# equal chance
            start_point[i+1,]<- df[s,]
          }
        }
      mp01[[i]]<- df
    } ## for loop for i ends
    mpdf<- do.call("rbind",mp01)
    mpdf$id<- t
    start_point$id<- t
    trajectory_list[[t]]<- start_point
    mp02[[t]]<- mpdf
  }  ## for loop for t ends
  tl<- do.call("rbind",trajectory_list)
  tl$sim<- m
  sim_list[[m]]<- tl
  mpdf01<- do.call("rbind",mp02)
  mpdf01$sim<- m
  mp03[[m]]<- mpdf01
}  
ul<- do.call("rbind",sim_list); mp_df<- do.call("rbind",mp03)
ul$br<- 0.8; mp_df$br<- 0.8
#ul$resident<- 0.2
ul$frag<- 0.9; mp_df$frag<- 0.9
ul$habitat<- 0.75; mp_df$habitat<- 0.75
al<-read.csv("Data_derived/new/ldd_fragments_05.csv",header = T)
gl<- rbind(al,ul)
write.csv(gl,"Data_derived/new/ldd_fragments_05.csv",row.names = F)
mp_al<-read.csv("Data_derived/new/mp_fragments_05.csv",header = T)
mp_gl<- rbind(mp_df,mp_al)
write.csv(mp_gl,"Data_derived/new/mp_fragments_05.csv",row.names = F)


x1<- data.frame(pr=gn_tvl_mat,group="Matrix")
x2<- data.frame(pr=gn_tvl_hab,group="Habitat")
x<- rbind(x1,x2)
kruskal.test(pr~group,data=x)
