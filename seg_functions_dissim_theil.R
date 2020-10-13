#Calculate Theil, Dissimilarity, and Gini segregation indices 
#for simulated neighborhoods with no pre-existing administrative divisions (eg census tracts)

#Preamble
library(ggplot2)
library(tidyr)
library(plyr);library(dplyr, warn.conflicts = F)
library(tictoc)
rm(list=ls())
s = function(x){summary(factor(x))}
####################################################################################################

##########function for spatial grid####################################################################
#These segregation indices rely on dividing the space into sub-units
#In the abstence of administrative divisions (eg census tracts) we can simply impose a grid
#This is useful for simulated neighborhoods
grid=function(lat,long,n){
  latbounds=c(min(lat),max(lat))
  longbounds=c(min(long),max(long))
  latgrid=seq(from=latbounds[1]-0.0005,to=latbounds[2]+0.0005,length.out=n)
  longgrid=seq(from=longbounds[1]-0.0005,to=longbounds[2]+0.0005,length.out=n)
  data.frame(lat=latgrid,long=longgrid)
}
####################################################################################################

####################################################################################################
#Function to count people in each grid space
count_ppl = function(data,grid,var,n){
  j=1
  count_matrix=matrix(nrow=(n-1)^2, ncol=length(levels(as.factor(data[,var]))))
  count_matrix=data.frame(count_matrix)
  names(count_matrix)=levels(as.factor(data[,var]))
  for(i in 1:(n-1)){ #loop over lat
    lat.min=grid$lat[i]
    lat.max=grid$lat[i+1]
    for(l in 1:(n-1)){ #loop over long
      long.min=grid$long[l]
      long.max=grid$long[l+1]
      inside=which(data$lat>lat.min&data$lat<lat.max&
                     data$long>long.min&data$long<long.max)
      C=data[inside,]
      
      for(k in 1:length(levels(as.factor(data[,var])))){ #loop over factor levels
        count_matrix[j,k]=sum(C[,var] == levels(as.factor(data[,var]))[k],na.rm=T) #na.rm=T added on 11-2; this isn't in earlier versions
      }
      j=j+1
    }
  }
  return(count_matrix)
}
####################################################################################################

####################################################################################################
#Function for Dissimilarity Index
dissimilarity_v2 = function(count_matrix){ #count_matrix is output from count_ppl
  #Choose majority, minority
  if(length(names(count_matrix)[which(apply(count_matrix,2,sum) > 0)]) == 1){ dis_index= NA #if only one group, code as NA
  } else {
    min=names(sort(apply(count_matrix,2,sum),decreasing=T))[2] #define minority group
    To=sum(count_matrix) #total population
    P=sum(count_matrix[,min])/To#total minority PROPORTION
    D=rep(0,nrow(count_matrix))
    for(i in 1:nrow(count_matrix)){ #loop over gridspaces
      ti=sum(count_matrix[i,]) #total pop in gridspace
      pi=count_matrix[i,min]/(max(ti,0)) #minority PROPORTION in gridspace
      D[i]= ti*abs(pi-P)/(2*To*P*(1-P))
    }
    D[is.na(D)]=0 #set na's to 0 (na's caused by 0/0, empty cell)
    dis_index=sum(D)
  }
}
####################################################################################################

####################################################################################################
#Function for Gini Index
gini=function(count_matrix,data,var){ #count_matrix is output from count_ppl
  if(nlevels(factor(data[,var])) == 1){ gini_index= NA #if only one group
  } else {
    maj=names(summary(factor(data[,var]))[order(summary(factor(data[,var])),decreasing=T)][1]) #define majority group
    min=names(summary(factor(data[,var]))[order(summary(factor(data[,var])),decreasing=T)][2]) #define minority group
    T=sum(count_matrix) #total population
    P=sum(count_matrix[,min])/T#total minority PROPORTION
    D=rep(0,nrow(count_matrix))
    t=rep(0,nrow(count_matrix))
    p=rep(0,nrow(count_matrix))
    for(i in 1:nrow(count_matrix)){ #loop over gridspaces
      t[i]=sum(count_matrix[i,]) #total pop in gridspace
      p[i]=count_matrix[i,min]/(max(t[i],0)) #minority PROPORTION in gridspace
    }
    p[which(is.na(p))]=0
    
    g=vector()
    for(i in 1:nrow(count_matrix)){
      gi=vector()
      for(j in 1:nrow(count_matrix)){
        gi[j]=t[i]*t[j]*abs(p[i]-p[j])/(2*T^2*P*(1-P))
      }
      g[i]=sum(gi)
    }
    gini_index=sum(g)
  }
}
####################################################################################################

####################################################################################################
#Function for Theil Index
theil = function(count_matrix, data, var) { #count_matrix is output from count_ppl
  if (nlevels(factor(data[,var])) == 1) {
    theil_index = NA #if only one group
  } else {
    #Define majority, minority
    maj = names(summary(factor(data[,var]))[order(summary(factor(data[,var])), decreasing =
                                                    T)][1])
    minn = names(summary(factor(data[,var]))[order(summary(factor(data[,var])), decreasing =
                                                     T)][2])
    Tot= sum(count_matrix) #total population
    P = sum(count_matrix[, minn]) / Tot#total minority PROPORTION
    #calculate entropy for neighborhood (for Theil index)
    E = P * log(1 / P) + (1 - P) * log(1 / (1 - P))
    e = rep(0, nrow(count_matrix))
    t = vector()
    for (i in 1:nrow(count_matrix)) {
      #loop over gridspaces
      t[i] = sum(count_matrix[i,]) #total pop in gridspace
      pi = count_matrix[i, minn] / (max(t[i], 0)) #minority PROPORTION in gridspace
      e[i] = pi * log(1 / pi) + (1 - pi) * log(1 / (1 - pi))
    }
    e[is.na(e)] = 0 #set na's to 0 (na's caused by 0/0, empty cell)
    theil_index = sum(t * (E - e) / (E * Tot))
  }
}
####################################################################################################

####################################################################################################
######Theil index for multiple groups 
theil_multi = function(count_matrix) { #count_matrix is output from count_ppl
  P = colSums(count_matrix) / sum(count_matrix) #"neighborhood" proportion of each group
  p = count_matrix / rowSums(count_matrix)
  p[is.na(p)] = 0
  p2 = p * log (1 / p)
  p2[is.na(p2)] = 0
  e = rowSums(p2) #entropy for each "neighborhood"
  E = sum(P * log (1 / P), na.rm = T)#entropy for "city"
  tm = sum(rowSums(count_matrix) * (E - e) / (E * sum(count_matrix)))
}
####################################################################################################

####################################################################################################
###########main function to calculate seg indices############################################
seg_calc = function (B, varname , gridsize = 11 ) {
  #invoke grid function
  B_grid=grid(B$lat,B$long,gridsize)
  B_grid=data.frame('lat'=B_grid[[1]],'long'=B_grid[[2]])
  #invoke count function
  count = count_ppl(B,B_grid,varname,gridsize)
  
  #calculate segregation indices
  #dis = dissimilarity(count,B,varname)
  dis = dissimilarity_v2(count)
  gi=gini(count,B,varname)
  th=theil(count,B,varname)
  th_m=theil_multi(count)
  #output list
  return(list("Dissimilarity"=dis,"Gini"=gi,"Theil"=th,"Theil_Multi"=th_m))
}
####################################################################################################

####################################################################################################
#function to break out results from list
br=function(out_religion){
  dis_relig=sapply(out_religion, "[[", 1)
  gini_relig=sapply(out_religion, "[[", 2)
  theil_relig=sapply(out_religion, "[[", 3)
  tm_relig=sapply(out_religion, "[[", 4)
  rel_seg=data.frame('Dissimilarity' = dis_relig, 'Gini'=gini_relig, 'Theil' = theil_relig, 'Theil_Multi' = tm_relig, 'Neigh' = factor(names(dis_relig)))
}
####################################################################################################