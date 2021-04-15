
PlotRoot2Leave<-function(temp,tips,edges,root,internal){
  edgenames <-apply(temp[,c(1,2)],1,function(X){paste0(X,collapse = '')})
  paths <- lapply(tips,function(tip){
    path<-GetPath(tip,edges,root,tips,internal)
    if(dim(path)[1]!=1){
      path_enames <- apply(path[,c(2,3)],1,function(X){paste0(X,collapse= '')})
    }else{
      path_enames <- paste0(path[1,c(2,3)],collapse= '')
    }
    temp <- temp[edgenames%in%path_enames,]
    plot(x=temp[,3],y=temp[,4],pch=16,cex=0.5,main=tip)    
  })
}

ImpulseEVFpertip <- function(phyla, edges,root,tips,internal, neutral, tip,Sigma,evf_center){
  beta <- runif(1,5,15)
  h1 <- runif(1,2,5)
  h2 <- runif(1,0,h1)
  t1 <- runif(1,0,0.9)
  t2 <- runif(1,t1,1)
  path<-GetPath(tip,edges,root,tips,internal)
  totaltime <- sum(path[,4])
  ImpulseParam <- c(beta,0,h1,h2,t1*totaltime,t2*totaltime)
  par(mfrow=c(1,3))
  plot(phyla,show.tip.label = F,lwd=2)
  tiplabels(cex=2)
  nodelabels(cex=2)
  neutral2 <- SampleSubtree(root,0,evf_center,edges,ncells,Sigma,neutral)
  feature_samples <- AddImpulse(path,edges,neutral2,ImpulseParam,T,tip,tips)
  return(feature_samples)
}


AddImpulse <- function(path,edges,neutral,ImpulseParam,plotting,main,tips){
  edgenames <-apply(neutral[,c(1,2)],1,function(X){paste0(X,collapse = '')})
  enames <- unique(edgenames)
  if(dim(path)[1]!=1){
    path_enames <- apply(path[,c(2,3)],1,function(X){paste0(X,collapse= '')})
  }else{
    path_enames <- paste0(path[1,c(2,3)],collapse= '')
  }
  path_neutral <- neutral[edgenames%in%path_enames,]
  timepoints <- path_neutral[,3]
  path_values <- evalImpulse(ImpulseParam,timepoints)
  path_neutral[,4] <- path_neutral[,4]+path_values
  #the bug here is that other_neutral looses some nodes? (didn't include children that are below other nodes)
  other_neutral <- neutral[!edgenames%in%path_enames,]
  node_impulse <- evalImpulse(ImpulseParam,path[,5])
  other_neutral<-lapply(c(1:length(path[,2])),function(i){
    par <- path[i,2]
    firstedge <- other_neutral[other_neutral[,1]==par,]
    firstedge[,4] <- firstedge[,4]+node_impulse[i]
    other_edges <- edges[!apply(edges[,c(2,3)],1,function(X){paste0(X,collapse = '')})%in%path_enames,]
    children <- getDescendants(other_edges, par,tips,NULL)
    other_neutral[other_neutral[,1]%in%children,4]=other_neutral[other_neutral[,1]%in%children,4]+node_impulse[i]
    return(rbind(other_neutral[other_neutral[,1]%in%children,],firstedge))
  })
  other_neutral <- do.call(rbind,other_neutral)
  newtree<-rbind(path_neutral,other_neutral)
  if(plotting==T){
    PlotPaths(path_neutral,main)
    PlotPaths(newtree,NA)
  }
  return(newtree)
}

#' @author Liam Revell
getDescendants<-function(edges,par,tips,curr=NULL){
  if(is.null(curr)) {curr<-vector()}
  children<-edges[edges[,2]==par,3]
  curr<-c(curr,children)
  w<-children[!children%in%tips]
  if(length(w)>0) for(x in w){
    curr<-getDescendants(edges,x,tips,curr)    
  }
  return(curr)
}

GetRandomTree<-function(ntips){
  phyla <- rtree(ntips)
  plot(phyla,show.tip.label = F,lwd=2)
  tiplabels(bg='yellow',col='black',cex=1)
  nodelabels() 
  return(phyla) 
}

GetPath <- function(tip,edges,root,tips,internal){
  path<-c()
  child<-tip
  par<- 0
  while(par!=root){
    path<- rbind(c(edges[edges[,3]==child,]),path)
    par <- edges[edges[,3]==child,2]
    child <- par
  }
  path <- rbind(edges[edges[,2]==root & edges[,3]==child,],path)
  depth <- c(0,cumsum(path[,4]))
  totaltime <- sum(path[,4])
  path <- cbind(path,depth[-length(depth)])
  colnames(path) <- c('id','par','child','edge_len','par_depth')
  return(path)
}

PlotPaths<-function(temp,main){
  edgenames <-apply(temp[,c(1,2)],1,function(X){paste0(X,collapse = '')})
  enames <- unique(edgenames)
  cols <- rainbow(length(enames))
  plot(x=temp[,3],y=temp[,4],pch=16,cex=0.2,main=main)
  for(i in c(1:length(enames))){
    points(x=temp[edgenames==enames[i],3],y=temp[edgenames==enames[i],4],pch=16,cex=0.5,col=cols[i])    
  }
}

SampleEdge <- function(edge,depth,anc_state,edges,ncells,t_sample=NA){
  if(is.na(t_sample[1])){
    #t_sample <- c(0,sort( runif(round(edge[4]*ncells/sum(edges[,4])),0,edge[4]) ))
    if (ceiling(edge[4]*ncells/sum(edges[,4]))-1 < 0) {stop("the total number of cells is too few.")}
    t_sample <- c(0,seq(0, edge[4], edge[4]/(ceiling(edge[4]*ncells/sum(edges[,4]))-1)))
    t_sample<-c(t_sample,edge[4])
  }else{
    t_sample<-sort(c(0,t_sample-depth))
  }
  t_interval<-diff(t_sample)
  x_change <- sapply(t_interval,function(sig){rnorm(1,0,sqrt(sig))})
  x_sample <- cumsum(x_change)
  result<-cbind(depth+t_sample[-1],x_sample+anc_state)
  result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
  #plot(result,type='l')
  return(result)
}

SampleSubtree <- function(par,depth,anc_state,edges,ncells,neutral=NA){
  children <- edges[edges[,2]==par,3] # get the children of the current node
  result<-lapply(c(1:length(children)),function(j){
    edge<-edges[edges[,2]==par & edges[,3]==children[j],] # given the parent and child, find the edge
    if(sum(edges[,2]==children[j])==0){ # this means the current node is a leaf
      if(is.na(neutral[1])){
        result <- SampleEdge(edge,depth,anc_state,edges,ncells)}else{
          t_sample <- neutral[neutral[,1]==edge[2] & neutral[,2]==edge[3],3]
          result <- SampleEdge(edge,depth,anc_state,edges,ncells,t_sample)
        }
      result <- result[c(1:(length(result[,1]-1))),]
    }else{
      if(is.na(neutral[1])){
        result <- SampleEdge(edge,depth,anc_state,edges,ncells)}else{
          t_sample <- neutral[neutral[,1]==edge[2] & neutral[,2]==edge[3],3]
          result <- SampleEdge(edge,depth,anc_state,edges,ncells,t_sample)
        }
      anc_state <- result[length(result[,1]),4]
      result <- result[c(1:(length(result[,1]-1))),]
      depth <- depth + edge[4]
      result1 <- SampleSubtree(children[j],depth,anc_state,edges,ncells,neutral)
      result <- rbind(result,result1)
    }
    return(result)
  })
  result<-do.call(rbind,result)
  return(result)
}

evalImpulse <- function(vecImpulseParam, vecTimepoints) {
  # beta is vecImpulseParam[1] h0 is vecImpulseParam[2] h1 is
  # vecImpulseParam[3] h2 is vecImpulseParam[4] t1 is vecImpulseParam[5]
  # t2 is vecImpulseParam[6]
  vecImpulseValue <- sapply(vecTimepoints, function(t) {
    (1/vecImpulseParam[3]) *
      (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
         (1/(1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
      (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
         (1/(1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
  })
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}