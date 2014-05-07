# plot graphs 4 at a time, 16 to a page
require(ggplot2)

multiplot <- function(plotlist) {
  require(grid)
  numPlots = length(plots)
  
    # Set up the page

    countx <- 1;
    county<-1;
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      if((i-1) %% 16 == 0) {
        countx<-1; county<-1;
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(4,4)))
      }
      
      print(plots[[i]], vp = viewport(layout.pos.row = countx,
                                      layout.pos.col = county))
      countx <- countx+1;
      if(countx>4) {
        county <- county + 1;
        countx<-1;
      }
    }
}

plots <- vector('list', length(80))
count<-1;
for( s in 1:4 ){    
  for(y in 1:d$nYears){  
    for( r in 1:(d$nRivers+1) ){
        plots[[count]] <- ggplot(data=as.data.frame(
                    qqplot(stanout$pBeta[s,y,r,], 
                           c(out$pBeta[s,y,r,,1],out$pBeta[s,y,r,,2],out$pBeta[s,y,r,,3]),plot=FALSE)), 
                    mapping=aes(x=x,y=y)) +
        geom_point() +
        geom_smooth(method="lm", se=FALSE) +
        geom_abline(slope = 1, intercept = 0) +
        xlab(paste("stan pBeta[",s,',',y,',',r,']')) +
        ylab(paste("jags pBeta[",s,',',y,',',r,']')) +
        ggtitle(paste("q-q plot of pBeta[",s,',',y,',',r,']'));
        count <- count + 1;
    }
  }
}
multiplot(plotlist=plots);

plots2 <- vector('list', length(80))
count<-1;
for( s in 1:4 ){    
  for(y in 1:d$nYears){  
    for( r in 1:(d$nRivers+1) ){
      plots2[[count]] <- ggplot(data=as.data.frame(
        qqplot(stanout$phiBeta[s,y,r,], 
               c(out$phiBeta[s,y,r,,1],out$phiBeta[s,y,r,,2],out$phiBeta[s,y,r,,3]),plot=FALSE)), 
        mapping=aes(x=x,y=y)) +
        geom_point() +
        geom_smooth(method="lm", se=FALSE) +
        geom_abline(slope = 1, intercept = 0) + 
        xlab(paste("stan phiBeta[",s,',',y,',',r,']')) +
        ylab(paste("jags phiBeta[",s,',',y,',',r,']')) +
        ggtitle(paste("q-q plot of phiBeta[",s,',',y,',',r,']'));
      count <- count + 1;
      
    }
  }
}
multiplot(plotlist=plots2);

plots3 <- vector('list', length(64))
count<-1;
for( s in 1:4 ){    
  for(r1 in 1:d$nRivers){  
    for( r2 in 1:(d$nRivers) ){
      plots3[[count]] <- ggplot(data=as.data.frame(
        qqplot(stanout$psiBeta[s,r1,r2,], 
               c(out$psiBeta[s,r1,r2,,1],out$psiBeta[s,r1,r2,,2],out$psiBeta[s,r1,r2,,3]),plot=FALSE)), 
        mapping=aes(x=x,y=y)) +
        geom_point() +
        geom_smooth(method="lm", se=FALSE) +
        geom_abline(slope = 1, intercept = 0) +
        xlab(paste("stan psiBeta[",s,',',r1,',',r2,']')) +
        ylab(paste("jags psiBeta[",s,',',r1,',',r2,']')) +
        ggtitle(paste("q-q plot of psiBeta[",s,',',r1,',',r2,']'));
      count <- count + 1;
      
    }
  }
}
multiplot(plotlist=plots3);
