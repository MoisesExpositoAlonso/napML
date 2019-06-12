ggpca<-function(pca,
                        dimensions=c(1,2),
                        arrownames=rownames(pca$rotation),
                        arrowscale=2.2){


labelx= paste0('standarized PC1 (', round(summary(pca)[6][[1]][2,1] ,digits=3)*100,'% explained var.)')
labely= paste0('standarized PC2 (', round(summary(pca)[6][[1]][2,2] ,digits=3)*100,'% explained var.)')

p<-ggplot() + xlab(labelx)+ylab(labely)

# add points
eigenvec=data.frame(pca$x)
p<-p+
  geom_point(data=eigenvec,aes(x=scale(PC1),y=scale(PC2)),color='black',size=3.5,shape=1) +
  theme_bw()

# add arrows
rot=data.frame(pca$rotation)
rot$group=row.names(rot)

p<-
p+
  geom_segment(data=rot,aes(x=0,y=0, xend=PC1*arrowscale,yend=PC2*arrowscale ),
               # color=exp2col(sort(rot$group)),
               arrow=arrow(length = unit(1/2, "picas")),lwd=1.5) +
   geom_text(data=rot,aes(x=PC1*arrowscale,y=PC2*arrowscale, label=arrownames),
            # color=exp2col(sort(rot$group)),
            size=3)

return(p)
}


ggpca_fitness<-function(pca,
                        dimensions=c(1,2),
                        dim3=dsumgrand$FloweringRank,
                        dim3col=brewer.pal(9,"Greys"),
                        dim3name='Flowering rank',
                        # arrowcol=experimentcolors(),
                        arrownames=rownames(pca$rotation),
                        arrowscale=2.2){


labelx= paste0('standarized PC1 (', round(summary(pca)[6][[1]][2,1] ,digits=3)*100,'% explained var.)')
labely= paste0('standarized PC2 (', round(summary(pca)[6][[1]][2,2] ,digits=3)*100,'% explained var.)')

p<-ggplot() + xlab(labelx)+ylab(labely)

# add points
eigenvec=data.frame(pca$x)
p<-p+
  geom_point(data=eigenvec,aes(x=scale(PC1),y=scale(PC2)),color='black',size=3.5,shape=1)+
  geom_point(data=eigenvec,aes(x=scale(PC1),y=scale(PC2), group=dim3,color=dim3),size=3)+
  scale_colour_gradientn(dim3name,colours = dim3col) + theme_bw()


# add arrows
rot=data.frame(pca$rotation)
rot$group=row.names(rot)

p<-
p+
  geom_segment(data=rot,aes(x=0,y=0, xend=PC1*arrowscale,yend=PC2*arrowscale ),
               color=exp2col(sort(rot$group)),
               arrow=arrow(length = unit(1/2, "picas")),lwd=1.5) +
   geom_text(data=rot,aes(x=PC1*arrowscale,y=PC2*arrowscale, label=arrownames),
            color=exp2col(sort(rot$group)),size=3)

return(p)
}



ggpca_env<-function(pcaenv,
                        dimensions=c(1,2),
                        dim3=validexp(T) ,
                        dim3col=experimentcolors(),
                        arrowcol='grey',
                        arrownames=rownames(pcaenv$rotation),
                        arrowscale=2.2){


labelx= paste0('standarized PC1 (', round(summary(pcaenv)[6][[1]][2,1] ,digits=3)*100,'% explained var.)')
labely= paste0('standarized PC2 (', round(summary(pcaenv)[6][[1]][2,2] ,digits=3)*100,'% explained var.)')

# add points
eigenvec=data.frame(pcaenv$x)
eigenvec$group=rownames(eigenvec)

p<-ggplot() + xlab(labelx)+ylab(labely)+ theme_bw()+
  geom_point(data=eigenvec,aes(x=scale(PC1),y=scale(PC2),color=group ),size=5) +
  # geom_point(data=eigenvec,aes(x=scale(PC1),y=scale(PC2)),color=exp2col( (eigenvec$group) ),size=5) +
  geom_text(data=eigenvec,aes(x=scale(PC1),y=scale(PC2)),label=eigenvec$group)
  # scale_color_manual('',values= exp2col( sort(eigenvec$group) ) )
p

# add arrows
rot=data.frame(pcaenv$rotation)

p<-p+geom_segment(data=rot,aes(x=0,y=0, xend=PC1*arrowscale,yend=PC2*arrowscale  ),
                  color='darkgrey',arrow=arrow(length = unit(1/2, "picas")),lwd=1.5) +
  geom_text(data=rot,aes(x=PC1*arrowscale,y=PC2*arrowscale, label=arrownames),color='darkgrey',size=5)
p

return(p)
}

