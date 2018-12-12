
#set this to path where files etc reside
setwd("D:/backup/PAPERS_COOKING/SCHWACKE/__DDD")

#read in old MapMan annotation
bins<-read.table(file="Ath_AGI_LOCUS_TAIR10_Aug2012.txt",sep="\t",header=T)
bins$IDENTIFIER<-sub("\\..+$","",bins$IDENTIFIER)

#read in new MapMan X4 annotation
binsX4<-read.table(file="ArabidopsisTha.results.txt",sep="\t",header=T)
binsX4$IDENTIFIER<-sub("\\..+$","",binsX4$IDENTIFIER)

#read in Expression matrix from Gene Cat where ambigous probesets were removed
xps<-read.table(file="ExpMatAra.EXP.sanitized.txt",sep="\t",header=F,row.names=1)
#convert to matrix
xps<-as.matrix(xps)

cnames<-rownames(xps)

rownum<-dim(xps)[1]

#remove unknown bins
bins$BINCODE<-as.vector(bins$BINCODE)


#only major bin
bins$BINCODE<-substr(bins$BINCODE,1,2)
bins$BINCODE<-sub("\\.","",bins$BINCODE)
bins=bins[bins$BINCODE!="26",]
bins=bins[bins$BINCODE!="35",]
bins=bins[bins$IDENTIFIER!="",]
bins$BINCODE<-as.vector(bins$BINCODE)

#remove unknown bins and bin 50
binsX4$BINCODE<-substr(binsX4$BINCODE,1,2)
binsX4$BINCODE<-sub("\\.","",binsX4$BINCODE)
binsX4=binsX4[binsX4$BINCODE!="50",]
binsX4=binsX4[binsX4$BINCODE!="35",]
binsX4=binsX4[binsX4$IDENTIFIER!="",]

binsX4$BINCODE<-as.vector(binsX4$BINCODE)


#we do have empty ""  identifiers like this
#test it


for (cormethod in c("pe","spe")){
	#generate the correalation matrix
	xx<-cor(t(xps),method=cormethod)
	
	#generate empty result matrix
	res<-matrix(nrow=4,ncol=9)
	for (TT in 0:8) {
		THRES<-0.7+TT/30 #set correlation threshold
		a=0;
		b=0;
		c=0;
		d=0;

		for (i in 1:(rownum-1)){
			for (j in (i+1):rownum){
				if (xx[i,j]>THRES){
					k= bins$BINCODE[bins$IDENTIFIER == rownames(xx)[i]]
					l= bins$BINCODE[bins$IDENTIFIER == colnames(xx)[j]]
					if ((length(k)>0) && (length (l)>0)){ #do we have a "meaningful" bincode in both cases
							a=a+1 #increase case number
							if (length(intersect(k,l))>0 ){ #is there an intersection between the two annotations?
								b=b+1 
							}
					}
					
					m= binsX4$BINCODE[binsX4$IDENTIFIER == rownames(xx)[i]]
					n= binsX4$BINCODE[binsX4$IDENTIFIER == colnames(xx)[j]]
					if ((length(m)>0) && (length (n)>0)){ #do we have a "meaningful" X4 bincode in both cases
							c=c+1  #increase case number
							if (length(intersect(m,n))>0 ){ #is there an intersection between the two X4 annotations?
								d=d+1 
							}
					}
					
				}
			}
		}
		res[1,TT+1]<-a
		res[2,TT+1]<-b
		res[3,TT+1]<-c
		res[4,TT+1]<-d
		print(".")
	}
	if (cormethod=="pe"){
		resP<-res  #Pearson
	}
}

plot(0.7+0:8/30, resP[2,]/resP[1,],type="l",col="red",main="Shared Major Bin",xlab="Correlation Threshold", ylab="Concordant Bin Pairs/Informative Bin Pairs")
lines(0.7+0:8/30, res[2,]/res[1,],col="red",lty=2)
lines(0.7+0:8/30, resP[4,]/resP[3,],col="blue")
lines(0.7+0:8/30, res[4,]/res[3,],col="blue",lty=2)



legend(0.7, 1.0, legend=c("MapMan 3 Pearson", "MapMan 4 Pearson","MapMan 3 Spearman", "MapMan 4 Spearman"),
       col=c("red", "blue","red","blue"), lty=c(1,1,2,2), cex=0.8)
