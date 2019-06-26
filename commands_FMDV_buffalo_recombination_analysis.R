# List of commands used to generate files, analyse data etc. 
# for the analysis of sequence data from artificial FMDV infection 
# experiments in African buffaloes.

# See papers:
# Maree et al, J. Virol. 2016
# Ferretti et al, Viruses 2018
# Cortey et al, J. Virol. 2019
# Ferretti et al, bioRxiv, https://doi.org/10.1101/271239

# Contains bash preprocessing commands as well.


######################################
# Analysis of NGS sequences from inoculum

## VP1 only
for POS1 in 1667 1670 1694 1724 1737 1760 1763 1805 1814 1831 1865 1980 2019 2037 2048 2069 2090 2129 2186 2221 2273 2302; do for POS2 in 1667 1670 1694 1724 1737 1760 1763 1805 1814 1831 1865 1980 2019 2037 2048 2069 2090 2129 2186 2221 2273 2302; do samtools view 88_S88_L001.vs.Inoculum.bam Inoculum:1640-2310 | gawk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{n=0; if($1~"^@"){print $0} else {n=split($6,a,"[A-Z]+"); if((n==2)&&($6~"^[0-9]+M$")){ mylen=a[1]; if (($4<'$POS1')&&($4+mylen-1>'$POS2')) {p1='$POS1'-$4+1; p2='$POS2'-$4+1; if(ord[substr($11,p1,1)]-33>20 && ord[substr($11,p2,1)]-33>20) print p1,p2,substr($10,p1,1),substr($10,p2,1) }}}}' | awk '{if (NF==4) print $3""$4}' | sort | uniq -c | sort -k 1nr | awk '{if(NR==1){n1=$1;a1=$2};if(NR==2){n2=$1;a2=$2; r1=substr(a1,1,1)""substr(a2,2,1); r2=substr(a2,1,1)""substr(a1,2,1);}; if(NR>2){if($2==r1){nr1=$1; ar1=$2} else { if($2==r2){nr2=$1; ar2=$2} else {no=no+$1}}}}END{print "'$POS1'","'$POS2'",a1,a2,"r"ar1,"r"ar2,n1,n2,nr1+0,nr2+0,no+0}'; done; done | awk '$1<$2 && NF==11' > table_LD
cat table_LD | awk 'BEGIN{p="0"}{if($1!=p){print $0; p=$1}}' |awk '{d=sqrt(($1-$2)^2);s=$7+$8+$9+$10; pa=($7+$9)/s; pb=($7+$10)/s; ld=$7/s-pa*pb;ldmax=pa*(1-pb); if(pb*(1-pa)<ldmax){ldmax=pb*(1-pa)}; if((pa-0.5)*(pb-0.5)<0){lds=-ld} else {lds=ld}; print $1,$2,d,lds,lds/d,-log(ld/ldmax)/d}' | awk '{for(i=$1;i<$2;i++){print i,$6}}'> table_recrate

d<-read.table("table_recrate",header=F)
plot(d,type="l",lwd=2,ylim=c(0,0.01),ylab="Effective recombination rate per base",xlab="Position along capside sequence")

## whole capsid
cat 88_S88_L001.vs.Inoculum.counts | grep Inoculum | awk '$15>0.95' | awk '{if(($4+$7>2000)&&(($7/($4+$7))^2>0.01^2)&&($15>0.95)) print $2,$4,$7,$7/($7+$4),$3,$6}' > 88_S88_L001.vs.Inoculum.haplotypeSNPS

d<-read.table("88_S88_L001.vs.Inoculum.haplotypeSNPS",header=F)
hist(d[,4],breaks=c(0:10)/20,ylab="Number of mutations",xlab="Minor frequency in the sample",main="Frequency distribution of SNPs in inoculum",col=c(rep("green",4),rep("grey40",6)))
mycols<-c("green","grey40")
plot(d[,c(1,4)],ylim=c(0,1),ylab="Frequency of the minor allele",xlab="Position along capsid sequence",pch=19,col=mycols[1+(d[,4]>0.2)])
lines(lowess(d[,c(1,4)],f=1/2),col="grey40",lwd=2)
dev.copy2pdf("fig_inoculum.pdf")


#PREVIOUSLY: for POS1 in `cat 88_S88_L001.vs.Inoculum.haplotypeSNPS | awk '{if($4>0.25) {print $1}}'`; do for POS2 in `cat 88_S88_L001.vs.Inoculum.haplotypeSNPS | awk '{if($4>0.25) {print $1}}'`; do if [ $(($POS2-$POS1)) -gt 0 ] && [ $(($POS2-$POS1)) -lt 400 ] ; then echo $POS1"-"$POS2 ; samtools view 88_S88_L001.vs.Inoculum.bam "Inoculum:"$(($POS1-1))"-"$(($POS2+1)) | gawk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{n=0; if($1~"^@"){print $0} else {n=split($6,a,"[A-Z]+"); if((n==2)&&($6~"^[0-9]+M$")){ mylen=a[1]; if (($4<'$POS1')&&($4+mylen-1>'$POS2')) {p1='$POS1'-$4+1; p2='$POS2'-$4+1; if(ord[substr($11,p1,1)]-33>20 && ord[substr($11,p2,1)]-33>20) print p1,p2,substr($10,p1,1),substr($10,p2,1) }}}}' | awk '{if (NF==4) print $3""$4}' | sort | uniq -c | sort -k 1nr | awk '{if(NR==1){n1=$1;a1=$2};if(NR==2){n2=$1;a2=$2; r1=substr(a1,1,1)""substr(a2,2,1); r2=substr(a2,1,1)""substr(a1,2,1);}; if(NR>2){if($2==r1){nr1=$1; ar1=$2} else { if($2==r2){nr2=$1; ar2=$2} else {no=no+$1}}}}END{print "'$POS1'","'$POS2'",a1,a2,"r"ar1,"r"ar2,n1,n2,nr1+0,nr2+0,no+0}'; fi ; done; done > tmp ; cat tmp | awk '$1<$2 && NF==11' > table_LD_capsid ; rm tmp
for POS1 in `cat 88_S88_L001.vs.Inoculum.haplotypeSNPS | awk '{if($4>0.25) {print $1}}'`; do for POS2 in `cat 88_S88_L001.vs.Inoculum.haplotypeSNPS | awk '{if($4>0.25) {print $1}}'`; do if [ $(($POS2-$POS1)) -gt 0 ] && [ $(($POS2-$POS1)) -lt 400 ] ; then echo $POS1"-"$POS2  > /dev/stderr ; samtools view 88_S88_L001.vs.Inoculum.bam "Inoculum:"$(($POS1-1))"-"$(($POS2+1)) | gawk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{n=0; if($1~"^@"){print $0} else {n=split($6,a,"[A-Z]+"); if((n==2)&&($6~"^[0-9]+M$")){ mylen=a[1]; if (($4<'$POS1')&&($4+mylen-1>'$POS2')) {p1='$POS1'-$4+1; p2='$POS2'-$4+1; if(ord[substr($11,p1,1)]-33>20 && ord[substr($11,p2,1)]-33>20) print p1,p2,substr($10,p1,1),substr($10,p2,1) }}}}' | awk '{if (NF==4) print $3""$4}' | sort | uniq -c | sort -k 1nr | awk 'BEGIN{while(getline myline < "88_S88_L001.vs.Inoculum.haplotypeSNPS") {split(myline,line," "); if (line[1]=='$POS1'){al1=line[5]; al1m=line[6]}; if (line[1]=='$POS2'){al2=line[5]; al2m=line[6]}}; print al1,al1m,al2,al2m; a1=al1""al2;a2=al1m""al2m;ar1=al1""al2m;ar2=al1m""al2}{if($2==a1){n1=$1} else { if($2==a2){n2=$1} else { if($2==ar1){nr1=$1} else { if($2==ar2){nr2=$1} else {no=no+$1}}}}}END{print "'$POS1'","'$POS2'",a1,a2,"r"ar1,"r"ar2,n1,n2,nr1+0,nr2+0,no+0}'; fi ; done; done > tmp ; cat tmp | awk '$1<$2 && NF==11' > table_LD_capsid ; rm tmp
cat table_LD_capsid | awk 'BEGIN{p="0"}{if($1!=p){print $0; p=$1}}' |awk '{d=sqrt(($1-$2)^2);s=$7+$8+$9+$10; pa=($7+$9)/s; pb=($7+$10)/s; ld=$7/s-pa*pb; if(ld>0){ldmax=pa*(1-pb); if(pb*(1-pa)<ldmax){ldmax=pb*(1-pa)}; print $1,$2,d,ld,ld/d,-log(ld/ldmax)/d} else {ldmax=pa*pb; if((1-pb)*(1-pa)<ldmax){ldmax=(1-pb)*(1-pa)}; print $1,$2,d,ld,ld/d,-log(-ld/ldmax)/d}}' | awk '{for(i=$1;i<$2;i++){print i,$6}}'> table_recrate_capsid
#if((pa-0.5)*(pb-0.5)<0){lds=-ld} else {lds=ld};


cat SAT1-Marti.vs.Inoculum.clustalo | sed '1,3d' | awk '{if (((NR-1) % 4)==0) {s=s""substr($0,22)}; if (((NR-1) % 4)==1) {s2=s2""substr($0,22)}}END{match(s2,/^-+/); s2=substr(s2,RLENGTH+1); s=substr(s,RLENGTH+1); match(s2,/-+$/); s2=substr(s2,1,RSTART-1); s=substr(s,1,RSTART-1); l=split(s,a,""); for(i=1;i<=l;i++){print i"\t"a[i]}}' > SAT1-Marti.vs.Inoculum.alleles



d<-read.table("table_LD_capsid",header=F)
dsnp<-read.table("88_S88_L001.vs.Inoculum.haplotypeSNPS",header=F)
dsat1<-read.table("SAT1-Marti.vs.Inoculum.alleles",header=F)
pos_flip<-dsnp[as.character(dsat1[dsnp[,1],2])==as.character(dsnp[,6]),1]
al1<-sapply(d[,1],function(x){as.character(dsnp[dsnp[,1]==x,5])})
al2<-sapply(d[,2],function(x){as.character(dsnp[dsnp[,1]==x,5])})
al<-paste(al1,al2,sep="")
s<-rowSums(d[,7:10]); pa<-(d[,7]+d[,9])/s; pb<-(d[,7]+d[,10])/s; ld<-d[,7]/s-pa*pb; ldmax<-((ld>0)*pmin(pa*(1-pb),pb*(1-pa))+(ld<0)*pmin(pa*pb,(1-pa)*(1-pb)))
#ld<-ld*(2*(sapply(strsplit(as.character(d[,3]),split=""),function(x){x[1]})==al1)-1)*(2*(sapply(strsplit(as.character(d[,3]),split=""),function(x){x[2]})==al2)-1)
nn<-which((d[,2]-d[,1]>0) & (d[,2]==sapply(d[,1],function(x){min(d[d[,2]>x,2])})))
plot((d[nn,1]+d[nn,2])/2,abs(ld[nn])/ldmax[nn],type="l",ylim=c(-1,1),ylab="Linkage disequilibrium D' among consecutive haplotype mutations",xlab="Position along capsid sequence",lwd=2)
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")

dev.new()
plot(dsnp[dsnp[,4]>0.25,1],0.5+cumprod(sign(c(1,ld[nn])))*abs(dsnp[dsnp[,4]>0.25,4]-0.5),pch=19,ylab="Frequency of the haplotype allele",xlab="Position along capsid sequence",,ylim=c(0,1))
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")

dev.new()
plot(dsnp[,1],dsnp[,4],pch=19,ylab="Frequency of the minor allele",xlab="Position along capsid sequence",ylim=c(0,0.5))
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="freqmajor.pdf")

dev.new()
plot((d[nn,1]+d[nn,2])/2,(ld[nn])/ldmax[nn],type="l",ylim=c(-1,1),ylab="Linkage disequilibrium D' among consecutive major mutations",xlab="Position along capsid sequence",lwd=2)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
lines((d[nn,1]+d[nn,2])/2,cumprod(sign(ld[nn]/ldmax[nn])),type="l",lty=3,col="green")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDconsmajor.pdf")

dev.new()
par(mfrow=c(1,2))
hist(d[,4],breaks=c(0:10)/20,ylab="Number of mutations",xlab="Minor frequency in the sample",main="Frequency distribution of SNPs in inoculum",col=c(rep("green",4),rep("grey40",6)))
par(mar=c(5,4,2,4))
plot(dsnp[dsnp[,4]>0.25,1],0.5+(-1)^(0+(-1)^(0+(dsnp[dsnp[,4]>0.25,1] %in% pos_flip)))*abs(dsnp[dsnp[,4]>0.25,4]-0.5),pch=19,ylab="Frequency of the derived allele",xlab="Position along capsid sequence",ylim=c(0,1),col="grey40")
points(dsnp[dsnp[,4]<=0.25,1],0.5+(-1)^(0+(-1)^(0+(dsnp[dsnp[,4]<=0.25,1] %in% pos_flip)))*abs(dsnp[dsnp[,4]<=0.25,4]-0.5),pch=19,col="green")
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891)-(3180-1640),lty=3,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
lines(lowess(dsnp[dsnp[,4]>0.25,1],0.5+(-1)^(0+(-1)^(0+(dsnp[dsnp[,4]>0.25,1] %in% pos_flip)))*abs(dsnp[dsnp[,4]>0.25,4]-0.5),f=1/2),col="grey40",lty=2)
lines((d[nn,1]+d[nn,2])/2,((ld[nn])/ldmax[nn]*(-1)^(0+(d[nn,1] %in% pos_flip)+(d[nn,2] %in% pos_flip))+1)/2,lwd=2,col="grey40")
axis(4,at=c(0:4)/4,labels=c((-2):2)/2)
mtext("Linkage disequilibrium |D'| among consecutive derived mutations", side=4,line=2.5)
dev.copy2pdf(file="fig_inoculum_snpfreqld.pdf")

dev.new()
plot((d[nn,1]+d[nn,2])/2,(ld[nn])/ldmax[nn]*(-1)^(0+(d[nn,1] %in% pos_flip)+(d[nn,2] %in% pos_flip)),type="l",ylim=c(-1,1),ylab="Linkage disequilibrium D' among consecutive derived mutations",xlab="Position along capsid sequence",lwd=2)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
#lines((d[nn,1]+d[nn,2])/2,cumprod(sign((ld[nn])/ldmax[nn]*(-1)^(0+(d[nn,1] %in% pos_flip)+(d[nn,2] %in% pos_flip)))),type="l",lty=3,col="green")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDconsderived.pdf")

dev.new()
mincov<-10000
plot(c(d[s>mincov & pa>0.25 & pb>0.25,1],d[s>mincov & pa>0.25 & pb>0.25,2]),rep(ld[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],2),ylim=c(-1,1),ylab="Linkage disequilibrium D' among major mutations",xlab="Position along capsid sequence",pch=19,cex=0.1)
arrows(d[s>mincov & pa>0.25 & pb>0.25,1],ld[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],d[s>mincov & pa>0.25 & pb>0.25,2],ld[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],length=0.01,angle=90,code=3,col="red")
points((d[s>mincov & pa>0.25 & pb>0.25,1]+d[s>mincov & pa>0.25 & pb>0.25,2])/2,ld[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],cex=0.5,pch=19)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDmajor.pdf")

dev.new()
mincov<-10000
ldd<-ld*(-1)^(0+(d[,1] %in% pos_flip)+(d[,2] %in% pos_flip))
plot(c(d[s>mincov & pa>0.25 & pb>0.25,1],d[s>mincov & pa>0.25 & pb>0.25,2]),rep(ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],2),ylim=c(-1,1),ylab="Linkage disequilibrium D' among derived mutations",xlab="Position along capsid sequence",pch=19,cex=0.1)
arrows(d[s>mincov & pa>0.25 & pb>0.25,1],ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],d[s>mincov & pa>0.25 & pb>0.25,2],ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],length=0.01,angle=90,code=3,col="red")
points((d[s>mincov & pa>0.25 & pb>0.25,1]+d[s>mincov & pa>0.25 & pb>0.25,2])/2,ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],cex=0.5,pch=19)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDderived.pdf")

dev.new()
mincov<-10000
ldd<-ld*(-1)^(0+(d[,1] %in% pos_flip)+(d[,2] %in% pos_flip))
plot(3180-1640+c(d[s>mincov & pa>0.25 & pb>0.25,1],d[s>mincov & pa>0.25 & pb>0.25,2]),rep(ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],2),ylim=c(-1,1),ylab="Linkage disequilibrium D' among derived mutations",xlab="Position along genome sequence",pch=19,cex=0.1)
arrows(3180-1640+d[s>mincov & pa>0.25 & pb>0.25,1],ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],3180-1640+d[s>mincov & pa>0.25 & pb>0.25,2],ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],length=0.01,angle=90,code=3,col="red",lwd=2)
points(3180-1640+(d[s>mincov & pa>0.25 & pb>0.25,1]+d[s>mincov & pa>0.25 & pb>0.25,2])/2,ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25],cex=0.5,pch=19)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353),lty=2,col="blue")
mtext(text="1A",side=3,at=(1605+1860)/2,col="blue")
mtext(text="1B",side=3,at=(1860+2517)/2,col="blue")
mtext(text="1C",side=3,at=(2517+3180)/2,col="blue")
mtext(text="1D",side=3,at=(3180+3837)/2,col="blue")
dev.copy2pdf(file="LDderived_final.pdf")

plot(abs(d[s>mincov & pa>0.25 & pb>0.25,1]-d[s>mincov & pa>0.25 & pb>0.25,2]),abs(ldd[s>mincov & pa>0.25 & pb>0.25]/ldmax[s>mincov & pa>0.25 & pb>0.25]),cex=0.5,pch=19,xlab="Distance between SNPs",ylab="|D'|")
dev.copy2pdf(file="LDdecay.pdf")

mincov<-10000
n2n<-sapply(c(1:(length(nn)-1)),function(x){which(d[,1]==d[nn[x],1] & d[,2]==d[nn[x+1],2])})
nn1<-nn[c(1:(length(nn)-1))]
nn2<-nn[c(2:(length(nn)))]
cond<-ldd[nn1]>0 & ldd[nn2]>0 & s[nn1]>mincov & s[nn2]>mincov & s[n2n]>mincov
plot(((ldd/ldmax)[nn1[cond]]*(ldd/ldmax)[nn2[cond]]),((ldd/ldmax)[n2n[cond]]),pch=19,xlab="Predicted D' among next-to-nearest mutations",ylab="D' among next-to-nearest mutations")
#points(pmin((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]]),((ldd/ldmax)[n2n[cond]]),col="red",pch=19)
abline(0,1,lty=2)
dev.copy2pdf(file="LDnn_compare.pdf")
dev.new()
plot(d[nn2[cond],1],((ldd/ldmax)[n2n[cond]])-((ldd/ldmax)[nn1[cond]]*(ldd/ldmax)[nn2[cond]]),xlab="Position along capsid sequence",ylab="D' - predicted D' among next-to-nearest mutations",pch=19)
abline(h=0,lty=2)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDnn_pos.pdf")
# test n2n sign
spp_n2n<-sum(pmin((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])<((ldd/ldmax)[n2n[cond]]) & pmax((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])<((ldd/ldmax)[n2n[cond]]))
spm_n2n<-sum(pmin((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])<((ldd/ldmax)[n2n[cond]]) & pmax((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])>((ldd/ldmax)[n2n[cond]]))
smm_n2n<-sum(pmin((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])>((ldd/ldmax)[n2n[cond]]) & pmax((ldd/ldmax)[nn1[cond]],(ldd/ldmax)[nn2[cond]])>((ldd/ldmax)[n2n[cond]]))
s_n2n<-spp_n2n+spm_n2n+smm_n2n
pv<-0;for(i in 0:s_n2n){for (j in 0:(s_n2n-i)){v<-dbinom(i,s_n2n,1/3)*dbinom(j,s_n2n-i,1/2); if(v<dbinom(spp_n2n,s_n2n,1/3)*(dbinom(spm_n2n,s_n2n-spp_n2n,1/2))){pv<-pv+v}}}
dev.new()
n3n<-sapply(c(1:(length(nn)-2)),function(x){which(d[,1]==d[nn[x],1] & d[,2]==d[nn[x+2],2])})
nn1a<-nn[c(1:(length(nn)-2))]
nn2a<-n2n[c(2:(length(n2n)))]
nn1b<-n2n[c(1:(length(n2n)-1))]
nn2b<-nn[c(3:(length(nn)))]
cond<-ldd[nn1a]>0 & ldd[nn1b]>0 & ldd[nn2a]>0 & ldd[nn2b]>0 & s[nn1a]>mincov & s[nn2a]>mincov & s[nn1b]>mincov & s[nn2b]>mincov & s[n3n]>mincov
plot(pmax(((ldd/ldmax)[nn1a[cond]]*(ldd/ldmax)[nn2a[cond]]),((ldd/ldmax)[nn1b[cond]]*(ldd/ldmax)[nn2b[cond]])),((ldd/ldmax)[n3n[cond]]),pch=19,xlab="Predicted D' among next-next-to-nearest mutations",ylab="D' among next-next-to-nearest mutations")
abline(0,1,lty=2)
dev.copy2pdf(file="LDnnn_compare.pdf")
dev.new()
plot((d[nn1a[cond],2]+d[nn2b[cond],1])/2,((ldd/ldmax)[n3n[cond]])-(pmax(((ldd/ldmax)[nn1a[cond]]*(ldd/ldmax)[nn2a[cond]]),((ldd/ldmax)[nn1b[cond]]*(ldd/ldmax)[nn2b[cond]]))),xlab="Position along capsid sequence",ylab="D' - predicted D' among next-next-to-nearest mutations",pch=19)
abline(h=0,lty=2)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353)-(3180-1640),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2-(3180-1640),col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2-(3180-1640),col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2-(3180-1640),col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2-(3180-1640),col="blue")
dev.copy2pdf(file="LDnnn_pos.pdf")
# test n3n sign
spp_n3n<-sum(pmin((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])<((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])<((ldd/ldmax)[n3n[cond]]))
#spp_n3n<-sum(pmin((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])<((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])<((ldd/ldmax)[n3n[cond]]))
spm_n3n<-sum(pmin((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])<((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])>((ldd/ldmax)[n3n[cond]]))
#spm_n3n<-sum(pmin((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])<((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])>((ldd/ldmax)[n3n[cond]]))
smm_n3n<-sum(pmin((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])>((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1a[cond]],(ldd/ldmax)[nn2a[cond]])>((ldd/ldmax)[n3n[cond]]))
#smm_n3n<-sum(pmin((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])>((ldd/ldmax)[n3n[cond]]) & pmax((ldd/ldmax)[nn1b[cond]],(ldd/ldmax)[nn2b[cond]])>((ldd/ldmax)[n3n[cond]]))
s_n3n<-spp_n3n+spm_n3n+smm_n3n
pv<-0;for(i in 0:s_n3n){for (j in 0:(s_n3n-i)){v<-dbinom(i,s_n3n,1/3)*dbinom(j,s_n3n-i,1/2); if(v<dbinom(spp_n3n,s_n3n,1/3)*(dbinom(spm_n3n,s_n3n-spp_n3n,1/2))){pv<-pv+v}}}




d<-read.table("table_recrate_capsid",header=F)
plot(d[,1]+(3180-1640),d[,2],type="l",ylab="Effective recombination rate per base",xlab="Position along genome sequence",ylim=c(0,0.02))
#ylim=c(0,1.1*max(d[,2])))
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353),lty=2,col="blue")
mtext(text="VP4",side=3,at=(1605+1860)/2,col="blue")
mtext(text="VP2",side=3,at=(1860+2517)/2,col="blue")
mtext(text="VP3",side=3,at=(2517+3180)/2,col="blue")
mtext(text="VP1",side=3,at=(3180+3837)/2,col="blue")
dev.copy2pdf(file="recrate_capsid.pdf")
#dev.new()
#plot(d[,1],convolve(d[,2],exp(-(c(0:200)-200)^2/2/50^2)/sum(exp(-(c(0:200)-200)^2/2/50^2)),type="open")[1:length(d[,1])],type="l",ylab="Effective recombination rate per base",xlab="Position along capsid sequence")
pdf(file="recrate_plot_capsid.pdf")
plot(d[,1]+(3180-1640),d[,2],type="l",ylab="Recombination rate per base in inoculum",xlab="Position along genome sequence",ylim=c(0,0.015),col="lightgray")
xx<-c(d[,1]+(3180-1640),rev(d[,1]+(3180-1640)))
yy<-c(rep(0,length(d[,1])),rev(d[,2]))
polygon(xx,yy,col="lightgray",border=NA)
lines((3180-1640)+d[,1],ma(d[,2]),col=1,lwd=2)
lines((3180-1640)+d1[1,1]-1+c(1:length(rectrace)),ma(rectrace),col=2,lwd=2)
abline(v=c(1008,1086,1605,1860,2517,3180,3837,3891,4353),lty=2,col="blue")
mtext(text="1A",side=3,at=(1605+1860)/2,col="blue")
mtext(text="1B",side=3,at=(1860+2517)/2,col="blue")
mtext(text="1C",side=3,at=(2517+3180)/2,col="blue")
mtext(text="1D",side=3,at=(3180+3837)/2,col="blue")
mtext(text="2A",side=3,at=(3837+3891)/2,col="blue")
mtext(text="2B",side=3,at=(3891+4353)/2,col="blue")
dev.off()
#dev.copy2pdf(file="recrate_capsid_averaged100.pdf")



# Model of suppression of recombination
pdf(file="suppressionofrecombination.pdf",height=4,width=8)
etlist<-c(); dclist<-c();  ic<-1; for(dc in c(0.05,0.1,0.2)){ie<-1; for (et in 1/rev(c(0.005,0.01,0.02,0.05))^2) { curve((1-exp(-x/dc))/(1+exp(et*x^2)),xlim=c(0.001,0.2),ylim=c(0.0001,1),add=((ie*ic)!=1),col=ie+1,lty=4-ic,ylab="Suppression factor",xlab="Genetic divergence between recombining strains",lwd=3,log="y",yaxt="n");dclist<-c(dclist,dc);etlist<-c(etlist,et);ie<-ie+1};ic<-ic+1}
axis(2,at=c(1,0.1,0.01,0.001,0.0001),labels=c(1,0.1,0.01,0.001,0.0001))
ic<-1; for(dc in c(0.05,0.1,0.2)){ie<-1; for (et in 1/c(0.005,0.01,0.02,0.05)^2) { curve(1-exp(-x/dc),xlim=c(0,0.2),ylim=c(0,1),add=T,col="gray",lty=4-ic,lwd=3);ie<-ie+1};ic<-ic+1}
#legend("right",c("cross-imm,epi",paste(dclist,1/sqrt(etlist),sep=",")),col=c(0:9),lwd=2)
legend("right",c("Cross-immunity:",paste(c(0.05,0.1,0.2),sep=","),"Epistasis:","none",paste(1/rev(c(0.005,0.01,0.02,0.05)),sep=",")),col=c(0,1,1,1,0,"gray",2,3,4,5),lty=c(1,3,2,1,1,1,1,1,1,1),lwd=3,bg="white")
dev.off()

######################################
# Analysis of sequences from buffalos

cat references+final.mafft-ginsi.fa | grep ">" | awk '{seq=substr($1,2); ls=split(seq,a,"|"); name=a[ls]; c=0; if(match(name,/(19|44|[Xx]4)/)){ind=toupper(substr(name,RSTART,RLENGTH)); c=c+1}; if(match(name,/([Pp][hH][Tt]|[Pp][tT][Tt]|[Dd][Ss][Pp])/)){tissue=toupper(substr(name,RSTART,RLENGTH)); c=c+1}; if(match(name,/(G[Cc][0-9]+|[Ee][pP][iI][0-9]+|[Cc][Rr][0-9]+)/)){cell=toupper(substr(name,RSTART,RLENGTH)); clstart=RSTART+RLENGTH; match(cell,/(G[Cc]|[Ee][pP][iI]|[Cc][Rr])/) ; celltype=toupper(substr(cell,RSTART,RLENGTH)); c=c+1;  clstr=toupper(substr(name,clstart)); if(match(clstr,/[0-9]+/)){cl=substr(clstr,RSTART,RLENGTH)} else {cl=0}}; if(c==3){print seq"\t"ind"\t"tissue"\t"celltype"\t"cell"\t"cl} else {print seq"\tNA\tNA\tNA\tNA\tNA"}}' > table.final




library(ape)
data<-read.table("table.final",header=F,stringsAsFactors=F)
data[1,]<-c(data[1,1],rep("INO",5))
data[2,]<-c(data[2,1],rep("REF",5))
rownames(data)<-as.vector(data[,1])
seqs<-read.dna("references+final.mafft-ginsi.fa",format="fasta")
seq<-as.character(seqs)
seq<-seq[,seq[1,]!="-"]
data[data[,4]=="CR" & !is.na(data[,4]),4]<-"EPI"
data[,7]<-apply(data,1,function(x){paste(x[c(2,3,5)],collapse="_")})
data[data[,7]=="NA_NA_NA",7]<-NA
data[,8]<-apply(data,1,function(x){paste(x[c(2,3)],collapse="_")})
data[data[,8]=="NA_NA",8]<-NA
data[,9]<-apply(data,1,function(x){paste(x[c(2,3,4)],collapse="_")})
data[data[,9]=="NA_NA_NA",9]<-NA
data[,10]<-apply(data,1,function(x){paste(x[c(2,4)],collapse="_")})
data[data[,10]=="NA_NA",10]<-NA
data[,11]<-apply(data,1,function(x){paste(x[c(3,4)],collapse="_")})
data[data[,11]=="NA_NA",11]<-NA
## begin: infer error rate
#et<-apply(seq[grep("\\|SAT1",names(seq[,1]))[-c(4,11)],-snps_inoc],2,function(x){y<-table(x); y<-y[order(y,decreasing=T)]; return(y[1])});
#er<-sum(9-et[et<9])/prod(dim(seq[grep("\\|SAT1",names(seq[,1]))[-c(4,11)],-snps_inoc]))
## end
seq<-seq[!is.na(data[,2]),]
d<-data[!is.na(data[,2]),]

list_seqs_buffalo<-unique(data[!is.na(data[,7]),7])[-c(1:2)]

l<-dim(seq)[2]
n<-dim(seq)[1]

#af<-apply(seq[-c(1:2),],2,function(x){y<-table(x); y<-y[names(y)!="-"]; y<-y[order(y,decreasing=T)]; return(sum(y)-y[1])})
af<-apply(seq[-c(1:2),],2,function(x){y<-table(x); y<-y[names(y)!="-"]; y<-y[order(y,decreasing=T)]; return(y[2])})
naf<-apply(seq[-c(1:2),],2,function(x){y<-table(x); y<-y[names(y)!="-"]; y<-y[order(y,decreasing=T)]; return(names(y[2]))})
af2<-apply(seq[-c(1:2),],2,function(x){y<-table(x); y<-y[names(y)!="-"]; y<-y[order(y,decreasing=T)]; return(y[3])})
sf<-apply(seq[-c(1:2),],2,function(x){y<-table(x); y<-y[names(y)!="-"]; y<-y[order(y,decreasing=T)]; return(sum(y))})
af[is.na(af)]<-0
naf[is.na(af)]<-"n"
af2[is.na(af2)]<-0
# Error model: fit exponential function to the first 10 frequency bins
plot(1:24,log10(hist(c(af,af2),breaks=-0.5+c(0:201),plot=F)$counts[2:25]),ylim=c(0,3),pch=c(rep(21,2),rep(19,8),rep(21,14)),ylab="Number of SNPs",xlab="Minor allele count",yaxt="n")
axis(2,at=c(0:3),labels=10^c(0:3))
abline(lm(log10(hist(c(af,af2),breaks=-0.5+c(0:201),plot=F)$counts[1+c(3:10)])~c(3:10)),lty=1)
dev.copy2pdf(file="Sanger_error_model.pdf")

# this fixes the minimum frequency that we trust at 21 (corresponding to <1 false SNP expected on average)
hist(af[af>20]/sf[af>20],breaks=20,col="brown",main="",ylab="Number of SNPs",xlab="minor allele frequency")
dev.copy2pdf(file="buffalo_freqdistr.pdf")





#first base 1640
startseq<-1639
snps<-read.table("/data/ribeca/Bryan/Inoculum/88_S88_L001.vs.Inoculum.counts",skip=3,header=F)
snps<-snps[(snps[,2]>startseq)&(snps[,2]<startseq+l+1),]
snps[,2]<-snps[,2]-startseq
base2num<-c(1:4);names(base2num)<-toupper(c("a","c","g","t"))
num2base<-toupper(c("a","c","g","t"))
snps_inoc<-which(pmin(snps[,7],snps[,4])/(snps[,7]+snps[,4])>0.01)




# Fst tests
samecl<-outer(colnames(dist.snps_hf),rownames(dist.snps_hf),Vectorize(function(x,y){data[x,7]==data[y,7]}))
samecl<-samecl+1
diag(samecl)<-0
samecl<-samecl-1
sameind<-outer(colnames(dist.snps_hf),rownames(dist.snps_hf),Vectorize(function(x,y){data[x,2]==data[y,2]}))
sameind<-sameind+1
diag(sameind)<-0
sameind<-sameind-1
sametis<-outer(colnames(dist.snps_hf),rownames(dist.snps_hf),Vectorize(function(x,y){data[x,3]==data[y,3]}))
sametis<-sametis+1
diag(sametis)<-0
sametis<-sametis-1
samecel<-outer(colnames(dist.snps_hf),rownames(dist.snps_hf),Vectorize(function(x,y){data[x,4]==data[y,4]}))
samecel<-samecel+1
diag(samecel)<-0
samecel<-samecel-1
#
#build weight vectors and matrices
wv_cl<-(1/table(d[-c(1,2),7]))[d[-c(1,2),7]]
wv_all<-(1/colSums(table(d[-c(1,2),c(7,9)])>0))[d[-c(1,2),9]]*wv_cl
wv_ind<-(1/colSums(table(d[-c(1,2),c(9,2)])>0))[d[-c(1,2),2]]*wv_all
wv_tis<-(1/colSums(table(d[-c(1,2),c(9,3)])>0))[d[-c(1,2),3]]*wv_all
wv_cel<-(1/colSums(table(d[-c(1,2),c(9,4)])>0))[d[-c(1,2),4]]*wv_all
wv_ind_tis<-(1/colSums(table(d[-c(1,2),c(9,8)])>0))[d[-c(1,2),8]]*wv_all
wv_ind_cel<-(1/colSums(table(d[-c(1,2),c(9,10)])>0))[d[-c(1,2),10]]*wv_all
wv_tis_cel<-(1/colSums(table(d[-c(1,2),c(9,11)])>0))[d[-c(1,2),11]]*wv_all
names(wv_cl)<-rownames(d[-c(1,2),])
names(wv_all)<-rownames(d[-c(1,2),])
names(wv_ind)<-rownames(d[-c(1,2),])
names(wv_tis)<-rownames(d[-c(1,2),])
names(wv_cel)<-rownames(d[-c(1,2),])
names(wv_ind_tis)<-rownames(d[-c(1,2),])
names(wv_ind_cel)<-rownames(d[-c(1,2),])
names(wv_tis_cel)<-rownames(d[-c(1,2),])
wm_cl<-outer(wv_cl,wv_cl)
wm_all<-outer(wv_all,wv_all)
wm_ind<-outer(wv_ind,wv_ind)
wm_tis<-outer(wv_tis,wv_tis)
wm_cel<-outer(wv_cel,wv_cel)
wm_ind_tis<-outer(wv_ind_tis,wv_ind_tis)
wm_ind_cel<-outer(wv_ind_cel,wv_ind_cel)
wm_tis_cel<-outer(wv_tis_cel,wv_tis_cel)



afall<-as.numeric(sapply(1:l,function(x){weighted.mean(seq[-c(1,2),x]==naf[x],wv_all)}))
afall[is.na(afall)]<-0
af19<-as.numeric(sapply(1:l,function(x){weighted.mean(seq[d[,2]=="19",x]==naf[x],wv_all[d[-c(1,2),2]=="19"])}))
af19[is.na(af19)]<-0
af44<-as.numeric(sapply(1:l,function(x){weighted.mean(seq[d[,2]=="44",x]==naf[x],wv_all[d[-c(1,2),2]=="44"])}))
af44[is.na(af44)]<-0
afX4<-as.numeric(sapply(1:l,function(x){weighted.mean(seq[d[,2]=="X4",x]==naf[x],wv_all[d[-c(1,2),2]=="X4"])}))
afX4[is.na(afX4)]<-0
pdf(file="buffalo_freqseq.pdf")
plot(which(af>20),afall[af>20],pch=19,ylab="Frequency",xlab="Position along VP1",ylim=c(0,0.5))
points(which(af>20),af44[af>20],pch=19,col="red")
points(which(af>20),af19[af>20],pch=19,col="blue")
points(which(af>20),afX4[af>20],pch=19,col="green")
legend("topright",c("all","19","44","X4"),col=c("black","blue","red","green"),pch=19)
dev.off()
pdf(file="buffalo_freqseq_inoc.pdf")
plot(which(af>20 & 1:l %in% snps_inoc),afall[af>20 & 1:l %in% snps_inoc],pch=19,ylab="Frequency",xlab="Position along VP1",ylim=c(0,0.5))
points(which(af>20 & 1:l %in% snps_inoc),af44[af>20 & 1:l %in% snps_inoc],pch=19,col="red")
points(which(af>20 & 1:l %in% snps_inoc),af19[af>20 & 1:l %in% snps_inoc],pch=19,col="blue")
points(which(af>20 & 1:l %in% snps_inoc),afX4[af>20 & 1:l %in% snps_inoc],pch=19,col="green")
legend("topright",c("all","19","44","X4"),col=c("black","blue","red","green"),pch=19)
dev.off()
mean(1-afall[af>20 & 1:l %in% snps_inoc])
mean(1-af19[af>20 & 1:l %in% snps_inoc])
mean(1-af44[af>20 & 1:l %in% snps_inoc])
mean(1-afX4[af>20 & 1:l %in% snps_inoc])
for(i in unique(d[-c(1:2),9])){afx<-as.numeric(sapply(1:l,function(x){weighted.mean(seq[d[,9]==i,x]==naf[x],wv_all[d[-c(1,2),9]==i])})) ; print(c(i,1-mean(afx[af>20 & 1:l %in% snps_inoc])))}
# internal nucleotide diversities:
dv<-c()
for(count in 1:3){
if(count==1){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)}
if(count==2){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
if(count==3){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="19")
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="44")
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="X4")
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("DSP"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PHT"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PTT"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_GC"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_EPI"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_GC"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_EPI"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_GC"))) )
subm<-(samecl<2 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
}
for(count in 1:3){
if(count==1){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)}
if(count==2){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
if(count==3){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="19")
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="44")
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="X4")
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("DSP"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PHT"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PTT"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_GC"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_EPI"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_GC"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_EPI"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_GC"))) )
subm<-(samecl==1 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
}
for(count in 1:3){
if(count==1){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)}
if(count==2){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
if(count==3){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="19")
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="44")
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="X4")
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("DSP"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PHT"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),3] %in% c("PTT"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("DSP_GC"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_EPI"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PHT_GC"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_EPI"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
subv<-(((d[-c(1,2),11] %in% c("PTT_GC"))) )
subm<-(samecl==0 & outer(subv,subv))
dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
}
diversity<-t(matrix(dv,ncol=9))
#rownames(diversity)<-c("all","inoc","new")
colnames(diversity)<-c("19","44","X4","DSP","PHT","PTT","DSP_EPI","DSP_GC","PHT_EPI","PHT_GC","PTT_EPI","PTT_GC")

pdf(file="diversity_within_tissues.pdf",height=4,width=7)
par(mar=c(5,4,4,4))
colsnd<-c("black","blue","red")
plot(rep(1:(dim(diversity)[2]),rep(9,(dim(diversity)[2]))),diversity/l,ylab="Nucleotide diversity",xlab="",xaxt="n",pch=rep(rep(c(19,20,21),c(3,3,3)),(dim(diversity)[2])),cex=2,col=colsnd[1+((c(1:length(diversity))-1) %% 3)],ylim=c(0,max(diversity)/l))
axis(1,at=1:(dim(diversity)[2]),labels=sub("_",",",colnames(diversity)),las=2)
abline(v=3.5,lty=2)
abline(v=6.5,lty=2)
abline(v=8.5,lty=3)
abline(v=10.5,lty=3)
dev.off()

# linear model
ind_t<-9 # 7
cond_t<- (af>20 & (c(1:length(af)) %in% snps_inoc)) # (af>20)
f_t <-function(x){x} # {log(x+0.1/l)}
dvi<-c(); dvin<-c(); for (i in unique(d[,ind_t])[order(unique(d[,ind_t]))]){if(sum(d[,ind_t]==i)>=2){dvin<-c(dvin,rep(i,sum(d[,ind_t]==i)*(sum(d[,ind_t]==i)-1)/2));dvi<-c(dvi,c(dist.dna(as.DNAbin(seq[d[,ind_t]==i,cond_t]),model="N")/l))}}
dvi_animal<-gsub("_.*","",dvin)
dvi_location<-(c(1:length(dvin)) %in% grep("EPI|CR",dvin) )
dvi_tissue<-sub("_[1-9A-Z]+","",sub("[1-9A-Z]+_","",dvin))
lmdvi<-lm(f_t(dvi) ~ dvi_animal + dvi_tissue + dvi_location)
coef(summary(lmdvi))


#......
#dvar<-c()
#for(count in 1:3){
#if(count==1){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)}
#if(count==2){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
#if(count==3){dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)}
##then
#subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="19")
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="44")
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT"))) & d[-c(1,2),2]=="X4")
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),3] %in% c("DSP"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),3] %in% c("PHT"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),3] %in% c("PTT"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("DSP_EPI"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("DSP_GC"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("PHT_EPI"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("PHT_GC"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("PTT_EPI"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#subv<-(((d[-c(1,2),11] %in% c("PTT_GC"))) )
#subm<-(samecl<2 & outer(subv,subv))
#dv<-c(dv,weighted.mean(dist.snps_hf[subm & sameind==1 & sametis==1 & samecel==1],wm_all[subm & sameind==1 & sametis==1 & samecel==1]))
#}

#Recombinants
seqs_snps_inoc<-seq[-c(1:2),af>20 & 1:l %in% snps_inoc]
seqs_snps_inoc<-seq[d[,2]=="X4",af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldp<-abs(ld/ldmax)
ldside<-ldp; ldside[row(ldside) != col(ldside)+1 ]<-0; ldpred<-solve(diag(length(ld[,1]))-ldside)
ldcomp<-lower.tri(ldp, diag = FALSE)*log(ldpred/ldp+10^-4)+upper.tri(ldp, diag = TRUE)*ldp
ldcomp2<-lower.tri(ldp, diag = FALSE)*( log10(abs(log(abs(ldp)+10^(-10))/log(abs(ldpred)+10^(-10)))+10^(-10)) )+upper.tri(ldp, diag = TRUE)*ldp
ldpredloc<-outer(1:20,1:20,Vectorize(function(i,j){if(abs(i-j)<2){return(ldp[i,j])} else{return(max(sapply((min(i,j)+1):(max(i,j)-1),function(k){ldp[i,k]*ldp[k,j]})))}}))
ldcomp3<-lower.tri(ldp, diag = FALSE)*( log10(abs(log(abs(ldp)+10^(-10))/log(abs(ldpredloc)+10^(-10)))+10^(-10)) )+upper.tri(ldp, diag = TRUE)*1000
q<-0.26
Rfun<-function(r,s){sr<-s/r;g<-sqrt((1+sr)^2-8*q*(1-q)*sr); return(-log(1-(1+sr-g/tanh(r*g/2+log((1+sr+g)/(1+sr-g))/2))/(4*q*(1-q)*sr)))}
infsel<-outer(1:20,1:20,Vectorize(function(i,j){if(abs(i-j)<2){return(0)} else{ if(-log(ldpred[max(i,j),min(i,j)])>0 & -log(ldpred[max(i,j),min(i,j)])>-log(ldp[max(i,j),min(i,j)])+0.1) {return(uniroot(function(x){Rfun(-log(ldpred[max(i,j),min(i,j)]),x)+log(ldp[i,j])},interval=pmax(c(0.001,10),c(0.01,100)*(-log(ldpred[max(i,j),min(i,j)]))))$root)} else {return(0)} }}))
infselloc<-outer(1:20,1:20,Vectorize(function(i,j){if(abs(i-j)<2){return(0)} else{ if(-log(ldpredloc[max(i,j),min(i,j)])>0 & -log(ldpredloc[max(i,j),min(i,j)])>-log(ldp[max(i,j),min(i,j)])+0.1) {return(uniroot(function(x){Rfun(-log(ldpredloc[max(i,j),min(i,j)]),x)+log(ldp[i,j])},interval=pmax(c(0.001,10),c(0.01,100)*(-log(ldpredloc[max(i,j),min(i,j)]))))$root)} else {return(0)} }}))
#ldcomp4<-lower.tri(ldp, diag = FALSE)*( log10(abs(log(abs(ldp)+10^(-10))/log(abs(ldpredloc)+10^(-10)))+10^(-10)) )+upper.tri(ldp, diag = TRUE)*ldp
infsellocX4<-infselloc
ldcompsel<-lower.tri(ldp, diag = FALSE)*( -infsel)+diag(dim(ldp)[1])*0+upper.tri(ldp, diag = FALSE)*( infselloc)
library(lattice)
levelplot(t(ldcompsel),column.values=snps_inoc[af[snps_inoc]>20],row.values=snps_inoc[af[snps_inoc]>20],at=c(-max(infsel)*c(100:1)/100,0,max(infselloc)*c(1:100)/100),panel = function(...){
panel.levelplot(...)
panel.abline(0,1)
},col.regions=c(
#colorRampPalette(c("violet","grey20","darkblue"))(401),
colorRampPalette(c("darkblue","blue","lightblue","white"))(100),"white",colorRampPalette(c("white","green","darkgreen","brown"))(100)),xlab="Position in 1D sequence - Epistatic coefficient (pairwise-corrected)",ylab="Position in 1D sequence - Epistatic coefficient")
dev.copy2pdf(file="buffalo_modelepistasis_44.pdf")
max(infselloc)
 max(infsel)
nn<-c(nn,(infselloc[ldind])[codesyn==2])
ns<-c(ns,(infselloc[ldind])[indcompnonsyn])
wilcox.test((infselloc[ldind])[codesyn==2],(infselloc[ldind])[indcompnonsyn],alternative="greater")
fisher.test(matrix(c(sum((infselloc[ldind])[codesyn==2]>0),sum((infselloc[ldind])[codesyn==2]==0),sum((infselloc[ldind])[indcompnonsyn]>0),sum(((infselloc[ldind])[indcompnonsyn]==0))),nrow=2),alternative="greater")
dev.copy2pdf(file="buffalo_modelepistasis_all.pdf")
#levelplot(t(ldcomp),column.values=snps_inoc[af[snps_inoc]>20],row.values=snps_inoc[af[snps_inoc]>20],at=c(-5*c(100:1)/100,c(1:100)/100),col.regions=c(colorRampPalette(c("darkblue","blue","lightblue","white","yellow","orange","red"))(200),"black"),xlab="",ylab="",colorkey=F)
dev.copy2pdf(file="buffalo_LDmatrixepistasis_19.pdf")
plot(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),type="l",ylab="Effective recombination rate per base",xlab="Position in capsid")
dev.copy2pdf(file="buffalo_recrate.pdf")


plot(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),type="l",ylab="Effective recombination rate per base",xlab="Position in capsid",ylim=c(0,0.05),lwd=2)
i<-1
for (inds in unique(d[-c(1:2),2])){
i<-i+1
seqs_snps_inoc<-seq[d[,2]==inds,af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
lines((snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),col=i)
print(inds)
print(mean(rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20]))))
}
seqs_snps_inoc<-seq[-c(1:2),af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
lines(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),col="black",lwd=2)
dtrec<-read.table("../Inoculum/table_recrate_capsid",header=F)
lines(dtrec,lwd=2,lty=2)
legend("topright",c("all inds",paste("ind",unique(d[-c(1:2),2]),sep=" "),"inoculum"),lty=c(1,1,1,1,2),lwd=c(2,1,1,1,2),col=c(1:4,1))
dev.copy2pdf(file="buffalo_recrate_compare.pdf")
# mean recombination rates:
#initial 0.0028
#19: 0.0067
#X4: 0.0072
#44: 0.0117

mycolind<-c("black","orange","darkred","blue")
plot(3180+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),type="l",ylab="Recombination rate per base",xlab="Position along genome sequence",ylim=c(0,0.05),lwd=2,col="white")
i<-1
for (inds in unique(d[-c(1:2),2])){
i<-i+1
seqs_snps_inoc<-seq[d[,2]==inds,af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
lines(3180+2*(i==4)+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),col=mycolind[i],lwd=2,lty=1+2*(i<4))
print(inds)
print(mean(rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20]))))
}
seqs_snps_inoc<-seq[d[,2] %in% c("X4","19"),af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
lines(3180+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20]))),rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20])),col="red",lwd=2)
dtrec<-read.table("../Inoculum/table_recrate_capsid",header=F)
lines(dtrec[,1]+3180-1640,dtrec[,2],lwd=2,lty=1,col="black")
legend("topright",c("inoculum","35 dpi", "400 dpi"),lty=1,lwd=2,col=c("black","red","blue"))
dev.copy2pdf(file="recrates_1D.pdf")

##for(tis in c("PHT_GC","PTT_GC","DSP_GC","PHT_EPI","PTT_EPI","DSP_EPI")){dtrec<-c()
##for(tis in c("PHT","PTT","DSP")){dtrec<-c()
# timepoints....
wl<-130
vals_rec_inoc<-c()
for(nw in 1:floor(650/wl)){
vals_rec_inoc<-c(vals_rec_inoc,mean(dtrec[dtrec[,1] %in% c((1669+wl*(nw-1)):(1669+wl*nw-1)),2]))
}
vals_rec_40d<-c()
seqs_snps_inoc<-seq[d[,2] %in% c("X4"),af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
for(nw in 1:floor(650/wl)){
vals_rec_40d<-c(vals_rec_40d,mean(rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20]))[(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20])))) %in% c((1669+wl*(nw-1)):(1669+wl*nw-1))]))
}
seqs_snps_inoc<-seq[d[,2] %in% c("19"),af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
for(nw in 1:floor(650/wl)){
vals_rec_40d[nw]<-vals_rec_40d[nw]+mean(rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20]))[(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20])))) %in% c((1669+wl*(nw-1)):(1669+wl*nw-1))])
}
vals_rec_40d<-vals_rec_40d/2
vals_rec_400d<-c()
seqs_snps_inoc<-seq[d[,2] %in% c("44"),af>20 & 1:l %in% snps_inoc]
ld<-outer(1:20,1:20,Vectorize(function(x,y){mean((seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]) & (seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
ldmax<-outer(1:20,1:20,Vectorize(function(x,y){max(mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x]),mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y]))-mean(seqs_snps_inoc[,x]==(seq[1,af>20 & 1:l %in% snps_inoc])[x])*mean(seqs_snps_inoc[,y]==(seq[1,af>20 & 1:l %in% snps_inoc])[y])}))
for(nw in 1:floor(650/wl)){
vals_rec_400d<-c(vals_rec_400d,mean(rep(-log(diag((ld/ldmax)[1:19,2:20]))/diff(snps_inoc[af[snps_inoc]>20]),diff(snps_inoc[af[snps_inoc]>20]))[(1640+(snps_inoc[af[snps_inoc]>20])[1]+c(1:sum(diff(snps_inoc[af[snps_inoc]>20])))) %in% c((1669+wl*(nw-1)):(1669+wl*nw-1))]))
}
vals_rec<-rbind(vals_rec_inoc,vals_rec_40d,vals_rec_400d)
vals_rec<-cbind(rowMeans(vals_rec),vals_rec)
##print(tis)
##print(vals_rec[,1])
##}
pdf(file="recrate_timepoints_final.pdf")
atv<-c(1,1+70/400,3)
colsw<-rainbow(5, s = 1, v = 1, start = 2/6, end = 4/6, alpha = 1)
plot(atv,vals_rec[,1],type="b",pch=19,lwd=3,ylim=c(0,max(vals_rec)),xlab="Days post infection",xaxt="n",ylab="Average recombination rate")
axis(side=1,at=atv,label=c("0","35","400"))
for (nw in 2:6){points(atv,vals_rec[,nw],pch=20,col=colsw[nw-1]);lines(atv,vals_rec[,nw],col=colsw[nw-1],lwd=2)}
points(atv,vals_rec[,1],pch=19);lines(atv,vals_rec[,1],lwd=3)
points(atv,c(0.33,4,5.5)/674,pch=19,col="red");lines(atv,c(0.33,4,5.5)/674,pch=19,col="red",lwd=2,lty=2)
legend("bottomright",c("Recombination:","whole sequence","5' end","3' end","","Heterozygosity"),col=c("white","black","green","blue","white","red"),lwd=c(3,3,2,2,2,2),lty=c(1,1,1,1,1,2))
dev.off()


haplo_snps_inoc<-apply(seq[-c(1:2),af>20 & 1:l %in% snps_inoc],1,function(x){paste(x,collapse="",sep="")})
haplotypes_inoc<-table(haplo_snps_inoc)
haplotypes_inoc<-haplotypes_inoc[order(haplotypes_inoc,decreasing=T)]
list_haplo<-names(haplotypes_inoc)
hap1<-strsplit(list_haplo[1],"")[[1]]
hap2<-strsplit(list_haplo[2],"")[[1]]
nrec<-rep(0,17)
for(myhap in list_haplo){smyhap<-strsplit(myhap,"")[[1]];eqt<-(smyhap==hap1)-(smyhap==hap2);for(p in 3:19){if(abs(eqt[p-1]+eqt[p-2]-eqt[p]+eqt[p+1])==4) {nrec[p]<-nrec[p]+haplotypes_inoc[myhap]}}}





#barplot frequencies
allfreqs<-c(
0.82,
0.56,
0.91,
0.99,
1,
0.62,
0.68,
0.66,
0.70,
0.56,
0.76,
0.72,
0.73,
0.75,
0.72,
0.78,
0.69,
0.65,
0.77
);
namesallfreqs<-c(paste("19",c(
"",
"PhT/GC",
"PhT/Epi",
"PtT/GC",
"PtT/Epi",
"DSP/Epi"
),sep=" "),paste("X4",c(
"",
"PhT/GC",
"PhT/Epi",
"PtT/GC",
"PtT/Epi",
"DSP/Epi"
),sep=" "),paste("44",c(
"",
"PhT/GC",
"PhT/Epi",
"PtT/GC",
"PtT/Epi",
"DSP/GC",
"DSP/Epi"
),sep=" "))
pdf(file="freqs_compare.pdf",height=5,width=7)
par(las=2)
par(mar=par("mar")+c(0.4,0,0,0))
barplot(t(matrix(c(0.54,rep(0,length(allfreqs)),0.02,1-allfreqs,0.44,allfreqs),ncol=3)),col=c("darkgreen","green","brown"),names.arg=c("Inoculum",namesallfreqs),space=c(0,0.2+0.8*(namesallfreqs %in% c("X4 ","44 "))+2.5*(namesallfreqs %in% c("19 "))))
dev.off()


#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################




# Weighted Fst approach (with AMOVA-style p-values by permutations)
nsim<-1000
library(permute)
res_perm<-c()
test_name<-c()
#
# Do samples have an effect?
# different clone in same sample vs different samples, but same ind/tissue/cell type
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,other")
#
#
# Do individuals have an effect?
# can test only on "DSP_EPI","PHT","PTT"
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_DSP,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_DSP,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_DSP,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PHT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),3] %in% c("PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-GC,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-GC,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-GC,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-EPI,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-EPI,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PHT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PHT-EPI,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-GC,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-GC,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_GC")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-GC,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-EPI,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-EPI,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("PTT_EPI")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds_PTT-EPI,other")
#
#
#
# Do tissues matter?
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_GC,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_GC,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_GC,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_EPI,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_EPI,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_EPI,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_19,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_19,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_19,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_44,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_44,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_44,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_X4,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_X4,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue_X4,other")
#
# Do cell types matter?
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_DSP,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_DSP,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_DSP,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PHT,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PHT,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PHT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PHT,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PTT,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PTT,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_PTT,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_19,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_19,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("19"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_19,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_44,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_44,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("44"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_44,other")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_X4,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_X4,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),2] %in% c("X4"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell_X4,other")
#

rownames(res_perm)<-test_name

par(mar=c(8,4,4,4))
cols_perm<-rep(c("black","blue","red"),((2+dim(res_perm)[1])%/%3))
res_perm<-res_perm[c(1:3,4:6,58:60,67:69,61:63,64:66,70:72,7:9,40:42,46:48,43:45,49:51,52:54,55:57,10:12,13:15,16:18,19:21,22:24,31:33,34:36,37:39),]
test_name<-test_name[c(1:3,4:6,58:60,67:69,61:63,64:66,70:72,7:9,40:42,46:48,43:45,49:51,52:54,55:57,10:12,13:15,16:18,19:21,22:24,31:33,34:36,37:39)]
plot((2+c(1:dim(res_perm)[1]))%/%3,res_perm[,1],ylab="Fst",xlab="",xaxt="n",pch=19,col=cols_perm,cex=
 0.75+0.75*(p.adjust(res_perm[,2]+1/nsim/2,method="BH")<0.05)+0.5*(p.adjust(res_perm[,2]+1/nsim/2,method="BH")<0.01)+0.5*(p.adjust(res_perm[,2]+1/nsim/2,method="BH")<0.001)
#((0.5-log2(p.adjust(res_perm[,2]+1/nsim/2,method="BH"))%/%4))/2
,ylim=c(-0.1,1));
axis(1,at=1:((2+dim(res_perm)[1])%/%3),labels=sub("_",",",sub("_"," (",sub("([A-Z])$","\\1)",unique(unlist(lapply(strsplit(test_name,","),function(x){as.character(x[1])})))))),las=2)
axis(1,at=1:((2+dim(res_perm)[1])%/%3),labels=sub("_",",",sub("_"," (",sub("([A-Z])$","\\1)",unique(unlist(lapply(strsplit(test_name,","),function(x){as.character(x[1])})))))),las=2)
abline(h=0,lty=2)
#abline(h=1,lty=2)
abline(v=1.5,lty=3)
abline(v=7.5,lty=3)
abline(v=14.5,lty=3)
legend("topright",c("all SNPS","SNPs from inoculum","new SNPs"),col=c("black","blue","red"),pch=19,bg="white")
legend("topleft",c("p<0.05","p<0.01"),col=c("black"),pch=19,pt.cex=c(1.5,2),bg="white")

dev.copy2pdf(file="fst_plot_global.pdf")

#
#
#








#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################




# Weighted Fst approach (with AMOVA-style p-values by permutations)
nsim<-1000
library(permute)
res_perm<-c()
test_name<-c()
#
# Do samples have an effect?
# different clone in same sample vs different samples, but same ind/tissue/cell type
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
vt1<-weighted.mean(dist.snps_hf[samecl==1],wm_cl[samecl==1])/weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1],wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])
vt2<-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==1]<x,wm_cl[samecl==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_hf[samecl==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_cl[samecl==0 & sameind==1 & sametis==1 & samecel==1])})
vt1<-1-vt1
vr1<-c();vr2<-c();
#vr01<-c();vr02<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),9], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl
wm_rand<-wm_cl[myperm,myperm]
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1])/weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==1]<x,wm_rand[samecl_rand==1])})-sapply(3*c(1:4),function(x){weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]<x,wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1])}))
#vr01<-c(vr01,weighted.mean(dist.snps_rand[samecl_rand==1],wm_rand[samecl_rand==1]))
#vr02<-c(vr02,weighted.mean(dist.snps_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1],wm_rand[samecl_rand==0 & sameind==1 & sametis==1 & samecel==1]))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"samples,other")
#
#
# Do individuals have an effect?
# can test only on "DSP_EPI","PHT","PTT"
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% c("DSP_EPI"))|(d[-c(1,2),3] %in% c("PHT","PTT")))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"inds,other")
#

# Do tissues matter?
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),4] %in% c("EPI","GC"))
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"tissue,other")
#

# Do cell types matter?
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,all")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,inoc")
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),3] %in% c("DSP","PHT","PTT"))
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,"cell,other")
#
#
#

#Fst between cell types
for(mysub in unique(d[-c(1,2),8])){
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("cell_",mysub,",all",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("cell_",mysub,",inoc",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-(samecl<2 & sameind==1 & sametis==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & samecel==1],wm_all[subm & samecel==1])/weighted.mean(dist.snps_hf[subm & samecel==0],wm_all[subm & samecel==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & samecel==1]<x,wm_all[subm & samecel==1]))-(weighted.mean(dist.snps_hf[subm & samecel==0]<x,wm_all[subm & samecel==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),8], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),8] %in% mysub)
subm<-( sameind==1 & sametis==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & samecel==1],wm_rand[subm & samecel==1])/weighted.mean(dist.snps_rand[subm & samecel==0],wm_rand[subm & samecel==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & samecel==1]<x,wm_rand[subm & samecel==1]))-(weighted.mean(dist.snps_rand[subm & samecel==0]<x,wm_rand[subm & samecel==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("cell_",mysub,",other",sep=""))
#
}
#
#
#


#Fst between tissues
for(mysub in unique(d[-c(1,2),10])){
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 ]),model="N",as.matrix=T)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("tissue_",mysub,",all",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("tissue_",mysub,",inoc",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-(samecl<2 & sameind==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sametis==1],wm_all[subm & sametis==1])/weighted.mean(dist.snps_hf[subm & sametis==0],wm_all[subm & sametis==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sametis==1]<x,wm_all[subm & sametis==1]))-(weighted.mean(dist.snps_hf[subm & sametis==0]<x,wm_all[subm & sametis==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),10], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-(d[-c(1,2),10] %in% mysub)
subm<-( sameind==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sametis==1],wm_rand[subm & sametis==1])/weighted.mean(dist.snps_rand[subm & sametis==0],wm_rand[subm & sametis==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sametis==1]<x,wm_rand[subm & sametis==1]))-(weighted.mean(dist.snps_rand[subm & sametis==0]<x,wm_rand[subm & sametis==0]))}))
}
sum(vr1>=vt1)/nsim
for(i in 1:4){print(sum(vr2[i,]>=vt2[i])/nsim)}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("tissue_",mysub,",other",sep=""))
#
}
#
#
#

#Fst between individuals
for(mysub in unique(d[-c(1,2),11])){
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("ind_",mysub,",all",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & (c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("ind_",mysub,",inoc",sep=""))
#
dist.snps_hf<-dist.dna(as.DNAbin(seq[-c(1,2),af>20 & !(c(1:length(af)) %in% snps_inoc)]),model="N",as.matrix=T)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-(samecl<2 & sametis==1 & samecel==1 & outer(subv,subv))
vt1<-weighted.mean(dist.snps_hf[subm & sameind==1],wm_all[subm & sameind==1])/weighted.mean(dist.snps_hf[subm & sameind==0],wm_all[subm & sameind==0])
vt2<-sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_hf[subm & sameind==1]<x,wm_all[subm & sameind==1]))-(weighted.mean(dist.snps_hf[subm & sameind==0]<x,wm_all[subm & sameind==0]))})
vt1<-1-vt1
vr1<-c();vr2<-c();
for(i in 1:nsim){
mycontrol<-how(within = Within(type = "free"), plots = Plots(strata = d[-c(1,2),11], type = "none"), blocks = NULL, nperm = nsim)
myperm<-permute(i, n-2, mycontrol)
dist.snps_rand<-dist.snps_hf[myperm,myperm]
samecl_rand<-samecl[myperm,myperm]
wtemp<-1/((sapply(unique(d[-c(1,2),9]),function(x){sum((wv_all[myperm])[d[-c(1:2),9]==x])}))[d[-c(1,2),9]]); wm_rand<-wm_all[myperm,myperm]*outer(wtemp,wtemp)
subv<-((d[-c(1,2),11] %in% mysub))
subm<-( sametis==1 & samecel==1 & outer(subv,subv))
vr1<-c(vr1,1-weighted.mean(dist.snps_rand[subm & sameind==1],wm_rand[subm & sameind==1])/weighted.mean(dist.snps_rand[subm & sameind==0],wm_rand[subm & sameind==0]))
vr2<-cbind(vr2,sapply(3*c(1:4),function(x){(weighted.mean(dist.snps_rand[subm & sameind==1]<x,wm_rand[subm & sameind==1]))-(weighted.mean(dist.snps_rand[subm & sameind==0]<x,wm_rand[subm & sameind==0]))}))
}
res1<-sum(vr1>=vt1)/nsim
res2<-c();for(i in 1:4){res2<-c(res2,sum(vr2[i,]>=vt2[i])/nsim)}
res_perm<-rbind(res_perm,c(vt1,res1,as.vector(t(matrix(c(vt2,res2),ncol=2)))))
test_name<-c(test_name,paste("ind_",mysub,",other",sep=""))
#
}



