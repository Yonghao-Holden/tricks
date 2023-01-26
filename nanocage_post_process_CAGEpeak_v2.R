args=commandArgs( trailingOnly=T )
library(stringr)
library(DESeq2)
ifile=args[1]
capfilter=args[2]
oprefix=args[3]

#ifile="/scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/peak/deepseq_20M/8cells_deepseq.peak.nrPassThreshold_1.paraclu.txt"
#capfilter=0.15
#oprefix="/scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/peak/deepseq_20M/8cells_deepseq.peak.nrPassThreshold_1.paraclu"

#ifile="nanocage_keratinocyte.CTSS"
#capfilter=0.15
#oprefix="nanocage_keratinocyte"

capsuffix=paste( paste(".capf", capfilter, sep=""), ".txt", sep="")
normsuffix=paste( paste(".capf", capfilter, sep=""), ".norm.txt", sep="")
normtpmsuffix=paste( paste(".capf", capfilter, sep=""), ".normtpm.txt", sep="")
allbedsuffix=paste( paste(".capf", capfilter, sep=""), ".bed", sep="" );

ofile_raw=paste(oprefix, ".raw.txt", sep="")
ofile_capfilter=paste(oprefix, capsuffix, sep="")
ofile_normal=paste(oprefix, normsuffix, sep="")
ofile_normal_tpm=paste(oprefix, normtpmsuffix, sep="")
ofile_bed=paste(oprefix, allbedsuffix, sep="")

df=read.table( ifile, header=T )
nSample=(dim(df)[2]-7)/4
#Input file column description
##Column 1-7:
##Column 8-N: TPM, TPM (UnencodedG), rawcount, rawcount(Unencoded G)

#Four tasks
#1. Print raw count
#2. CapFilter
#3. Normalization
#4. TXT to bed Note txt is 1-based coordinate

#1. Print raw count only
print("Extract raw count")
nStartIndex=7+nSample*2+1
nEndIndex=nStartIndex+nSample-1
df.raw=df[, c(1:7, nStartIndex:nEndIndex)]
#Rename the column name
vecSampleName=colnames(df.raw[, 8:(dim(df.raw)[2])])
#df.index=str_locate(vecSampleName, ".nanoCAGE")
df.index=str_locate(vecSampleName, "_R1.fastq.gz")
vecSampleName=substr( vecSampleName, 9, df.index[, 1]-1)
colnames( df.raw ) = c( colnames(df.raw[, 1:7]), vecSampleName )
write.table( df.raw, ofile_raw, quote=F, row.names=F, col.names=T, sep="\t")


#2. CapFilter
print("CapFilter")
df.capf=df.raw[df.raw$maxUGPercentage>=capfilter, ] 
write.table( df.capf, ofile_capfilter, quote=F, row.names=F, col.names=T, sep="\t")


#3. Normalization
print("Normalization")
df.capf.data=df.capf[, 8:(8+nSample-1)]

#metadata
#vecCell=c("COV413a", "DMS454", "DMS53", "H2110", "H2228", "H716", "SKMEL2", "SNU475" )
#df.metadata=as.data.frame( matrix( data=c("OESC", "SCLC", "SCLC", "NSCLC", "NSCLC", "COL", "MEL", "HEPC"), nrow=length(vecCell)), row.names=vecCell)
#colnames(df.metadata)=c("cancer")

df.metadata=as.data.frame(matrix(data = rep(c("sample"), times=nSample)), nrow=nSample, row.names=vecSampleName)
colnames(df.metadata)=c("sample")

#Double check the metadata and the column of data match
if( FALSE==all(rownames(df.metadata) == colnames(df.capf.data)) ){
  print("Metadata is not valid")
}

#Normalization
#Let's say two set of cells: treated, untreated. Dispersion in two conditions may differ then add the "treatment" to the design and let the DESEQ2 consider "design" when calculating dispersion
dds_blind <- DESeqDataSetFromMatrix(countData = df.capf.data,
                                    colData = df.metadata,
                                    design = ~ 1) #design = ~1 indiciate not considering the experimental design.
dds_blind=estimateSizeFactors( dds_blind ) 

counts.norm=counts( dds_blind, normalized=T )
df.capf.norm=cbind(df.capf[, 1:7], counts.norm)
write.table( df.capf.norm, ofile_normal, quote=F, row.names=F, col.names=T, sep="\t")


#4. Normalized count to TPM
#Normalized count depends on the library size. Depending on the library size, the noise may have a several tags
#To call peaks out of noise, TPM is widely used.
#Using the median size of normalized library size, I will calculate the scaling factor. 
#For each sample, I will divide the normalized tag count by the scaling factor to make it normalized count per million
counts.normtpm=fpm( dds_blind, robust= T )
df.capf.norm.tpm=cbind(df.capf[, 1:7], counts.normtpm)
write.table( df.capf.norm.tpm, ofile_normal_tpm, quote=F, row.names=F, col.names=T, sep="\t")


#5. TXT to BED conversion
#chr start end conscluster norm_tpm strand maxUGpercentage norm_rawcount
#
print("BED conversion")
odir=dirname( oprefix )
df.capf.norm.tpm.bed=df.capf.norm.tpm;
df.capf.norm.tpm.bed$start=df.capf.norm.tpm.bed$start-1
df.capf.norm.tpm.bed$peakname=paste( "conscluster", df.capf.norm.tpm.bed$consensus.cluster, sep="" )
write.table( df.capf.norm.tpm.bed[, c(2, 3, 4, 7+nSample+1, 6, 5, 7:(7+nSample) )], ofile_bed, quote=F, row.names=F, col.names=T, sep="\t")


#for(i in vecSampleName){
#  obedname=paste( i, paste(".peak.capf", capfilter, sep=""), sep="" );
#  obedname=paste( obedname, ".bed", sep="");
#  obedname=paste( odir, obedname, sep="/");
#  
#  df.persample=cbind( df.capf.norm.tpm.bed[, c("chr", "start", "end", "peakname", i, "strand", "maxUGPercentage")], df.capf.norm[, c(i)] )
#  write.table( df.persample, obedname, quote=F, row.names=F, col.names=F, sep="\t")
#} 


#0. CAGEr TPM
#This value is not normlized across samples.
#ofile_cagertpm=paste(oprefix, "capf0.15.cager_tpm.txt", sep="")
#df.cager_tpm=df[, c(1:(8+nSample-1))]
#df.cager_tpm.capf=df.cager_tpm[df.cager_tpm$maxUGPercentage>=capfilter, ]
#colnames( df.cager_tpm.capf ) = c( colnames(df.cager_tpm[, 1:7]), vecSampleName )
#write.table( df.cager_tpm.capf, ofile_cagertpm, quote=F, row.names=F, col.names=T, sep="\t")



#############################################
##Let's see whether the division by the normalized library size will change the set of genes with comparable expression level##
##                                         ##
##                                         ##
#############################################
#tdf.norm.ratio=sweep( counts.norm, 1, rowMedians( as.matrix( counts.norm) ), FUN='/')
#tdf.normtpm.ratio=sweep( counts.normtpm, 1, rowMedians( as.matrix( counts.normtpm) ), FUN='/')

#tdf.std=data.frame( clusterid=df.capf$consensus.cluster, std=rowSds(as.matrix( tdf.norm.ratio) ), stdtpm=rowSds( as.matrix(  tdf.normtpm.ratio) ) )
#tdf.std=tdf.std[ !is.na(tdf.std$std) & !is.na(tdf.std$stdtpm), ]
#ggplot()+
#  geom_point( data=tdf.std, aes(x=std, y=stdtpm))+xlab("Standard deviation of the ratio (Normalized expression count)")+ylab("Standard deviation of the ratio (Normalized TPM)")
#############################################
##                                         ##
#############################################
#colMins( as.matrix(vst_blind.tpm))
#colMaxs( as.matrix(vst_blind.tpm) )


##Let's see whether the normalized count in edgeR the same as DEseq2
##                                         ##
#############################################
#library(edgeR)
#df.edger=DGEList(assay( dds_blind ) )
#df.edger=calcNormFactors( df.edger, method="RLE")
#df.edger$samples$norm.factors
#counts.normtpm.edger=cpm( df.edger, normalized.lib.sizes = T, log=F)

#df.tdf1=as.data.frame( counts.normtpm); df.tdf1$cluster=rownames( counts.normtpm); df.tdf1$type="DEseq2"; df.tdf1=gather(df.tdf1, "cell", "tpm", COV413a:SNU475)
#df.tdf2=as.data.frame( counts.normtpm.edger); df.tdf2$cluster=rownames( counts.normtpm.edger); df.tdf2$type="edgeR"; df.tdf2=gather(df.tdf2, "cell", "tpm", COV413a:SNU475)
#df.plot=df.tdf1; df.plot$tpm_edger=df.tdf2$tpm
#df.plot$ratio=df.plot$tpm_edger/df.plot$tpm
#ggplot()+
#  geom_point( data=df.plot, aes(x=tpm, y=tpm_edger))+
#  facet_grid(.~cell)+xlab("Normalized TPM by DEseq2")+ylab("Normalized TPM by edgeR")
#ggplot()+
#  geom_point( data=df.plot, aes(x=tpm, y=ratio))+
#  facet_grid(.~cell)+xlab("Normalized TPM by DEseq2")+ylab("Ratio of normalized TPM (edgeR/DEseq2)")


#counts.normtpm2=fpm( dds_blind, robust= T )
#df.tdf1=as.data.frame( counts.normtpm2); df.tdf1$cluster=rownames( counts.normtpm2); df.tdf1$type="DEseq2"; df.tdf1=gather(df.tdf1, "cell", "tpm", COV413a:SNU475)
#df.tdf2=as.data.frame( counts.normtpm.edger); df.tdf2$cluster=rownames( counts.normtpm.edger); df.tdf2$type="edgeR"; df.tdf2=gather(df.tdf2, "cell", "tpm", COV413a:SNU475)
#df.plot=df.tdf1; df.plot$tpm_edger=df.tdf2$tpm
#df.plot$ratio=df.plot$tpm_edger/df.plot$tpm
#ggplot()+
#  geom_point( data=df.plot, aes(x=tpm, y=ratio))+
#  facet_grid(.~cell)
#############################################
##                                         ##
#############################################


