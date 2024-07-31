library("dplyr")

args<-commandArgs(T)




oligosumm<-data.table::fread(paste(args[1],"/02.APT/count/Oligo_APT_ALL_barcode_count.csv",sep=""))

#barcode<-data.table::fread(paste(args[1],"/01.scRNA/",args[2],"/02.count/",args[2],"_combined_list.txt",sep=""),header=F)
barcode<-data.table::fread(paste(args[1],"/01.scRNA/",args[2],"/02.count/","barcodeTranslate.txt",sep=""),header=F)

colnames(barcode)<-c("V1","cell")


APT<-inner_join(barcode,oligosumm,by="V1")

APT$type<-"oligo"

cdna<-data.table::fread(paste(args[1],"/02.APT/count/CDNA_APT_ALL_barcode_count.csv",sep=""))


APTcdna<-inner_join(barcode,cdna,by="V1")


APTcdna$type<-"cdna"


#colnam<-colnames(APT)[which(stringr::str_detect(colnames(APT),"APT"))]
colnam<-c("PTK7_APT", "Contr", "split_APT")

col<-c("V1","cell","type",colnam)

#col<-c("V1","cell","type","PTK7-APT","MET-APT","CD71-APT","PD-L1-APT","CD49C-APT","PTPRF-APT","CD318-APT","EPCAM-APT","ZNRF3-APT","Contr")

t<-setdiff(col,colnames(APT))

df <- data.frame()

for (i in 1:nrow(APT)) {
    vectr <- c(rep(0,length(t)))
    df <- rbind(df,vectr) # fill data
}

colnames(df)<-t

APT<-cbind(APT,df)

APT<-as.data.frame(APT)

id<-match(col,colnames(APT))

APT<-APT[,id]


t<-setdiff(col,colnames(APTcdna))

df <- data.frame()

for (i in 1:nrow(APTcdna)) {
    vectr <- c(rep(0,length(t)))
    df <- rbind(df,vectr) # fill data
}

colnames(df)<-t

APTcdna<-cbind(APTcdna,df)

APTcdna<-as.data.frame(APTcdna)

id<-match(col,colnames(APTcdna))

APTcdna<-APTcdna[,id]





allAPT<-rbind(APT,APTcdna)


#all<-aggregate(cbind(`PTK7-APT`,`MET-APT`,`CD71-APT`,`PD-L1-APT`,`CD49C-APT`,`PTPRF-APT`,`CD318-APT`,`EPCAM-APT`,`ZNRF3-APT`,`Contr`)~cell,allAPT,sum)

validCell <- as.data.frame(data.table::fread(paste(args[1],"/01.scRNA/",args[2],"/02.count/","beads_barcodes.txt",sep=""),header=F))

allAPT2<-allAPT[,-c(1,3)]
allAPT2 <- allAPT2[which(allAPT2$cell %in% validCell$V1),]

all<-aggregate(.~cell,allAPT2,sum)

write.csv(all,file=paste(args[1],"/02.APT/count/ALL_APT_cell.csv",sep=""),quote=F,row.names=F)


mean<-round(colMeans(all[,-1]),0)

sum<-colSums(all[,-1])

cellnumber<-length(unique(all$cell))

count<-cbind(mean,sum)
count<-as.data.frame(count)
count$APT<-rownames(count)
count$cellnumber<-cellnumber

count$sample<-args[2]
count<-count[,c(5,3,4,2,1)]

write.table(count,file=paste(args[1],"/02.APT/count/Count_mean.txt",sep=""),sep="\t",row.names=F,quote=F)
