library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(data.table)
library(stringr)

args<-commandArgs(T)



seurat_data <- Read10X(args[1],gene.column = 1)
#seurat_obj <- CreateSeuratObject(counts = seurat_data)
#filter<-seurat_obj[c("PTK7-APT","MET-APT","CD71-APT","PD-L1-APT","CD49C-APT","PTPRF-APT","CD318-APT","EPCAM-APT","ZNRF3-APT","Contr"),]


#expr<-filter[["RNA"]]@counts

#APT<-rownames(seurat_data)[which(str_detect(rownames(seurat_data),"APT_"))]
#expr<-seurat_data[APT,]
expr<-seurat_data[which(rownames(seurat_data) %in% c("PTK7_APT","Contr","split_APT")),]
expr2<-t(as.data.frame(expr))
write.table(expr2,file =paste(args[2],args[3],"_count.csv",sep=""),row.names=T,sep="\t",quote=F)

#expr<-expr%>%filter(rownames(expr)%in%c("Contr","CD318_APT","CD49C_APT","CD71_APT","EPCAM_APT","MET_APT","PD_L1_APT","PRNP_APT","PTK7_APT","PTPRF_APT","ZNRF3_APT"))
