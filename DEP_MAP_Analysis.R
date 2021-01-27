#Sherri Smart MERTK analysis 

library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(tidyverse)
library(patchwork)
library(plyr)

setwd("C:/Users/12298/Desktop/Data_Analytics/Sherrie_MERTK/")


Genes <- c("MERTK (10461)", "TYRO3 (7301)", "AXL (558)", "PROS1 (5627)", "GAS6 (2621)", "LGALS3 (3958)", "TUB (7275)", "TULP1 (7287)")
#Bringing in files from DEPMAP and creating character vectors -----------
#CRISPR_Score <- fread("https://ndownloader.figshare.com/files/25770029")
#Expression_Set <- fread("https://ndownloader.figshare.com/files/25797011")

Genes <- c("MERTK (10461)", "TYRO3 (7301)", "AXL (558)", "PROS1 (5627)", "GAS6 (2621)", "LGALS3 (3958)", "TUB (7275)", "TULP1 (7287)")
Cell_Line <- c("ACH-001033",	"ACH-001038",	"ACH-001192", "	ACH-001193",	"ACH-001205",	"ACH-001283",	"ACH-001430",	"ACH-001431",	"ACH-000279",	"ACH-000087",	"ACH-000041",	"ACH-000039",	"ACH-000424",	"ACH-000391",	"ACH-000499",	"ACH-000052",	"ACH-001001",	"ACH-000748",	"ACH-001715",	"ACH-001526",	"ACH-001814",	"ACH-000364",	"ACH-000082",	"ACH-000410")
EWS_Cell_Lines <- c("ACH-001033",	"ACH-001038",	"ACH-001192", "	ACH-001193",	"ACH-001205",	"ACH-001283",	"ACH-001430",	"ACH-001431",	"ACH-000279",	"ACH-000087",	"ACH-000041",	"ACH-000039",	"ACH-000424",	"ACH-000391",	"ACH-000499",	"ACH-000052")
OS_Cell_Lines <- c("ACH-001001",	"ACH-000748",	"ACH-001715",	"ACH-001526",	"ACH-001814",	"ACH-000364",	"ACH-000082",	"ACH-000410")
Cell_Name <- c("CHLA57",	"COGE352",	"SKNEP1",	"SKPNDW",	"TC32",	"TC106",	"TC138",	"TC205",	"EWS502",	"SKES1",	"RDES",	"SKNMC",	"TC71",	"MHHES1",	"EW8",	"A673",	"143B",	"SJSA1",	"CAL72",	"HUO9",	"OS252",	"U2OS",	"G292CLONEA141B1",	"SAOS2")

Cell_Matrix <- cbind(Cell_Line, Cell_Name) %>% as.data.frame() %>% `colnames<-`(c("Cell_Line", "Cell_Name"))

#Extracting relevant cell lines and then genes - Doing separately so no information loss--------------
CRISPR_Score_Ext <- CRISPR_Score %>% t() 
colnames(CRISPR_Score_Ext) = CRISPR_Score_Ext[1,]
CRISPR_Score_Ext <- subset((setDT(as.data.frame(CRISPR_Score_Ext), keep.rownames = TRUE)[]), rn %in% Genes) %>% pivot_longer(names_to = "Cell_Line", values_to = "CRISPR_Score", -rn) %>% `colnames<-`(c("Gene_Name", "Cell_Line", "CRISPR_Score"))

Expression_Set_Ext <- Expression_Set %>% t() 
colnames(Expression_Set_Ext) = Expression_Set_Ext[1,]
Expression_Set_Ext <- subset((setDT(as.data.frame(Expression_Set_Ext), keep.rownames = TRUE)[]), rn %in% Genes) %>% pivot_longer(names_to = "Cell_Line", values_to = "Expression", -rn) %>% `colnames<-`(c("Gene_Name", "Cell_Line", "Expression"))

#Saving RDS files for these two data frames - no merging due to potential loss of information
saveRDS(CRISPR_Score_Ext, "All_Cell_Lines_CRISPR.rds")
saveRDS(Expression_Set_Ext, "ALL_Cell_Lines_Expr.rds")
ALL_Cell_Lines_CRISPR <- readRDS("All_Cell_Lines_CRISPR.rds")
ALL_Cell_Lines_Expr <- readRDS("ALL_Cell_Lines_Expr.rds")
#----------------------------------------------------------------------


for (i in 1:length(Genes)) {
  
  Specific_Gene = Genes[i]

#Creating individual data frames ---------
All_Cell_Lines_CRISPR_Spec <- subset(ALL_Cell_Lines_CRISPR, Gene_Name == Specific_Gene)
ALL_Cell_Lines_Expr_Spec<- subset(ALL_Cell_Lines_Expr, Gene_Name == Specific_Gene)

#Expression
EWS_OS_Cell_Expr <- ALL_Cell_Lines_Expr_Spec[ALL_Cell_Lines_Expr_Spec$Cell_Line %in% EWS_Cell_Lines,]
EWS_OS_Cell_Expr <- merge(EWS_OS_Cell_Expr, Cell_Matrix, by="Cell_Line")
EWS_OS_Cell_Expr <- cbind(EWS_OS_Cell_Expr, rep("EWS", length(1:nrow(EWS_OS_Cell_Expr)))) %>% `colnames<-`(c("Cell_Line", "Gene_Name", "Expression", "Cell_Name", "Cell_Type"))
#EWS_OS_Cell_Expr$Cell_Type <- ifelse(EWS_OS_Cell_Expr$Cell_Line %in% EWS_Cell_Lines, "EWS_Cells", "OS_Cells")

All_Others_Exp <- ALL_Cell_Lines_Expr_Spec[!ALL_Cell_Lines_Expr_Spec$Cell_Line %in% Cell_Line,]
All_Others_Exp <- cbind(All_Others_Exp, rep("Others", length(1:nrow(All_Others_Exp)))) %>% `colnames<-`(c("Gene_Name", "Cell_Line", "Expression", "Cell_Type"))
All_Others_Exp <- All_Others_Exp[c("Cell_Line", "Gene_Name", "Expression", "Cell_Type")]

EXPR_DF <- bind_rows(EWS_OS_Cell_Expr[,-4], All_Others_Exp)


#CERES Scores
EWS_OS_CRISPR <- All_Cell_Lines_CRISPR_Spec[All_Cell_Lines_CRISPR_Spec$Cell_Line %in% EWS_Cell_Lines,]
EWS_OS_CRISPR <- merge(EWS_OS_CRISPR, Cell_Matrix, by="Cell_Line")
EWS_OS_CRISPR <- cbind(EWS_OS_CRISPR, rep("EWS", length(1:nrow(EWS_OS_CRISPR)))) %>% `colnames<-`(c("Cell_Line", "Gene_Name", "CRISPR_Score", "Cell_Name", "Cell_Type"))
#EWS_OS_CRISPR$Cell_Type <- ifelse(EWS_OS_CRISPR$Cell_Line %in% EWS_Cell_Lines, "EWS_Cells", "OS_Cells")

All_Others_CRISPR <-All_Cell_Lines_CRISPR_Spec[!All_Cell_Lines_CRISPR_Spec$Cell_Line %in% Cell_Line,]
All_Others_CRISPR <- cbind(All_Others_CRISPR, rep("Others", length(1:nrow(All_Others_CRISPR)))) %>% `colnames<-`(c("Gene_Name", "Cell_Line", "CRISPR_Score",  "Cell_Type"))
All_Others_CRISPR <- All_Others_CRISPR[c("Cell_Line", "Gene_Name", "CRISPR_Score", "Cell_Type")]

CRISPR_DF <- bind_rows(EWS_OS_CRISPR[,-4], All_Others_CRISPR)

#MERTK plots across Select cancer cell lines---------
#EWS/OS Cell Lines CERES score
p1 <- ggplot(as.data.frame(EWS_OS_CRISPR), aes(x=fct_reorder(Cell_Name, as.numeric(CRISPR_Score), .desc = TRUE), y=as.numeric(CRISPR_Score))) + geom_bar(stat="identity") + scale_y_continuous() + ylab("CERES Dependency Score") + xlab("Cell Line") + ylim(-1,1) + coord_flip() + ggtitle(paste0(stringr::str_split(Specific_Gene, pattern = " ")[[1]][1]))
#EWS/OS Cell Lines Expression
p2 <- ggplot(as.data.frame(EWS_OS_Cell_Expr), aes(x=fct_reorder(Cell_Name, as.numeric(Expression)), y=as.numeric(Expression))) + geom_bar(stat="identity") +  ylab("2Log Expression") + xlab("Cell Line") + ylim(0,10) + coord_flip() + ggtitle(paste0(stringr::str_split(Specific_Gene, pattern = " ")[[1]][1]))

#MERTK averages between 22 cell lines vs the rest - Set x axis from -1:1
#All cells - Cell type CERES score
p3 <- ggplot(CRISPR_DF, aes(x=Cell_Type, y=as.numeric(CRISPR_Score))) + scale_y_continuous() + geom_boxplot() + ylab("CERES Dependency Score") + xlab("Cell Type") + geom_hline(aes(yintercept=0.0), linetype="dashed", color="red") + ylim(-1,1) + coord_flip() + ggtitle(paste0(stringr::str_split(Specific_Gene, pattern = " ")[[1]][1]))

#All cells - Cell type Expression
p4 <- ggplot(as.data.frame(EXPR_DF), aes(x=Cell_Type, y=as.numeric(Expression))) + geom_boxplot() + ylab("2Log Expression") + xlab("Cell Type") + geom_hline(aes(yintercept=0.0), linetype="dashed", color="red") + ylim(0, 10)+ coord_flip() + ggtitle(paste0(stringr::str_split(Specific_Gene, pattern = " ")[[1]][1]))

p1+p2+p3+p4 
ggsave(paste0(Specific_Gene, ".pdf"), width = 10.5, height = 7.2)

#Mann-Whitney-wilcox Test - assumes data are not normally distributed which we cannot assume based on how few data points are present for the EWS group ------ 
wilcox_test <- (wilcox.test(as.numeric(EWS_OS_CRISPR[,3]), as.numeric(All_Others_CRISPR[,3]), paired=FALSE))
chars <- capture.output(print(wilcox_test), print(Specific_Gene))
writeLines(chars, con = file(paste0(Specific_Gene, "pValues.txt")))


}


#END---------
