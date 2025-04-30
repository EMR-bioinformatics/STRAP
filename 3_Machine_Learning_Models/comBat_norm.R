load("merge_R4RA.RData")
load("merge_STRAP.RData")

####Merge####

colIDs <- c("SeqID","Sample.ID","Visit","Pain","CRP","Tender","Swollen","DAS2c","Pathotype", "Treatment")

R4RA_meta_ALL2 <- R4RA_meta_ALL[ ,c("SeqID","Patient.I.D.","Visit","Arthritis.Activity","CRP","Number.of.Tender.Joints","Number.of.Swollen.Joints",
                                    "DAS2c","Pathotype","Treatment")]

colnames(R4RA_meta_ALL2) <- colIDs
R4RA_meta_ALL2$Cohort <- "R4RA"

STRAP_meta_ALL2 <- STRAP_meta[ ,c("SeqID","Patient.I.D.","Visit","Pain","CRP","Number.of.Tender.Joints","Number.of.Swollen.Joints",
                                  "DAS2c","Pathotype","Treatment")]


colnames(STRAP_meta_ALL2) <- colIDs
STRAP_meta_ALL2$Cohort <- "STRAP"

meta_ALL <- rbind(R4RA_meta_ALL2,STRAP_meta_ALL2)

meta_RNA_samples <- meta_ALL[which(is.na(meta_ALL$SeqID)==FALSE), ]

genes <- intersect(rownames(R4RA_vst_ALL),rownames(STRAP_vst_all))

R4RA_vst_sub <- R4RA_vst_ALL[which(rownames(R4RA_vst_ALL)%in%genes), ]
R4RA_vst_sub2 <- R4RA_vst_sub[genes, ]

STRAP_vst_sub <- STRAP_vst_all[which(rownames(STRAP_vst_all)%in%genes), ]
STRAP_vst_sub2 <- STRAP_vst_sub[genes, ]

vst_all <- cbind(R4RA_vst_sub2, STRAP_vst_sub2)

#all(colnames(vst_all)==meta_RNA_samples$SeqID) #TRUE

####comBat normalisation####

library(sva)

vst_all_comBat <- ComBat(vst_all, mean.only = TRUE,
                         batch = meta_RNA_samples$Cohort,
                         ref.batch = "STRAP")

save.image("merge_R4RA_STRAP.RData")
