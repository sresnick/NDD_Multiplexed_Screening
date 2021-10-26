library(dplyr)
library(Rsamtools)
library(reshape2)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
library(metap)
library(matrixStats)
library(stringr)

#######################################
####PRELIMINARY BAM FILE PROCESSING####
#######################################
library(plyr)
bam <- scanBam() #file path
bam <- as.data.frame(bam)
bam.0 <- filter(bam, bam$flag == 0)
bam.0numdup <- ddply(bam.0, .(qname,rname), nrow)
colnames(bam.0numdup) <- c("well", "barcode", "count")
melted <- melt(bam.0numdup, id.vars = c("well", "barcode"))
casted <- dcast(melted, barcode~well)
casted[is.na(casted)] <- 0
rownames(casted) <- casted$barcode
casted[1] <- NULL
write.csv(casted, file = "",quote = F, row.names = T) # file path
rm(bam, bam.0, bam.0numdup)

####################################################
####READING IN CASTED TABLE - INITIAL FILTERING ####
####################################################
NBR <- 12#Number of bioreps to help with column spacing and division
casted <- read.csv("", header = T, row.names = 1) # file path to casted table you just output
directory <- "" #an id for a sub directory that you want to output files to eg. "Plate1_"

# Define the column names of controls # do this manually, example columns used in mauscript
#control.match.chap.6st <- c("Chap1_6H", "Chap1_7H", "Chap1_8H", "Chap2_8H", "Chap2_9A", "Chap2_9B", "Chap1_9H", "Chap1_10A", "Chap1_10B", "Chap2_8B", "Chap2_8C", "Chap2_8D", "Chap1_3H", "Chap1_4H", "Chap1_5H", "Chap2_5H", "Chap2_6H", "Chap2_7H")
control.match <- "control.match.chap.6st"

# filter out columns that have low number of reads
a <- colSums(casted) < 15000  # Change number as appropriate ~25% of what you hope to get is a good cut off
low.col.read.match <- vector(mode = "character", length = sum(a))
for (i in 1:sum(a)) { 
  #low.col.read.match[i] <- strsplit(colnames(casted)[a][i], "[.]")[[1]][1]
  low.col.read.match[i] <- colnames(casted)[a][i]
}
casted_filt <- casted[, !grepl(paste(low.col.read.match, collapse = "|"), colnames(casted))]
# reorganize casted table so that controls are at the front of the matrix
control.col.loc <- which(grepl(paste(control.match, collapse="|"), colnames(casted_filt)))
condition.col.loc <- which(!grepl(paste(control.match, collapse="|"), colnames(casted_filt)))
casted_filt_counts <- casted_filt[, c(control.col.loc,condition.col.loc) ]

#relative abundances
rel_abundance <- as.data.frame(scale(casted_filt_counts, center=F, scale = colSums(casted_filt_counts)))

# counts per million
cpm <- rel_abundance*1000000

#get od table
od <- read.csv(sprintf("~/%s_od.csv", directory), header = T, row.names = 1) # diretory to file. format is a header row with all names and a row below with all values

#subsetting od table to only include OD measurements of tests that survived filtering
#1 filter OD so that it only includes names from the post filtered cpm_averaged
#might need to change back to cpm_averaged depending on where you are
#column_list <- colnames(cpm)
# relevant_ods <- od[, colnames(od) %in% column_list] # use if you have master list

#reordering od columns names to match column names of cpm_averaged
relevant_ods_reordered <- relevant_ods[names(cpm)]
#correcting for OD measurements by multiplying ODs by the cpm average value
cpm_averaged_od_corrected <- sweep(cpm, 2, as.numeric(relevant_ods_reordered[1,]), "*")


#remove bio-rep information (remove ".B1" or ".B2") from column names
str_sub(colnames(cpm_averaged_od_corrected), -3) <- ""

#averaging bio reps
if (NBR > 1) {
  cpm_averaged_od_corrected <- as.data.frame(sapply(unique(names(cpm_averaged_od_corrected)), function(col) rowMeans(cpm_averaged_od_corrected[names(cpm_averaged_od_corrected) == col])))
} else {
  cpm_averaged_od_corrected <- cpm_averaged_od_corrected
}

# REmove other identifiers from plate name
str_sub(colnames(cpm_averaged_od_corrected), 1,7) <- ""

# getting gene names and storing them in dataframe genes in column Model
genes <- as.data.frame(word(rownames(cpm_averaged_od_corrected),1,sep = "_"))
colnames(genes) <- "Model"
genes$Model <- as.character(genes$Model)

# defining numbers of controls and tests
num_controls <- length(control.col.loc)/NBR
num_rescuers <- ncol(cpm_averaged_od_corrected)- num_controls

#make control only df
control_casted <- cpm_averaged_od_corrected[,1:num_controls]

###### write out input file #####
write.table(cpm_averaged_od_corrected, sprintf("~/%s.csv", directory), quote = F, row.names = T, col.names = NA, sep = ",")

##########################################################
#### Mean-Variance Modeling and individual BC testing ####
##########################################################

all_mean_var <- NULL
for (i in 1:num_rescuers) {
  chap <- colnames(cpm_averaged_od_corrected[num_controls+i])
  test_column <- which(grepl(chap, colnames(cpm_averaged_od_corrected)))
  test <- cpm_averaged_od_corrected[,test_column]
  tbl <- cbind(control_casted,test)
  colnames(tbl)[num_controls+1] <- chap
  tbl$gene <- genes$Model
  tbl[tbl == 0] <- NA
  
  tbl <- tbl[rowSums(is.na(tbl[, 1:num_controls+1]))<3,]
  tbl[is.na(tbl)] <- 0
  genes_filtered <- as.data.frame(tbl$gene)
  colnames(genes_filtered) <- "Model"
  tbl$gene <- NULL
  
  controls <- tbl[,1:num_controls]
  mean_var <- data.frame(matrix(0, nrow = nrow(tbl), ncol = 2))
  colnames(mean_var) <- c("Mean", "Variance")
  rownames(mean_var) <- rownames(controls)
  mean_var$Mean <- rowMeans(controls)
  mean_var$Variance <- rowVars(as.matrix(controls))
  mean_var$testMean <- tbl[,num_controls+1]
  
  ### mageck style variance modeling where variance = mean + k(mean)^b 
  regres <- glm(log(Variance-Mean) ~ log(Mean), data = mean_var)
  adj_int <- (exp(1))^(regres$coefficients[1])
  mean_var$mag_adj_var <- (mean_var$Mean) + (adj_int*((mean_var$Mean)^regres$coefficients[2]))
  
  mean_var$mag_phigh <- apply(mean_var, 1, function(x)  pnorm(((x[3]-x[1])/sqrt(x[7])), lower.tail = FALSE))
  mean_var$mag_plow <- apply(mean_var, 1, function(x) pnorm(((x[3]-x[1])/sqrt(x[7])), lower.tail = TRUE))
  mean_var$mag_ptwosided <- apply(mean_var, 1, function(x) pnorm(-abs((x[3]-x[1])/(sqrt(x[7]))))*2)
  
  ### add necessary stuff to mean var table
  mean_var$Barcode <- rownames(controls)
  mean_var$Gene <- genes_filtered$Model
  mean_var$Well <- chap
  mean_var$Log2FC <- log2(mean_var$testMean/mean_var$Mean)
  ## do some pvalue adjustment for comparisons within each well
  mean_var$mag_phigh_fdr <- p.adjust(mean_var$mag_phigh, method = "BH")
  mean_var$mag_plow_fdr <- p.adjust(mean_var$mag_plow, method = "BH")
  all_mean_var <- rbind(all_mean_var, mean_var)
  
}

##########################################
#### Mean-Variance Metrics for Output ####
##########################################
###### computing mean, variance, adjusted variance, and coefficient of variation for each barcode in the controls
mean_var_metrics <-  data.frame(matrix(0, nrow = nrow(control_casted), ncol = 2))
colnames(mean_var_metrics) <- c("Mean", "Variance")
rownames(mean_var_metrics) <- rownames(control_casted)
mean_var_metrics$Mean <- rowMeans(control_casted)
mean_var_metrics$Variance <- rowVars(as.matrix(control_casted))
mean_var_metrics$Std_dev <- sqrt(mean_var_metrics$Variance)
mean_var_metrics$coef_of_variation <- mean_var_metrics$Std_dev/mean_var_metrics$Mean
mean_var_metrics$mageck_adjusted_variance <- (mean_var_metrics$Mean) + (adj_int*((mean_var_metrics$Mean)^regres$coefficients[2]))
write.table(mean_var_metrics, sprintf("~/%s.csv", directory), quote = F, row.names = T, col.names = NA, sep = ",")

#####################################################################
#### FDR corrections and Data correction for P-value combination ####
#####################################################################

#converting to numeric and dealing with p values = 1 or 0 since sumz doesnt like them
all_mean_var$mag_phigh_fdr <- as.numeric(all_mean_var$mag_phigh_fdr)
all_mean_var$mag_phigh_fdr[all_mean_var$mag_phigh_fdr == 0] <- 1e-200
all_mean_var$mag_phigh_fdr[all_mean_var$mag_phigh_fdr == 1] <- 0.99999999999999
all_mean_var$mag_plow_fdr <- as.numeric(all_mean_var$mag_plow_fdr)
all_mean_var$mag_plow_fdr[all_mean_var$mag_plow_fdr == 0] <- 1e-200
all_mean_var$mag_plow_fdr[all_mean_var$mag_plow_fdr == 1] <- 0.99999999999999

## write out all_mean_var before fitlering out anything
write.table(all_mean_var, sprintf("/~/barcode_tests/%s.csv", directory), quote = F, row.names = T, col.names = NA, sep = ",")

# adding a count for number of barcoded strains in each test
# filtering out models with less than 3 barcodes per strain
all_mean_var <- add_count(all_mean_var, Gene, Well)
all_mean_var <- filter(all_mean_var, n > 2)

sumtable <- all_mean_var %>% group_by(Well,Gene) %>% dplyr::summarize(avglfc = mean(Log2FC), num_bc = n(), num_sig_mag_phigh_raw_bc = sum(mag_phigh < 0.1), num_sig_mag_plow_raw_bc = sum(mag_plow < 0.1))
sumtable$stouf_phigh_mag <- 0
sumtable$stouf_plow_mag <- 0

#combining p value information
for (j in 1:nrow(sumtable)){
  chap <- sumtable$Well[j]
  geneid <- sumtable$Gene[j]
  
  blahsubset <- all_mean_var %>% ungroup() %>% filter(Well == chap, Gene == geneid) 
  # individual test adjusted - mageck modeling
  sumtable$stouf_phigh_mag[j] <- sumz(blahsubset$mag_phigh_fdr)$p[1,1]
  sumtable$stouf_plow_mag[j] <- sumz(blahsubset$mag_plow_fdr)$p[1,1]
  
}
  
sumtable_adjusted <- sumtable
sumtable_adjusted <- sumtable_adjusted %>% mutate(stouf_phigh_mag_fdr = p.adjust(stouf_phigh_mag, method = "BH"))
sumtable_adjusted <- sumtable_adjusted %>% mutate(stouf_plow_mag_fdr = p.adjust(stouf_plow_mag, method = "BH"))
  
# for plotting
sumtable_adjusted$stouf_phigh_mag_fdr_CLASS <- cut(sumtable_adjusted$stouf_phigh_mag_fdr,  breaks = c(-1,1E-60,1E-45,1E-30,1E-15,0.05,1), labels = c("[0]-[1E^-60]","[1E^-60]-[1E^-45]","[1E^-45]-[1E^-30]","[1E^-30]-[1E^-15]", "[1E^-15]-[0.05]", "[>0.05]"))
sumtable_adjusted$stouf_phigh_mag_fdr_CLASS <- as.factor(sumtable_adjusted$stouf_phigh_mag_fdr_CLASS)
sumtable_adjusted$stouf_plow_mag_fdr_CLASS <- cut(sumtable_adjusted$stouf_plow_mag_fdr,  breaks = c(-1,1E-60,1E-45,1E-30,1E-15,0.05,1), labels = c("[0]-[1E^-60]","[1E^-60]-[1E^-45]","[1E^-45]-[1E^-30]","[1E^-30]-[1E^-15]", "[1E^-15]-[0.05]", "[>0.05]"))
sumtable_adjusted$stouf_plow_mag_fdr_CLASS <- as.factor(sumtable_adjusted$stouf_plow_mag_fdr_CLASS)

#bring in name index
name_index <- read.csv("~/chaperone_key.csv", header = F)
colnames(name_index) <-c("Well", "Chaperone")
sumtable_adjusted_chap <- left_join(sumtable_adjusted, name_index, by="Well")

# write out sumtable_adjusted_chap!!
write.table(sumtable_adjusted_chap, sprintf("~/output_tables/%s.csv", directory), quote = F, row.names = T, col.names = NA, sep = ",")

######################################
############### Plotting #############
######################################
st <- ggplot(sumtable_adjusted_chap, aes(Gene, Well))
#st <- ggplot(sumtable_adjusted_chap, aes(Gene, Chaperone))

st + geom_raster(aes(fill=stouf_phigh_mag_fdr_CLASS) ) + theme(axis.text.x = element_text(angle=90, hjust=1) ) + scale_fill_manual(values = c("gray1","red4","orangered" ,"orange","gold", "white"), drop = F) + geom_vline(xintercept=seq(0.5,90.5,by=1), color="lightgray") + geom_hline(yintercept=seq(0.5,185.5,by=1), color="lightgray") + labs(fill = "FDR") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Model") + ggtitle("Rescuers")
st + geom_raster(aes(fill=avglfc) ) + scale_fill_gradient2(na.value = "white", limits = c(0,2), breaks = c(-2,-1,0,1,2), oob = scales::squish) +  theme(axis.text.x = element_text(angle=90, hjust=1) ) + geom_vline(xintercept=seq(0.5,90.5,by=1), color="lightgray") + geom_hline(yintercept=seq(0.5,185.5,by=1), color="lightgray") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Model") + labs(fill = "Log2 Fold Change")

















