setwd('~/Desktop')
allele_freq <- read.csv(file = 'pop_allele_freq_gnomAD.csv') #import allele freq csv into R
View(allele_freq)
allele_freq_filter <- allele_freq[ -c(2:8,12,14,15) ] #filter out unneeded columns from allele freq
allele_freq_filter <- transform(allele_freq_filter, AF_Afr_Am = Allele.Count.African.African.American / Allele.Number.African.African.American)
allele_freq_filter <- transform(allele_freq_filter, AF_Lat_Am = Allele.Count.Latino.Admixed.American / Allele.Number.Latino.Admixed.American)
allele_freq_filter <- transform(allele_freq_filter, AF_Ashkenazi_Jewish = Allele.Count.Ashkenazi.Jewish / Allele.Number.Ashkenazi.Jewish)
allele_freq_filter <- transform(allele_freq_filter, AF_Eas_As = Allele.Count.East.Asian / Allele.Number.East.Asian)
allele_freq_filter <- transform(allele_freq_filter, AF_Eur_Fin = Allele.Count.European..Finnish. / Allele.Number.European..Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AF_Eur_non_Fin = Allele.Count.European..non.Finnish. / Allele.Number.European..non.Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AF_South_As = Allele.Count.South.Asian / Allele.Number.South.Asian)
# above are calculations for allele frequencies for each population
allele_freq_filter <- allele_freq_filter[ -c(6:37) ] #filter out unneeded columns from allele freq
allele_freq_filter[allele_freq_filter==0] <- NA #changed 0's to NA in order to visually see better
library(dplyr) #load dplyr package
allele_freq_filter <- allele_freq_filter %>% filter(Annotation == "missense_variant" | Annotation == "stop_gained") #filtered out rows that weren't equal to missense mutation or stop gained
HBB_SNPs <- read.delim2('HBB_SNP_clinvar_tabular.txt') #import HBB SNPs txt file into R
View(HBB_SNPs)
HBB_SNP_filter <- HBB_SNPs[ -c(1:3,6:8,10:15) ] #filter out unneeded columns from HBB SNPs
View(HBB_SNP_filter)
colnames(HBB_SNP_filter)[3] <- 'Position' #change col name to match allele_freq_filter
colnames(HBB_SNP_filter)[1] <- 'Condition(s)' #change col name 
colnames(HBB_SNP_filter)[2] <- 'Clinical_Significance' #change col name 
HBB_SNP_filter <- HBB_SNP_filter %>%filter(!grepl('Benign', Clinical_Significance)) #filters out rows containing "Benign"
HBB_SNP_filter <- HBB_SNP_filter %>%filter(!grepl('benign', Clinical_Significance))
HBB_SNP_filter <- HBB_SNP_filter[ -c(4) ] #filter out unneeded columns from HBB SNPs
HBB_SNP_with_Allele_freq <- full_join(allele_freq_filter,HBB_SNP_filter,by=c("Position" = "Position")) #merge two tables based on matching position of SNP
View(HBB_SNP_with_Allele_freq)
HBB_SNP_with_Allele_freq <- HBB_SNP_with_Allele_freq %>%filter(!grepl('lc_lof', Flags)) #filters out rows that were flagged (low confidence)
HBB_SNP_with_Allele_freq <- HBB_SNP_with_Allele_freq[ -c(3) ] #filter out unneeded flag column after filtering out flagged variants
HBB_SNP_with_Allele_freq <- HBB_SNP_with_Allele_freq %>% filter(Annotation == "missense_variant" | Annotation == "stop_gained") #filtered out rows that weren't equal to missense mutation or stop gained
colnames(HBB_SNP_with_Allele_freq)[5] <- 'Allele_Frequency_African_American'
colnames(HBB_SNP_with_Allele_freq)[6] <- 'Allele_Frequency_Latino_Admixed_American'
colnames(HBB_SNP_with_Allele_freq)[7] <- 'Allele_Frequency_Ashkenazi_Jewish'
colnames(HBB_SNP_with_Allele_freq)[8] <- 'Allele_Frequency_East_Asian'
colnames(HBB_SNP_with_Allele_freq)[9] <- 'Allele_Frequency_European_Finnish'
colnames(HBB_SNP_with_Allele_freq)[10] <- 'Allele_Frequency_European_non_Finnish'
colnames(HBB_SNP_with_Allele_freq)[11] <- 'Allele_Frequency_South_Asian'
HBB_SNP_with_Allele_freq[150,5:11] <- colSums(HBB_SNP_with_Allele_freq[,5:11], na.rm=TRUE) # added row to show sum of allele frequencies for each SNP for each population
Total_freq <- c( 3.271881e-01, 1.704877e-02, 3.569365e-02, 4.461999e-03, 7.476229e-03 ,4.585113e-03 ,4.416007e-02)
par(mar=c(5, 5 , 3, 5))
barplot(Total_freq,
      main = "Total Allele Frequencies for Deleteroius SNPs in the HBB Gene for a Given Population ",
      xlab = "Populations",
      cex.axis = 0.8,
      cex.lab = 0.8,
      cex.names = 0.8,
      cex.main = 0.8,
      ylab = "Total Allele Frequencies",
      names.arg = c("African_American", "Latino Admixed American", "Ashkenazi Jewish", "East Asian", "European Finnish", "European non Finnish", "South Asian"),
      col = "blue",)

HBB_SNP_with_Allele_freq[151,4:11] <- colSums(HBB_SNP_with_Allele_freq[,4:11], na.rm=TRUE) # added row to show sum of allele frequencies for each SNP for each population

