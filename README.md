# Project 2

## Introduction
I chose to study the Hemoglobin Subunit Beta (HBB) gene, which encodes for the Hemoglobin subunit beta protein.  Hemoglobin subunit beta, along with Hemoglobin Alpha 1, forms Hemoglobin A in humans, and serves the vital function of carrying oxygen from the lungs to different tissues. Mutations in HBB can lead to genetic blood disorders including sickle cell anemia and beta thalassemias. 

![](https://upload.wikimedia.org/wikipedia/commons/thumb/1/1f/Protein_HBB_PDB_1a00.png/300px-Protein_HBB_PDB_1a00.png)
*source: https://upload.wikimedia.org/wikipedia/commons/thumb/1/1f/Protein_HBB_PDB_1a00.png/300px-Protein_HBB_PDB_1a00.png

 
For this project, I was interested in seeing how affected  different ethnic populations are by HBB mutations and the frequency at which these mutations occur. I focused on Single Nucleotide Polymorphisms (SNPs) as the mutation since deleterious SNPs are what lead to the genetic blood disorders.

![](https://upload.wikimedia.org/wikipedia/commons/a/a8/Thalassemia_beta.jpg)
*source:
https://upload.wikimedia.org/wikipedia/commons/a/a8/Thalassemia_beta.jpg

## Methods
To explore how HBB mutations affect different ethnic populations, I chose to use data that listed out all of the possible mutations in the HBB gene, as well as data that showed the allele freqency for HBB in different populations.

### Data Acquisition 
For the data relating to mutations in the HBB gene, I download a .txt file from ClinVar which lists out all of the possible variations in the HBB gene. In terms of population allele frequency, I used gnomAD to download a .csv file. The .csv shows the allele frequency for each HBB variant along with lots of other information including: allele count, allele number, and annotation.


### Data Integration
### Filtering and Calculations for Data Set 1 (Population Allle Frequency)
In order to integrate the two data sets, I used R studio and began filtering the population allele frequency data. 
I imported the file using read.csv, while at the same time assigning it the name "allele_freq". I then filtered out columns with unnecessary information using the code shown: 
`allele_freq_filter <- allele_freq[ -c(2:8,12,14,15) ] #filter out unneeded columns from allele freq`
Columns that were filtered out include: Reference (nucleotide), Alternate (nucleotide), Consequence, and Homozygote and Hemizygote count for each population. The columns I did keep in the data frame were "Position", "Annotation", "Flags",and "Allele Count" and "Allele Number" for each population. I then calculated the allele frequency for each individual population, in order to use these frequencies later as a basis for comparison. I also assigned these individual frequencies to a variable named "allele_freq_filter"
```
allele_freq_filter <- transform(allele_freq_filter, AF_Afr_Am = Allele.Count.African.African.American / Allele.Number.African.African.American)
allele_freq_filter <- transform(allele_freq_filter, AF_Lat_Am = Allele.Count.Latino.Admixed.American / Allele.Number.Latino.Admixed.American)
allele_freq_filter <- transform(allele_freq_filter, AF_Ashkenazi_Jewish = Allele.Count.Ashkenazi.Jewish / Allele.Number.Ashkenazi.Jewish)
allele_freq_filter <- transform(allele_freq_filter, AF_Eas_As = Allele.Count.East.Asian / Allele.Number.East.Asian)
allele_freq_filter <- transform(allele_freq_filter, AF_Eur_Fin = Allele.Count.European..Finnish. / Allele.Number.European..Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AF_Eur_non_Fin = Allele.Count.European..non.Finnish. / Allele.Number.European..non.Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AF_South_As = Allele.Count.South.Asian / Allele.Number.South.Asian)
```
I then changed all of the zeroes in the table into "NA" in order to better visualize which population groups had an allele frequency for each mutation. Following this, I filtered out any annotation that was not a missense variant or stop gained variant, because these were the only annotations that stemmed from deleterious SNPs.
```
allele_freq_filter <- allele_freq_filter %>% filter(Annotation == "missense_variant" | Annotation == "stop_gained") #filtered out rows that weren't equal to missense mutation or stop gained
```

### Filtering and Calculations for Data Set 2 (HBB SNP data)
For data set 2, I again imported the data, this time using read.delim2 in order to use the .txt file. I then started filtering out columns, only keeping "Conditions", "Clinical Significance", and  "GRCh37 Location". I renamed these three columns into names that were easier to use and read. For example "GRCh37 Location" became "Position" and "Clinical_Significance". Under the "Clinical_Significance" column, I then filtered out any rows of variants marked either as "benign" since these variations do not lead to the genetic blood disorders.
```
HBB_SNP_filter <- HBB_SNP_filter %>%filter(!grepl('Benign', Clinical_Significance)) #filters out rows containing "Benign"
HBB_SNP_filter <- HBB_SNP_filter %>%filter(!grepl('benign', Clinical_Significance))
```
At this point, the two sets of data were ready to be merged together into one table.

### Merging, Further Filtering/Organzing, and Visualization of Data
In order to merge the two sets of data together, I full joined them together, under the condition that values under each of their respective "Position" columns were the same. I named the new table of data "HBB_SNP_with_Allele_freq".
```
HBB_SNP_with_Allele_freq <- full_join(allele_freq_filter,HBB_SNP_filter,by=c("Position" = "Position")) #merge two tables based on matching position of SNP
View(HBB_SNP_with_Allele_freq)
```
I then filtered the data more, by removing any rows that were flagged as "lc_lof" (for low confidence) and then removing the "Flags" column after this.
```
HBB_SNP_with_Allele_freq <- HBB_SNP_with_Allele_freq %>%filter(!grepl('lc_lof', Flags)) #filters out rows that were flagged (low confidence)
HBB_SNP_with_Allele_freq <- HBB_SNP_with_Allele_freq[ -c(3) ] #filter out unneeded flag column after filtering out flagged variants
```
I then renamed all of the column names pertaining to each population because I realized that what the names mean would not be clear in a visual representation of the data. Below is an exmapled of the 3 (of 7) columns whose names were changed:
```
colnames(HBB_SNP_with_Allele_freq)[5] <- 'Allele_Frequency_African_American'
colnames(HBB_SNP_with_Allele_freq)[6] <- 'Allele_Frequency_Latino_Admixed_American'
colnames(HBB_SNP_with_Allele_freq)[7] <- 'Allele_Frequency_Ashkenazi_Jewish'
```
In order to compare the frequency with which each population was affected by mutations in the HBB gene, I summed up allele frequency of each variant for each population. I added this new data to the table by creating a new row, with the bottom value of each "Population" column representing the total frquency for all deleterious mutations.
```
HBB_SNP_with_Allele_freq[150,(5:11)] <- colSums(HBB_SNP_with_Allele_freq[,5:11], na.rm=TRUE) # added row to show sum of allele frequencies for each SNP for each population
```
Following this, I plotted the values for each population on a bar graph, in order to visualize just how prevalent mutations are in each group. (Graph included in next section)
```
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

```

## Results

The output of the data integration was succesful and I was able to see how differently each population is affected by deleterious SNPs in their HBB gene. Below is the graph created from the data merge:

![](https://i.imgur.com/WDV3EDG.png)

Based on the graph, it is evident that African Americans is the group most affected by harmful HBB mutations. Following them, the next most affected group are South Asians, with Ashkenazi Jewish showing only a slightly lower frequnecy. The least affected groups are East Asians, Finnish Europeans, and non- Finnish Europeans.


## Discussion
### Successes
There were several successes in this project. The first one was finding enough database entries for both of the data types I was searching for. Another was filtering out each data set enough to only keep the pertinent information and then successfully merging them together. Uses of different tools and commands in R were successful, many of which I looked up online while trying to figure out how to perform a specific function/action.

### Failures and Areas Needing Improvement
In terms of failures, I found that when trying to clean up each data set there were many lines of code that were repetitive. I feel that there are better ways to excute what I was trying to achieve, but I was not quite sure how to. 

There were also several areas of the project I feel that I could have improved on. One would be using more sources of data. For example, the data regarding population allele frequency only had several groups. I feel that it would have been more interesting and informative to include more populations in the data analysis. Another area which could be improved is altering the final output of data in the graph so that it can be better understood, such as changing the total allele freqency for each population into a percentage to better convey how one population is more/less affected than another. I feel that the graph could be improved as well, I changed many of the parameters and font sizes to fit everything together, but the text is a bit small still.


### Future Applications
In terms of future applications, I think this project/data could be used again by incorporating more populations affected by HBB mutations. Also, a breakdown of different genetic disease could be included in the graph in order to show how one population might be more affected by one type of disease than another. Additionally, the framework of this project could be used to study how other genetic disease affect different populations.

