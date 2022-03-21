library(coin)
library(qvalue)
library(vegan)
library(ape)
set.seed(500)

#In this script... 
#Read in tables
#Metadata group assingnments
#Abundance filtering
#Get taxonomy for filtered OTU tables
#Sum across each phylotype to get abundances
#Alpha diversity - Shannon index
#Kruskal for Shannon
#Ordination - PCoA
#Adonis/PERMANOVA
#Phylotype comparisons (with kruskal and qvalue)


#fresh start
rm(list = ls())

setwd("~/Documents/Oregon State University/Sharpton Lab/Dog Project")

##################
# Read in tables #
##################
#otu table
even <- read.table("table_even40000.biom.txt",
                   sep="\t", header = TRUE, check.names = FALSE, row.names = 1)
#metadata
# agg and non-agg in order so vegan stuff will work
meta <- read.csv("meta.3.csv", header = TRUE, check.names = FALSE, row.names = 1)

w.uni.filt <- read.table("weighted_unifrac_test.df4_biom.txt",
                         header=TRUE, sep="\t", check.names = FALSE, row.names = 1)
unw.uni.filt <- read.table("unweighted_unifrac_test.df4_biom.txt",
                           header=TRUE, sep="\t", check.names = FALSE, row.names = 1)
gsig <- read.table("group_significance_Group.txt",
                   header=TRUE, sep="\t")
##############################
# Metadata Group Assignments #
##############################


Agg.ids <- rownames(meta[which(meta[,4] == "aggressive"),])

Non.Agg.ids <- rownames(meta[which(meta[,4] == "not.aggressive"),])

All.ids <- c(Agg.ids, Non.Agg.ids)

Agg.df  <- even[, Agg.ids]
Non.Agg.df  <- even[, Non.Agg.ids]

test.df <- cbind( Agg.df, Non.Agg.df )
classes <- c( rep(0, length(Agg.ids)), rep(1, length(Non.Agg.ids)))

#############
# Filtering #
#############

#drop rows where there are no occurrances (check for previous qc)
test.df2 <- subset( test.df, apply( test.df, 1, function(x) sum(x) > 0 ))

#filter OTUs only found in a few samples
occ.thresh <- 3 #CHANGE num of samples you want to be the cut-off
test.df3 <- subset( test.df2, apply( test.df2, 1, function(x) length( x[x!=0]) > occ.thresh ))

#sum of otu count for each line
test.df3.sum <- as.data.frame((apply(test.df3, 1, function(x) sum(x))))
colnames(test.df3.sum) <- "Sum"
test.df3.sum.order <- test.df3.sum[order(test.df3.sum), , drop = FALSE]

#filter out OTUs that were observed fewer then x amount of times across all samples
freq <- 20
test.df3.sum.order.20 <- subset(test.df3.sum.order, test.df3.sum.order$Sum <= freq)
test.df3.20 <- test.df3[-which(rownames(test.df3) %in% rownames(test.df3.sum.order.20)),]
dim(test.df3.20)
test.df4 <- test.df3.20

#beta_diversity.py on test.df4 to get weighted and unweighted unifrac values(w.uni.filt and unw.uni.filt)

#######################################
# Get taxonomy for filtered OTU table #
#######################################

#match up the OTUs in the filtered OTU table with their taxonomy found in gsig
filt.tax <- gsig[which(gsig$OTU %in% rownames(test.df4)),]
otu.gsig.order <- test.df4[match(filt.tax[,1], rownames(test.df4)),] #subsets test.df4 OTUs in the order they are found in gsig 
clean.filt.tax <- filt.tax[,c(1,8)] #column 1 are OTUs, column2 is the taxonomy of those OTUs
tax.no.bracs <- as.matrix(gsub("\\[|\\]", "", clean.filt.tax$taxonomy))

#############################
# Taxonomy string splitter #
############################

#come back to string splitter and make it shorter (function in a function?)
#function for each individual unsplit group - and can then have a vector of QIIME notation or something
string_splitter  <- function(tax_col){
  #takes a one column matrix of qiime taxonomy output and splits indiviual taxonomic levels into their own columns
  species <- tax_col
  genus <- tax_col
  family <- tax_col
  order <- tax_col
  class <- tax_col
  phylum <- tax_col
  kingdom <- tax_col
  for(i in 1:length(tax_col)){
    #Split before kingdom name
    kingdom[i]     <- unlist(strsplit(kingdom[i], split='k__',   fixed=T))[2]   #Split after 'k__' in taxa string
    #Split after kingdom name
    kingdom[i]     <- unlist(strsplit(kingdom[i], split=';p__', fixed=T))[1]   #Split before ';p__' in taxa string
    #Split before phylum name
    phylum[i]     <- unlist(strsplit(phylum[i], split='p__',   fixed=T))[2]   #Split after 'p__' in taxa string
    #Split after phylum name
    phylum[i]     <- unlist(strsplit(phylum[i], split=';c__', fixed=T))[1]   #Split before ';c__' in taxa string
    #Split before class name
    class[i]     <- unlist(strsplit(class[i], split='c__',   fixed=T))[2]   #Split after 'c__' in taxa string
    #Split after class name
    class[i]     <- unlist(strsplit(class[i], split=';o__', fixed=T))[1]   #Split before ';f__' in taxa string
    #Split before class name
    order[i]     <- unlist(strsplit(order[i], split='o__',   fixed=T))[2]   #Split after 'o__' in taxa string
    #Split after class name
    order[i]     <- unlist(strsplit(order[i], split=';f__', fixed=T))[1]   #Split before ';f__' in taxa string
    #Split before family name
    family[i]     <- unlist(strsplit(family[i], split='f__',   fixed=T))[2]   #Split after 'f__' in taxa string
    #Split after family name
    family[i]     <- unlist(strsplit(family[i], split=';g__', fixed=T))[1]   #Split before ';g__' in taxa string
    #Split before genus name
    genus[i]     <- unlist(strsplit(genus[i], split='g__',   fixed=T))[2]   #Split after 'g__' in taxa string
    #Split after genus name
    genus[i]     <- unlist(strsplit(genus[i], split=';s__', fixed=T))[1]   #Split before ';s__' in taxa string
    #Split before species name
    species[i]     <- unlist(strsplit(genus[i], split='s__',   fixed=T))[2]   #Split after 'g__' in taxa string
    
  }
  all_tax <- cbind(kingdom, phylum, class, order, family, genus, species)
  colnames(all_tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  return(all_tax)
}
split.tax <- string_splitter(tax.no.bracs)
split.tax[which(split.tax == "")] <- NA #setting all the blanks to NA
which(split.tax == "") # testing to see if the blanks are gone

otu.tax.df <- cbind(otu.gsig.order, split.tax) #matches the OTU abundances with the split taxonmy 

###################
# Phylotype Sums #
##################

phylo_sums <- function (otu.gsig.order, otu.df.tax, phylo.name) {
  #sums the OTU counts for each phylotype
  #
  #otu.df.tax = data frame with including OTU counts and taxonomy; rownames are OTU names, colnames are sample IDS until
  #you get to taxonomies and then the last colnames are phylotype levels e.g. "Kingdom" to "Species"
  #phylo.names = string to specify the phylotpe level you wanted summed e.g. "Phyla","Class",..."Species"
  unq.phyla <- NULL
  level.index <- which(colnames(otu.df.tax) == phylo.name)
  unq.phyla <- as.matrix(unique(otu.df.tax[,level.index]))
  unq.phyla <- as.matrix(unq.phyla[-which(unq.phyla %in% NA)])
  print(unq.phyla)
  all.taxa.sums <- NULL
  taxa.sums <- NULL
  for (i in 1:nrow(unq.phyla)) {
    toi.index <- which(otu.df.tax[level.index] == unq.phyla[i,])
    toi <- otu.gsig.order[toi.index,]
    taxa.sums <- as.matrix(colSums(toi))
    colnames(taxa.sums) <- unq.phyla[i,]
    all.taxa.sums <- cbind(all.taxa.sums, taxa.sums)
  }
  return(t(all.taxa.sums))
}

phyla.df <- phylo_sums(otu.gsig.order, otu.tax.df, "Phylum")
class.df <- phylo_sums(otu.gsig.order, otu.tax.df, "Class")
order.df <- phylo_sums(otu.gsig.order, otu.tax.df, "Order")
family.df <- phylo_sums(otu.gsig.order, otu.tax.df, "Family")
genus.df <- phylo_sums(otu.gsig.order, otu.tax.df, "Genus")





###################################
# Alpha diversity - Shannon Index #
###################################
#######################
# Kruskal for Shannon #
#######################
get_kruskal_pvalue<- function( dataframe.row, classes){
  df = data.frame(abundance=dataframe.row, feature = classes )
  p <- pvalue( 
    kruskal_test(
      abundance~factor(feature), 
      data = df,
      distribution= approximate(B = 5000)
    )
  )
  return(p)
}

main.df <- test.df4
ids.1 <- Agg.ids
ids.2 <- Non.Agg.ids
group.1 <- "Aggressive"
group.2 <- "Non.Aggressive"

#get_shannon function will separate out the two groups of interest, calc their shannon indeces, compare the shannon values from
#the 2 groups with kruskall-wallis (and give you a p-value), and output boxplots so you can visually compare the alpha diversity
get_shannon <- function(main.df, ids.1, ids.2, group.1, group.2){
  #sample IDs should be the colnames
  
  classes <- c( rep(group.1, length(ids.1)), rep(group.2, length(ids.2)) )
  #classes <- c( rep("Aggressive", length(Agg.ids)), rep("Non-Aggressive", length(Non.Agg.ids)) )
  
  #df1 <- test.df4[,Agg.ids]
  #df2 <- test.df4[,Non.Agg.ids]
  
  df1 <- main.df[,ids.1]
  df2 <- main.df[,ids.2]
  
  div.1 <- diversity(t(df1), index = "shannon", MARGIN = 1, base = exp(1))
  div.1 <- as.data.frame(div.1)
  #div.1 <- df1[order(-div.1$div.1), , drop = FALSE]
  div.1 <- as.data.frame(div.1)
  colnames(div.1) <- "shannon"
  
  div.2 <- diversity(t(df2), index = "shannon", MARGIN = 1, base = exp(1))
  div.2 <- as.data.frame(div.2)
  #div.2 <- df2[order(-div.2$div.2), , drop = FALSE]
  div.2 <- as.data.frame(div.2)
  colnames(div.2) <- "shannon"
  
  all.the.divs <- rbind(div.1, div.2)
  
  print(all.the.divs)
  
  rawp <-NULL
  rawp <- apply( t(all.the.divs), 1, get_kruskal_pvalue, classes ) #using get_kruskal_pvalue function for shannon
  rawp.df <- as.data.frame(rawp)
  rawp.sort <- as.data.frame((rawp.df[order(rawp.df), , drop =FALSE]))
  
  cat("P-Value is: ", rawp.sort[1,])
  
  #plotting by samples in exp vs control
  div.com <- cbind(all.the.divs, classes)
  boxplot(shannon ~ classes, data = div.com, col = (c("red", "purple")), ylab = "Shannon Index")
  
  return(rawp.sort[1,])
  
}

#get phyla.df, phyla.class, etc from the Dog_OTU_Phylo_Count_Sums.R script
shan.otu <- get_shannon(test.df4, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")
shan.phyla <- get_shannon(phyla.df, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")
shan.class <- get_shannon(class.df, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")
shan.order <- get_shannon(order.df, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")
shan.family <- get_shannon(family.df, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")
shan.genus <- get_shannon(genus.df, Agg.ids, Non.Agg.ids, "Aggressive", "Non-Aggressive")

#####################
# Ordination - PCoA #
#####################
classes <- c( rep("Aggressive", length(Agg.ids)), rep("Not_Aggressive", length(Non.Agg.ids)))
df <- test.df4
dist.meth <- "bray"   
dist <- vegdist(t(df), method=dist.meth)
k=2
dims <- c(1,2)
groups <- classes
names(groups) <- colnames(df)
ellp.kind <- "se"
#PCOA
#choose unweighted UniFrac, weighted UniFrac, or Bray-Curtis only to proceed
object <- cmdscale(w.uni.filt, k=k, eig = TRUE, add = TRUE) #use this for unifrac
#object <- cmdscale(unw.uni.filt, k=k, eig = TRUE, add = TRUE) #use this for unifrac
#object <- cmdscale(dist, k=k) #Bray-Curtis

#NMDS
#object <- metaMDS( t(df), distance = dist.meth, k=k)

#Finding the % variation principal coordinate axes account for
eigsum <- sum(object$eig)
object$eig[1]/eigsum #pco 1 - 46%
object$eig[2]/eigsum #pco 2 - 12%

#build the plot
##aggression
ids.1 <- Agg.ids
ids.2 <- Non.Agg.ids
mds.fig <- ordiplot(object, display="sites", type = "none", choices=dims, xlab = "PC1  (46%)", ylab = "PC2  (12%)" )
points(mds.fig, "sites", pch = 19, col = "green4", select = ids.1)
points(mds.fig, "sites", pch = 19, col = "purple3", select = ids.2 )
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("green4"), lwd = 2, show.groups = "Aggressive")
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("purple3"), lwd = 2, show.groups = "Not_Aggressive")
#ordispider(object, groups, label = TRUE, choices=dims )
levels(classes) <- c("Aggressive","Not Aggressive")
legend("bottomleft", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("green4", "purple3"), y.intersp = 2)


ef <- envfit(object, as.data.frame(meta[,4], perm = 9999)) #aggresive vs. non
#ef <- envfit( object, (meta[,c(2,4)]), perm = 4999)
ef
####################
# Adonis/PERMANOVA #
####################
#removing unknown from ages in metadata
unk <- meta[which(meta$Age == "Unknown"),]
no.unk <- meta[-which(meta$Age == "Unknown"),]
w.uni.filt.unk <- w.uni.filt[-which(rownames(w.uni.filt) %in% rownames(unk)),]
w.uni.filt.unk <- w.uni.filt.unk[,-which(colnames(w.uni.filt.unk) %in% rownames(unk))]
unw.uni.filt.unk <- unw.uni.filt[-which(rownames(unw.uni.filt) %in% rownames(unk)),]
unw.uni.filt.unk <- unw.uni.filt.unk[,-which(colnames(unw.uni.filt.unk) %in% rownames(unk))]
no.unk <- no.unk[match(rownames(unw.uni.filt.unk), rownames(no.unk)),]

set.seed(123)
Agg <- factor(meta$Dog.Aggression)
adonis(w.uni.filt ~ Sex, data = meta[match(rownames(w.uni.filt), rownames(meta)),], permutations = 9999)
adonis(unw.uni.filt ~ Sex, data = meta[match(rownames(w.uni.filt), rownames(meta)),], permutations = 9999)

adonis(w.uni.filt.unk ~ Age, data = no.unk, permutations = 9999)
adonis(unw.uni.filt.unk ~ Age, data = no.unk, permutations = 9999)

adonis(dist ~ Dog.Aggression, data = meta, permutations = 9999)#use bray

#########################
# Phylotype Comparisons #
#########################

krusk_phylo_qval <- function(group.1, group.2, phylo.df){
  #function that will prepare bacterial phylotypes of two groups for the get_kruskal_pvalue function
  #and then mutiple value correct with qvalue and then rearrange those results so the output is easy to read
  #group.1 = the sample ids for the first group you are interested in comparing
  #group.2 = the samples ids for the second group you are interested in comparing
  #phylo.df = the data frame of phylotypes you want to run the use kruskall-wallis tests on
  df.1      <- NULL
  df.2      <- NULL
  rawp      <- NULL
  rawp.sort <- NULL
  qval      <- NULL
  qvalp     <- NULL
  df.1      <- phylo.df[, group.1]
  df.2      <- phylo.df[, group.2]
  com.df    <- cbind(df.1, df.2)
  classes   <- c(rep(0, length(group.1)), rep(1, length(group.2)))
  rawp      <- apply(phylo.df, 1, get_kruskal_pvalue, classes)
  rawp.df   <- as.data.frame(rawp)
  rawp.sort <- as.data.frame((rawp.df[order(rawp.df), , drop =FALSE]))
  print(rawp.sort)
  qval      <- qvalue(rawp.sort$rawp)
  hist(rawp)
  qvalp     <- cbind(rawp.sort, qval$qvalues)
  print(qvalp)
  return(qvalp)
}

phyla.an.q  <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, phyla.df)
class.an.q  <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, class.df)
order.an.q  <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, order.df)
family.an.q  <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, family.df)
genera.an.q  <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, genus.df)

otu.an.q <- krusk_phylo_qval(Agg.ids, Non.Agg.ids, test.df4)


