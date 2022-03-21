library(qvalue)
library(ggplot2)
library(ape)
library(phylobase)
library(phytools)
library(gplots)
library(ggtree)
library(metacoder)

setwd("/Users/sharptot/Dropbox/sharpton/projects/dogs")

data.dir <- "kirchoff_2017/rarefied/"

p.test <- read.table(file=paste(data.dir, "p_test_results.tab_stats.txt",sep=""),
                     stringsAsFactors = FALSE)
tax    <- read.table(file=paste(data.dir, "dog_claatu_tax_dict.tab",sep=""),
                     stringsAsFactors = FALSE)
counts <- read.table( file=paste(data.dir, "claatu_counts_rarefied.tab",sep=""),
                      stringsAsFactors = FALSE, 
                      header = TRUE,
                      row.names = 1)
t.counts <- t(counts)
map    <- read.table( file=paste(data.dir, "mapping_file.txt",sep=""),
                      stringsAsFactors = FALSE, 
                      header = TRUE,
                      row.names = 1)
otus   <- read.table( file = paste(data.dir, "table_even40000.biom.txt",sep=""),
                      stringsAsFactors = FALSE, 
                      header = TRUE,
                      row.names = 1)
otu.tax <- read.table( file=paste(data.dir, "rep_set_tax_assignments.txt",sep=""),
                       stringsAsFactors = FALSE,
                       sep = "\t")
tree   <- read.tree( file = paste(data.dir, "new_prepped_tree.tre",sep=""))
#write.csv(map, "Supplemental_Table_1.csv", row.names=FALSE)

## this takes a minute....
tree.2 <- midpoint.root( tree )                    
# are there significantly conserved clades?
## hist the pvalues
hist( p.test[,6])
q <- qvalue( p.test[,6] )
hist( q$qvalues)
## we have a mass over small pvalues
idx    <- which(q$qvalues < 0.05)
nnames <- p.test[idx,1]
sig.mat <- matrix(unlist( strsplit( nnames, c("_")) ), ncol=2, byrow=T)

# are any of these clades > in abundance in agg v. noagg
## reorder counts, dim will be smaller because some nodes sig in add & noagg
sub.tcnt <- subset(t.counts, rownames(t.counts) %in% sig.mat[,1])

### testing between the two data frames....
#test.tab <- t.counts #no power...
test.tab <- sub.tcnt

## match on samples
sub.cnt  <- t(test.tab)
sub.ord  <- sub.cnt[ match(rownames(map), rownames(sub.cnt)),]
## run kruskal and correct p values
agg.sub <- sub.ord[ map$Dog.Aggression == "aggressive", ]
nog.sub <- sub.ord[ map$Dog.Aggression == "not.aggressive", ]
merged  <- rbind( agg.sub, nog.sub )
#aggression factors
facs <- c( rep( 1, dim(agg.sub)[1]), rep( 2, dim(nog.sub)[1]) )
res  <- apply( merged, 2, function(x) kruskal.test( x, facs  )$p.value)
res.q <- qvalue(res)$qvalues
hist(res.q)
## What are the set of significant hits?
hits <- res.q[which(res.q < 0.2) ]
hit.counts <- merged[ , colnames(merged) %in% names(hits) ]
#greater in aggressive or non-aggressive?
med.agg  <- apply( agg.sub, 2, mean )
med.nog  <- apply( nog.sub, 2, mean )
med.rat  <- ( med.agg + 1) / ( med.nog + 1 ) #add small to avoid 0 or Inf
hits.rat <- med.rat[ which( names( med.rat ) %in% names(hits) ) ]
agg.rat  <- hits.rat[which(hits.rat > 1)]
## dataviz: heatmap
heatmap.2( log10(hit.counts+1), Rowv = rownames(hit.counts), trace = "none")
#build si_4 here
hit.tdf   <- subset( tax, tax$V1 %in% names(hits))
s.hit.tdf <- hit.tdf[ order(hit.tdf[,2]), ]
m.hit     <- match( names(hits), s.hit.tdf[,1])
s.hits    <- hits[order(m.hit)]
s.med.agg   <- med.agg[match( names(hits), names(med.agg))]
s.med.rat   <- med.rat[match( names(s.hits), names(med.rat))]

si.4 <- data.frame( node = names(s.hits), q.value = s.hits, taxon.string = s.hit.tdf$V2, 
                    fold.diff.aggressive.over.non = s.med.rat)
write.csv( si.4, file="Supplemental_Table_4.csv", row.names = FALSE)
## dataviz: boxplots
pdf( file = "dog_clade_strats.pdf" )
for( a in 1:length(s.hits)){
  hit <- s.hits[a]
  vec <- merged[,names(hit)]
  name <- tax[tax$V1 == names(s.hits[a]),]
  boxplot( vec~facs, names=c("Aggressive", "Not aggressive"), main=name )
}
dev.off()

# let's do some tree stuff. Something interesting with Lacto....
phy       <- as(tree.2, "phylo4")

get_tax_colors <- function( phy, taxon, tax.tab ){
  #phy is a phylo4 object
  #taxon is a taxon label associated with a node
  #tax.tab is a node to taxon lookup (df)
  ## get ordered node names
  labs  <- labels(phy)
  new   <- tax.tab$V2[match(unlist(labs), tax.tab$V1)]
  bin   <- new == taxon
  bin[which( bin == FALSE )] = "black"
  bin[which( is.na(bin))]    = "black"
  bin[which( bin == TRUE )]  = "red"
  return(bin)
} 

get_node_colors <- function( subtree, sig.nodes){
  labs <- labels(subtree)  
  to.color <- sig.nodes
  for( a in 1:length(sig.nodes)){
    #might throw warnings, but we literally could not care about them
    node.set <- unlist(descendants( subtree, sig.nodes[a], type=c("all") ))
    to.color <- c(to.color, names(node.set))
  }
  bin  <- labs %in% to.color 
  bin[which( bin == FALSE )] = "black"
  bin[which( bin == TRUE )]  = "red"
  return(bin)
}

taxon <- "g__Lactobacillus"
taxon <- "f__Lachnospiraceae"
taxon <- "g__Turicibacter"
d  <- data.frame( node=1:length(nodeId(phy)), color=get_tax_colors(phy, taxon, tax))
ggtree(phy) %<+% d + aes(color=I(color))

## grab a smaller set of the tree that enables visualization of 
## tree trends around nodes of interest
#clade    <- "node7141"
#clade   <- "node6841"

find.common.ancestor <- function( clade1, clade2, phy ){
  tmp1 <- subset( phy, node.subtree = clade1)
  my.node1 <- getNode(phy, node=c(clade1))
  anc1 <- ancestors( phy, my.node1 ) # use this to select nodes
  tmp2 <- subset( phy, node.subtree = clade2)
  my.node2 <- getNode(phy, node=c(clade2))
  anc2 <- ancestors( phy, my.node2 )
  idx  <- min(which(names(anc2) %in% names(anc1)))
  common <- names(anc2[idx])
  return(common)
}


## Turicibacter
clade1 <- "node1504"
clade2 <- "node1573"
common <- find.common.ancestor(clade1, clade2, phy)
taxon <- "g__Turicibacter"

## Blautia
clade1 <- "node6698"
clade2 <- "node5903"
common <- find.common.ancestor(clade1, clade2, phy)
taxon  <- "g__Blautia"

## Lactobacillus
clade1 <- "node7144"
clade2 <- "node3276"
common <- find.common.ancestor(clade1, clade2, phy)
taxon  <- "g__Lactobacillus"


#sub.node <- "node7140" #can also use the node value (e.g., 13929)
sub.node <- common
subtree <- subset( phy, node.subtree = sub.node)
d  <- data.frame( node=1:length(nodeId(subtree)), color=get_tax_colors(subtree, taxon, tax))
ggtree(subtree) %<+% d + aes(color=I(color)) #tips aren't red because we don't have their taxonomy

## Now let's plot OTU abundance across the samples for the otus within the clade? 
## Not actually used any more
otu.set <- descendants( subtree, sub.node, type=c("tips") )
tax.set <- otu.tax$V2[match(names(otu.set), otu.tax$V1)]

### heatmap...
build.node.map.full <- function( nodeval, subtree, otus, sig.nodes, hit.clade ){
  # make phylo object
  sub.phylo <- as(subtree, "phylo")
  # convert into ultrametric
  sub.clade <-  multi2di.phylo( sub.phylo )
  sub.clade$edge.length[which(sub.clade$edge.length == 0)] <- 0.00001
  # get otu data
  sub.otus <- subset(otus, rownames(otus) %in% sub.clade$tip.label)
  # rename and order animals
  colnames(sub.otus) <- gsub(x = colnames(sub.otus), pattern = "X", replacement = "")
  i        <- map$Dog.Aggression
  names(i) <- rownames(map)
  s        <- sort(i)
  order.map       <- map[ match( names(s), rownames(map)),]
  order.otus      <- sub.otus[, match( rownames(order.map), colnames( sub.otus) )]
  #prep clade dendro
  if( length( sub.clade$tip.label) > 2 ){
    clade.tree_um  <- chronopl(sub.clade,
                               lambda = 0.1,
                               tol = 0)
    clade.dend     <- as.dendrogram(as.hclust.phylo(clade.tree_um))
    # build otu dendo & order otus
    #otuclade.name  <- clade$tip.label
    otuclade.name  <- labels(clade.dend)
    otus.new.order <- match( otuclade.name,
                             row.names(order.otus))
    order.otus <- order.otus[otus.new.order,]
  } 
  else {
    otuclade.name  <- sub.clade$tip.label
    otus.new.order <- match( otuclade.name,
                             row.names(order.otus))
    order.otus <- order.otus[otus.new.order,]
  }  
  t.otus <- t( order.otus )
  run.heatmap = 0

  ## ggtree method
  d <- data.frame( node=1:length(nodeId(subtree)), color=get_node_colors( subtree, hit.clade))
  p <- ggtree( subtree, branch.length = "none" ) %<+% d + aes(color=I(color))
  p2<- gheatmap( p, log10(order.otus+0.1), low="white", high="black", colnames_angle = 90) 
  p2
}

clades <- c( clade1, clade2)

p <- build.node.map.full( sub.node, subtree, otus, sig.nodes, clades )
ggsave(  "Turicibacter_clades.pdf", p)

for( a in 1:length(sig.nodes)){
  clade   <- sig.nodes[a]
  my.node <- getNode(phy, node=c(clade), missing="fail")
  anc     <- ancestors( phy, my.node ) # use this to select nodes
  anc.node <- names(anc[2])
  #sub.node <- "node7140" #can also use the node value (e.g., 13929)
  #sub.node <- "node6833"
  subtree <- subset( phy, node.subtree = anc.node)
  taxon   <- subset( tax, tax$V1 == clade )
  file= paste(clade, ".heat_tree.pdf", sep="")
  p<- build.node.map( anc.node, subtree, otus, sig.nodes, clade )
  p <- p + ggtitle( taxon$V2 )
  ggsave( file, p )
}

# Three stories seem to emerge:
# 1. There are distinct clades w/in a genus enriched in no-agg (Blautia)
# 2. There are distinct clades w/in a genus enriched in agg (Lactobacillus)
# 3. There are clades within the genus Turicibacter that differ in association
# 
# This is in addition to clade-level analyses
# Let's build metacoder results and call it good.

## METACODER ANALYSIS
fmt  <- as.data.frame(apply( otu.tax[,1:2], 2, function(x) gsub( " ", "", x ) ))
fmt2 <- paste(fmt$V1, fmt$V2, sep=" ")
data <- extract_taxonomy( fmt2, 
                  regex = "^(.*) (.*)", 
                  key = c(id= "obs_info", "class"), 
                  class_sep = ";" )
heat_tree(data, node_size = n_obs, node_label = name, node_color = n_obs)

# need to get the values for each sample and group accordingly
otu.dat  <- t(otus)
rownames(otu.dat) <- gsub( "^X", "", rownames(otu.dat))
otu.ord  <- otu.dat[ match(rownames(map), rownames(otu.dat)),]
a.otu.tab <- as.data.frame(t(otu.ord))
a.otu.tax <- otu.tax[match( rownames( a.otu.tab ), otu.tax$V1 ),  2 ]
a.otu.tax <- gsub( " ", "",  a.otu.tax )
a.otu.tax <- gsub( ";", "|", a.otu.tax)
a.otu.df  <- cbind( rownames(a.otu.tab), a.otu.tab, a.otu.tax )
colnames(a.otu.df) <- noquote( colnames(a.otu.df) ) #gsub( "\"", "", colnames(a.otu.df))
write.table( a.otu.df, file="a.otu.df.tab", sep="\t", row.names = FALSE, quote = FALSE)
a.otu.data <- parse_taxonomy_table(file.path("a.otu.df.tab"), 
                                 taxon_col = c("class" = -1), 
                                 header = TRUE,
                                 class_sep = "\\|", 
                                 sep = "\t", 
                                 row.names = 1,
                                 quote = "")
map[,"sample_id"] <- rownames(map)
a.otu.data$mapping = map
data2 <- data
data  <- a.otu.data
colnames(data$obs_data) <- gsub( "\\\"", "", colnames(data$obs_data) )

calculate_abundance <- function(data, col = colnames(data$obs_data) %in% data$mapping$sample_id, col_name = "abundance") {
  data$taxon_data[[col_name]] <- vapply(obs(data), 
                                        function(i) sum(data$obs_data[i, col]), numeric(1))
  return(data)
}
calculate_prop <- function(data, col = colnames(data$obs_data) %in% data$mapping$sample_id, col_name = "abundance") {
  data$taxon_data[[col_name]] <- vapply(obs(data), 
                                        function(i) sum(data$obs_data[i, col]) , numeric(1)) / sum(data$obs_data[, col])
  return(data)
}

a.otu.data <- calculate_abundance(a.otu.data)


plot_all <- function(data, output_name, seed = 1) {
  set.seed(seed)
  data %>%
    filter_taxa(abundance >= min_read_count) %>%
    filter_taxa(name != "") %>% # Some taxonomic levels are not named
    heat_tree(node_size = n_obs,
              node_size_axis_label = "Number of OTUs",
              node_color = abundance,
              node_color_trans = "area",
              node_color_axis_label = "Number of reads",
              node_label = name,
              node_label_max = 100,
              overlap_avoidance = 1,
              output_file = output_name)
}

#use site for aggressive and not.aggressive

calculate_prop <- function(data, site) {
  sample_cols <- as.character(unlist(data$mapping[data$mapping$Dog.Aggression == site, "sample_id"]))
  sample_cols <- sample_cols[sample_cols %in% colnames(data$obs_data)]
  obs_indexes <- obs(data)
  total_counts <- vapply(sample_cols, function(s) sum(data$obs_data[, s]), numeric(1))
  col_name <- paste0(site, "_median_prop")
  lapply(obs_indexes,
         function(i) {
           vapply(sample_cols, 
                  function(s) sum(data$obs_data[i, s]) / total_counts[s],
                  numeric(1))})
}

calculate_prop_diff <- function(data, site_1, site_2) {
  
  remove_inf <- function(values) {
    values[values == Inf] <- 10000000000000000000000000
    values[values == -Inf] <- -10000000000000000000000000
    values[is.nan(values)] <- 0
    return(values)
  }
  
  props_1 <- calculate_prop(data, site_1)
  med_1   <- lapply( props_1, median )
  props_2 <- calculate_prop(data, site_2)
  med_2   <- lapply( props_2, median )
   p_value <- mapply(function(x, y) wilcox.test(x, y)$p.value,
                    props_1, props_2)
  p_value <- p.adjust(p_value, method = "fdr" )
  diff    <- mapply( "-", med_1, med_2, SIMPLIFY = FALSE )
  result <- ifelse( diff == 0, 0,
                   log2(vapply(props_1, median, numeric(1)) / vapply(props_2, median, numeric(1))))
  remove_inf(result)
}

plot_body_site_diff <- function(data, site_1, site_2, output_name, seed = 1) {
  set.seed(seed)
  data %>%
    mutate_taxa(median_prop_diff = calculate_prop_diff(data, site_1, site_2)) %>%
    filter_taxa( name == "k__Bacteria", subtaxa = TRUE ) %>% 
    filter_taxa( abs(median_prop_diff) > 0 ) %>%
    #filter_taxa(abundance >= min_read_count) %>%
    #filter_taxa(name != "") %>% # Some taxonomic levels are not named
    heat_tree(node_size_axis_label = "Number of OTUs",
              node_size = n_obs,
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color = median_prop_diff,
              node_color_range = diverging_palette(),
              node_color_trans = "linear",
              node_color_interval = color_interval,
              edge_color_interval = color_interval,
              node_label = name,
              node_label_max = 200,
              node_label_size_range = c(0.015, 0.02),
              #overlap_avoidance = 20,
              output_file = "Figure_2.pdf"
              #output_file = result_path(paste0(output_name, "--", site_1, "_vs_", site_2))
              )
  return( data )
}

color_interval <- c(-4, 4)
min_read_count <- 1
plot_all(a.otu.data, "test.pdf")
#rename to get rid of s__ and g__ strings:
a.otu.data$taxon_data[['name']][which( a.otu.data$taxon_data[['name']] == "g__" )] <- ""
a.otu.data$taxon_data[['name']][which( a.otu.data$taxon_data[['name']] == "s__" )] <- ""
a <- plot_body_site_diff( a.otu.data, "aggressive",  "not.aggressive" )



########

agg.otu <- otu.ord[ map$Dog.Aggression == "aggressive", ]
agg.tab <- as.data.frame(t(agg.otu))
agg.tax <- otu.tax[match( rownames( agg.tab ), otu.tax$V1 ),  2 ]
agg.tax <- gsub( " ", "", agg.tax )
agg.tax <- gsub( ";", "|", agg.tax)
agg.df  <- cbind( rownames(agg.tab), agg.tab, agg.tax )
write.table( agg.df, file="agg.df.tab", sep="\t", row.names = FALSE)
agg.data <- parse_taxonomy_table(file.path("agg.df.tab"), 
                             taxon_col = c("class" = -1), 
                             header = TRUE,
                             class_sep = "\\|", 
                             sep = "\t", 
                             row.names = 1 )

calculate_abundance <- function(data, col = colnames(data$obs_data) %in% data$mapping$sample_id, col_name = "abundance") {
  data$taxon_data[[col_name]] <- vapply(obs(data), 
                                        function(i) sum(data$obs_data[i, col]), numeric(1))
  return(data)
}

nog.otu <- otu.ord[ map$Dog.Aggression == "not.aggressive", ]
nog.tab <- as.data.frame(t(nog.otu))
nog.tax <- otu.tax[match( rownames( nog.tab ), otu.tax$V1 ),  2 ]
nog.tax <- gsub( " ", "", nog.tax )
nog.tax <- gsub( ";", "|", nog.tax)
nog.df  <- cbind( rownames(nog.tab), nog.tab, nog.tax )
write.table( nog.df, file="nog.df.tab", sep="\t", row.names = FALSE)
nog.data <- parse_taxonomy_table(file.path("nog.df.tab"), 
                                 taxon_col = c("class" = -1), 
                                 header = TRUE,
                                 class_sep = "\\|", 
                                 sep = "\t", 
                                 row.names = 1 )





agg.vec <- colSums(agg.otu)
nog.vec <- colSums(nog.otu)

agg.tax <- otu.tax[match( names( agg.vec ), otu.tax$V1 ),  2 ] 
agg.tax <- agg.tax[ complete.cases( agg.tax ), ]
agg.df  <- cbind( agg.tax, agg.vec  )
agg.fmt <- paste(agg.df$V1, agg.df$V2, agg.df$V3, sep=" ")
agg.data <- extract_taxonomy( agg.fmt, 
                          regex = "^(.*) (.*) (.*)", 
                          key = c(id= "obs_info", "class", "taxon_info"), 
                          class_sep = ";" )

nog.tax <- otu.tax[match( otu.tax$V1, names( nog.vec ) ),  1:2 ]
nog.tax <- nog.tax[ complete.cases( nog.tax ), ]
nog.df  <- cbind( nog.tax, nog.vec  )


# Things I need from Nicole
#1. OTU taxonomy - might have this on the server = have it
#2. How was this tree assembled?
#3. Random note: we don't have complete taxonomy strings, so this does weird viz stuff in some analyses
#4. Need tool that takes hit node, walks up tree to a reasonable ancestor, builds node map, and on
#   on that map indicates the node in question

build.node.map.full <- function( nodeval, subtree, otus, sig.nodes, hit.clade ){
  # make phylo object
  sub.phylo <- as(subtree, "phylo")
  # convert into ultrametric
  sub.clade <-  multi2di.phylo( sub.phylo )
  sub.clade$edge.length[which(sub.clade$edge.length == 0)] <- 0.00001
  # get otu data
  sub.otus <- subset(otus, rownames(otus) %in% sub.clade$tip.label)
  # rename and order animals
  colnames(sub.otus) <- gsub(x = colnames(sub.otus), pattern = "X", replacement = "")
  i        <- map$Dog.Aggression
  names(i) <- rownames(map)
  s        <- sort(i)
  order.map       <- map[ match( names(s), rownames(map)),]
  order.otus      <- sub.otus[, match( rownames(order.map), colnames( sub.otus) )]
  #prep clade dendro
  if( length( sub.clade$tip.label) > 2 ){
    clade.tree_um  <- chronopl(sub.clade,
                               lambda = 0.1,
                               tol = 0)
    clade.dend     <- as.dendrogram(as.hclust.phylo(clade.tree_um))
    # build otu dendo & order otus
    #otuclade.name  <- clade$tip.label
    otuclade.name  <- labels(clade.dend)
    otus.new.order <- match( otuclade.name,
                             row.names(order.otus))
    order.otus <- order.otus[otus.new.order,]
  } 
  else {
    otuclade.name  <- sub.clade$tip.label
    otus.new.order <- match( otuclade.name,
                             row.names(order.otus))
    order.otus <- order.otus[otus.new.order,]
  }  
  t.otus <- t( order.otus )
  run.heatmap = 0
  if( run.heatmap ){
    pdf( file=paste( nodeval, ".pdf", sep="" ) )
    if( length( sub.clade$tip.label) > 2 ){
      x <- as.matrix(log10(order.otus+0.01))
      
      heatmap.2(#t.otus.bin, 
        x,
        dendrogram="row",
        Rowv=clade.dend,
        Colv=FALSE,
        #Rowv=host.dend, 
        #RowSideColors = diet.cols,
        #margins=c(10,10),
        trace="none",
        #col=c(gray.colors(2)[2],gray.colors(2)[1]),
        #col=c("grey88","black"),
        #rowsep=c(1:dim(t.otus)[1]),
        #colsep=c(1:dim(t.otus)[2]),
        #sepcolor="white",
        labCol = colnames(x),
        key = TRUE
      )
    } 
    else{
      heatmap.2(t.otus.bin, 
                dendrogram="both",
                #Colv=clade.dend,
                Rowv=host.dend, 
                RowSideColors = diet.cols,
                margins=c(10,10),
                trace="none",
                #col=c(gray.colors(2)[2],gray.colors(2)[1]),
                col=c("grey88","black"),
                rowsep=c(1:dim(t.otus.bin)[1]),
                colsep=c(1:dim(t.otus.bin)[2]),
                sepcolor="white",
                labCol = FALSE,
                key = FALSE
      )
    }
    dev.off()
  }
  
  ## ggtree method
  d <- data.frame( node=1:length(nodeId(subtree)), color=get_node_colors( subtree, hit.clade))
  p <- ggtree( subtree, branch.length = "none" ) %<+% d + aes(color=I(color))
  p2<- gheatmap( p, log2(order.otus+0.01), low="blue", high="orange", colnames_angle = 90) 
  p2
}
