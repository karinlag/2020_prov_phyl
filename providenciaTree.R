library(ggtree)
library(phangorn)
library(ape)
library(dplyr)
library(ggplot2)
library(tidytree)

setwd("C:/Users/vi1452/GitProjects/2020_prov_phyl")
linelist = "2020_10_12_curated_metadata.tsv"
linelist_df = read.table(linelist, sep="\t", header = TRUE)

# read tree, strip names to make them more readable
treefile = "gubbins.filtered_polymorphic_sites.phylip.contree"
tree = read.tree(treefile)
tree$tip.label <- sub("_pilon_spades.fasta", "", tree$tip.label)


# this section does the matching to get the tree tiplabels
# into the linelist. The purpose of the old_values is to figure
# out which ones that aren't in the tree, since the linelist contains
# two isolates that aren't in the tree
match_vector <- c()
old_values <- c()

for (i in linelist_df$Saksnr) {
  for (j in tree$tip.label){
    match <- grepl(i, j)
    if (match == TRUE) {
      match_vector <- c(match_vector, j)
      old_values <- c(old_values, i)
    }
  }
}

matching_df <- data.frame("Saksnr" = old_values, "tiplabel" = match_vector)

linelist_df <- left_join(linelist_df, matching_df, by="Saksnr") %>% 
  select(tiplabel, everything())


# so, now midpointing the tree
mid_tree = midpoint(tree, node.labels = "support")

# this next part is where I actually make the tree
test <- tidytree::as_tibble(mid_tree) %>%
  select(node, label) %>%
  mutate(test = if_else(as.numeric(label) >= 95, TRUE, FALSE)) %>%
  filter(test == TRUE)

nodes <- test$node

prettytree <- ggtree(mid_tree) %<+% linelist_df +
  geom_treescale() +
  
  geom_nodepoint(aes(subset = node %in% nodes, fill=))
  geom_tiplab(aes(label=Saksnr), align=T, linetype=NA, 
              size=2, offset=0.1, hjust=0.5) +
  geom_tiplab(aes(label=metadata.ID.for.phylogeny), align=T, linetype=NA, 
              size=2, offset=0.2, hjust=0.5)
prettytree


prettytree <- ggtree(mid_tree) %<+% linelist_df +
  geom_treescale() +
  geom_tiplab(aes(label=Saksnr), align=T, linetype=NA, 
            size=2, offset=0.1, hjust=0.5) +
  geom_tiplab(aes(label=metadata.ID.for.phylogeny), align=T, linetype=NA, 
              size=2, offset=0.2, hjust=0.5) +
  # size of the entire tree
  xlim(0, 0.3) +
  geom_nodepoint(aes(fill = ifelse(as.numeric(label) >= 95, ">= 95", "< 95")),
                 pch = 21,
                 size = 2) +
  labs(fill="STring")
   geom_tiplab(size = 1.6, offset=0.00005)

prettytree