library(ape)
phytophthora.tree <- read.nexus("Fig3.7locus.mrbayes.tre")
# 7 locus phylogeny for Phytophthora
# Bourret, Tyler; Abad, Gloria; Burgess, Treena (2023), 
# “Supporting data for taxonomic and phylogenetic revision of the genus 
# Phytophthora based on the types”, Mendeley Data, V1, doi: 10.17632/yg2rmfzstw.1
# https://data.mendeley.com/datasets/yg2rmfzstw/1


phytophthora.tree$tip.label
# some cleaning of phytophthora names

# clean the tip labels to match the trait database
phytophthora.tree$tip.label <- sapply(str_split(phytophthora.tree$tip.label, pattern = "_ET|_EE|_EN|_EpT|_NT|_RS|_SE"),
                                      "[[", 1)  

phytophthora.tree$tip.label <- gsub("P_", "Phytophthora_", phytophthora.tree$tip.label)
phytophthora.tree$tip.label <- gsub("_", " ", phytophthora.tree$tip.label)
phytophthora.tree$tip.label <- trimws(phytophthora.tree$tip.label)

setdiff(traits$phyto, phytophthora.tree$tip.label)

phytophthora.tree$tip.label <- gsub("Phytophthora x cambivora",
                                    "Phytophthora cambivora",
                                    phytophthora.tree$tip.label)
phytophthora.tree$tip.label <- gsub("gloveri", 
                                    "glovera", 
                                    phytophthora.tree$tip.label)
phytophthora.tree$tip.label <- gsub("Phytophthora andina", 
                                    "Phytophthora x andina", 
                                    phytophthora.tree$tip.label)
phytophthora.tree$tip.label <- gsub("Phytophthora marrasii", 
                                    "Phytophthora marrassi", 
                                    phytophthora.tree$tip.label)




# check all species in the trait database are represented in the phylogeny
setdiff(traits$phyto, phytophthora.tree$tip.label) 
# these are not represented on IDPhy
# https://idtools.org/phytophthora/index.cfm?packageID=1131

# create a covariance matrix
PhytophthoraA <- ape::vcv.phylo(phytophthora.tree)

# prune the phylogenies to only represent species in the dataset
# normally not important as brms only requires all species in the data to be represented in the phylogeny
# but here the matrix calculations for the cophylogenetic effects aren't possible in R
# and reducing the size of the phylogeny might help to calculate the kronecker product of the two covariance matrices

#phytophthora.tree <- keep.tip(phytophthora.tree, intersect(host_trees$phytophthora, phytophthora.tree$tip.label))
#length(phytophthora.tree$tip.label) # 157

# check if ultrametric
# is.ultrametric(phytophthora.tree)
# plot(phytophthora.tree)
# # and probably not just a numerical issue: the distnaces from root 
# # to tip are definitely quite different
# 
# #https://brettkvo.com/how-to-create-a-ultrametric-dichotomous-phylogeny-in-r-using-vertlife-org-or-birdtree-org/
# phytophthora.tree <- chronos(phytophthora.tree, lambda = 0)
# is.ultrametric(phytophthora.tree)
# 
# # check if bifurcating 
# is.binary(phytophthora.tree)
# 
# phytophthora.tree <- multi2di(phytophthora.tree)
# is.binary(phytophthora.tree)

#check_cophenentic <- cophenetic.phylo(phytophthora.tree)
#diag(check_cophenentic) <- NA
check_cophenentic <- reshape2::melt(cophenetic.phylo(phytophthora.tree))
check_cophenentic <- check_cophenentic[check_cophenentic$value>0,]
check_cophenentic[order(check_cophenentic$value),][1:20,]


# create covarinace matrix for modelling
A1 <- ape::vcv(phytophthora.tree, corr = TRUE) # # correlation matrix
A2 <- ape::vcv(phytophthora.tree, corr = FALSE) # vcv matrix


# these give quite different matrices 
ggplot() +
  geom_abline(slope = 1, color = "blue") +
  geom_point(aes(A1, A2))

# according to Shinichi Nakagawa, to estimate phylogenetic signal, we must use the 
# correlation matrix
# # https://discourse.mc-stan.org/t/covariance-matrix-phylogenetic-models/20477/3

PhytophthoraA <- ape::vcv.phylo(phytophthora.tree, corr = TRUE)
saveRDS(PhytophthoraA, file = "PhytophthoraA.rds")




