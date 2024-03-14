install.packages("devtools")
library(devtools)
install_github("roblanf/sangeranalyseR",ref="master")
install_github("roblanf/sangeranalyseR",ref="develop")
library(sangeranalyseR)
BiocManager::install("DECIPHER")
BiocManager::install("sangerseqR")
BiocManager::install("BiocStyle")
rawDataDir <- system.file("C:/Users/LENOVO/Desktop/BIOINFORMÁTICA/Secuencias ATPasa/Placa_TORATP/PRÁCTICA/ABI", package = "sangeranalyseR")
my_aligned_contigs <- SangerAlignment(ABIF_Directory     = "C:/Users/LENOVO/Desktop/BIOINFORMÁTICA/Secuencias ATPasa/Placa_TORATP/PRÁCTICA/ABI",processMethod = "REGEX",
                                 REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                 REGEX_SuffixReverse = "_[0-9]*_R.ab1$")
launchApp(SangerAlignment)
library("Biostrings")
library("BiocManager")

dna<-readDNAStringSet(fas)

#####TAREA PRÁCTICA####
##https://www.ncbi.nlm.nih.gov/popset/2450646312 ###base de datos

#://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html ###pagina saque parte del script

library(DECIPHER)
fas <- "sequence.fasta"

library(Biostrings)

dna <- readDNAStringSet(fas)



###clustal
library(msa)

mySequences <- readAAStringSet(fas)

myFirstAlignment <- msa(mySequences)

###ClustalW

myClustalWAlignment <- msa(mySequences, "ClustalW")

myClustalWAlignment


###ClustalOmega
myClustalOmegaAlignment <- msa(mySequences, "ClustalOmega")
myClustalOmegaAlignment

###Muscle
myMuscleAlignment <- msa(mySequences, "Muscle")
myMuscleAlignment


###para salvar el alineamiento

writeXStringSet(unmasked(myClustalWAlignment), file="testclustal3.fasta")
##TERMINO DE TAREA PENDIENTE DEL JUEVES 29 DE FEBRERO#

#####analisis phylogeneticos

library(ape)
# Carga las secuencias desde un archivo FASTA
my_sequences <- read.dna("testclustal3.fasta", format = "fasta")


# Estima el modelo de evolución
install.packages("phangorn")
library(phangorn)



###obtener el modelo de evolucion molecular###

library(phangorn)
library(ape)


canis <- read.dna("testclustal3.fasta", format="fasta")
canis_phyDat <- phyDat(canis, type = "DNA", levels = NULL)



mt <- modelTest(canis_phyDat)

mt

write.csv(mt, "mt2.csv")


HKY+I## fue el mejor modelo con BIC


##para correr Mrbayes necesitamos nuestra base en formato nexus, para convertir el formato fasta to nexus###
##EN ESTE PUNTO NO ME SALE BIEN PORQUE:#

library(seqinr)

data = read.fasta("testclustal3.fasta") #ME DICE QUE NO ES POSIBLE ABRIR EL ARCHIVO PORQUE NO EXISTE EL ARCHIVO O EL DIRECTORIO#
##13Marzo: Ya lo resvolví, mi documento era testclustal3 no testclustal2
library(ape)


write.nexus.data(data, file = "my_sequences3.nexus", interleave = FALSE)



####maximumlikelihood in R

library(ape)
library(phangorn)
library(msa)
library(Biostrings)
library(seqinr)

fit <- as.pml(mt, "HKY+I")
fit <- as.pml(mt, "BIC")

fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt


bs <- bootstrap.pml(fit_mt, bs=100, optNni=TRUE,
                    control = pml.control(trace = 0))



plotBS(midpoint(fit_mt$tree), p = .5, type="p", digits=2, main="Ultrafast bootstrap")

plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", main="Standard bootstrap")

plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", digits=0, method = "TBE", main="Transfer bootstrap")


###cortar nombres de las secuencias

fasta <-readAAStringSet("testclustal3.fasta")
print(fasta) 


# Loop through the sequences and trim the pattern
for (i in 1:length(fasta)) {sequence <- as.character(fasta[[i]])}  # convert the sequence to a character string
  sequence <- gsub(".*(Canis.*lupus).*", "\\1", sequence)  # remove everything that comes before the pattern
  fasta[[i]] <- AAString(sequence)  # convert the trimmed sequence back to a AAString object
  

  # Trim the name of the sequence
  header <- names(fasta)[i]
  header <- gsub("^>\\s*", "", header) # remove the ">" and any leading whitespace
  header <- gsub("\\s.*", "", header) # remove any additional information after the name
  names(fasta)[i] <- header # set the new name for the sequence
  
  # Write the trimmed sequences to a new FASTA file
  writeXStringSet(fasta, "file2.fasta")
  
################Repetir todo desde escoger el mejor modelo mutacional
  
  
  mammals <- read.dna("file2.fasta", format="fasta")
  mammals_phyDat <- phyDat(mammals, type = "DNA", levels = NULL)
  
  
  
  mt <- modelTest(mammals_phyDat)
  
  print(mt)
  
  write.csv(mt, "mt.csv")
  
  
  HKY+I## fue el mejor modelo con BIC
  
  
  ####maximumlikelihood in R
  
  library(ape)
  library(phangorn)
  
  
  data = read.fasta("file2.fasta")
  
  
  fit <- as.pml(mt, "HKY+I")
  fit <- as.pml(mt, "BIC")
  
  fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
  fit_mt
  
  
  bs <- bootstrap.pml(fit_mt, bs=100, optNni=TRUE,
                      control = pml.control(trace = 0))
  
  #for fasta2.fasta
  
  plotBS(midpoint(fit_mt$tree), p = .5, type="p", digits=2, main="Ultrafast bootstrap")
  
  plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", main="Standard bootstrap")
  
  plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", digits=0, method = "TBE", main="Transfer bootstrap")
  
  
  # assigning standard bootstrap values to our tree; this is the default method
  tree_stdbs <- plotBS(fit_mt$tree, bs, type = "n")
  
  # assigning transfer bootstrap values to our tree
  tree_tfbs <- plotBS(fit_mt$tree, bs, type = "n", method = "TBE")
  
  
  
  cnet <- consensusNet(bs, p=0.2)
  plot(cnet, show.edge.label=TRUE)
  

#####################Exporting a tree

Now that we have our tree with bootstrap values, we can easily write it to a file in Newick-format:

# tree with ultrafast bootstrap values
write.tree(fit_mt$tree, "canis1.tree")

# tree with standard bootstrap values
write.tree(tree_stdbs, "canis2.tree")

# tree with transfer bootstrap values
write.tree(tree_tfbs, "canis3.tree")

write.tree(tree_stdbs, file = "canis2.nwk")


#####visualizaciones mas bonitos para los arboles########


##https://yulab-smu.top/treedata-book/chapter4.html

library (BiocManager)
BiocManager::install("treeio")
a
BiocManager::install("ggtree")
a
library("treeio")
library("ggtree")



nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()


mytree <- read.tree("canis2.tree")

ggtree(mytree) + 
  geom_point(aes(shape=isTip, color=isTip), size=3)

p <- ggtree(mytree) + 
  geom_nodepoint(color="lightpink2", alpha=1/4, size=10) 
p + geom_tippoint(color="steelblue", shape=8, size=3)



library(ape)

# Read in .tree file
mytree <- read.tree("canis2.tree")

# Check the class of the object
class(mytree)

##If the file is in a different format, such as a Newick tree file, you can still use read.tree() to read in the file and convert it to a "phylo" object. For example:
  
  r


# Convert to "phylo" object
mytree <- as.phylo(mytree)

# Check the class of the object
class(mytree)


# Read in Newick file
mytree <- read.tree("canis2.nwk")



####hacer una red de haplotipos#####

library("ape")
install.packages ("pegas")
library("pegas")
naso<-read.dna("file2.fasta", format="fasta") 
naso

# list of indivduals
table(rownames(naso))
nasoHaps1 <- haplotype(naso)
nasoHaps1

ind.hap<-with(
  stack(setNames(attr(nasoHaps1, "index"), rownames(nasoHaps1))), 
  table(hap=ind, pop=rownames(naso)[values])
)
ind.hap[1:6, 1:3]  #print just a chunk


net <- haploNet(nasoHaps1)
plot(net, size = attr(net, "freq"))

plot(nasoNet, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8, pie=ind.hap)
legend(50,50, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
#AQUI ME SALIO ERROR Y E NOMBRE NO ES nasoNet ENTONCES LO CAMBIE SOLO A naso#

####YA NO PUEDO PASAR DE AQUI, EL my_dist no esta
my_dist <- dist.haplo(naso)
my_network <- haploNet(naso, dist = my_dist)
plot(my_network, size = attr(my_network, "freq"), cex = 2)
