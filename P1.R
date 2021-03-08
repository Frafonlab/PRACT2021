install.packages("UniprotR")
install.packages("RCurl")
install.packages("bitops")
install.packages("devtools")
install.packages("dplyr", dependencies = TRUE)
install.packages("tidyverse")
install.packages("BiocManager")
install.packages("writexl")
install.packages("ALL")

BiocManager::install("UniProt.ws")
BiocManager::install("hipathia")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("topGO")
BiocManager::install("GO.db")
BiocManager::install("ALL")

library(hipathia)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(BiocManager)
library(UniProt.ws)
library(UniprotR)
library(devtools)
library(dplyr)
library(tidyverse)
library(writexl)
library(topGO)
library(GO.db)
library(ALL)

pathways <- load_pathways("hsa") ## Cargamos la base de datos de pathways en humanos = hsa
genes <- pathways$all.genes ## Cargamos los genes de la base de datos de pathways de Hipathia
genes
genes_ids <- data.frame( entrez = genes,
                         symbol = mapIds(x = org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENTREZID"),
                         stringsAsFactors = F)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
view(org.Hs.eg.db)
uniprots <- data.frame( entrez = genes,
                        uniprot = mapIds(x = org.Hs.eg.db, keys = genes, column = "UNIPROT", keytype = "ENTREZID"),
                        stringsAsFactors = F)
#GETPROTEINGOINFO

Obj_ProteinGOinfo <- GetProteinGOInfo(uniprots$uniprot)
BP_ProteinGOinfo <- select(Obj_ProteinGOinfo, Gene.ontology..biological.process.)

names(BP_ProteinGOinfo) #NOMBRE DE DATAFRAME
BP_ProteinGOinfo<-rename(BP_ProteinGOinfo, BiologicalProcess = Gene.ontology..biological.process.) #RENOMBRAR COLUMNA DE DATAFRAME

BPNA<-filter(BP_ProteinGOinfo, is.na(BiologicalProcess)) #VECTOR DE NULOS
BPVALOR<-filter(BP_ProteinGOinfo, !is.na(BiologicalProcess)) #VECTOR DE VALORES

uniprotsfinal <- uniprots
uniprotsfinal$BiologicalProcess <- "NA" ## Creamos una columna que va a contener la info de anotaciones que está en
## BP_ProteinGOinfo$BiologicalProcess, y la llenamos de "NA" para que así las proteínas que no tengan anotación
##se queden simplemente con el "NA". Es la forma de evitar campos vacíos y diferente número de filas.
idx_sub <- which (uniprotsfinal$uniprot %in% rownames(BP_ProteinGOinfo)) ## A %in% B: "%in%" lo que hace es buscar
## A en B y te da un vector de TRUE/FALSE, con "which" lo pasas a índice (ej. 1,5,9,7,19,30)
## en este caso una vez hagas el "which" el índice te va a indicar la posición de elementos
## en el vector uniprotsfinal$uniprot, primer argumento.                                                                 
idx_replace <- match(uniprotsfinal$uniprot, rownames(BP_ProteinGOinfo)) ## El match te da directamente el índice
## (ej. 1,4,8,7,12,34) con el orden en el que encuentra los identificadores uniprot$uniprot en 
## los rownames(BP_ProteinGOinfo). En este caso funciona al revés que "%in%", porque el índice corresponde con
## la posición en el vector rownames(BP_ProteinGOinfo), que es el segundo argumento de la función. 
## Utilizamos "match" porque en este caso quieres acceder a los elementos de BP_ProteinGOinfo$BiologicalProcess,
## pero te interesa que estén (que son las anotaciones de las proteinas) en el orden del data.frame uniprotsfinal.
v1= uniprotsfinal$BiologicalProcess[idx_sub]
v2= BP_ProteinGOinfo$BiologicalProcess[idx_replace]
uniprotsfinal$BiologicalProcess <- str_replace(v1, v2) 
write_xlsx(uniprotsfinal, "uniprotsfinal.xlsx")

#GETPROTEINFUNCTION

Obj_ProteinFunction <- GetProteinFunction(uniprots$uniprot)
FCC_ProteinFunction <- select(Obj_ProteinFunction, Function..CC.)

names(FCC_ProteinFunction) #NOMBRE DE DATAFRAME
FCC_ProteinFunction<-rename(FCC_ProteinFunction, FunctionCC = Function..CC.) #RENOMBRAR COLUMNA DE DATAFRAME

FCCNA<-filter(FCC_ProteinFunction, is.na(FunctionCC)) #VECTOR DE NULOS
FCCVALOR<-filter(FCC_ProteinFunction, !is.na(FunctionCC)) #VECTOR DE VALORES

uniprotsfinal2 <- uniprots
uniprotsfinal2$FuncionCC <- "NA" ## Creamos una columna que va a contener la info de anotaciones que está en
## BP_ProteinGOinfo$FuncionCC, y la llenamos de "NA" para que así las proteínas que no tengan anotación
##se queden simplemente con el "NA". Es la forma de evitar campos vacíos y diferente número de filas.
idx_sub2 <- which (uniprotsfinal$uniprot %in% rownames(FCC_ProteinFunction)) ## A %in% B: "%in%" lo que hace es buscar
## A en B y te da un vector de TRUE/FALSE, con "which" lo pasas a índice (ej. 1,5,9,7,19,30)
## en este caso una vez hagas el "which" el índice te va a indicar la posición de elementos
## en el vector uniprotsfinal$uniprot, primer argumento.                                                                 
idx_replace2 <- match(uniprotsfinal$uniprot, rownames(FCC_ProteinFunction)) ## El match te da directamente el índice
## (ej. 1,4,8,7,12,34) con el orden en el que encuentra los identificadores uniprot$uniprot en 
## los rownames(FCC_ProteinFunction). En este caso funciona al revés que "%in%", porque el índice corresponde con
## la posición en el vector rownames(FCC_ProteinFunction), que es el segundo argumento de la función. 
## Utilizamos "match" porque en este caso quieres acceder a los elementos de FCC_ProteinFunction$FuncionCC,
## pero te interesa que estén (que son las anotaciones de las proteinas) en el orden del data.frame uniprotsfinal2.
v3= uniprotsfinal2$FuncionCC[idx_sub2]
v4= FCC_ProteinFunction$FuncionCC[idx_replace2]
uniprotsfinal2$FuncionCC <- str_replace(v3, v4)
write_xlsx(uniprotsfinal2, "uniprotsfinal2.xlsx")

#saveRDS(Obj_ProteinFunction, file="Obj_ProteinFunction")

#################################################################################

#Paquete uniprot.ws

up <- UniProt.ws(taxId=9606)
uni<-up@taxIdUniprots
keytypes(up)
columns(up)

res <- UniProt.ws::select(up, keys = genes, columns = "UNIPROTKB", keytype =  "ENTREZ_GENE")
res2 <- UniProt.ws::select(up, keys = genes, columns = "KEYWORDS", keytype =  "ENTREZ_GENE")
vectorvalores<-filter(res, !is.na(res$FUNCTION))

#################################################################################

#Paquete go.db

columns(GO.db)
keytypes(GO.db)
data(ALL)
data

#BPtermsGO <- ls(GOTERM)
#head(BPtermsGO)

res2 <- GO.db::select(up, keys = genes, columns = "GO.db", keytype =  "ENTREZ_GENE")

# xy <- as.list(GOTERM)
# xxy <- as.list(GOBPCHILDREN)
# xyy <-as.list(GOBPPARENTS)
# xyyy <- as.list(GOBPANCESTOR)
# xxxy <- as.list(GOBPOFFSPRING)
 
vec1<-list()
for(i in 1:length(xy)){
  if(xy[[i]]@Ontology=="BP"){
    vec1<-c(vec1, xy[[i]])
  }
}

GOID<-list()
for(i in 1:length(vec1)){
  GOID<-c(GOID,vec1[[i]]@GOID)
}

TERM<-list()
for(i in 1:length(vec1)){
  TERM<-c(TERM,vec1[[i]]@Term)
}

unlist(GOID)
unlist(TERMS)
goterms<-data.frame(GOID=unlist(GOID), TERM=unlist(TERMS))

################################################################################

#Paquete topGO

BPterms <- ls(GOBPTerm)
head(BPterms)

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "entrez")
a<-annFUN.GO2genes("BP", feasibleGenes = genes, xx)
b<-stack(a)
b<-data.frame(values=b$ind, ind=b$values)
d<-unstack(b)


