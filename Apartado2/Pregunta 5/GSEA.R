#En primer lugar, instalamos los paquetes necesarios:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("tidyverse")
BiocManager::install("apeglm")


#Ahora cargamos las librerías necesarias:

library("tidyverse")
library("DESeq2")
library("apeglm")

# Para esta práctica, he tenido en cuenta la Vignette de 
# DESeq2, en relación con el input desde la matriz de cuentas.
# Se puede consultar en el siguiente enlace:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 


# Establecemos nuestro directorio de trabajo:
setwd(dir = "/home/vant/Escritorio/transcriptomic-final-exercise/Apartado2/")



# Cargamos los datos obtenidos de la pregunta anterior:
dds2 <- readRDS("dds2.rds")

# Revisamos los metadatos:
colData(dds2)
dds$patient
dds$group

# Revisamos el diseño:
design(dds2)

# Hacemos la DEG de DPN vs Control a las 24h.
res <- results(dds2, alpha = 0.05, contrast = c("group", "DPN_24h", "Control_24h"))
summary(res)
res

# Como vamos a realizar un GSEA preranked, necesitamos crear el
# archivo .rnk
# Para ello, en primer lugar, reducimos el log fold change para
# encoger los genes con bajos recuentos con sus LFCs, ya que
# generan imprecisión. De esta manera, nos quedamos con los genes
# de mayor recuento para mejorar la visualización del ranking de genes.
# Usamos la función lfcShrink de DESeq2 con el método apeglm:
resultsNames(dds2)
res.ape <- lfcShrink(dds2, coef = "group_DPN_24h_vs_Control_24h", type = "apeglm",
                     res = res)
summary(res.ape) # Mismo número de genes up/down padj < 0.05

# Explicación visual del encogimiento de los genes con bajos recuentos:
par(mfrow = c(1, 2))
plotMA(res, ylim = c(-3, 3))
plotMA(res.ape, ylim = c(-3, 3))


# Creamos el archivo .rnk
rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)
head(rnk)

# Guardamos el archivo .rnk (sin la cabecera y separado por tabuladores)
write.table(rnk, file = "DPN_Control_24h.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

# Pasamos a ejecutar GSEA en el terminal de bash (Ver indicaciones en el informe)


#Una vez obtenidos los resultados de GSEA, voy a filtrar aquellos
# genes enriquecidos:

# Muestras tratadas con DPN

DPN_perturbed <- read.table(file = "GSEA_output/DPN_Control_24h.GseaPreranked.1745243641608/DPN_PERTURBED.tsv",
                            header = T, sep = "\t")
#Filtro para quedarme con los que aparezcan en el Leading Edge:
DPN_perturbed_refined <- DPN_perturbed[DPN_perturbed[,6] == "Yes", c(2,3:5)]
DPN_perturbed_refined

# Guardamos los resultados en un tsv:
write.table(DPN_perturbed_refined,
            file = "GSEA_output/DPN_perturbed_refined_24h.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Muestras control

DPN_unperturbed <- read.table(file = "GSEA_output/DPN_Control_24h.GseaPreranked.1745243641608/DNP_UNPERTURBED.tsv",
                            header = T, sep = "\t")
#Filtro para quedarme con los que aparezcan en el Leading Edge:
DPN_unperturbed_refined <- DPN_unperturbed[DPN_unperturbed[,6] == "Yes", c(2,3:5)]
DPN_unperturbed_refined

# Guardamos los resultados en un tsv:
write.table(DPN_unperturbed_refined,
            file = "GSEA_output/DPN_unperturbed_refined_24h.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
