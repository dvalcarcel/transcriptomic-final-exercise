
#En primer lugar, instalamos los paquetes necesarios:
install.packages("readr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")
install.packages("RColorBrewer")
BiocManager::install("EnhancedVolcano")
BiocManager::install("mygene")

#Ahora cargamos las librerías necesarias:

library("readr")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library("mygene")

# Para esta práctica, he tenido en cuenta la Vignette de 
# DESeq2, en relación con el input desde la matriz de cuentas.
# Se puede consultar en el siguiente enlace:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 


# Establecemos nuestro directorio de trabajo:
setwd(dir = "/home/vant/Escritorio/transcriptomic-final-exercise/Apartado2/")


#Cargamos los datasets


# Leemos el archivo de datos crudos, como un csv, indicándole 
# que use como separadores los tabuladores.
# Dejamos el nombre de las filas como el nombre de los genes.
# Luego convertimos el data frame en una matriz, para que 
# pueda usarla DESeq2.
cts <- as.matrix(read.csv(file = "input/rawcounts.tsv", 
                             sep = "\t", row.names = 1))

# Para los metadatos del experimento, leemos el tsv como csv 
#y dejamos las filas con el nombre de los genes.
coldata <- read.csv(file = "input/metadata.tsv", 
                    sep = "\t", row.names = 1)

#Convertimos en factores las columnas del experimento, para 
# que lo pueda manejar DESeq2:

coldata$patient <- factor(coldata$patient)
coldata$agent <- factor(coldata$agent)
coldata$time <- factor(coldata$time)

# Examinamos la matriz de conteo y los datos de coldata 
# para ver si son coherentes en cuanto al orden de la muestra.

head(cts, 2)
coldata

# Se observa de manera visual que las muestran siguen 
# el mismo orden y el mismo nombre, por lo que no hay 
# que modificarlas.
# Sin embargo, vamos a comprobar, por un lado, si contamos
# con las mismas muestras en ambas:

all(rownames(coldata) %in% colnames(cts))

# y si se encuentran en el mismo orden:
all(rownames(coldata) == colnames(cts))

# Al devolvernos dos TRUE, verificamos estos dos criterios.

# Antes de crear el objeto DESeq, voy a unir las variables
# agent y time en una sola:

coldata$group <- as.factor(paste(coldata$agent, coldata$time, sep = "_"))

# Creamos el objeto DESeq,ajustando el modelo por el tipo de
# paciente, poniendo como factor de interés el grupo, es decir,
# el tratamiento o control a lo largo del tiempo.
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ patient + group)
dds

# Prefiltrado 
# Vamos a eliminar genes que contengan menos de 10 cuentas,
# para ahorrar memoria.

keep <- rowSums(counts(dds)) >= 10 

# Para eliminar los genes con esas cuentas bajas, aplicamos
# el siguiente código: 
dds <- dds[keep, ]

dds

# Análisis exploratorio
# Vamos a usar vst (variance stabilizing transformation)
# para estabilizar y normalizar la varianza de las cuentas.
# Esto nos sirve para la visualización con PCA.


vsd <- vst(dds, blind = TRUE)

# Hacemos el PCA con la variable patient:

plotPCA(vsd, intgroup = "patient")

# Se observa que el PC1 explica el 59% de la varianza y 
# el PC2 el 17%. En el gráfico se puede observar que los
# pacientes se separan unos de otros, salvo una muestra del
# paciente 4. 
# El PC1 diferencia bien los pacientes 1 y 2 del 3 y 4; mientras
# que el PC2 diferencia el paciente 1 y 3 del 2 y 4.

plotPCA(vsd, intgroup = "group")

# En este plot, no vemos una clara diferencia entre tratamientos.
# Vemos la misma distribución por pacientes.

# Ahora voy a generar un plot en el que combine el grupo con
# el paciente, para verlo de manera más clara:

pcaData <- plotPCA(vsd, intgroup=c("patient", "group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group, shape=patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# A continuación, vamos a calcular la matriz de distancias
# a partir de las cuentas normalizadas (vst) para analizar
# cuanto se separa cada muestra (paciente - grupo)

#Transponemos los valores, calculamos las distancias y lo convertimos
# a matriz:
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Ponemos los nombres del grupo en las filas y eliminamos los nombres
# de las columnas.
rownames(sampleDistMatrix) <- paste(vsd$patient, vsd$group, sep="-")
colnames(sampleDistMatrix) <- NULL

# Creamos el heatmap:
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# En el heatmap volvemos a ver que las muestras de un mismo 
# paciente tienden a agruparse juntas.

# Vamos a correr la función DESeq para realizar todos los pasos
# de DESeq2 y por consiguiente la expresión diferencial. 
# Además, realizaremos el modelo lineal teniendo en cuenta el test de Wald 

dds2 <- DESeq(dds, test = "Wald")

# Realizamos la estimación de la dispersión mediante un gráfico.
plotDispEsts(dds2)

# Ahora vamos a realizar el plotMa. 
# En DESeq2, la función plotMA muestra los cambios de log2 fold
# atribuibles a una variable dada sobre la media de recuentos normalizados 
# para todas las muestras del DESeqDataSet. 
# Los puntos se colorearán de azul si el valor p ajustado es inferior a 0,1.
# Los puntos que quedan fuera de la ventana se representan como triángulos 
# abiertos que apuntan hacia arriba o hacia abajo.

plotMA(dds2)
# Observamos que cuanto mayor es la media, mayor es el número
# de cambios significativos (azul).

## DEG de genes para el tratamiento con OHT vs Control 24h ##


# Obtención de nuestra lista de genes DEG para OHT vs Control a las 24h:
results_OHT_vs_Control_24h <- results(object = dds2,
                      contrast = c("group","OHT_24h", "Control_24h"), 
                      alpha = 0.05,
                      pAdjustMethod = "BH", 
                      tidy = TRUE
)

head(results_OHT_vs_Control_24h)

# Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
genesID <-mygene::queryMany(results_OHT_vs_Control_24h$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID <- genesID[!duplicated(genesID$query),]
results_OHT_vs_Control_24h$row <- ifelse(is.na(genesID$symbol),genesID$query,genesID$symbol)

head(results_OHT_vs_Control_24h)


# Realizamos el Heatmap de los 30 genes TOP DGE por p-valor ajustado
# para OHT vs Control 24h con los genes anotados por los símbolos de ensembl.

gene_symbols <- ifelse(is.na(genesID$symbol), genesID$query, genesID$symbol)

rownames(vsd) <- gene_symbols

top30_OHT_vs_Control_24h <- head(order(results_OHT_vs_Control_24h$padj), 30)
mat_OHT_Control <- assay(vsd)[top30_OHT_vs_Control_24h, ]

pheatmap(mat_OHT_Control, 
         scale = "row", 
         main = "Top 30 DEG OHT vs Control 24h",
         fontsize_row = 8)

# Ahora vamos a filtrar la DEG por los genes que son estadísticamente
# significativos en el p ajustados y los vamos a ordenar según su
# log2FoldChange:

deg_filtered_OHT_vs_Control_24h <- results_OHT_vs_Control_24h[results_OHT_vs_Control_24h$padj < 0.05, ]
deg_filtered_OHT_vs_Control_24h_clean <- na.omit(deg_filtered_OHT_vs_Control_24h)

print(deg_filtered_OHT_vs_Control_24h_clean$row)
nrow(deg_filtered_OHT_vs_Control_24h_clean)

deg_filtered_OHT_vs_Control_24h_clean <- 
  deg_filtered_OHT_vs_Control_24h_clean[order(deg_filtered_OHT_vs_Control_24h_clean$log2FoldChange, 
                                              decreasing = TRUE), ]
head(deg_filtered_OHT_vs_Control_24h_clean)

# Observamos que hay 48 genes que se encuentran diferencialmente 
# expresados entre las muestras pertenecientes al grupo tratado 
# con OHT con respecto al control tras 24h.

# Los genes CGA, CBFA2T3 y HLA-DRB5 se encuentran más expresados
# en las muestras con OHT frente al control a las 24h.

# Los genes PRB4, GDA y NEFH se encuentran menos expresados
# en las muestras con OHT frente al control a las 24h.

# Guardamos la tabla para añadirla al informe:

write.csv(deg_filtered_OHT_vs_Control_24h_clean, 
          "DEG_filtered_OHT_vs_Control_24h.csv", row.names = FALSE)

# Creamos el Volcano plot para visualizarlo.

EnhancedVolcano(results_OHT_vs_Control_24h,
                lab = results_OHT_vs_Control_24h$row,
                x = "log2FoldChange",
                y = "padj",
                title = "OHT vs Control 24h",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 4.0)



## DEG de genes para el tratamiento con DPN vs Control 24h ##


# Obtención de nuestra lista de genes DEG para DPN vs Control a las 24h:
results_DPN_vs_Control_24h <- results(object = dds2,
                                      contrast = c("group","DPN_24h", "Control_24h"), 
                                      alpha = 0.05,
                                      pAdjustMethod = "BH", 
                                      tidy = TRUE
)

head(results_DPN_vs_Control_24h)

# Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
genesID_DPN <-mygene::queryMany(results_DPN_vs_Control_24h$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID_DPN <- genesID_DPN[!duplicated(genesID_DPN$query),]
results_DPN_vs_Control_24h$row <- ifelse(is.na(genesID_DPN$symbol),genesID_DPN$query,genesID_DPN$symbol)

head(results_DPN_vs_Control_24h)


# Realizamos el Heatmap de los 30 genes TOP DGE por p-valor ajustado
# para DPN vs Control 24h con los genes anotados por los símbolos de ensembl.

gene_symbols_DPN <- ifelse(is.na(genesID_DPN$symbol), genesID_DPN$query, genesID_DPN$symbol)

rownames(vsd) <- gene_symbols_DPN

top30_DPN_vs_Control_24h <- head(order(results_DPN_vs_Control_24h$padj), 30)
mat_DPN_Control <- assay(vsd)[top30_DPN_vs_Control_24h, ]

pheatmap(mat_DPN_Control, 
         scale = "row", 
         main = "Top 30 DEG DPN vs Control 24h",
         fontsize_row = 8)

# Ahora vamos a filtrar la DEG por los genes que son estadísticamente
# significativos en el p ajustados y los vamos a ordenar según su
# log2FoldChange:

deg_filtered_DPN_vs_Control_24h <- results_DPN_vs_Control_24h[results_DPN_vs_Control_24h$padj < 0.05, ]
deg_filtered_DPN_vs_Control_24h_clean <- na.omit(deg_filtered_DPN_vs_Control_24h)

print(deg_filtered_DPN_vs_Control_24h_clean$row)
nrow(deg_filtered_DPN_vs_Control_24h_clean)

deg_filtered_DPN_vs_Control_24h_clean <- 
  deg_filtered_DPN_vs_Control_24h_clean[order(deg_filtered_DPN_vs_Control_24h_clean$log2FoldChange, 
                                              decreasing = TRUE), ]
head(deg_filtered_DPN_vs_Control_24h_clean)

# Observamos que hay 87 genes que se encuentran diferencialmente 
# expresados entre las muestras pertenecientes al grupo tratado 
# con DPN con respecto al control tras 24h.

# Los genes THAP9-AS1, CTTNBP2 y SCD se encuentran más expresados
# en las muestras con DPN frente al control a las 24h.

# Los genes PRB4, GDA y ATP1A3 se encuentran menos expresados
# en las muestras con DPN frente al control a las 24h.

# Guardamos la tabla para añadirla al informe:

write.csv(deg_filtered_DPN_vs_Control_24h_clean, 
          "DEG_filtered_DPN_vs_Control_24h.csv", row.names = FALSE)

# Creamos el Volcano plot para visualizarlo.

EnhancedVolcano(results_DPN_vs_Control_24h,
                lab = results_DPN_vs_Control_24h$row,
                x = "log2FoldChange",
                y = "padj",
                title = "DPN vs Control 24h",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 4.0)

#Por último, guardo el objeto dds2 para poder utilizarlo en la
# siguiente pregunta:

saveRDS(dds2, file = "dds2.rds")
