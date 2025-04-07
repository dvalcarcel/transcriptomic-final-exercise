#!/bin/bash


# Ruta al archivo de adaptadores descargado de https://raw.githubusercontent.com/cnio-bu/cluster_rnaseq/refs/heads/master/resources/trimming/adapters.fa
ADAPTERS="/home/vant/Escritorio/transcriptomic-final-exercise/Apartado1/input/trimming/adapters.fa"

# Loop para cada muestra
for sid in $(ls data/*.fastq | cut -d"/" -f2 | cut -d"_" -f1 | sort | uniq) #Defino el sample id de los archivos fastq que tengo en la carpeta data.
do
    echo "Procesando trimming de sample: $sid..."

    bbduk.sh \
      in1=data/${sid}_1.fastq \ 
      in2=data/${sid}_2.fastq \
      out1=trimming/${sid}_1.trimmed.fastq \
      out2=trimming/${sid}_2.trimmed.fastq \
      ref=$ADAPTERS \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1 \
      tbo \
      tpe \
      qtrim=rl \
      trimq=20 \
      stats=trimming/${sid}_bbduk_stats.txt

      #Defino los dos fastq de entrada (in1 e in2) y los dos fastq de salida (out1 e out2).
      #Defino el archivo fasta con los adaptadores (ref).
      #ktrim=r: recorta adaptadores de 3' a 5'.
      #k=23: longitud de la secuencia de búsqueda mediante el kmer.
      #mink=11: longitud mínima del kmer.
      #hdist=1: distancia máxima permitida entre el kmer y la secuencia de búsqueda.
      #tbo: especifica que también recorte los adaptadores en función de la detección de superposición de pares utilizando BBMerge.
      #tpe: recorta ambos reads a la misma longitud.
      #qtrim=rl: recortará la calidad a Q10 utilizando el algoritmo Phred.
      #trimq=20: recorta la calidad mínima después del recorte a 20.
      #stats: genera un archivo de estadísticas de recorte.

    echo "Trimming de sample: $sid terminado."
done
