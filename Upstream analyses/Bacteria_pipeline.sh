[TOC]


    wd=/e/2025Crop_Rotation_Sequencing_Analysis_Overall/Bacteria
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}


    csvtk -t stat result/metadata_raw.txt
    cat -A result/metadata_raw.txt | head -n3
    sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
    cat -A result/metadata.txt | head -n3

    cat temp/*.merged.fq > temp/all.fq
    ls -lsh temp/all.fq
    head -n 6 temp/all.fq|cut -c1-60
    awk 'NR%4==1{print}' temp/all.fq | cut -d'.' -f1 | sort | uniq | wc -l


time vsearch --fastx_filter temp/all.fq \
  --fastq_stripleft 0 --fastq_stripright 0 \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/filtered.fa \
  2> fastx_filter.log
  
  cat fastx_filter.log
  
    head temp/filtered.fa
    
    head -n 20 temp/filtered.fa
    
    awk 'NR%4==1{print substr($1,2)}' temp/all.fq | cut -d'.' -f1 | sort | uniq > all_samples.txt
    grep "^>" temp/filtered.fa | sed 's/^>//' | cut -d'.' -f1 | sort | uniq > filtered_samples.txt
 
    comm -23 all_samples.txt filtered_samples.txt > lost_samples.txt
    wc -l lost_samples.txt
    head lost_samples.txt    

    head temp/filtered.fa


    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 10 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 

    ls -lsh temp/uniques.fa

    head -n 2 temp/uniques.fa


    vsearch --cluster_unoise temp/uniques.fa \
  --minsize 150 \
  --centroids temp/otus.fa \
  --relabel ASV_

head -n 2 temp/otus.fa
      
    mkdir -p result/raw
    cp -f temp/otus.fa result/raw/otus.fa


    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 24 \
    	--otutabout result/raw/otutab.txt 
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    csvtk -t stat result/raw/otutab.txt

    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/rdp_16s_v18.fa \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax | cat -A
    sed -i 's/\r//' result/raw/otus.sintax

    wc -l result/raw/otutab.txt
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    wc -l result/otutab.txt
cut -f 1 result/otutab.txt | tail -n +2 > result/otutab.id
vsearch --fastx_getseqs result/raw/otus.fa \
  --labels result/otutab.id \
  --fastaout result/otus.fa
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax
    cp result/raw/otu* result/

Rscript -e "x<-read.delim('result/otutab.txt',check.names=FALSE); \
s<-colSums(x[,-1,drop=FALSE]); \
out<-data.frame(Samples=ncol(x)-1,Features=nrow(x)-1, \
Min=min(s),Q1=as.numeric(quantile(s,0.25)),Median=median(s), \
Mean=mean(s),Q3=as.numeric(quantile(s,0.75)),Max=max(s)); \
write.table(out,'result/otutab_sample_depth.txt',sep='\t',row.names=FALSE,quote=FALSE)"
cat result/otutab.stat

Rscript -e "x<-read.delim('result/otutab.txt',check.names=FALSE); \
s<-colSums(x[,-1,drop=FALSE]); \
out<-data.frame(Sample=names(s),Depth=as.numeric(s)); \
write.table(out,'result/otutab_sample_depth.txt',sep='\t',row.names=FALSE,quote=FALSE)"
cat result/otutab_sample_depth.txt

mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
  --depth 10000 --seed 1 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt
Rscript -e "x<-read.delim('result/otutab_rare.txt',check.names=FALSE); \
s<-colSums(x[,-1,drop=FALSE]); \
out<-data.frame(Samples=ncol(x)-1,Features=nrow(x)-1, \
Min=min(s),Q1=as.numeric(quantile(s,0.25)),Median=median(s), \
Mean=mean(s),Q3=as.numeric(quantile(s,0.75)),Max=max(s)); \
write.table(out,'result/otutab_rare.stat',sep='\t',row.names=FALSE,quote=FALSE)"
cat result/otutab_rare.stat

    Rscript ${db}/script/otu_mean.R -h
    Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 200 --all TRUE --type mean \
      --output result/otutab_mean.txt
    head -n3 result/otutab_mean.txt

    mkdir -p result/beta/
Rscript - <<'EOF'
suppressPackageStartupMessages({ library(vegan) })

tab <- read.delim("result/otutab_rare.txt", check.names = FALSE, quote = "", comment.char = "")
X <- as.matrix(tab[,-1, drop = FALSE]); rownames(X) <- tab[,1]
M <- t(X)  

bray <- as.matrix(vegdist(M, method = "bray"))
jac  <- as.matrix(vegdist(M > 0, method = "jaccard", binary = TRUE))
# Euclidean & Manhattan
euc  <- as.matrix(dist(M, method = "euclidean"))
man  <- as.matrix(dist(M, method = "manhattan"))

w <- function(mat, file) {
  out <- data.frame(Sample = rownames(mat), as.data.frame(mat), check.names = FALSE)
  write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
}
dir.create("result/beta", showWarnings = FALSE, recursive = TRUE)
w(bray, "result/beta/bray_curtis.txt")
w(jac,  "result/beta/jaccard.txt")
w(euc,  "result/beta/euclidean.txt")
w(man,  "result/beta/manhattan.txt")
EOF

    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    mkdir -p result/tax
#
Rscript - <<'EOF'
tab <- read.delim("result/otutab_rare.txt", check.names = FALSE, quote = "", comment.char = "")
X   <- as.matrix(tab[,-1, drop = FALSE]); rownames(X) <- tab[,1]
cs  <- colSums(X); cs[cs==0] <- 1
Xp  <- sweep(X, 2, cs, "/") * 100  # 百分比
out <- data.frame(OTUID = rownames(Xp), round(Xp, 4), check.names = FALSE)
write.table(out, "result/otutab_rare_pct.txt", sep="\t", quote=FALSE, row.names=FALSE)
EOF

#
Rscript - <<'EOF'
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
otu <- read.delim("result/otutab_rare_pct.txt", check.names = FALSE, quote = "", comment.char = "")
feat <- otu[[1]]; mat <- as.matrix(otu[,-1, drop = FALSE]); rownames(mat) <- feat
sx <- read.delim("result/otus.sintax", header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
sx <- sx[,1:2]; colnames(sx) <- c("Feature","Tax")
get_rank <- function(tx, r){
  m <- regexpr(paste0("(^|,)", r, ":[^,]+"), tx, perl = TRUE)
  s <- ifelse(m > 0, regmatches(tx, m), NA)
  s <- gsub("^,?", "", s); s <- sub(paste0("^", r, ":"), "", s); s <- sub("\\([^)]*\\)$", "", s)
  s[is.na(s)] <- "Unassigned"; s <- gsub("[()\"#]", "", s); s <- gsub("/chloroplast", "", s, fixed = TRUE); s
}
taxdf <- data.frame(Feature = sx$Feature,
                    p = get_rank(sx$Tax,"p"), c = get_rank(sx$Tax,"c"),
                    o = get_rank(sx$Tax,"o"), f = get_rank(sx$Tax,"f"),
                    g = get_rank(sx$Tax,"g"), check.names = FALSE)
common <- intersect(rownames(mat), taxdf$Feature)
mat <- mat[common,,drop=FALSE]; taxdf <- taxdf[match(common, taxdf$Feature), ]
write_sum <- function(rank_char, out_file){
  lab <- taxdf[[rank_char]]
  agg <- rowsum(mat, group = lab, reorder = FALSE)  # 
  out <- data.frame(Taxon = rownames(agg), round(agg,4), check.names = FALSE)
  out <- out[order(rowMeans(out[,-1, drop=FALSE]), decreasing = TRUE), ]
  write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
dir.create("result/tax", showWarnings = FALSE, recursive = TRUE)
write_sum("p", "result/tax/sum_p_pct.txt")
write_sum("c", "result/tax/sum_c_pct.txt")
write_sum("o", "result/tax/sum_o_pct.txt")
write_sum("f", "result/tax/sum_f_pct.txt")
write_sum("g", "result/tax/sum_g_pct.txt")
EOF


wc -l result/tax/sum_*.txt
head -n 3 result/tax/sum_g_pct.txt

