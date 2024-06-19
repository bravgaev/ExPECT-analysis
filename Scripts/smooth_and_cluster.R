
# R version 3.5.2

# -----
# Parameters
# -----

options(scipen=999)
# no scientific notation
options(digits=5)
# rounds numbers to save diskspace

sample.dirs <- c('') # multiple allowed
output.dir <- ''

#prefix <- c()  # multiple prefixes allowed; can also be a list of samples
prefix <- c("DNA")
                                              
chr.in.seqnames <- T                           
# state whether 'chr1' or '1' is used          
                                               
chrs <- as.character(c(1:22))                  
# or e.g. c(1:22, 'X', 'Y'); which chrs should be included?

min.CpG.depth <- 1
# makes results more robust (& saves discspace); default = 1 or more (no filtering)

fraction.NA.allowed <- .5 # used to be .75 !
# removes rare observations (& saves discspace); default = more than 75% of samples are NA --> CpG removed

n.cores <- 18
# cores used during smoothening

do.cluster <- T
# whether CpGs should be clustered

cluster.max.len.between.CpG <- 50
# max length between CpGs of same clusters

cluster.min.CpG <- 2
# min number of CpGs in one cluster

cluster.min.length <- 0
# min length of cluster

cluster.min.coverage <- 5
# min mean coverage of cluster

cluster.min.cpg.coverage <- 1
# min mean coverage per CpG in cluster

cluster.fraction.NA.allowed <- .5
# removes rare observations (& saves discspace); default = more than 50% of samples are NA --> cluster removed

mean.of.reads <- T # only applies to non-smoothened data; alternative: (default: mean of sites; F)

# use `tail -f smooth.and.cluster.o` to check progress

# -----
# Library
# -----

suppressMessages(library("bsseq"))
suppressMessages(library("data.table"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("HDF5Array"))
# -----
# Main
# -----

# 1. Get file & sample names

files <- list.files(sample.dirs, full.names = T)
files.nodir <- sapply(strsplit(files, '/'), function(x) x[length(x)])
mask <- rowSums(sapply(prefix, function(x) grepl(paste0('^', x), files.nodir))) > 0 & grepl('cov.gz', files.nodir) & !grepl('dedupl', files.nodir) & !grepl('lambda', files.nodir) & !grepl('trimmed', files.nodir)
files <- files[mask]
samples <- sapply(strsplit(files.nodir[mask], '.', fixed = T), function(x) x[1])
samples <- sapply(strsplit(samples, '_', fixed = T), function(x) x[1])
rm(files.nodir, mask)

dups <- !duplicated(samples)
files <- files[dups]
samples <- samples[dups]

print(samples)
print(length(samples))

# 2. Load files & apply filters

load.profile <- function(file){
  profile <- fread(cmd=paste0('gzip -dc ', file), stringsAsFactors = F)
  profile[,c(3,4)] <- NULL
  colnames(profile) <- c('chr', 'pos', 'M', 'm')
  profile$cov <- profile$M + profile$m
  profile <- profile[profile$cov >= min.CpG.depth,]
  if (chr.in.seqnames) profile$chr = gsub('chr', '', profile$chr)
  profile <- profile[profile$chr %in% chrs, ]
  return(profile)
}

cat(paste0('Defining CpG set ...\n'))

for (i in 1:length(samples)){
  sample = samples[i]
  cat(paste0(' -- loading sample ', sample, ' ...\n'))
  profile <- load.profile(files[i])
  profile[,3:5] <- NULL
  profile$count.new <- 1
  if (i == 1){
    me.matrix <- profile
    colnames(me.matrix)[3] <- 'count'
    cat(paste0(c('  --> current number of CpG islands = ', nrow(me.matrix), '\n')))
    next
  }
  me.matrix <- merge(me.matrix, profile, by = c('chr', 'pos'), all = TRUE)
  me.matrix <- cbind(me.matrix[,1:2], rowSums(me.matrix[,3:4], na.rm = T))
  colnames(me.matrix)[3] <- 'count'
  cat(paste0(c('  --> current number of CpG islands = ', nrow(me.matrix), '\n')))
}

cat(c(' -- Removing CpG islands with too many NAs ...\n'))
me.matrix <- me.matrix[which(me.matrix$count / length(samples) >= 1 - fraction.NA.allowed),]
cat(c('  --> final number of CpG islands = ', nrow(me.matrix), '\n'))

me.matrix$count <- NULL
CpG.annot <- me.matrix

cat(paste0('Creating raw CpG matrix ...\n'))
for (i in 1:length(samples)){
  sample = samples[i]
  cat(paste0(' -- reloading sample ', sample, ' ...\n'))
  profile <- load.profile(files[i])
  profile$m <- NULL
  colnames(profile)[3:4] <- c(sample, paste0(sample, '.coverage'))
  me.matrix <- merge(me.matrix, profile, by = c('chr', 'pos'), all.x = TRUE)
}

rm(sample, profile, sample.dirs, min.CpG.depth, fraction.NA.allowed, i, files, chr.in.seqnames, load.profile)

# 3. Write raw data

me.matrix[is.na(me.matrix)] <- 0
me.matrix <- as.data.frame(me.matrix, stringsAsFactors = F)

me.matrix.M <- me.matrix[, which(!grepl('coverage', colnames(me.matrix)))[-c(1,2)]]
me.matrix.coverage <- me.matrix[, grepl('coverage', colnames(me.matrix))]
colnames(me.matrix.coverage) <- samples

cat(c('Saving raw data as raw.Rds ...\n'))
bs.raw <- me.matrix.M / me.matrix.coverage
ids <- paste0(CpG.annot$chr, ':', CpG.annot$pos)
rownames(bs.raw) <- ids
saveRDS(round(bs.raw, 3), paste0(output.dir, '/raw.Rds'))

# 4. Create cluster

merge.sum <- function(sample, vec, cluster.gr){
  gr.object <- GRanges(cbind(CpG.annot, vec))
  agg <- aggregate(gr.object, hits, vec = sum(vec))
  elementMetadata(cluster.gr)[[sample]][countQueryHits(hits) > 0L] <- agg$vec
  return(cluster.gr)
}

merge.mean <- function(sample, vec, cluster.gr){
  gr.object <- GRanges(cbind(CpG.annot, vec))
  agg <- aggregate(gr.object, hits, vec = mean(vec))
  elementMetadata(cluster.gr)[[sample]][countQueryHits(hits) > 0L] <- agg$vec
  return(cluster.gr)
}

if(do.cluster){
  
  cat(c('Creating cluster ...\n'))
  
  cluster <- matrix(nrow = 0, ncol = 4)
  
  for (chr in chrs){
    pos <- CpG.annot$pos[CpG.annot$chr == chr]
    
    starts <- c(pos[1], pos[which(pos[-1] - pos[-length(pos)] > cluster.max.len.between.CpG) + 1])
    ends <- c(pos[which(pos[-1] - pos[-length(pos)] > cluster.max.len.between.CpG)], pos[length(pos)])
    amount <- diff(c(0, which(pos[-1] - pos[-length(pos)] > cluster.max.len.between.CpG), length(pos)))
    
    chr.cluster <- matrix(nrow = length(ends), ncol = 4)
    chr.cluster[,1] <- rep(chr, length(starts))
    chr.cluster[,2] <- starts
    chr.cluster[,3] <- ends
    chr.cluster[,4] <- amount
    
    cluster <- rbind(cluster, chr.cluster)
  }
  
  cat(paste0(c('  --> initial number of clusters = ', nrow(cluster), '\n')))
  
  cluster <- as.data.frame(cluster, stringsAsFactors = F)
  for (i in 2:4) cluster[,i] <- as.numeric(cluster[,i])
  colnames(cluster) <- c('chr', 'start', 'end', 'CpG')
  
  cat(paste0(' -- filtering CpG content ...\n'))
  cluster <- cluster[cluster$CpG >= cluster.min.CpG,]
  cat(paste0(c('  --> current number of clusters = ', nrow(cluster), '\n')))
  
  cat(paste0(' -- filtering CpG length ...\n'))
  cluster <- cluster[cluster$end - cluster$start + 1 >= cluster.min.length,]
  cat(paste0(c('  --> current number of clusters = ', nrow(cluster), '\n')))
  
  rownames(cluster) = paste0(cluster$chr, ':', cluster$start, '-', cluster$end)
  cluster.annot <- cluster
  cluster$CpG <- NULL
  
  # 4b. Map CpG to cluster
  
  cat(c('Mapping raw data to cluster ...\n'))
  
  cluster.gr <- GRanges(cluster)
  CpG.annot$start = CpG.annot$pos ; CpG.annot$end = CpG.annot$pos
  CpG.annot$pos <- NULL
  
  hits <- findOverlaps(cluster.gr, GRanges(CpG.annot))
  
  cluster.raw.coverage <- cluster.gr
  cluster.raw.M <- cluster.gr
  
  for (sample in samples){
    cat(paste0(' -- working on sample ', sample, ' ...\n'))
    i = which(samples == sample)
    cluster.raw.coverage <- merge.sum(sample, me.matrix.coverage[,i], cluster.raw.coverage)
    cluster.raw.M <- merge.sum(sample, me.matrix.M[,i], cluster.raw.M)
  }
  
  cluster.raw.M <- as.data.frame(elementMetadata(cluster.raw.M))
  rownames(cluster.raw.M) <- rownames(cluster)
  
  cluster.raw.coverage <- as.data.frame(elementMetadata(cluster.raw.coverage))
  rownames(cluster.raw.coverage) <- rownames(cluster)
  
  cluster.annot$mean.coverage <- rowMeans(cluster.raw.coverage)
  cluster.annot$fraction.not.missing <- rowSums(cluster.raw.coverage != 0) / length(samples)
  cluster.annot$mean.coverage.per.cpg <- cluster.annot$mean.coverage / cluster.annot$CpG
  
  cat(paste0(' -- filtering by minimum mean coverage of cluster ...\n'))
  mask1 <- cluster.annot$mean.coverage >= cluster.min.coverage
  cat(paste0(c('  --> current number of clusters = ', length(which(mask1)), '\n')))
  
  cat(paste0(' -- filtering by minimum mean coverage per CpG in cluster ...\n'))
  mask2 <- cluster.annot$mean.coverage.per.cpg >= cluster.min.cpg.coverage
  cat(paste0(c('  --> current number of clusters = ', length(which(mask1 & mask2)), '\n')))
  
  cat(paste0(' -- Removing clusters with too many NAs ...\n'))
  mask3 <- cluster.annot$fraction.not.missing >= 1 - cluster.fraction.NA.allowed
  mask <- mask1 & mask2 & mask3
  cat(c('  --> final number of clusters = ', length(which(mask)), '\n'))
  
  # 4c. Save raw cluster & CpG annotations
  
  saveRDS(cluster.annot[mask,], paste0(output.dir, '/cluster.annot.Rds'))
  if (mean.of.reads){
      saveRDS(cluster.raw.M[mask,] / cluster.raw.coverage[mask,], paste0(output.dir, '/cluster.raw.Rds'))
  }
  else {
      cluster.large.frame <- cluster.gr
      cat(c('Mapping raw data to cluster again, since --mean.of.reads F ...\n'))
      for (sample in samples){
        cat(paste0(' -- working on sample ', sample, ' ...\n'))
        i = which(samples == sample)
        cluster.large.frame <- merge.mean(sample, bs.raw[,i], cluster.large.frame)
      }
      cluster.large.frame <- as.data.frame(elementMetadata(cluster.large.frame))
      rownames(cluster.large.frame) <- rownames(cluster)
      saveRDS(cluster.large.frame[mask,], paste0(output.dir, '/cluster.raw.Rds'))
  }
  rm(bs.raw)
}

# 5. Smooth data

bs.object <- BSseq(chr = me.matrix$chr, pos = me.matrix$pos,
                   M = as.matrix(me.matrix.M), Cov = as.matrix(me.matrix.coverage))
rm(me.matrix)

be.smooth <- function(name, ns = 70, h = 1000){
  cat(paste0('Smoothening ', name, '...\n'))
  bs.smooth <- BSmooth(bs.object, ns = ns, h = h,
                       BPPARAM = BiocParallel::MulticoreParam(workers = n.cores, progressbar = T))
  cat(paste0('Saving as ', name,'.Rds ...\n'))
  large.frame <- as.data.frame(getMeth(bs.smooth))
  rownames(large.frame) <- ids
  saveRDS(round(large.frame, 3), paste0(output.dir, '/', name, '.Rds'))
  #Save the bsseq object as a HDF5Array 
  saveHDF5SummarizedExperiment(bs.smooth, dir = output.dir, replace = T) 
 
  if(do.cluster){
    # 5b. Map CpG to cluster
    
    cat(c('Mapping smoothened data to cluster ...\n'))
    
    cluster.large.frame <- cluster.gr
    
    for (sample in samples){
      cat(paste0(' -- working on sample ', sample, ' ...\n'))
      i = which(samples == sample)
      cluster.large.frame <- merge.mean(sample, large.frame[,i], cluster.large.frame)
    }
    
    cluster.large.frame <- as.data.frame(elementMetadata(cluster.large.frame))
    rownames(cluster.large.frame) <- rownames(cluster)
    
    saveRDS(cluster.large.frame[mask,], paste0(output.dir, '/cluster.', name,'.Rds'))
  }
}

be.smooth('light') # default parameters
#be.smooth('heavy', ns = 500, h = 20000) # 'block analysis', Wulfridge et al, 2019

cat('Done!\n')
