#all source functions for pre-processing (adapted to PCAWG)
# /data/beroukhim1/Shared_Data/PCAWG/TCGA/SVs/Delly

## using vcfs, processes, logs, scales, and then edits scaled file to "remove" breakpoint specific features 
### (homlen, insertionlen, hom_gc, insertion_gc)
## by shuffling the average values, then classifies with GaTSV

## DELLY detection, parse vcf to dt
library(data.table)

vcf_to_dt <- function(vcf_path, sample_name) {
  cat(paste0(vcf_path, "\n"))
  if (!file.exists(vcf_path)) {
    stop(paste("File does not exist", vcf_path))
  }
  
  cat("Reading file...\n")
  lines <- readLines(vcf_path)
  
  # ONLY skip lines starting with "##" so that data.table natively grabs the "#CHROM" header
  skip_lines <- sum(startsWith(lines, "##")) 
  vcf_dt <- data.table::fread(vcf_path, skip = skip_lines)
  
  if (nrow(vcf_dt) == 0) {
    return(vcf_dt)
  }
  
  # Safely rename ONLY the first 9 standard VCF columns, leaving the rest alone
  setnames(vcf_dt, 1:9, c("seqnames", "start", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "GENO"))
  
  # Handle the sample columns dynamically based on how many there are (10, 11, 13, etc.)
  sample_cols <- colnames(vcf_dt)[10:ncol(vcf_dt)]
  tumor_idx <- grep("tumor|somatic|primary", sample_cols, ignore.case = TRUE)
  normal_idx <- grep("control|normal|blood|germline", sample_cols, ignore.case = TRUE)
  
  if (length(tumor_idx) > 0) {
    vcf_dt$TUMOR <- vcf_dt[[ sample_cols[tumor_idx[1]] ]]
    if (length(normal_idx) > 0) {
      vcf_dt$NORMAL <- vcf_dt[[ sample_cols[normal_idx[1]] ]]
    }
  } else {
    # Fallback: if columns don't have recognizable names, assume standard positional layout
    if (length(sample_cols) >= 2) {
      vcf_dt$NORMAL <- vcf_dt[[ sample_cols[1] ]]
      vcf_dt$TUMOR <- vcf_dt[[ sample_cols[2] ]]
    } else if (length(sample_cols) == 1) {
      vcf_dt$TUMOR <- vcf_dt[[ sample_cols[1] ]]
    }
  }
  
  cat("Gathering Metadata...\n")
  vcf_dt$sample <- sample_name
  vcf_dt[, sid := sample_name]
  vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
  vcf_dt[, end := start] # Initialize end to start
  
  # ---------------------------------------------------------
  # DELLY PARSING LOGIC
  # ---------------------------------------------------------
  vcf_dt[, uid := ID]
  vcf_dt[, altchr := gsub(".*?;CHR2=([^;]+).*", "\\1", INFO)]
  vcf_dt[, altpos := as.integer(gsub(".*?;END=([0-9]+).*", "\\1", INFO))]
  vcf_dt[, SPAN := ifelse(seqnames == altchr, abs(altpos - start), -1)]
  
  # Delly connection types (CT): 3to5 (DEL), 5to3 (DUP), 3to3 (INV), 5to5 (INV)
  vcf_dt[, CT := gsub(".*?;CT=([^;]+).*", "\\1", INFO)]
  vcf_dt[, strand := fcase(CT == "3to5", "+", CT == "5to3", "-", CT == "3to3", "+", CT == "5to5", "-", default="*")]
  vcf_dt[, altstrand := fcase(CT == "3to5", "-", CT == "5to3", "+", CT == "3to3", "+", CT == "5to5", "-", default="*")]
  vcf_dt[, inv := strand == altstrand]
  
  # Read counts: FORMAT is GT:GL:GQ:FT:RC:DR:DV:RR:RV
  if ("TUMOR" %in% colnames(vcf_dt)) {
    tumor_splits <- strsplit(vcf_dt$TUMOR, ":", fixed = TRUE)
    vcf_dt[, TUMALT := as.integer(sapply(tumor_splits, `[`, 7)) + as.integer(sapply(tumor_splits, `[`, 9))]
    vcf_dt[, TUMCOV := as.integer(sapply(tumor_splits, `[`, 6)) + as.integer(sapply(tumor_splits, `[`, 8)) + TUMALT]
    vcf_dt[, TUMLOD := NA_real_] 
  }
  
  if ("NORMAL" %in% colnames(vcf_dt)) {
    normal_splits <- strsplit(vcf_dt$NORMAL, ":", fixed = TRUE)
    vcf_dt[, NORMALT := as.integer(sapply(normal_splits, `[`, 7)) + as.integer(sapply(normal_splits, `[`, 9))]
    vcf_dt[, NORMCOV := as.integer(sapply(normal_splits, `[`, 6)) + as.integer(sapply(normal_splits, `[`, 8)) + NORMALT]
    vcf_dt[, NORMLOD := NA_real_]
  }
  
  cat("Cleaning up...\n")
  canonical_contigs <- c(c(1:24), c('X','Y'), paste0('chr', c(1:24)), paste0('chr', c('X','Y')))
  bad.ix <- vcf_dt[!(seqnames %in% canonical_contigs), uid]
  vcf_dt <- vcf_dt[!uid %in% bad.ix]
  vcf_dt[, seqnames := ifelse(grepl('chr', seqnames), seqnames, paste0("chr", seqnames))]
  vcf_dt[, altchr := ifelse(grepl('chr', altchr), altchr, paste0("chr", altchr))]
  
  return(vcf_dt)
}

############################################################

build_bedpe_with_metadata <- function(merged_dt) {
  cat("Building bedpe...\n")
  if (nrow(merged_dt) == 0) return(data.table()) 
  
  bedpe <- data.table(
    chrom1 = merged_dt$seqnames,
    start1 = as.numeric(merged_dt$start),
    end1 = as.numeric(merged_dt$end),
    chrom2 = merged_dt$altchr,
    start2 = as.numeric(merged_dt$altpos),
    end2 = as.numeric(merged_dt$altpos),
    name = paste0(merged_dt$sample, "_", merged_dt$uid),
    score = merged_dt$QUAL,
    strand1 = merged_dt$strand,
    strand2 = merged_dt$altstrand,
    REF_1 = merged_dt$REF,
    ALT_1 = merged_dt$ALT,
    REF_2 = "N",          
    ALT_2 = merged_dt$ALT,
    SPAN = merged_dt$SPAN,
    FILTER = merged_dt$FILTER,
    sample = merged_dt$sample,
    truth = merged_dt$truth,        # <--- Added Truth Column
    TUMALT = merged_dt$TUMALT,
    TUMCOV = merged_dt$TUMCOV,
    GENO = merged_dt$GENO,
    TUMOR = merged_dt$TUMOR,
    MAPQ_1 = merged_dt$MAPQ,
    MAPQ_2 = merged_dt$MAPQ   
  )
  return(bedpe)
}

filter_paired <- function(germline_pth, somatic_pth, sample) {
  
  # Process Germline
  if (file.exists(germline_pth)) {
    germ_dt <- vcf_to_dt(germline_pth, sample)
    if (nrow(germ_dt) > 0) germ_dt[, truth := "GERMLINE"]
  } else { germ_dt <- data.table() }
  
  # Process Somatic
  if (file.exists(somatic_pth)) {
    som_dt <- vcf_to_dt(somatic_pth, sample)
    if (nrow(som_dt) > 0) som_dt[, truth := "SOMATIC"]
  } else { som_dt <- data.table() }
  
  # Combine
  vcf_dt <- rbind(germ_dt, som_dt, fill = TRUE)
  vcf_bedpe <- build_bedpe_with_metadata(vcf_dt)
  
  if (nrow(vcf_bedpe) == 0) return(vcf_bedpe)
  
  tmp <- vcf_bedpe[MAPQ_1 >= 60 | MAPQ_2 >= 60] 
  tmp <- tmp[TUMCOV > 1] 
  tmp <- tmp[SPAN > 49 | SPAN == -1] 
  
  return(tmp)
}

#ANNOTATION SCRIPTS
check_reformat = function(i, df){
  
  row <- df[i,]
  
  if (row$chrom1 > row$chrom2) {
    return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[,11:ncol(df)]))
  } else if (row$chrom1==row$chrom2 & row$start1>row$start2){
    return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[, 11:ncol(df)]))
  } else {
    return(row)
  }
}

fuzzy_filter_germline = function(itter = NULL, bed = NULL, g = NULL) {
  sub <- bed[itter,]
  ## reorder for filtering
  if(sub$chrom1 > sub$chrom2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else if (sub$chrom1 == sub$chrom2 & sub$start1 > sub$start2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else {
    sub_ord <- sub
  }
  ref_sub <- g[chrom1 == sub_ord$chrom1 & chrom2 == sub_ord$chrom2]
  ### calculate distances
  ref_sub[,str_dist := abs(start - as.numeric(as.character(sub_ord$start1)))]
  ref_sub[,end_dist := abs(end - as.numeric(as.character(sub_ord$start2)))]
  ref_sub[,gnomad_dist := (str_dist + end_dist)]
  ### choose closest match
  ref_min <- ref_sub[which.min(ref_sub$gnomad_dist)]
  if(nrow(ref_min) == 0) {
    sub <- cbind(sub, gnomad_dist = 2*1e9, gnomad_dist_avg = (1e9)/2 ,gnomad_d_bkpt1 = 1e9,gnomad_d_bkpt2 = 1e9) #remove the old gnomad_dist, add the new one,decouple bkpt distances
    return(sub)
  } else {
    sub <- cbind(sub, gnomad_dist = ref_min$gnomad_dist, gnomad_dist_avg = (ref_min$gnomad_dist)/2,gnomad_d_bkpt1 = ref_min$str_dist,gnomad_d_bkpt2 = ref_min$end_dist)
    return(sub)
  }
}

closest_germline = function(bp = NULL, cores = 1, genome = NULL) {
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  if(as.character(genome) == 'hg19') {
    gnomad_germline = gnomad_hg19
  } else if (as.character(genome) == 'hg38') {
    gnomad_germline = gnomad_hg38
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  
  cat("Comparing against known germline...")
  annotated_bedpe <- rbindlist(mclapply(1:nrow(bp), fuzzy_filter_germline, bp, g = gnomad_germline, mc.cores = cores))
  cat("done.\n")
  return(annotated_bedpe)
}

find_closest_match_line = function(i, bedpe_l, LINE_dt = NULL, LINE_dt_ranges = NULL){
  
  row_l <- bedpe_l[i,] #the particular row of the bedpe file we are working with
  
  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- LINE_dt[ seqnames == row_l$chrom1 | seqnames == row_l$chrom2]
  ref_sub_ranges <- LINE_dt_ranges[seqnames(LINE_dt_ranges) == row_l$chrom1 | seqnames(LINE_dt_ranges) == row_l$chrom2]
  
  row_l_str1 <- GRanges(row_l$chrom1, IRanges(as.numeric(as.character(row_l$start1)),width=1))
  row_l_str2 <- GRanges(row_l$chrom2, IRanges(as.numeric(as.character(row_l$start2)),width=1))
  
  check_str1_overlap <- ref_sub_ranges %&% row_l_str1
  check_str2_overlap <- ref_sub_ranges %&% row_l_str2
  
  ref_sub[,str1_dist_l := ifelse(length(check_str1_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_l$start1))), abs(end - as.numeric(as.character(row_l$start1)))))]
  ref_sub[,str2_dist_l := ifelse(length(check_str2_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_l$start2))), abs(end - as.numeric(as.character(row_l$start2)))))]
  
  str1_min <- min(ref_sub$str1_dist_l)
  str2_min <- min(ref_sub$str2_dist_l)
  line_dist <- str1_min+str2_min
  if(line_dist == Inf){ #modeling the case where there is no minimum;both str1_dist and str2_dist columns will be empty
    return(cbind(row_l, line_dist=""))
  } else {
    return (cbind(row_l, line_dist))
  }
}

find_closest_match_sine = function(i, bedpe_s, SINE_dt=NULL, SINE_dt_ranges = NULL){
  row_s <- bedpe_s[i,] 
  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- SINE_dt[ seqnames == row_s$chrom1 | seqnames == row_s$chrom2]
  ref_sub_ranges <- SINE_dt_ranges[seqnames(SINE_dt_ranges) == row_s$chrom1 | seqnames(SINE_dt_ranges) == row_s$chrom2]
  
  row_s_str1 <- GRanges(row_s$chrom1, IRanges(as.numeric(as.character(row_s$start1)),width=1))
  row_s_str2 <- GRanges(row_s$chrom2, IRanges(as.numeric(as.character(row_s$start2)),width=1))
  
  check_str1_overlap <- ref_sub_ranges %&% row_s_str1
  check_str2_overlap <- ref_sub_ranges %&% row_s_str2
  
  ref_sub[,str1_dist_s := ifelse(length(check_str1_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_s$start1))), abs(end - as.numeric(as.character(row_s$start1)))))]
  ref_sub[,str2_dist_s := ifelse(length(check_str2_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_s$start2))), abs(end - as.numeric(as.character(row_s$start2)))))]
  
  str1_min <- min(ref_sub$str1_dist_s)
  str2_min <- min(ref_sub$str2_dist_s)
  sine_dist <- str1_min+str2_min
  if(sine_dist == Inf){ #modeling the case where there is no minimum;both str1_dist and str2_dist columns will be empty
    return(cbind(row_s, sine_dist=""))
  }
  else {
    return (cbind(row_s, sine_dist))
  }
}

closest_line_sine = function(bp = NULL, genome = NULL, cores = 1) {
  
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  if(as.character(genome) == 'hg19') {
    LINE_dt = LINE_dt_hg19
    colnames(LINE_dt)[1:3] <- c('seqnames','start','end')
    LINE_dt_ranges = makeGRangesFromDataFrame(LINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
    SINE_dt = SINE_dt_hg19
    colnames(SINE_dt)[1:3] <- c('seqnames','start','end')
    SINE_dt_ranges = makeGRangesFromDataFrame(SINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
  } else if(as.character(genome) == 'hg38') {
    LINE_dt = LINE_dt_hg38
    LINE_dt_ranges = makeGRangesFromDataFrame(LINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
    SINE_dt = SINE_dt_hg38
    SINE_dt_ranges = makeGRangesFromDataFrame(SINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  
  cat("Checking format...")
  bp_ord <- rbindlist(mclapply(1:nrow(bp), check_reformat, bp, mc.cores = cores))
  cat("done.\n")
  bp_ord[chrom1 == 23, chrom1 := "X"]
  bp_ord[chrom2 == 23, chrom2 := "X"]
  bp_ord[chrom1 == 24, chrom1 := "Y"]
  bp_ord[chrom2 == 24, chrom2 := "Y"]
  cat("Comparing against LINE elements...")
  line_annotated <- rbindlist(mclapply(1:nrow(bp_ord), find_closest_match_line, bp_ord, LINE_dt = LINE_dt, LINE_dt_ranges = LINE_dt_ranges, mc.cores = cores))
  cat("done.\n")
  
  cat("Comparing against SINE elements...")
  sine_line_annotated <- rbindlist(mclapply(1:nrow(line_annotated), find_closest_match_sine, line_annotated, SINE_dt = SINE_dt, SINE_dt_ranges = SINE_dt_ranges, mc.cores = cores))
  cat("done.\n")
  
  return(sine_line_annotated)
}

find_closest_sv = function(i, bedpe_l){
  
  row_l <- bedpe_l[i,]
  
  ###create the reference dt here with the row in question removed
  bedpe_l_wo_row <- bedpe_l[-i,]
  sample_intra_events <- bedpe_l_wo_row[chrom1 == chrom2,]
  sample_inter_events <- bedpe_l_wo_row[chrom1 != chrom2,]
  ref_intra <- data.table(sample_intra_events$chrom1,sample_intra_events$start1,sample_intra_events$start2)
  ref_inter1 <-  data.table(sample_inter_events$chrom1,sample_inter_events$start1,sample_inter_events$end1)
  ref_inter2 <- data.table(sample_inter_events$chrom2,sample_inter_events$start2,sample_inter_events$end2)
  sample_ref_dt <- rbind(ref_intra,ref_inter1,ref_inter2) #this is a reference table of SVs that lacks the row in question
  colnames(sample_ref_dt) <- c('seqnames','start','end')
  
  row_l_str1 <- GRanges(row_l$chrom1, IRanges(as.numeric(as.character(row_l$start1)),width=1))
  row_l_str2 <- GRanges(row_l$chrom2, IRanges(as.numeric(as.character(row_l$start2)),width=1))
  
  if (row_l$chrom1 == row_l$chrom2){  
    ### both bedpe and ref should be sorted so lower bkpt comes first 
    ref_sub <- sample_ref_dt[seqnames == row_l$chrom1] #includes all events affecting that chromsome
    if (nrow(ref_sub) == 0){
      return(cbind(row_l, sv_dist=0))
    } else{
      ref_sub_ranges <- makeGRangesFromDataFrame(ref_sub, seqnames.field="seqnames", start.field="start", end.field ="end")
      
      check_str1_overlap <- ref_sub_ranges %&% row_l_str1
      check_str2_overlap <- ref_sub_ranges %&% row_l_str2
      
      
      ref_sub[,str1_dist_l := ifelse(length(check_str1_overlap)>=1, 0, min(abs(as.numeric(start) - as.numeric(as.character(row_l$start1))), abs(as.numeric(end) - as.numeric(as.character(row_l$start1)))))] #automatically returns the min value
      ref_sub[,str2_dist_l := ifelse(length(check_str2_overlap)>=1, 0, min(abs(as.numeric(start) - as.numeric(as.character(row_l$start2))), abs(as.numeric(end) - as.numeric(as.character(row_l$start2)))))]
    }
    
  } else { #if we are dealing with an interchromosomal event
    #work on the first chrom
    ref_sub_chr1 <- sample_ref_dt[seqnames == row_l$chrom1]
    if (nrow(ref_sub_chr1)==0){
      ref_sub_chr1 <- data.table(cbind("dummy","dummy","dummy")) #create a dummy table with at least 1 observation
      ref_sub_chr1[,str1_dist_l := 0]
    } else{
      ref_sub_ranges_chr1 <- makeGRangesFromDataFrame(ref_sub_chr1, seqnames.field="seqnames", start.field="start", end.field ="end")
      check_str1_overlap <- ref_sub_ranges_chr1 %&% row_l_str1
      ref_sub_chr1[,str1_dist_l := ifelse(length(check_str1_overlap)>=1, 0, min(abs(as.numeric(start) - as.numeric(as.character(row_l$start1))), abs(as.numeric(end) - as.numeric(as.character(row_l$start1)))))] #automatically returns the min value
    }
    #work on the second chromosome
    ref_sub_chr2 <- sample_ref_dt[seqnames == row_l$chrom2] #includes all events affecting that chromsome
    if (nrow(ref_sub_chr2)==0){
      ref_sub_chr2<- data.table(cbind("dummy","dummy","dummy"))
      ref_sub_chr2[,str2_dist_l := 0]
    }else{
      ref_sub_ranges_chr2 <- makeGRangesFromDataFrame(ref_sub_chr2, seqnames.field="seqnames", start.field="start", end.field ="end")
      check_str2_overlap <- ref_sub_ranges_chr2 %&% row_l_str2
      ref_sub_chr2[,str2_dist_l := ifelse(length(check_str2_overlap)>=1, 0, min(abs(as.numeric(start) - as.numeric(as.character(row_l$start2))), abs(as.numeric(end) - as.numeric(as.character(row_l$start2)))))]
    }
    ref_sub <- cbind(ref_sub_chr1[1:min(nrow(ref_sub_chr1),nrow(ref_sub_chr2))],ref_sub_chr2$str2_dist_l[1:min(nrow(ref_sub_chr1),nrow(ref_sub_chr2))])
    colnames(ref_sub) <- c('seqnames','start','end','str1_dist_l','str2_dist_l')
  }
  
  
  str1_min <- min(ref_sub$str1_dist_l)
  str2_min <- min(ref_sub$str2_dist_l)
  sv_dist <- str1_min+str2_min
  if(sv_dist == Inf){ #modeling the case where there is no minimum;both str1_dist and str2_dist columns will be empty
    return(cbind(row_l, sv_dist=0))
  } else {
    return (cbind(row_l, sv_dist))
  }
}

count_sv_5mbp = function (i,bedpe_l) {
  row_l <- bedpe_l[i,]
  bedpe_l_wo_row <- bedpe_l[-i,] #the rest of the data table to compare against
  row_dt <- rbind(data.table(row_l$chrom1,row_l$start1),data.table(row_l$chrom2,row_l$start2))
  colnames(row_dt) <- c('chrom','bkpt') #a datatable containing the breakpoints for that sv
  row_dt$bkpt <- as.numeric(row_dt$bkpt)
  bedpe_dt <- rbind(data.table(bedpe_l_wo_row$chrom1,bedpe_l_wo_row$start1),data.table(bedpe_l_wo_row$chrom2,bedpe_l_wo_row$start2))
  colnames(bedpe_dt) <- c('chrom','bkpt') #a table containing the search space of other SVs.
  bedpe_dt$bkpt <- as.numeric(bedpe_dt$bkpt)
  bedpe_dt_chrom1 <- subset(bedpe_dt,chrom == row_dt$chrom[1]) #subset the reference space by the chrom of interest
  find_sv_left <- bedpe_dt_chrom1[which(bedpe_dt_chrom1$bkpt > (row_dt$bkpt[1] - 2.5*10^6) & bedpe_dt_chrom1$bkpt < (row_dt$bkpt[1] + 2.5*10^6)),]
  bedpe_dt <- bedpe_dt[!(bedpe_dt$bkpt %in%find_sv_left$bkpt & bedpe_dt$chrom %in%find_sv_left$chrom),] #remove those SVs/breakpoints that have been counted already
  bedpe_dt_chrom2 <- subset(bedpe_dt,chrom == row_dt$chrom[2])
  find_sv_right <- bedpe_dt_chrom2[which(bedpe_dt_chrom2$bkpt > (as.numeric(row_dt$bkpt[2]) - 2.5*10^6) & bedpe_dt_chrom2$bkpt < (as.numeric(row_dt$bkpt[2]) + 2.5*10^6)),]
  sv_count_5Mbp <- sum(nrow(find_sv_left),nrow(find_sv_right))
  
  return (cbind(row_l, sv_count_5Mbp))
}

rep_time <- function(itter,bedpe,genome=NULL){
  row_l <- bedpe[itter,]
  row_l[,chrom1:=paste0('chr',chrom1)];row_l[,chrom2:=paste0('chr',chrom2)]
  
  if(as.character(genome)=='hg38'){ #if hg38, liftover SV breakpoints to hg19
    left_GRanges <- GRanges(row_l$chrom1,IRanges(as.numeric(as.character(row_l$start1)),width = 1))
    right_GRanges <- GRanges(row_l$chrom2,IRanges(as.numeric(as.character(row_l$start2)),width = 1))
    row_l_left = unlist(liftOver(left_GRanges, hg38ToHg19.chain))
    row_l_right = unlist(liftOver(right_GRanges, hg38ToHg19.chain))
  }
  else if (as.character(genome)=='hg19'){
    row_l_left <- GRanges(row_l$chrom1, IRanges(as.numeric(as.character(row_l$start1)),width=1))
    row_l_right <- GRanges(row_l$chrom2, IRanges(as.numeric(as.character(row_l$start2)),width=1))
  }else{
    stop('Please state hg19 or hg38 as genome')
  }
  
  sub_ranges_reptime_left <- reptimedata_hg19[seqnames(reptimedata_hg19) == row_l$chrom1]
  if (is_empty(sub_ranges_reptime_left)){
    sv_reptime_left <- 0
  }else{
    check_left_overlap <-  gr2dt(sub_ranges_reptime_left %&% row_l_left)
    check_left_overlap <- na.omit(check_left_overlap) #incase there are overlapping reference regions
    if (nrow(check_left_overlap)==0){ #if there are only NA values in the check_overlap datatable
      sv_reptime_left <- 0
    }else{
      if (nrow(check_left_overlap)>1){
        sv_reptime_left <- check_left_overlap[which.min(width),]$score[1] #if there are intervals that are tied for width, choose 1
      }else{
        sv_reptime_left <- check_left_overlap$score
      }
    }
  }
  
  sub_ranges_reptime_right <- reptimedata_hg19[seqnames(reptimedata_hg19) == row_l$chrom2]
  if (is_empty(sub_ranges_reptime_right)){
    sv_reptime_right <- 0
  }else{
    check_right_overlap <-  gr2dt(sub_ranges_reptime_right %&% row_l_right)
    check_right_overlap <- na.omit(check_right_overlap) #incase there are overlapping reference regions
    if (nrow(check_right_overlap)==0){ #if there are only NA values in the check_overlap datatable
      sv_reptime_right <- 0
    }else{
      if (nrow(check_right_overlap)>1){
        
        sv_reptime_right <- check_right_overlap[which.min(width),]$score[1] #if there are intervals that are tied for width, choose 1
      }else{
        sv_reptime_right <- check_right_overlap$score
      }
    }
    
  }
  row_l[,chrom1:=gsub('chr','',chrom1)];row_l[,chrom2:=gsub('chr','',chrom2)]
  return (cbind(row_l, sv_reptime_left, sv_reptime_right))    
}

annot_geneexon <- function(itter, bp, genome=NULL) {
  if(as.character(genome) == 'hg19') {
    hg_genes = hg19_genes
    hg_exons = hg19_exons
  } else if (as.character(genome) == 'hg38') {
    hg_genes = hg38_genes
    hg_exons = hg38_exons
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  sub <- bp[itter,]
  
  if(sub$chrom1 == sub$chrom2) {
    sub_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), as.numeric(as.character(sub$start2))))
    
    genes_overlapped <- hg_genes %&% sub_gr
    exons_overlapped <- hg_exons %&% sub_gr
    
    genes_fract <- genes_overlapped %O% sub_gr
    
    if(any(genes_fract == 1)) {
      sub[, CN_annot := 1]
    } else {
      sub[, CN_annot := 0]
    }
    
    if(length(exons_overlapped) > 0) {
      sub[,exon_annot := 1]
    } else {
      sub[,exon_annot := 0]
    }
    
  } else { #we expect that interchromosomal events do not overlap genes
    bp1_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), width = 1))
    bp2_gr <- GRanges(sub$chrom2, IRanges(as.numeric(as.character(sub$start2)), width = 1))
    
    exons_bp1 <- hg_exons %&% bp1_gr
    exons_bp2 <- hg_exons %&% bp2_gr
    
    exons_impacted <- length(exons_bp1) + length(exons_bp2)
    
    if(exons_impacted > 0) {
      sub[, CN_annot := 0]
      sub[,exon_annot := 1]
    } else {
      sub[, CN_annot := 0]
      sub[,exon_annot := 0]
    }
  }
  sub[,svtype := ifelse(chrom1 == chrom2, ifelse(strand1 == strand2, ifelse(strand1 == '+', "h2hINV", "t2tINV"), ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]
  return(sub)
}

check_tp53 <- function(bedpe){
  filename <- bedpe$sample[1] #unique identifier for the sample
  mut_status <- metadata$tp53_mutation_status[which(metadata$sample == filename)]
  return(bedpe[,tp53_status:=mut_status])
}

process_file <- function(germline_pth, somatic_pth, sample, n_cores, genome, output_path='./') {
  filt_file <- suppressWarnings(filter_paired(germline_pth, somatic_pth, sample))
  
  bedpe <- filt_file
  bedpe[, chrom1 := gsub('chr','',chrom1)]
  bedpe[, chrom2 := gsub('chr','',chrom2)]
  bedpe$start1 <- as.numeric(bedpe$start1)
  bedpe$start2 <- as.numeric(bedpe$start2)
  bedpe$end1 <- as.numeric(bedpe$end1)
  bedpe$end2 <- as.numeric(bedpe$end2)
  
  fuzzy <- closest_germline(bp = bedpe, cores = n_cores, genome = genome)
  
  line_sine <- closest_line_sine(bp = fuzzy, genome = genome, cores = n_cores)
  
  cat("Finding distance to closest SV... \n")
  nearest_sv_dist <- rbindlist(mclapply(1:nrow(line_sine), find_closest_sv, line_sine,mc.cores=n_cores))
  cat("done. \n")
  
  cat("Finding no of SVs in 5Mbp window... \n")
  sv_annotated <- rbindlist(mclapply(1:nrow(nearest_sv_dist), count_sv_5mbp, nearest_sv_dist, mc.cores = n_cores))
  cat("done. \n")
  
  cat("Adding replication timing info... \n")
  reptime_added <- rbindlist(mclapply(1:nrow(sv_annotated), rep_time, sv_annotated, genome=genome, mc.cores = n_cores))
  cat("done. \n")
  
  chroms <-  c(as.character(c(1:22)),'X','Y')
  reptime_added$chrom1 <- as.character(reptime_added$chrom1)
  reptime_added$chrom2 <- as.character(reptime_added$chrom2)
  store_indices = which(!reptime_added$chrom1 %in% chroms| !reptime_added$chrom2 %in% chroms) #done to address the observation that some SVs were mapping to hpv
  bedpe_clean <- reptime_added[!store_indices,]
  bedpe_clean$start1 <- as.numeric(bedpe_clean$start1)
  bedpe_clean$start2 <- as.numeric(bedpe_clean$start2)
  
  cat("Performing gene/exon annotation...")
  bedpe_annot <- rbindlist(mclapply(1:nrow(bedpe_clean), annot_geneexon, bedpe_clean, genome=genome,mc.cores = n_cores)) 
  cat("done. \n")
  
  cat("Checking TP53 status...")
  tp53_added <- check_tp53(bedpe_annot)
  cat("done. \n")
  
  filename <- tp53_added$sample[1]
  write.table(tp53_added, paste0(output_path,filename,'_processed.bedpe'), row.names = F, col.names = T, sep = "\t", quote = F)
  
  return(tp53_added)
}

add_last_feat <- function(df) {
  df <- as.data.table(df)
  
  df[, del:=ifelse(svtype=='DEL', 1, 0)]
  df[, dup:=ifelse(svtype=='DUP', 1, 0)]
  df[, inv:=ifelse(grepl(paste0('h2hINV', '|' ,'t2tINV'), svtype), 1, 0)]
  df[, inter:=ifelse(svtype=='INTER', 1, 0)]
  
  sample_svs <- as.data.table(table(df$sample)) #number of SVs in each sample 
  colnames(sample_svs) <- c('sample','num_sv_sample')
  df <- merge(df,sample_svs,by='sample')
  
  return(df)
}

features_tolog <- c('gnomad_d_bkpt1','gnomad_d_bkpt2','line_dist','sine_dist','num_sv_sample',
                    'sv_dist','sv_count_5Mbp')

log_feat <- function(x){
  log_x <- log(x, base = 10)
  log_x[which(is.na(log_x))] <- min(subset(log_x, !is.na(log_x)))
  log_x[which(is.infinite(log_x))] <- 0
  return((log_x))
}

features_toscale <- c('log_SPAN', 'log_gnomad_d_bkpt1', 'log_gnomad_d_bkpt2', 'del', 'dup', 'inv', 'inter',
                      'log_line_dist', 'log_sine_dist', 'log_num_sv_sample', 'CN_annot', 'exon_annot',
                      'log_sv_dist','log_sv_count_5Mbp', 'sv_reptime_left','sv_reptime_right','tp53_status')

run_GaTSV <- function(germline_pth, somatic_pth, sample, n_cores=1, genome='hg19', output_path='./'){
  cat(paste0('Reference genome: ',genome,'\n'))
  cat(paste0('Writing outputs to: ', output_path,'\n'))
  
  tmp <- process_file(germline_pth, somatic_pth, sample=sample, n_cores=n_cores, genome=genome, output_path=output_path)
  cutoff_prob <- 0.2684 #optimal tpr+ppv cutoff
  
  test <- add_last_feat(tmp)
  test <- test[SPAN>=1e3 | SPAN==-1,]
  test[, log_SPAN := log(SPAN, base = 10)] 
  test$log_SPAN[which(is.na(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.na(test$log_SPAN)))
  test$log_SPAN[which(is.infinite(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.infinite(test$log_SPAN)))
  
  test_sub_prelog <- test[, .SD, .SDcols = features_tolog]
  test_sub_log <- data.table(apply(test_sub_prelog, 2, FUN = log_feat))
  colnames(test_sub_log) <- paste0('log_', colnames(test_sub_log))
  
  # preserves column names
  test <- cbind(test, test_sub_log)
  test_df <- as.data.frame(test[, .SD, .SDcols=features_toscale]) #select the features to scale
  
  test_scaled <- data.table()
  for (i in colnames(test_df)){
    row_l <- scaling_mat[which(scaling_mat$feature == i),]
    feature_col <- test_df[, grepl(i, colnames(test_df))]
    scaled_feature <- (feature_col - (row_l$mean)) / row_l$sd
    test_scaled <- cbind(test_scaled, scaled_feature)
  }
  
  # na values are making the svm mad...
  test_scaled[is.na(test_scaled)] <- 0
  cat("Performing classification...\n")
  y_pred_radial <- predict(GaTSV, newdata = test_scaled, decision.values = T, probability = T)
  
  # Extract probabilities (this is a standard R matrix)
  prob_mat <- attr(y_pred_radial, 'probabilities')
  
  # Dynamically grab the column name that contains "1" (catches "1", "X1", etc.)
  pos_col <- colnames(prob_mat)[grepl("1", colnames(prob_mat))]
  
  # Fallback just in case: if it can't find a '1', use the second column
  if(length(pos_col) == 0) {
    pos_col <- colnames(prob_mat)[2]
  } else {
    pos_col <- pos_col[1] 
  }
  
  # Extract the raw numeric vector of probabilities for the SOMATIC class
  pos_probs <- prob_mat[, pos_col]
  
  # Calculate predictions directly onto your final 'test' data.table
  test[, predicted_class := ifelse(pos_probs >= cutoff_prob, 'SOMATIC', 'GERMLINE')]
  
  filename <- test$sample[1]
  write.table(test, paste0(output_path, filename, '_classified.bedpe'), row.names = F, col.names = T, sep = "\t", quote = F)
  cat('done.\n')
  
  return(test)
}