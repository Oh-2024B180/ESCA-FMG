###############################################################################
## Stage 1 + Stage 2: Somatic Variant QC & FMG Identification
##
## Inputs:
## data/variants_eso_final_masked.csv (CODA tumor-only Mutect2, AF>=0.05 applied)
## data/recurrent_site_rescue_list.csv (optional, auto-template if absent)
##
## Outputs (tables/, each <name>_v10_1.csv; OUT_PREFIX below):
## [1] FMG_list_primary (Stage 2 FMG list, primary), also written as Tab1_2_FMG_list_primary.csv (read by stages 03-07)
## [2] FMG_list_sensitivity (Stage 2 FMG list, sensitivity)
## [3] mutation_matrix_primary (sample x gene binary)
## [4] mutation_matrix_sensitivity
## [5] variants_somatic_primary (annotated variants, primary)
## [6] variants_somatic_sensitivity
## [7] gene_prevalence_primary (all genes, primary)
## [8] gene_prevalence_sensitivity
## [9] gene_support_qc (gene-level TxDb dependency)
## [10] variant_qc_flags (per-variant flags + reasons)
## [11] recurrent_site_review (review-target recurrent sites)
## [12] reviewworthy_variants_primary
##
## Cache:
## cache/cache_myvariant_stage01_02.rds
###############################################################################

# 0. Working directory + output directory
# run from the repository root (see README.md)

INPUT_FILE <- "data/variants_eso_final_masked.csv"
OUTDIR <- "tables"
LOGDIR <- "logs"
CACHE_DIR <- "cache"
RECURRENT_RESCUE_FILE <- "data/recurrent_site_rescue_list.csv"

for (d in c(OUTDIR, LOGDIR, CACHE_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 0.5. Reproducibility
set.seed(1262118388)

# 0.6. Logging start
log_file <- file.path(LOGDIR, "01_02_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE) # screen + file parallel record
sink(log_con, type = "message")

# 1.
required_pkgs <- c(
  "data.table", "httr", "jsonlite", "stringr",
  "GenomicRanges", "IRanges", "AnnotationDbi",
  "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    paste0(
      "following installbecome present nothave: ",
      paste(missing_pkgs, collapse = ", "),
      "\nfirst installagain executeplease \n",
      " if (!require('BiocManager')) install.packages('BiocManager')\n",
      " BiocManager::install(c('GenomicRanges','IRanges','AnnotationDbi',",
      "'org.Hs.eg.db','TxDb.Hsapiens.UCSC.hg19.knownGene'))\n",
      " install.packages(c('data.table','httr','jsonlite','stringr'))"
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  library(stringr)
  library(GenomicRanges)
  library(IRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})

options(timeout = 300)
cat("=== Stage 1+2: FMG Pipeline v10.1 start ===\n")
cat(sprintf("Run start: %s\n\n", Sys.time()))

# 2. Parameters
OUT_PREFIX <- "v10_1"
GENOME_BUILD <- "hg19"

# Stage 1D - population AF
GNOMAD_COMMON_CUTOFF <- 0.01 # >= 0.01: common germline -> exclude
GNOMAD_POSSIBLE_CUTOFF <- 0.001 # >= 0.001: possible germline -> exclude (primary only)

# Stage 1E - recurrent site
RECURRENCE_REVIEW_CUTOFF <- 5L
PRIMARY_DROP_UNRESCUED_RECURRENT_SITE <- FALSE

# Stage 2A - FMG cutoff
PREV_CUTOFF <- 0.05

# Stage 2D - MyVariant
MYVARIANT_BATCH <- 400L
MYVARIANT_FALLBACK <- 150L
MYVARIANT_MAX_RETRY <- 3L
MYVARIANT_CACHE_RDS <- file.path(CACHE_DIR, "cache_myvariant_stage01_02.rds")

# TxDb coding-only handling
KEEP_TXDB_CODING_ONLY <- FALSE

# gene-level txdb dependency QC
TXDB_ONLY_GENE_FRAC_CUTOFF <- 0.50
TXDB_ONLY_GENE_MIN_SAMPLES <- 5L
WHITELIST_TXDB_GENES <- character(0)

# 3. Helper functions
ensure_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

ensure_chr_prefix <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("M", "MT", "Mt")] <- "M"
  paste0("chr", x)
}

collapse_atomic_strings <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x) & x != "." & x != "NA" & x != "NULL"]
  if (!length(x)) return(NA_character_)
  paste(unique(x), collapse = ";")
}

collapse_group_strings <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  vals <- unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  vals <- vals[nzchar(vals) & vals != "." & vals != "NA" & vals != "NULL"]
  if (!length(vals)) return(NA_character_)
  paste(unique(vals), collapse = ";")
}

max_atomic_numeric <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_first_char <- function(x) {
  if (is.null(x)) return(character())
  if (methods::is(x, "AtomicList") || methods::is(x, "List")) {
    return(vapply(x, function(el) {
      if (length(el) == 0L) NA_character_ else as.character(el[[1]])
    }, character(1)))
  }
  as.character(x)
}

first_token <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  z <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  z <- trimws(z)
  z <- z[nzchar(z)]
  if (!length(z)) return(NA_character_)
  z[1]
}

extract_col_as_char <- function(dt, col) {
  if (!col %in% names(dt)) return(rep(NA_character_, nrow(dt)))
  x <- dt[[col]]
  if (is.list(x)) {
    vapply(x, collapse_atomic_strings, character(1))
  } else {
    vapply(as.list(x), collapse_atomic_strings, character(1))
  }
}

extract_col_as_num <- function(dt, col) {
  if (!col %in% names(dt)) return(rep(NA_real_, nrow(dt)))
  x <- dt[[col]]
  if (is.list(x)) {
    vapply(x, max_atomic_numeric, numeric(1))
  } else {
    vapply(as.list(x), max_atomic_numeric, numeric(1))
  }
}

split_tokens <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  x <- gsub(",", ";", x, fixed = TRUE)
  z <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  z <- trimws(z)
  unique(z[nzchar(z)])
}

normalize_consequence_token <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub(" ", "_", x0)
  
  map <- c(
    "missense" = "missense_variant",
    "nonsense" = "stop_gained",
    "frameshift" = "frameshift_variant",
    "splice_donor" = "splice_donor_variant",
    "splice_acceptor" = "splice_acceptor_variant",
    "splice_site" = "splice_region_variant",
    "splicesite" = "canonical_splice_variant_txdb",
    "protein_altering" = "protein_altering_variant",
    "initiator_codon" = "initiator_codon_variant",
    "non_synonymous" = "missense_variant",
    "non_synonymous_coding" = "missense_variant",
    "frame_shift" = "frameshift_variant",
    "canonical_splice" = "canonical_splice_variant",
    "coding" = "coding_sequence_variant_txdb",
    "coding_region" = "coding_sequence_variant_txdb",
    "gene_body" = "gene_body_variant_txdb"
  )
  
  if (x0 %in% names(map)) return(unname(map[[x0]]))
  x0
}

pick_best_consequence <- function(x) {
  toks <- split_tokens(x)
  if (!length(toks)) return(NA_character_)
  toks <- vapply(toks, normalize_consequence_token, character(1))
  toks <- unique(toks)
  
  rank_tbl <- c(
    "splice_acceptor_variant" = 1L,
    "splice_donor_variant" = 1L,
    "canonical_splice_variant" = 1L,
    "canonical_splice_variant_txdb" = 1L,
    "stop_gained" = 2L,
    "frameshift_variant" = 3L,
    "stop_lost" = 4L,
    "start_lost" = 4L,
    "initiator_codon_variant" = 4L,
    "inframe_insertion" = 5L,
    "inframe_deletion" = 5L,
    "inframe_indel" = 5L,
    "missense_variant" = 6L,
    "protein_altering_variant" = 7L,
    "splice_region_variant" = 8L,
    "synonymous_variant" = 20L,
    "coding_sequence_variant_txdb" = 30L,
    "gene_body_variant_txdb" = 40L
  )
  
  rk <- rank_tbl[toks]
  rk[is.na(rk)] <- 999L
  toks[order(rk, toks)][1]
}

is_functional_primary <- function(x) {
  if (is.na(x) || !nzchar(x)) return(FALSE)
  x <- normalize_consequence_token(x)
  x %in% c(
    "splice_acceptor_variant",
    "splice_donor_variant",
    "canonical_splice_variant",
    "canonical_splice_variant_txdb",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "initiator_codon_variant",
    "inframe_insertion",
    "inframe_deletion",
    "inframe_indel",
    "missense_variant",
    "protein_altering_variant"
  )
}

build_hgvs <- function(chrom, pos, ref, alt) {
  n <- length(chrom)
  hgvs <- rep(NA_character_, n)
  
  for (i in seq_len(n)) {
    ch <- ensure_chr_prefix(chrom[i])
    p <- suppressWarnings(as.integer(pos[i]))
    r <- as.character(ref[i])
    a <- as.character(alt[i])
    
    if (is.na(ch) || is.na(p) || is.na(r) || is.na(a)) next
    if (grepl(",", a, fixed = TRUE)) next
    
    lr <- nchar(r)
    la <- nchar(a)
    
    if (lr == 1L && la == 1L) {
      hgvs[i] <- sprintf("%s:g.%d%s>%s", ch, p, r, a)
    } else if (lr == 1L && la > 1L) {
      ins <- substring(a, 2L)
      if (nchar(ins) <= 100L) {
        hgvs[i] <- sprintf("%s:g.%d_%dins%s", ch, p, p + 1L, ins)
      }
    } else if (lr > 1L && la == 1L) {
      ds <- p + 1L
      de <- p + lr - 1L
      dseq <- substring(r, 2L)
      if (nchar(dseq) <= 500L) {
        if (ds == de) {
          hgvs[i] <- sprintf("%s:g.%ddel%s", ch, ds, dseq)
        } else {
          hgvs[i] <- sprintf("%s:g.%d_%ddel%s", ch, ds, de, dseq)
        }
      }
    } else {
      if (max(lr, la) <= 100L) {
        hgvs[i] <- sprintf("%s:g.%d_%ddelins%s", ch, p, p + lr - 1L, a)
      }
    }
  }
  
  hgvs
}

query_myvariant <- function(hgvs_ids, fields,
                            batch_size = 400L,
                            fallback_size = 150L,
                            max_retry = 3L) {
  if (length(hgvs_ids) == 0L) return(data.table())
  
  results <- list()
  ridx <- 1L
  n_batch <- ceiling(length(hgvs_ids) / batch_size)
  
  for (i in seq_len(n_batch)) {
    s <- (i - 1L) * batch_size + 1L
    e <- min(i * batch_size, length(hgvs_ids))
    batch <- hgvs_ids[s:e]
    
    cat(sprintf(" Batch %d/%d (%d ids)...", i, n_batch, length(batch)))
    ok <- FALSE
    
    for (att in seq_len(max_retry)) {
      resp <- tryCatch(
        POST(
          "https://myvariant.info/v1/variant",
          body = list(ids = paste(batch, collapse = ","), fields = fields),
          encode = "form",
          timeout(180)
        ),
        error = function(e) NULL
      )
      
      if (!is.null(resp) && httr::status_code(resp) == 200L) {
        txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        obj <- fromJSON(txt, flatten = TRUE)
        results[[ridx]] <- as.data.table(obj)
        ridx <- ridx + 1L
        ok <- TRUE
        cat(" OK\n")
        break
      }
      
      cat(sprintf(" [retry %d]", att))
      Sys.sleep(att * 2)
    }
    
    if (!ok) {
      cat(sprintf(" -> fallback (%d)\n", fallback_size))
      n_sub <- ceiling(length(batch) / fallback_size)
      
      for (j in seq_len(n_sub)) {
        ss <- (j - 1L) * fallback_size + 1L
        se <- min(j * fallback_size, length(batch))
        sub_batch <- batch[ss:se]
        
        cat(sprintf(" Sub %d/%d (%d ids)...", j, n_sub, length(sub_batch)))
        sub_ok <- FALSE
        
        for (att in seq_len(max_retry)) {
          resp2 <- tryCatch(
            POST(
              "https://myvariant.info/v1/variant",
              body = list(ids = paste(sub_batch, collapse = ","), fields = fields),
              encode = "form",
              timeout(180)
            ),
            error = function(e) NULL
          )
          
          if (!is.null(resp2) && httr::status_code(resp2) == 200L) {
            txt2 <- httr::content(resp2, as = "text", encoding = "UTF-8")
            obj2 <- fromJSON(txt2, flatten = TRUE)
            results[[ridx]] <- as.data.table(obj2)
            ridx <- ridx + 1L
            sub_ok <- TRUE
            cat(" OK\n")
            break
          }
          
          cat(sprintf(" [retry %d]", att))
          Sys.sleep(att * 2)
        }
        
        if (!sub_ok) cat(" FAIL\n")
        Sys.sleep(0.3)
      }
    }
    
    Sys.sleep(0.7)
  }
  
  if (!length(results)) return(data.table())
  rbindlist(results, fill = TRUE, use.names = TRUE)
}

standardize_myvariant_result <- function(dt) {
  if (nrow(dt) == 0L) return(data.table())
  
  dt <- as.data.table(dt)
  
  if ("query" %in% names(dt) && !"hgvs_id" %in% names(dt)) {
    setnames(dt, "query", "hgvs_id")
  }
  if ("_id" %in% names(dt) && !"hgvs_id" %in% names(dt)) {
    setnames(dt, "_id", "hgvs_id")
  }
  if ("notfound" %in% names(dt)) {
    dt <- dt[is.na(notfound) | notfound != TRUE]
  }
  if (!"hgvs_id" %in% names(dt)) {
    stop("MyVariant resulthgvs_id/query columnhave.")
  }
  
  unique(dt, by = "hgvs_id")
}

get_myvariant_with_cache <- function(hgvs_ids, fields, cache_file,
                                     batch_size = 400L,
                                     fallback_size = 150L,
                                     max_retry = 3L) {
  ensure_dir(dirname(cache_file))
  
  cache_dt <- if (file.exists(cache_file)) {
    readRDS(cache_file)
  } else {
    data.table(hgvs_id = character())
  }
  cache_dt <- as.data.table(cache_dt)
  
  if (!"hgvs_id" %in% names(cache_dt)) {
    cache_dt <- data.table(hgvs_id = character())
  }
  
  missing_ids <- setdiff(hgvs_ids, cache_dt$hgvs_id)
  
  cat(sprintf(" Cached annotations: %d\n", nrow(cache_dt)))
  cat(sprintf(" Missing annotations to query: %d\n", length(missing_ids)))
  
  if (length(missing_ids) > 0L) {
    new_dt <- query_myvariant(
      hgvs_ids = missing_ids,
      fields = fields,
      batch_size = batch_size,
      fallback_size = fallback_size,
      max_retry = max_retry
    )
    new_dt <- standardize_myvariant_result(new_dt)
    
    if (nrow(new_dt) > 0L) {
      cache_dt <- rbindlist(list(cache_dt, new_dt), fill = TRUE, use.names = TRUE)
      cache_dt <- unique(cache_dt, by = "hgvs_id")
      saveRDS(cache_dt, cache_file)
    }
  }
  
  cache_dt[hgvs_id %in% hgvs_ids]
}

make_splice_windows <- function(exons_by_gene, flank = 2L) {
  out <- list()
  idx <- 1L
  
  for (gid in names(exons_by_gene)) {
    gr <- exons_by_gene[[gid]]
    if (length(gr) == 0L) next
    
    left <- GRanges(seqnames(gr), IRanges(pmax(start(gr) - flank, 1L), start(gr) + flank - 1L))
    right <- GRanges(seqnames(gr), IRanges(pmax(end(gr) - flank + 1L, 1L), end(gr) + flank))
    
    s <- c(left, right)
    mcols(s)$gene_id <- gid
    out[[idx]] <- s
    idx <- idx + 1L
  }
  
  if (!length(out)) return(GRanges())
  suppressWarnings(do.call(c, out))
}

fallback_txdb_overlap <- function(miss_dt, txdb) {
  if (nrow(miss_dt) == 0L) return(data.table())
  
  seqinf <- seqinfo(txdb)
  valid_seq <- names(seqlengths(seqinf))
  
  dt <- copy(miss_dt)
  dt[, chrom := ensure_chr_prefix(chrom)]
  dt <- dt[chrom %in% valid_seq]
  if (nrow(dt) == 0L) return(data.table())
  
  maxlen <- seqlengths(seqinf)[dt$chrom]
  rw <- pmax(nchar(dt$ref), 1L)
  keep <- !is.na(maxlen) & !is.na(dt$pos) & dt$pos >= 1L & (dt$pos + rw - 1L) <= maxlen
  dt <- dt[keep]
  rw <- rw[keep]
  if (nrow(dt) == 0L) return(data.table())
  
  gr_var <- GRanges(dt$chrom, IRanges(dt$pos, width = rw))
  mcols(gr_var)$hgvs_id <- dt$hgvs_id
  
  cds_gr <- cds(txdb, columns = c("gene_id"))
  cds_hits <- findOverlaps(gr_var, cds_gr, ignore.strand = TRUE)
  
  cds_dt <- if (length(cds_hits) > 0L) {
    data.table(
      hgvs_id = as.character(mcols(gr_var)$hgvs_id[queryHits(cds_hits)]),
      gene_id = safe_first_char(mcols(cds_gr)$gene_id[subjectHits(cds_hits)]),
      txdb_region = "coding_region"
    )
  } else {
    data.table(hgvs_id = character(), gene_id = character(), txdb_region = character())
  }
  
  exons_by_gene <- exonsBy(txdb, by = "gene")
  splice_gr <- make_splice_windows(exons_by_gene, flank = 2L)
  splice_dt <- data.table(hgvs_id = character(), gene_id = character(), txdb_region = character())
  
  if (length(splice_gr) > 0L) {
    splice_hits <- findOverlaps(gr_var, splice_gr, ignore.strand = TRUE)
    if (length(splice_hits) > 0L) {
      splice_dt <- data.table(
        hgvs_id = as.character(mcols(gr_var)$hgvs_id[queryHits(splice_hits)]),
        gene_id = safe_first_char(mcols(splice_gr)$gene_id[subjectHits(splice_hits)]),
        txdb_region = "spliceSite"
      )
    }
  }
  
  gene_gr <- suppressWarnings(genes(txdb))
  gene_hits <- findOverlaps(gr_var, gene_gr, ignore.strand = TRUE)
  gene_dt <- if (length(gene_hits) > 0L) {
    data.table(
      hgvs_id = mcols(gr_var)$hgvs_id[queryHits(gene_hits)],
      gene_id = names(gene_gr)[subjectHits(gene_hits)],
      txdb_region = "gene_body"
    )
  } else {
    data.table(hgvs_id = character(), gene_id = character(), txdb_region = character())
  }
  
  ann_dt <- unique(rbindlist(list(splice_dt, cds_dt, gene_dt), fill = TRUE))
  if (nrow(ann_dt) == 0L) return(data.table())
  
  ann_dt[, loc_pri := fifelse(txdb_region == "spliceSite", 1L,
                              fifelse(txdb_region == "coding_region", 2L, 9L))]
  
  gids <- unique(ann_dt$gene_id[!is.na(ann_dt$gene_id) & ann_dt$gene_id != ""])
  sym <- if (length(gids) > 0L) {
    as.data.table(AnnotationDbi::select(
      org.Hs.eg.db,
      keys = as.character(gids),
      columns = "SYMBOL",
      keytype = "ENTREZID"
    ))
  } else {
    data.table(ENTREZID = character(), SYMBOL = character())
  }
  
  ann_dt <- merge(ann_dt, sym, by.x = "gene_id", by.y = "ENTREZID", all.x = TRUE)
  
  ann_dt[order(hgvs_id, loc_pri),
         .SD[1],
         by = hgvs_id][
           , txdb_consequence := fifelse(
             txdb_region == "spliceSite", "canonical_splice_variant_txdb",
             fifelse(txdb_region == "coding_region", "coding_sequence_variant_txdb",
                     "gene_body_variant_txdb")
           )][
             , .(hgvs_id, txdb_gene = SYMBOL, txdb_region, txdb_consequence)]
}

build_matrix_and_prevalence <- function(var_dt, all_samples, prevalence_cutoff = 0.05) {
  if (nrow(var_dt) == 0L) {
    mat_df <- data.frame(row.names = all_samples)
    gprev <- data.table(gene = character(), n_mutated = integer(), prevalence = numeric())
    fmg <- copy(gprev)
    return(list(mat = mat_df, prevalence = gprev, fmg = fmg))
  }
  
  sg <- unique(var_dt[, .(sample_id_masked, genename)], by = c("sample_id_masked", "genename"))
  sg[, mutated := 1L]
  
  mat <- dcast(sg, sample_id_masked ~ genename, value.var = "mutated", fill = 0L)
  mat_df <- as.data.frame(mat)
  rownames(mat_df) <- mat_df$sample_id_masked
  mat_df$sample_id_masked <- NULL
  
  missing_samples <- setdiff(all_samples, rownames(mat_df))
  if (length(missing_samples) > 0L) {
    zero_block <- matrix(
      0L,
      nrow = length(missing_samples),
      ncol = ncol(mat_df),
      dimnames = list(missing_samples, colnames(mat_df))
    )
    mat_df <- rbind(mat_df, as.data.frame(zero_block))
  }
  
  mat_df <- mat_df[order(rownames(mat_df)), , drop = FALSE]
  
  if (ncol(mat_df) == 0L) {
    gprev <- data.table(gene = character(), n_mutated = integer(), prevalence = numeric())
    fmg <- copy(gprev)
    return(list(mat = mat_df, prevalence = gprev, fmg = fmg))
  }
  
  gprev <- data.table(
    gene = colnames(mat_df),
    n_mutated = as.integer(colSums(mat_df)),
    prevalence = as.numeric(colSums(mat_df) / length(all_samples))
  )
  setorder(gprev, -prevalence, -n_mutated, gene)
  
  fmg <- gprev[prevalence >= prevalence_cutoff]
  
  list(mat = mat_df, prevalence = gprev, fmg = fmg)
}

write_matrix_csv <- function(mat_df, file) {
  out <- cbind(sample_id_masked = rownames(mat_df), mat_df)
  write.csv(out, file, row.names = FALSE)
}

prepare_recurrent_rescue_file <- function(file) {
  if (!file.exists(file)) {
    tmpl <- data.table(
      hgvs_id = character(),
      rescue_reason = character()
    )
    fwrite(tmpl, file)
    cat(sprintf(" Rescue template created: %s\n", file))
  }
}

read_recurrent_rescue <- function(file) {
  if (!file.exists(file)) {
    return(data.table(hgvs_id = character(), rescue_reason = character()))
  }
  
  x <- fread(file)
  if (nrow(x) == 0L) {
    return(data.table(hgvs_id = character(), rescue_reason = character()))
  }
  
  if (!"rescue_reason" %in% names(x)) {
    x[, rescue_reason := "manual_rescue"]
  }
  
  if (!"hgvs_id" %in% names(x)) {
    req <- c("chrom", "pos", "ref", "alt")
    if (!all(req %in% names(x))) {
      stop("recurrent_site_rescue_list.csvhgvs_id chrom,pos,ref,alt ")
    }
    x[, hgvs_id := build_hgvs(chrom, pos, ref, alt)]
  }
  
  x <- x[!is.na(hgvs_id) & nzchar(hgvs_id)]
  unique(x[, .(hgvs_id, rescue_reason)], by = "hgvs_id")
}

# 4. Data load
cat(">>> STEP 1: Data load\n")

if (!file.exists(INPUT_FILE)) {
  stop(sprintf("input file have: %s\ncurrent Working directory: %s",
               INPUT_FILE, getwd()))
}

raw <- fread(INPUT_FILE)

if ("filter" %in% names(raw)) {
  raw <- raw[is.na(filter) | filter %in% c("PASS", ".")]
}
if ("FILTER" %in% names(raw)) {
  raw <- raw[is.na(FILTER) | FILTER %in% c("PASS", ".")]
}

req_cols <- c("sample_id_masked", "chrom", "pos", "ref", "alt")
miss_cols <- setdiff(req_cols, names(raw))
if (length(miss_cols) > 0L) {
  stop(sprintf("required column missing: %s", paste(miss_cols, collapse = ", ")))
}

raw[, chrom := ensure_chr_prefix(chrom)]
raw[, pos := as.integer(pos)]

if ("AF" %in% names(raw)) {
  raw[, AF := suppressWarnings(as.numeric(AF))]
} else {
  raw[, AF := NA_real_]
}

raw_sv <- unique(
  raw[, .(sample_id_masked, chrom, pos, ref, alt, AF)],
  by = c("sample_id_masked", "chrom", "pos", "ref", "alt")
)

N_TOTAL <- uniqueN(raw_sv$sample_id_masked)
ALL_SAMPLES <- sort(unique(raw_sv$sample_id_masked))

cat(sprintf(" %d sample-variant rows, %d samples\n", nrow(raw_sv), N_TOTAL))

pos_stats <- raw_sv[, .(
  n_samples = uniqueN(sample_id_masked),
  mean_tumor_af = suppressWarnings(mean(AF, na.rm = TRUE)),
  median_tumor_af = suppressWarnings(median(AF, na.rm = TRUE))
), by = .(chrom, pos, ref, alt)]

pos_stats[!is.finite(mean_tumor_af), mean_tumor_af := NA_real_]
pos_stats[!is.finite(median_tumor_af), median_tumor_af := NA_real_]

uv_all <- unique(raw_sv[, .(chrom, pos, ref, alt)], by = c("chrom", "pos", "ref", "alt"))
uv_all[, hgvs_id := build_hgvs(chrom, pos, ref, alt)]
n_hgvs_fail <- sum(is.na(uv_all$hgvs_id))

uv <- uv_all[!is.na(hgvs_id)]
cat(sprintf(" Unique variants: %d (HGVS OK), %d HGVS build fail\n", nrow(uv), n_hgvs_fail))

# 5. MyVariant annotation with cache
cat("\n>>> STEP 2: MyVariant.info annotation (+ cache)\n")

fields <- paste(
  c(
    "dbnsfp.genename",
    "cadd.gene.genename",
    "cadd.consequence",
    "dbnsfp.annovar.effect",
    "snpeff.ann.effect",
    "gnomad_genome.af.af",
    "gnomad_genome.af.af_eas",
    "gnomad_exome.af.af",
    "gnomad_exome.af.af_eas"
  ),
  collapse = ","
)

annot_raw <- get_myvariant_with_cache(
  hgvs_ids = uv$hgvs_id,
  fields = fields,
  cache_file = MYVARIANT_CACHE_RDS,
  batch_size = MYVARIANT_BATCH,
  fallback_size = MYVARIANT_FALLBACK,
  max_retry = MYVARIANT_MAX_RETRY
)

annot_sum <- data.table(
  hgvs_id = character(),
  all_genename = character(),
  all_consequence_raw = character(),
  gnomad_genome_af = numeric(),
  gnomad_genome_af_eas = numeric(),
  gnomad_exome_af = numeric(),
  gnomad_exome_af_eas = numeric()
)

if (nrow(annot_raw) > 0L) {
  annot_raw <- as.data.table(annot_raw)
  
  annot_raw[, gene_mv1 := extract_col_as_char(.SD, "dbnsfp.genename")]
  annot_raw[, gene_mv2 := extract_col_as_char(.SD, "cadd.gene.genename")]
  annot_raw[, csq_mv1 := extract_col_as_char(.SD, "cadd.consequence")]
  annot_raw[, csq_mv2 := extract_col_as_char(.SD, "dbnsfp.annovar.effect")]
  annot_raw[, csq_mv3 := extract_col_as_char(.SD, "snpeff.ann.effect")]
  
  annot_raw[, gnomad_genome_af := extract_col_as_num(.SD, "gnomad_genome.af.af")]
  annot_raw[, gnomad_genome_af_eas := extract_col_as_num(.SD, "gnomad_genome.af.af_eas")]
  annot_raw[, gnomad_exome_af := extract_col_as_num(.SD, "gnomad_exome.af.af")]
  annot_raw[, gnomad_exome_af_eas := extract_col_as_num(.SD, "gnomad_exome.af.af_eas")]
  
  annot_sum <- annot_raw[, .(
    all_genename = collapse_group_strings(c(gene_mv1, gene_mv2)),
    all_consequence_raw = collapse_group_strings(c(csq_mv1, csq_mv2, csq_mv3)),
    gnomad_genome_af = max_or_na(gnomad_genome_af),
    gnomad_genome_af_eas = max_or_na(gnomad_genome_af_eas),
    gnomad_exome_af = max_or_na(gnomad_exome_af),
    gnomad_exome_af_eas = max_or_na(gnomad_exome_af_eas)
  ), by = hgvs_id]
}

cat(sprintf(" MyVariant annotated variants: %d / %d\n", nrow(annot_sum), nrow(uv)))

# 6. TxDb fallback
cat("\n>>> STEP 2-b: TxDb fallback\n")

uv_annot <- merge(uv, annot_sum, by = "hgvs_id", all.x = TRUE)
missing_gene <- uv_annot[is.na(all_genename)]
cat(sprintf(" Gene-missing variants for TxDb fallback: %d\n", nrow(missing_gene)))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txr <- fallback_txdb_overlap(missing_gene[, .(hgvs_id, chrom, pos, ref, alt)], txdb)

if (nrow(txr) > 0L) {
  uv_annot <- merge(uv_annot, txr, by = "hgvs_id", all.x = TRUE)
} else {
  uv_annot[, `:=`(
    txdb_gene = NA_character_,
    txdb_region = NA_character_,
    txdb_consequence = NA_character_
  )]
}

uv_annot[is.na(all_genename) & !is.na(txdb_gene), all_genename := txdb_gene]

uv_annot[, genename := vapply(all_genename, first_token, character(1))]
uv_annot[, flag_multigene_annotation := !is.na(all_genename) & grepl(";", all_genename, fixed = TRUE)]

uv_annot[, all_consequence := collapse_group_strings(c(all_consequence_raw, txdb_consequence)), by = hgvs_id]
uv_annot[, picked_consequence := vapply(all_consequence, pick_best_consequence, character(1))]

uv_annot[, annot_source := fifelse(
  !is.na(all_consequence_raw) & !is.na(txdb_consequence), "myvariant+txdb",
  fifelse(!is.na(all_consequence_raw), "myvariant",
          fifelse(!is.na(txdb_consequence), "txdb", NA_character_))
)]

uv_annot <- merge(uv_annot, pos_stats, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)

uv_annot[, gnomad_af_max := apply(
  cbind(gnomad_genome_af, gnomad_genome_af_eas, gnomad_exome_af, gnomad_exome_af_eas),
  1,
  function(z) {
    z <- suppressWarnings(as.numeric(z))
    z <- z[is.finite(z)]
    if (!length(z)) NA_real_ else max(z)
  }
)]

cat(sprintf(" gene annotation: %d / %d (%.1f%%)\n",
            sum(!is.na(uv_annot$genename)),
            nrow(uv_annot),
            100 * sum(!is.na(uv_annot$genename)) / nrow(uv_annot)))

# 7. Rescue file
cat("\n>>> STEP 2-c: recurrent rescue file n")

prepare_recurrent_rescue_file(RECURRENT_RESCUE_FILE)
rescue_dt <- read_recurrent_rescue(RECURRENT_RESCUE_FILE)
cat(sprintf(" Rescue-listed recurrent sites: %d\n", nrow(rescue_dt)))

uv_annot <- merge(uv_annot, rescue_dt, by = "hgvs_id", all.x = TRUE)
uv_annot[, flag_recurrent_site_rescued := !is.na(rescue_reason) & nzchar(rescue_reason)]

# 8. Variant-level QC flags
cat("\n>>> STEP 3: Variant-level QC flags\n")

uv_annot[, flag_common_germline := !is.na(gnomad_af_max) & gnomad_af_max >= GNOMAD_COMMON_CUTOFF]
uv_annot[, flag_possible_germline := !is.na(gnomad_af_max) &
           gnomad_af_max >= GNOMAD_POSSIBLE_CUTOFF &
           gnomad_af_max < GNOMAD_COMMON_CUTOFF]

uv_annot[, flag_recurrent_site_review := !flag_common_germline &
           is.na(gnomad_af_max) &
           !is.na(n_samples) &
           n_samples >= RECURRENCE_REVIEW_CUTOFF]

uv_annot[, flag_consequence_missing := is.na(picked_consequence)]
uv_annot[, flag_txdb_coding_only := picked_consequence == "coding_sequence_variant_txdb"]
uv_annot[, flag_txdb_only := annot_source == "txdb"]
uv_annot[, flag_functional_primary := vapply(picked_consequence, is_functional_primary, logical(1))]

if (KEEP_TXDB_CODING_ONLY) {
  uv_annot[flag_txdb_coding_only == TRUE, flag_functional_primary := TRUE]
}

uv_annot[, flag_primary_review_only := flag_recurrent_site_review & !flag_recurrent_site_rescued]

uv_annot[, flag_sensitivity_drop := flag_common_germline |
           is.na(genename) |
           flag_multigene_annotation |
           !flag_functional_primary]

uv_annot[, flag_primary_drop_variant := flag_common_germline |
           flag_possible_germline |
           is.na(genename) |
           flag_multigene_annotation |
           !flag_functional_primary]

if (PRIMARY_DROP_UNRESCUED_RECURRENT_SITE) {
  uv_annot[flag_recurrent_site_review == TRUE & !flag_recurrent_site_rescued,
           flag_primary_drop_variant := TRUE]
}

if (!KEEP_TXDB_CODING_ONLY) {
  uv_annot[flag_txdb_coding_only == TRUE, flag_primary_drop_variant := TRUE]
  uv_annot[flag_txdb_coding_only == TRUE, flag_sensitivity_drop := TRUE]
}

uv_annot[, exclusion_reason_primary := NA_character_]
uv_annot[flag_common_germline == TRUE, exclusion_reason_primary := "common_germline_af"]
uv_annot[is.na(exclusion_reason_primary) & flag_possible_germline == TRUE,
         exclusion_reason_primary := "possible_germline_af"]
if (PRIMARY_DROP_UNRESCUED_RECURRENT_SITE) {
  uv_annot[is.na(exclusion_reason_primary) &
             flag_recurrent_site_review == TRUE & !flag_recurrent_site_rescued,
           exclusion_reason_primary := "recurrent_site_review_unrescued"]
}
uv_annot[is.na(exclusion_reason_primary) & is.na(genename),
         exclusion_reason_primary := "no_gene_annotation"]
uv_annot[is.na(exclusion_reason_primary) & flag_multigene_annotation == TRUE,
         exclusion_reason_primary := "multigene_annotation"]
uv_annot[is.na(exclusion_reason_primary) & flag_txdb_coding_only == TRUE & !KEEP_TXDB_CODING_ONLY,
         exclusion_reason_primary := "txdb_coding_region_only"]
uv_annot[is.na(exclusion_reason_primary) & flag_consequence_missing == TRUE,
         exclusion_reason_primary := "missing_consequence"]
uv_annot[is.na(exclusion_reason_primary) & flag_functional_primary == FALSE,
         exclusion_reason_primary := "non_primary_functional_consequence"]

uv_annot[, exclusion_reason_sensitivity := NA_character_]
uv_annot[flag_common_germline == TRUE, exclusion_reason_sensitivity := "common_germline_af"]
uv_annot[is.na(exclusion_reason_sensitivity) & is.na(genename),
         exclusion_reason_sensitivity := "no_gene_annotation"]
uv_annot[is.na(exclusion_reason_sensitivity) & flag_multigene_annotation == TRUE,
         exclusion_reason_sensitivity := "multigene_annotation"]
uv_annot[is.na(exclusion_reason_sensitivity) & flag_txdb_coding_only == TRUE & !KEEP_TXDB_CODING_ONLY,
         exclusion_reason_sensitivity := "txdb_coding_region_only"]
uv_annot[is.na(exclusion_reason_sensitivity) & flag_consequence_missing == TRUE,
         exclusion_reason_sensitivity := "missing_consequence"]
uv_annot[is.na(exclusion_reason_sensitivity) & flag_functional_primary == FALSE,
         exclusion_reason_sensitivity := "non_primary_functional_consequence"]

cat(sprintf(" Common germline: %5d\n", sum(uv_annot$flag_common_germline, na.rm = TRUE)))
cat(sprintf(" Possible germline: %5d\n", sum(uv_annot$flag_possible_germline, na.rm = TRUE)))
cat(sprintf(" Recurrent site review: %5d\n", sum(uv_annot$flag_recurrent_site_review, na.rm = TRUE)))
cat(sprintf(" Recurrent site rescued: %5d\n", sum(uv_annot$flag_recurrent_site_rescued, na.rm = TRUE)))
cat(sprintf(" Multi-gene annotation: %5d\n", sum(uv_annot$flag_multigene_annotation, na.rm = TRUE)))
cat(sprintf(" Functional primary: %5d\n", sum(uv_annot$flag_functional_primary, na.rm = TRUE)))
cat(sprintf(" Primary review-only variants: %5d\n", sum(uv_annot$flag_primary_review_only, na.rm = TRUE)))

uv_keep_sens <- uv_annot[flag_sensitivity_drop == FALSE]
uv_keep_primary0 <- uv_annot[flag_primary_drop_variant == FALSE]

cat(sprintf(" Sensitivity-kept variants: %d\n", nrow(uv_keep_sens)))
cat(sprintf(" Primary pre-geneQC kept: %d\n", nrow(uv_keep_primary0)))

if (nrow(uv_keep_primary0) == 0L) {
  stop("Primary variant set = 0. criterionreconsider.")
}

recurrent_review_tbl <- uv_annot[flag_recurrent_site_review == TRUE, .(
  hgvs_id, chrom, pos, ref, alt,
  genename, all_genename,
  picked_consequence, annot_source,
  gnomad_af_max, n_samples, mean_tumor_af, median_tumor_af,
  flag_recurrent_site_rescued, rescue_reason,
  flag_primary_review_only
)]
setorder(recurrent_review_tbl, -flag_primary_review_only, -n_samples, genename, chrom, pos)

# 9. Variant-to-sample merge
cat("\n>>> STEP 4: Sample-variant merge\n")

var_sens <- merge(
  raw_sv[, .(sample_id_masked, chrom, pos, ref, alt, AF)],
  uv_keep_sens[, .(
    chrom, pos, ref, alt, hgvs_id, genename,
    all_genename, all_consequence, picked_consequence,
    gnomad_af_max, n_samples, mean_tumor_af, median_tumor_af,
    annot_source, flag_possible_germline, flag_recurrent_site_review,
    flag_recurrent_site_rescued, rescue_reason,
    flag_multigene_annotation, flag_txdb_only
  )],
  by = c("chrom", "pos", "ref", "alt"),
  all = FALSE
)

var_primary0 <- merge(
  raw_sv[, .(sample_id_masked, chrom, pos, ref, alt, AF)],
  uv_keep_primary0[, .(
    chrom, pos, ref, alt, hgvs_id, genename,
    all_genename, all_consequence, picked_consequence,
    gnomad_af_max, n_samples, mean_tumor_af, median_tumor_af,
    annot_source, flag_possible_germline, flag_recurrent_site_review,
    flag_recurrent_site_rescued, rescue_reason,
    flag_multigene_annotation, flag_txdb_only, flag_primary_review_only
  )],
  by = c("chrom", "pos", "ref", "alt"),
  all = FALSE
)

cat(sprintf(" Sensitivity sample-variant pairs: %d\n", nrow(var_sens)))
cat(sprintf(" Primary pre-geneQC pairs: %d\n", nrow(var_primary0)))

# 10. Gene-level txdb dependency QC
cat("\n>>> STEP 5: Gene-level support QC\n")

if (nrow(var_primary0) > 0L) {
  gene_support <- var_primary0[, .(
    total_mutated_samples = uniqueN(sample_id_masked),
    txdb_only_samples = uniqueN(sample_id_masked[annot_source == "txdb"]),
    myvariant_supported_samples = uniqueN(sample_id_masked[annot_source %in% c("myvariant", "myvariant+txdb")]),
    recurrent_review_samples = uniqueN(sample_id_masked[flag_recurrent_site_review == TRUE]),
    n_variants = uniqueN(hgvs_id)
  ), by = genename]
  
  gene_support[, txdb_only_fraction := fifelse(
    total_mutated_samples > 0L,
    txdb_only_samples / total_mutated_samples,
    NA_real_
  )]
  
  gene_support[, recurrent_review_fraction := fifelse(
    total_mutated_samples > 0L,
    recurrent_review_samples / total_mutated_samples,
    NA_real_
  )]
  
  gene_support[, flag_gene_txdb_dependency := (
    total_mutated_samples >= TXDB_ONLY_GENE_MIN_SAMPLES &
      txdb_only_fraction > TXDB_ONLY_GENE_FRAC_CUTOFF &
      !(genename %in% WHITELIST_TXDB_GENES)
  )]
  
  setorder(gene_support, -flag_gene_txdb_dependency, -txdb_only_fraction, -total_mutated_samples, genename)
} else {
  gene_support <- data.table(
    genename = character(),
    total_mutated_samples = integer(),
    txdb_only_samples = integer(),
    myvariant_supported_samples = integer(),
    recurrent_review_samples = integer(),
    n_variants = integer(),
    txdb_only_fraction = numeric(),
    recurrent_review_fraction = numeric(),
    flag_gene_txdb_dependency = logical()
  )
}

cat(sprintf(" TxDb-dependent genes flagged: %d\n",
            sum(gene_support$flag_gene_txdb_dependency, na.rm = TRUE)))

var_primary <- copy(var_primary0)
if (nrow(gene_support) > 0L) {
  bad_genes <- gene_support[flag_gene_txdb_dependency == TRUE, genename]
  var_primary <- var_primary[!genename %in% bad_genes]
} else {
  bad_genes <- character(0)
}

cat(sprintf(" Primary post-geneQC pairs: %d\n", nrow(var_primary)))

if (nrow(gene_support) > 0L) {
  uv_annot <- merge(
    uv_annot,
    gene_support[, .(
      genename,
      flag_gene_txdb_dependency,
      txdb_only_fraction,
      recurrent_review_fraction
    )],
    by = "genename",
    all.x = TRUE
  )
} else {
  uv_annot[, `:=`(
    flag_gene_txdb_dependency = FALSE,
    txdb_only_fraction = NA_real_,
    recurrent_review_fraction = NA_real_
  )]
}
uv_annot[is.na(flag_gene_txdb_dependency), flag_gene_txdb_dependency := FALSE]

reviewworthy_primary_tbl <- uv_annot[
  flag_primary_drop_variant == FALSE &
    (flag_primary_review_only == TRUE | flag_gene_txdb_dependency == TRUE),
  .(
    hgvs_id, chrom, pos, ref, alt,
    genename, all_genename,
    picked_consequence, annot_source,
    n_samples, mean_tumor_af, median_tumor_af,
    flag_recurrent_site_review, flag_recurrent_site_rescued, rescue_reason,
    flag_primary_review_only,
    flag_gene_txdb_dependency,
    txdb_only_fraction,
    recurrent_review_fraction
  )
]
setorder(reviewworthy_primary_tbl, -flag_primary_review_only, -flag_gene_txdb_dependency, -n_samples, genename)

# 11. Matrix / FMG
cat("\n>>> STEP 6: Mutation matrix & FMG\n")

res_sens <- build_matrix_and_prevalence(
  var_dt = var_sens,
  all_samples = ALL_SAMPLES,
  prevalence_cutoff = PREV_CUTOFF
)

res_primary <- build_matrix_and_prevalence(
  var_dt = var_primary,
  all_samples = ALL_SAMPLES,
  prevalence_cutoff = PREV_CUTOFF
)

fmg_sens <- merge(
  res_sens$fmg,
  gene_support[, .(
    genename, total_mutated_samples, txdb_only_samples,
    myvariant_supported_samples, recurrent_review_samples,
    txdb_only_fraction, recurrent_review_fraction,
    flag_gene_txdb_dependency
  )],
  by.x = "gene", by.y = "genename", all.x = TRUE
)

fmg_primary <- merge(
  res_primary$fmg,
  gene_support[, .(
    genename, total_mutated_samples, txdb_only_samples,
    myvariant_supported_samples, recurrent_review_samples,
    txdb_only_fraction, recurrent_review_fraction,
    flag_gene_txdb_dependency
  )],
  by.x = "gene", by.y = "genename", all.x = TRUE
)

setorder(fmg_sens, -prevalence, -n_mutated, gene)
setorder(fmg_primary, -prevalence, -n_mutated, gene)

cat(sprintf(" Sensitivity matrix: %d samples x %d genes\n",
            nrow(res_sens$mat), ncol(res_sens$mat)))
cat(sprintf(" Primary matrix: %d samples x %d genes\n",
            nrow(res_primary$mat), ncol(res_primary$mat)))

cat(sprintf(" Sensitivity FMG (>= %.0f%%): %d genes\n",
            PREV_CUTOFF * 100, nrow(fmg_sens)))
cat(sprintf(" Primary FMG (>= %.0f%%): %d genes\n",
            PREV_CUTOFF * 100, nrow(fmg_primary)))

# 12.
cat("\n>>> STEP 7: result n")

f <- function(name) file.path(OUTDIR, sprintf("%s_%s.csv", name, OUT_PREFIX))

write.csv(fmg_primary, f("FMG_list_primary"), row.names = FALSE)
write.csv(fmg_primary, file.path(OUTDIR, "Tab1_2_FMG_list_primary.csv"), row.names = FALSE)   # the name stages 03-07 read
cat(sprintf(" [1] %s\n", f("FMG_list_primary")))

write.csv(fmg_sens, f("FMG_list_sensitivity"), row.names = FALSE)
cat(sprintf(" [2] %s\n", f("FMG_list_sensitivity")))

write_matrix_csv(res_primary$mat, f("mutation_matrix_primary"))
cat(sprintf(" [3] %s\n", f("mutation_matrix_primary")))

write_matrix_csv(res_sens$mat, f("mutation_matrix_sensitivity"))
cat(sprintf(" [4] %s\n", f("mutation_matrix_sensitivity")))

var_primary_out <- copy(var_primary)
list_cols1 <- names(var_primary_out)[vapply(var_primary_out, is.list, logical(1))]
if (length(list_cols1) > 0L) var_primary_out[, (list_cols1) := NULL]
write.csv(var_primary_out, f("variants_somatic_primary"), row.names = FALSE)
cat(sprintf(" [5] %s\n", f("variants_somatic_primary")))

var_sens_out <- copy(var_sens)
list_cols2 <- names(var_sens_out)[vapply(var_sens_out, is.list, logical(1))]
if (length(list_cols2) > 0L) var_sens_out[, (list_cols2) := NULL]
write.csv(var_sens_out, f("variants_somatic_sensitivity"), row.names = FALSE)
cat(sprintf(" [6] %s\n", f("variants_somatic_sensitivity")))

write.csv(res_primary$prevalence, f("gene_prevalence_primary"), row.names = FALSE)
cat(sprintf(" [7] %s\n", f("gene_prevalence_primary")))

write.csv(res_sens$prevalence, f("gene_prevalence_sensitivity"), row.names = FALSE)
cat(sprintf(" [8] %s\n", f("gene_prevalence_sensitivity")))

write.csv(gene_support, f("gene_support_qc"), row.names = FALSE)
cat(sprintf(" [9] %s\n", f("gene_support_qc")))

qc_out <- copy(uv_annot)[order(flag_primary_drop_variant, -n_samples, genename, chrom, pos)]
list_cols3 <- names(qc_out)[vapply(qc_out, is.list, logical(1))]
if (length(list_cols3) > 0L) qc_out[, (list_cols3) := NULL]
write.csv(qc_out, f("variant_qc_flags"), row.names = FALSE)
cat(sprintf(" [10] %s\n", f("variant_qc_flags")))

write.csv(recurrent_review_tbl, f("recurrent_site_review"), row.names = FALSE)
cat(sprintf(" [11] %s\n", f("recurrent_site_review")))

write.csv(reviewworthy_primary_tbl, f("reviewworthy_variants_primary"), row.names = FALSE)
cat(sprintf(" [12] %s\n", f("reviewworthy_variants_primary")))

# 13. Summary
cat("\n")
cat("========================================================================\n")
cat(" Stage 1+2 Pipeline v10.1 Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" Input sample-variant rows : %8d\n", nrow(raw_sv)))
cat(sprintf(" Total samples : %8d\n", N_TOTAL))
cat(sprintf(" Unique variants all : %8d\n", nrow(uv_all)))
cat(sprintf(" HGVS build fail : %8d\n", n_hgvs_fail))
cat(sprintf(" HGVS OK variants : %8d\n", nrow(uv)))
cat(sprintf(" Gene annotated : %8d\n", sum(!is.na(uv_annot$genename))))
cat(sprintf(" Common germline : %8d\n", sum(uv_annot$flag_common_germline, na.rm = TRUE)))
cat(sprintf(" Possible germline : %8d\n", sum(uv_annot$flag_possible_germline, na.rm = TRUE)))
cat(sprintf(" Recurrent review : %8d\n", sum(uv_annot$flag_recurrent_site_review, na.rm = TRUE)))
cat(sprintf(" Recurrent rescued : %8d\n", sum(uv_annot$flag_recurrent_site_rescued, na.rm = TRUE)))
cat(sprintf(" Multigene annotation : %8d\n", sum(uv_annot$flag_multigene_annotation, na.rm = TRUE)))
cat(sprintf(" Sensitivity kept vars : %8d\n", nrow(uv_keep_sens)))
cat(sprintf(" Primary kept vars preQC : %8d\n", nrow(uv_keep_primary0)))
cat(sprintf(" TxDb-dependent genes : %8d\n", sum(gene_support$flag_gene_txdb_dependency, na.rm = TRUE)))
cat(sprintf(" Primary kept pairs postQC: %8d\n", nrow(var_primary)))
cat(sprintf(" Sensitivity matrix : %4d samples x %4d genes\n",
            nrow(res_sens$mat), ncol(res_sens$mat)))
cat(sprintf(" Primary matrix : %4d samples x %4d genes\n",
            nrow(res_primary$mat), ncol(res_primary$mat)))
cat(sprintf(" FMG threshold : %7.2f%%\n", PREV_CUTOFF * 100))
cat(sprintf(" Sensitivity FMG genes : %8d\n", nrow(fmg_sens)))
cat(sprintf(" Primary FMG genes : %8d\n", nrow(fmg_primary)))
cat("========================================================================\n")

cat("\n>>> Top 20 Primary FMG\n")
print(head(fmg_primary, 20))

cat("\n>>> Top 20 Sensitivity FMG\n")
print(head(fmg_sens, 20))

cat("\n>>> Primary exclusion reason distribution\n")
print(qc_out[, .N, by = exclusion_reason_primary][order(-N)])

cat("\n>>> Sensitivity exclusion reason distribution\n")
print(qc_out[, .N, by = exclusion_reason_sensitivity][order(-N)])

cat("\n>>> Recurrent review table top 20\n")
print(head(recurrent_review_tbl, 20))

cat("\n>>> TxDb-dependent genes\n")
print(gene_support[flag_gene_txdb_dependency == TRUE][order(-txdb_only_fraction, -total_mutated_samples)])

# 14. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage01_02.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 1+2 Pipeline v10.1 done (%s) ===\n", Sys.time()))

# end: log file close
sink(type = "message")
sink()
close(log_con)
