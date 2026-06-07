#### HPRC R2 - MEI call comparison (Shadi 2026 - Anna-Sophie team)
# libraries ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(VariantAnnotation)

sys_dir <- "/home/shadi/Desktop/PanTE_Human_R2/MEI_call_comp"

# VCF annotation v2 - loop -------------------------------------------------------

## define your mapping df of VCF ID and RM headers

id_mapping <- read.table(file.path(sys_dir, "id_mapping.bed"), header=FALSE, sep="\t", col.names=c("full_id", "short_id"))

id_mapping <- id_mapping %>%
  mutate(
    base_id = gsub("_(REF|ALT[0-9]*)$", "", short_id),
    vcf_id  = gsub("^chr[^_]+_[0-9]+_(.+)_(REF|ALT[0-9]*)$", "\\1", full_id))

nrow(id_mapping)

## define the header INFO for the new RM in vcf

RM_fields <- DataFrame(
  Number = rep(".", 10),
  Type = rep("String", 10),
  Description = c(
    "ALT RM percent divergence",
    "ALT RM begin positions in variant",
    "ALT RM end positions in variant",
    "ALT RM matching repeat names",
    "ALT RM repeat class/family",
    "REF RM percent divergence",
    "REF RM begin positions",
    "REF RM end positions",
    "REF RM matching repeat names",
    "REF RM repeat class/family"),
  row.names = c("ALT_PERC_DIV","ALT_BEGIN","ALT_END","ALT_MATCHING_REPEAT","ALT_REPEAT_CLASSFAMILY",
                "REF_PERC_DIV","REF_BEGIN","REF_END","REF_MATCHING_REPEAT","REF_REPEAT_CLASSFAMILY")) %>% as.data.frame()

chr_vec <- c(paste0("chr", 1:22), "chrX", "chrY")

for (chr in chr_vec) {
  
  ## define paths and load files
  
  vcf_path <- file.path(sys_dir, paste0("vcf_chm13_chr/hprc-v2.0-mc-chm13_norm_indx_mod_len_", chr, ".wave.vcf.gz"))
  bed_path <- file.path(sys_dir, paste0("repeatmasker_output_vcf2fasta_chm13_chr/", chr, "/TE_only_", chr, ".bed"))
  id_mapping_path <- file.path(sys_dir, "id_mapping.bed")
  
  vcf        <- readVcf(vcf_path)
  bed        <- read.table(bed_path, header=TRUE, sep="\t")

  ## tag each row as REF or ALT and extract base ID (chr_pos_shortID with no ALT or REF)
  
  bed_mod <- bed %>%
    mutate(
      allele   = case_when(
        grepl("_REF$", query) ~ "REF",
        grepl("_ALT",  query) ~ "ALT",
        TRUE ~ NA_character_),
      base_id = gsub("_(REF|ALT[0-9]*)$", "", query))
  
  ## collapse RM hits per seq query
  
  bed_collapsed <- bed_mod %>%
    filter(!is.na(allele), query != "", base_id != "") %>%
    arrange(base_id, allele, begin) %>%
    group_by(query, base_id, allele) %>%
    summarise(
      perc_div        = paste(perc_div, collapse = ","),
      begin           = paste(begin, collapse = ","),
      end             = paste(end, collapse = ","),
      matching_repeat = paste(matching_repeat, collapse = ","),
      repeat_class    = paste(repeat_classfamily, collapse = ","),
      .groups = "drop") %>%
    arrange(base_id, allele, query) %>%
    ungroup()
  
  ## attach VCF IDs to bed via id_mapping
  
  bed_final <- bed_collapsed %>%
    left_join(
      id_mapping %>% dplyr::select(short_id, vcf_id),
      by = c("query" = "short_id"))
  
  # length(unique(bed_final$vcf_id))
  # length(bed_final[bed_final$allele == "REF",]$vcf_id)
  # length(unique(bed_collapsed$query))
  # nrow(bed_collapsed)
  # nrow(bed_final)
  
  ## split bed_final into ALT and REF projections
  
  bed_final_alt <- bed_final %>% filter(allele == "ALT")
  bed_final_ref <- bed_final %>% filter(allele == "REF")
  
  vcf_ids <- names(vcf)
  alt_idx <- match(vcf_ids, bed_final_alt$vcf_id)
  ref_idx <- match(vcf_ids, bed_final_ref$vcf_id)
  
  ## annotate VCF and store it
  
  info(header(vcf)) <- rbind(info(header(vcf)), RM_fields)
  
  info(vcf)$ALT_PERC_DIV <- bed_final_alt$perc_div[alt_idx]
  info(vcf)$ALT_BEGIN <- bed_final_alt$begin[alt_idx]
  info(vcf)$ALT_END <- bed_final_alt$end[alt_idx]
  info(vcf)$ALT_MATCHING_REPEAT <- bed_final_alt$matching_repeat[alt_idx]
  info(vcf)$ALT_REPEAT_CLASSFAMILY <- bed_final_alt$repeat_class[alt_idx]
  
  info(vcf)$REF_PERC_DIV <- bed_final_ref$perc_div[ref_idx]
  info(vcf)$REF_BEGIN <- bed_final_ref$begin[ref_idx]
  info(vcf)$REF_END <- bed_final_ref$end[ref_idx]
  info(vcf)$REF_MATCHING_REPEAT <- bed_final_ref$matching_repeat[ref_idx]
  info(vcf)$REF_REPEAT_CLASSFAMILY <- bed_final_ref$repeat_class[ref_idx]
  
  names(vcf) <- sub("_[ACGTN>].*$", "", names(vcf))
  
  writeVcf(vcf, file.path(sys_dir, "annotated_vcf_RM_Rscript", paste0("hprc_v2_mc_chm13_norm_mod_len_annotated_RM_", chr, ".vcf")))
  
}


