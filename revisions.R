library(data.table)
#dt <- fread("~/aggRSquare/release-build/ASW/240304isecASW_withsingletondiploid2x_EUR.RSquare") # I see this is with singleton! 
#maybe that's where I've gotten it wrong!!! And that's where the extra variants are coming from?
setwd("~/aggRSquare/release-build")

##### Reviewer 1 Question 3 ############

EUR_R2 <- fread("~/aggRSquare/release-build/ASW/240304isecASW_no1kgsingletondiploid2x_EUR.RSquare")
EUR_AF <- fread("unrelatedphase3EURAF.txt")
AFR_AF <- fread("unrelatedphase3AFRAF.txt")

EUR_R2
setnames(EUR_AF, "V2", "EUR_AF")
setnames(AFR_AF, "V2", "AFR_AF")
setkey(EUR_R2, `#SNP.ID`)
setkey(EUR_AF, V1)
setkey(AFR_AF, V1)
tmp <- EUR_R2[EUR_AF, nomatch = 0]
final <- AFR_AF[tmp]
final$AF_diff <- abs(final$AFR_AF - final$EUR_AF)
cor(final$Imputation.R2, final$AF_diff)
cor.test(final$Imputation.R2, final$AF_diff)
cor.test(final$EUR_AF, final$Imputation.R2)

######okay we now have the correct # of polymorphic variants###
t <- fread("ASW/240304isecASW_no1kgsingletondiploid8x_plain_zerodosage.aggRSquare")
sum(t$`#Variants`)

### Reviewer 2 Question 2 ################
true_r2 <- fread("ASW/240304isecASW_no1kgsingletondiploid0.1x_mixed_zerodosage.RSquare")
est_r2 <- fread("/net/fantasia/home/kiranhk/HMM/est_r2_ASW_0.1x.txt")

true_r2 <- fread("251118LOO_2x_UKBB.RSquare")
est_r2 <- fread("/net/fantasia/home/kiranhk/HMM/est_r2_2x_EAS_UKBB.txt")
setnames(est_r2, "V2", "estR2")

setkey(true_r2, `#SNP.ID`)
setkey(est_r2, V1)
r2_comp <- true_r2[est_r2, nomatch =0]
#setnames(r2_comp, "V2", "estR2")
r2_comp[, MAF:= ifelse(Allele.Frequency > 0.5, 1 - Allele.Frequency, Allele.Frequency)]
#I might need to recompute the R2s afterwards because the Minor Allele flipped...maybe not because
#both imputed and estimated R2 will have the allele flip, so the difference should be the same

#plot(r2_comp$MAF, r2_comp$estR2)
#plot(r2_comp$MAF, r2_comp$Imputation.R2)
r2_comp[estR2 > 1, estR2:= 1] # set R2 greater than 1 to be 1

cor.test(r2_comp$Imputation.R2, r2_comp$estR2, method = "spearman")
cor.test(r2_comp$Imputation.R2, r2_comp$estR2)
cor(r2_comp$MAF, r2_comp$Imputation.R2) #MAF correlation

r2_comp[, diffr2 := Imputation.R2 - estR2]
common <- r2_comp[MAF > 0.1]
#plot(common$MAF, common$diffr2)
#loess_fit <- loess(common$diffr2 ~ common$MAF, span = 0.7)
# Add the smoothed line to the plot
#lines(common$MAF, predict(loess_fit), col = "red", lwd = 2)
cor(common$Imputation.R2, common$estR2)
ggplot(common, aes(x = MAF, y = diffr2)) +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Density Plot of Difference in R2 (True - Estimated) VS MAF for Common Variants") + 
  xlab("Minor Allele Frequency(MAF)") + 
  ylab("True R2 - Estimated R2")

ggplot(r2_comp, aes(x = log10(MAF), y = diffr2)) +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Density Plot of Difference in R2 (True - Estimated) VS MAF ") + 
  xlab("Minor Allele Frequency(MAF)") + 
  ylab("True R2 - Estimated R2")

mean(r2_comp$diffr2)
sd(r2_comp$diffr2)


######## Reviewer 1 Question 5######

aggRSquare2Dt <- function(file_path, panel, target_ancestry, coverage) {
  # Compose filename based on your convention
  full_file <- paste0("~/aggRSquare/release-build/", file_path, coverage, "x_", panel, ".aggRSquare")
  
  # Read tab-delimited file, no header
  dt <- fread(full_file, sep = "\t", header = FALSE)
  dt <- dt[!grepl("^#", dt[[1]])]
  
  # Adjust columns as needed
  setnames(dt, c("binned_MAF", "Avg_MAF", "N", "aggregated_R2", "gold_MAF", "imputed_MAF"))
  
  # Add metadata columns
  dt[, Panel := panel]
  dt[, TargetAncestry := target_ancestry]
  dt[, Coverage := coverage]
  return(dt)
}

cov = c("0.1", "0.5", "1", "2", "4", "6", "8")

# Define all your configs with multiple coverages per ancestry/panel
panel_configs <- list(
  list(
    file_path = "250901isec10aDNA_no1kgsingletondiploid", 
    target_ancestry = "aDNA Targets", 
    panels = c("AFR_EUR", "AFR", "EUR", "mixed_fb_zerodosage", "mixed_phased_zerodosage", "plain_fb_zerodosage"), 
    coverages = cov
  ),
  list(
    file_path = "ASW/240304isecASW_no1kgsingletondiploid", 
    target_ancestry = "African-American Targets", 
    panels = c("mixed_zerodosage", "plain_zerodosage", "AFR", "EUR", "AFR_EUR" "mixed_zerodosagephased"), 
    coverages = cov
  ),
  list(
    file_path = "240308isec2panelEASdiploid_no1kgsingleton_", 
    target_ancestry = "East Asian Targets", 
    panels = c("AFR_EUR", "mixed", "AFR", "EUR"), 
    coverages = cov
  ),
  list(
    file_path = "EUR/240703LOO_", 
    target_ancestry = "European Targets", 
    panels = c("EURmega", "EURA", "mixedzerodosage", "EURB", "plainzerodosage"), 
    coverages = cov
  ),
  list(
    file_path = "251118LOO_", 
    target_ancestry = "East Asian Targets", 
    panels = c("EAS", "UKBB", "mixed"), 
    coverages = "2"
  )
)

# Expand all configurations into a data.table of jobs to run
all_configs <- rbindlist(
  lapply(panel_configs, function(cfg) {
    # For every combination of panel and coverage
    expand.grid(
      file_path = cfg$file_path,
      panel = cfg$panels,
      target_ancestry = cfg$target_ancestry,
      coverage = cfg$coverages,
      stringsAsFactors = FALSE
    )
  }), fill = TRUE
)

# Read all files, storing results
results_list <- lapply(1:nrow(all_configs), function(i) {
  cfg <- all_configs[i]
  # Try to read the file, if fails return NULL
  tryCatch({
    aggRSquare2Dt(cfg$file_path, cfg$panel, cfg$target_ancestry, cfg$coverage)
  }, error = function(e) {
    message(sprintf("Could not read file for %s %s %s %s", 
                    cfg$file_path, cfg$panel, cfg$target_ancestry, cfg$coverage))
    NULL
  })
})

# Combine all successful reads
full_data_table <- rbindlist(
  results_list[!sapply(results_list, is.null)], 
  use.names = TRUE, fill = TRUE
)

full_data_table[, gold_MAF:= NULL]
full_data_table[, imputed_MAF := NULL]
full_data_table[Panel %in% c("mixed", "mixed_zerodosage", "mixed_fb_zerodosage", "mixedzerodosage"), Panel := "MetaGLIMPSE"]
full_data_table[Panel %in% c("plain", "plain_zerodosage", "plain_fb_zerodosage", "plainzerodosage"), Panel := "MetaGLIMPSE-plain"]
full_data_table[Panel %in% c("mixed_phased_zerodosage", "mixed_zerodosagephased"), Panel := "MetaGLIMPSE-phased"]
# Save to CSV

full_data_table
fwrite(full_data_table, "~/aggRSquare/release-build/ImputationResults_Aggregated.csv")

cat("Aggregated imputation results table written to ~/aggRSquare/ImputationResults_Aggregated.csv\n")
