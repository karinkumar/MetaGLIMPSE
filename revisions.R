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
#read in INFO scores from Simone


setkey(true_r2, `#SNP.ID`)
setkey(est_r2, V1)
r2_comp <- true_r2[est_r2, nomatch =0]
setnames(r2_comp, "V2", "estR2")
r2_comp[, MAF:= ifelse(Allele.Frequency > 0.5, 1 - Allele.Frequency, Allele.Frequency)]
#I might need to recompute the R2s afterwards because the Minor Allele flipped...maybe not because
#both imputed and estimated R2 will have the allele flip, so the difference should be the same
cor.test(r2_comp$Imputation.R2, r2_comp$estR2)
cor(r2_comp$MAF, r2_comp$Imputation.R2)
plot(r2_comp$MAF, r2_comp$estR2)
plot(r2_comp$MAF, r2_comp$Imputation.R2)
r2_comp[estR2 > 1, estR2:= 1] # set R2 greater than 1 to be 1

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



####################Reviewer 1 Question1############
setwd("~/1kg30xASW/scripts")
files <- list.files(pattern = "6*indv.switch")

load_file <- function(cov, anc) {
  file_name <- sprintf("%d_asw_switch_%s.diff.indv.switch", anc, cov)
  if (file.exists(file_name)) {
    df <- fread(file_name)
    df[, cov := cov]
    df[, panel := anc]
    return(df)
  } else {
    print(file_name)
    return(NULL)
  }
}
prefix = ""
fs <- list()
for (p in prefix) {
  for (cov in 6) {
    # AMR
    dfs[[paste0(cov, p)]] <- load_file(p, chr, "1KG_AMR")
    # AFR
    dfs[[paste0(cov, p)]] <- load_file(p, chr, "1KG_AFR")
    # EUR
    dfs[[paste0("EUR_chr", chr, p)]] <- load_file(p, chr, "1KG_EUR")
  }
}
# Remove NULL entries (files that weren't loaded)
dfs <- dfs[!sapply(dfs, is.null)]

# Combine all data
all <- rbindlist(dfs, fill = TRUE)