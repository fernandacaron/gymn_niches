
cyc_matk <- read.FASTA("data/cyc_matk.fasta")
cyc_rbcl <- read.FASTA("data/cyc_rbcl.fasta")

# pegar nomes de cada spp
names_cyc_matk <- names(cyc_matk)
bin_cyc_matk <- character()
for (j in 1:length(names_cyc_matk)) {
  bin_cyc_matk[j] <- paste(strsplit(names_cyc_matk[j], split = " ")[[1]][c(2, 3)], collapse = "_")
}

names_cyc_rbcl <- names(cyc_rbcl)
bin_cyc_rbcl <- character()
for (j in 1:length(names_cyc_rbcl)) {
  bin_cyc_rbcl[j] <- paste(strsplit(names_cyc_rbcl[j], split = " ")[[1]][c(2, 3)], collapse = "_")
}

#contar número de bases em cada sequências
counts  <-  numeric()
for (i in 1:length(names_cyc_matk)) counts[i]  <-  sum(base.freq(cyc_matk[i], freq = TRUE))

counts  <-  numeric()
for (i in 1:length(names_cyc_rbcl)) counts[i]  <-  sum(base.freq(cyc_rbcl[i], freq = TRUE))

#escolher a maior sequencia de cada spp
dat_cyc_matk <- cbind(1:length(bin_cyc_matk), bin_cyc_matk, counts)
spp_cyc_matk <- names(table(bin_cyc_matk))
del_cyc_matk <- numeric()
for (i in 1:length(spp_cyc_matk)) {
  dat_cyc_matk_ <- subset(dat_cyc_matk, bin_cyc_matk == spp_cyc_matk[i])
  if(length(dat_cyc_matk_) > 3) {
    dat_cyc_matk_ <- dat_cyc_matk_[order(as.numeric(dat_cyc_matk_[, 3]), decreasing = T), ]
    dat_cyc_matk_ <- dat_cyc_matk_[-1,]
    ifelse(length(dat_cyc_matk_) == 3, 
           del_cyc_matk <- c(del_cyc_matk, as.numeric(dat_cyc_matk_[1])), 
           del_cyc_matk <- c(del_cyc_matk, as.numeric(dat_cyc_matk_[, 1])))
  }
}
cyc_matk2 <- cyc_matk[- del_cyc_matk]
write.FASTA(cyc_matk2, "data/cyc_matk2.fasta")

dat_cyc_rbcl <- cbind(1:length(bin_cyc_rbcl), bin_cyc_rbcl, counts)
spp_cyc_rbcl <- names(table(bin_cyc_rbcl))
del_cyc_rbcl <- numeric()
for (i in 1:length(spp_cyc_rbcl)) {
  dat_cyc_rbcl_ <- subset(dat_cyc_rbcl, bin_cyc_rbcl == spp_cyc_rbcl[i])
  if(length(dat_cyc_rbcl_) > 3) {
    dat_cyc_rbcl_ <- dat_cyc_rbcl_[order(as.numeric(dat_cyc_rbcl_[, 3]), decreasing = T), ]
    dat_cyc_rbcl_ <- dat_cyc_rbcl_[-1,]
    ifelse(length(dat_cyc_rbcl_) == 3, 
           del_cyc_rbcl <- c(del_cyc_rbcl, as.numeric(dat_cyc_rbcl_[1])), 
           del_cyc_rbcl <- c(del_cyc_rbcl, as.numeric(dat_cyc_rbcl_[, 1])))
  }
}
cyc_rbcl2 <- cyc_rbcl[- del_cyc_rbcl]
write.FASTA(cyc_rbcl2, "data/cyc_rbcl2.fasta")

# edit by hand
cyc_matk <- read.FASTA("data/cyc_matk3.fasta")
cyc_rbcl <- read.FASTA("data/cyc_rbcl3.fasta")

names_cyc_matk <- names(cyc_matk)
bin_cyc_matk <- character()
for (j in 1:length(names_cyc_matk)) {
  bin_cyc_matk[j] <- paste(strsplit(names_cyc_matk[j], split = " ")[[1]][c(2, 3)], collapse = "_")
}
names(cyc_matk) <- bin_cyc_matk

names_cyc_rbcl <- names(cyc_rbcl)
bin_cyc_rbcl <- character()
for (j in 1:length(names_cyc_rbcl)) {
  bin_cyc_rbcl[j] <- paste(strsplit(names_cyc_rbcl[j], split = " ")[[1]][c(2, 3)], collapse = "_")
}
names(cyc_rbcl) <- bin_cyc_rbcl

names_cyc <- names(cyc_rbcl)[names(cyc_rbcl) %in% names(cyc_matk)]

cyc_matk <- cyc_matk[names(cyc_matk) %in% names_cyc]
cyc_rbcl <- cyc_rbcl[names(cyc_rbcl) %in% names_cyc]

write.FASTA(cyc_matk, "data/cyc_matk_final.fasta")
write.FASTA(cyc_rbcl, "data/cyc_rbcl_final.fasta")

# align

cyc_files <- list.files(path = "data", pattern = "align.fasta$", 
                        recursive = TRUE, full.names = TRUE)

cyc_conc <- read.multiFASTA(cyc_files)
cyc_final <- concatenate(cyc_conc)

write.FASTA(cyc_final, "data/cyc_seq.fasta")

cyc_seq <- read.FASTA("data/cyc_seq.fasta")