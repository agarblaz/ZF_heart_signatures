# Clean environment
rm(list = ls())

# Load required libraries
library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(data.table)
library(DoubletFinder)

# Activate garbage collector
gc(verbose = FALSE)

# Function to process a single dataset
process_dataset <- function(file_path, age) {
  
  print(paste("Processing dataset for age:", age))
  print(paste("File path:", file_path))
  # Load the dataset
  data <- Read10X(data.dir = file_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = paste0("zebrafish_", age)
  )
  
  # Add age metadata
  seurat_obj$age <- age
  
  ## Basic preprocessing
  # Identify mitochondrial content
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^mt-"
  )
  
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^hb[ab]"
  )
  
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 2500 &
      percent.mt < 30 &
      percent.hb < 10
  )
  
  # Estimate doublet formation rate (e.g., 0.075 for 5,000 cells)
  doublet_rate <- ncol(seurat_obj) / 1000 * 0.01  # Adjust if needed
  nExp <- round(doublet_rate * ncol(seurat_obj))
  
  # Find optimal pK (optional but recommended)
  sweep.res <- paramSweep(seurat_obj, PCs = 1:30, sct = TRUE)
  gc()
  
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  print(paste("Best pK:", best_pK, " | nExp:", nExp))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = 1:30,
    pN = 0.25,       # Default
    pK = best_pK,    # Or manually set (e.g., 0.09)
    nExp = nExp,
    sct = TRUE      # Set TRUE if using SCTransform
  )
  
  DF_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(DF_col) != 1) {
    strop("Could not uniquely identify the doublet classification column. Please check the metadata columns.")
  }
  print(paste("Using doublets column:", DF_col))
  # Añadir la columna como una nueva variable explícita para filtrar
  seurat_obj$DoubletStatus <- seurat_obj@meta.data[[DF_col]]
  seurat_obj@meta.data[[DF_col]] <- NULL
  seurat_obj_clean <- subset(
    seurat_obj,
    subset = DoubletStatus == "Singlet"
  )
  
  print(paste("Number of cells after doublet removal:", ncol(seurat_obj_clean)))
  
  # Normalization
  seurat_obj_clean <- SCTransform(
    seurat_obj_clean,
    vst.flavor = "v2",
    method = "glmGamPoi"
  )
  gc()
  
  # Run PCA
  seurat_obj_clean <- RunPCA(seurat_obj_clean)
  
  # Run UMAP
  seurat_obj_clean <- RunUMAP(seurat_obj_clean, dims = 1:30)
  
  # Find neighbors
  seurat_obj_clean <- FindNeighbors(seurat_obj_clean, dims = 1:30)
  
  # Find clusters
  seurat_obj_clean <- FindClusters(seurat_obj_clean, resolution = 0.5)
  
  gc(verbose = FALSE)
  
  return(seurat_obj_clean)
}

# List of datasets and their corresponding ages
datasets <- list(
  list(path = "data/adult_rep1", age = "adult"),
  list(path = "data/adult_rep2", age = "adult"),
  list(path = "data/adult_rep3", age = "adult"),
  list(path = "data/adult_rep4", age = "adult")
)

# Process all datasets
objs_list <- map(datasets, ~process_dataset(.x$path, .x$age))

saveRDS("output/adult_raw_data.rds")

gc(verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = objs_list, nfeatures = 5000)
objs_list <- PrepSCTIntegration(object.list = objs_list, anchor.features = features)

# Integrate data using Seurat's integration method
print("Integrating datasets...")
print("Finding integration anchors...")
anchors <- FindIntegrationAnchors(
  object.list = objs_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:50,
  k.anchor = 5)
print("Integration anchors found.")
gc(verbose = FALSE)

# Integrate the data
print("Integrating data...")
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT")
print("Data integrated.")
gc(verbose = FALSE)

# Switch back to the RNA assay
DefaultAssay(integrated) <- "integrated"
gc(verbose = FALSE)

# Run dimensionality reduction and clustering on integrated data
print("Running PCA, UMAP, and clustering on integrated data...")
print("Running PCA...")
integrated <- RunPCA(integrated, npcs = 50)
print("PCA complete.")

print("Running UMAP...")
integrated <- RunUMAP(integrated, dims = 1:30)
print("UMAP complete.")

print("Finding neighbors...")
integrated <- FindNeighbors(integrated, dims = 1:30, k.param = 20)
print("Neighbors found.")

print("Finding clusters...")
integrated <- FindClusters(integrated, resolution = 0.4)
print("Clusters found.")


DefaultAssay(integrated) <- "SCT"
# Save the seurat object as rds file
saveRDS(
  integrated,
  file = "output/adult_integrated.rds"
)

# Identify markers for each cluster
print("Finding markers for each cluster...")
markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE)
print("Markers identified.")

writexl::write_xlsx(
  markers,
  path = "output/adult_markers.xlsx"
)

gc(verbose = FALSE)

# Annotate cell types based on marker genes (this step requires manual curation)
cell_type_markers <- list(
  ## Markers from Winata's AVC paper:
  ## https://link.springer.com/article/10.1007/s00018-021-03939-y#Sec22
  "AVC" = c(
    "tnfrsf11a", "dact1", "fgfr3", "zfpm1", "dag1", "ctnnb1", "lamb2l", "fn1a",
    "hectd1", "cx30.3", "fbn2b", "fgf3", "sema3b", "cacna1c", "nfkb2", "slit2",
    "mdm4", "nfkb1", "has1", "tbx5a", "crip2", "vcana", "tln1", "myoc", "gja3",
    "cxadr", "adamts9", "agr2", "pdgfra", "hspb1", "myl7", "snai1b", "flna",
    "eya1", "acana", "tnfrsf1a", "tbx20", "nkx2.5"),
  ## Markers from Winatas' SAN papers:
  ## https://doi.org/10.1186/s12864-021-08016-z 
  ## https://www.sciencedirect.com/science/article/pii/S2589004224013087
  "CM_SINOATRIAL" = c(
    "wnt1", "bmp4", "tbx5a", "tbx18", "isl1", "lyve1a", "rsg6", "cpne5b",
    "chrm2a", "vsnl1a", "fgf13a", "tbx3a", "shox2", "hcn4", "hcn1", "smtnl1",
    "hmga1a", "jun", "nr2f1a", "atf3", "sox11b", "tbx2a", "mef2aa", "xbp1",
    "cebpd", "irf9", "tob1b", "gata6", "fosl2", "fosb", "atp2a2a", "fhl27",
    "vdac3", "slc20a1b", "atp2b1a", "sec61g", "clic2", "rrad", "atp1b3b",
    "slc43a3b", "tmem38a", "kcnq1.1", "lama5", "hapln1a", "ctnnb1", "bsg", 
    "tenm3"),
  ## markers from Junker's Fibroblast paper:
  ## https://www.nature.com/articles/s41588-022-01129-5
  "FIBROBLAST" = c(
    "col1a2", "mdka", "col1a1a", "col1a1b", "col5a1", "dpt", "fn1b", "ccl25b",
    "sparc", "clu", "dcn", "mfap5", "postnb", "c4", "fn1a", "col12a1a", "c6.1",
    "si:ch1073-459j12.1", "htra1b", "pmp22a", "hapln1a", "ckba", "col6a3",
    "gstm.3", "col6a1", "tnfaip6", "tfa", "mmp2", "f3b", "rbp4", "col5a2a",
    "CR936442.1", "ifitm1", "sostdc1a", "anxa2a", "pcolcea", "olfml3b", "tmsb2",
    "si:ch211-105c13.3", "si:ch211-137i24.10", "c1qtnf5", "ppib",
    "si:ch211-106h4.12", "col6a2", "mgp", "cxcl12a", "edil3a", "dkk3b",
    "si:ch211-198c19.3", "fstl1b"),
  ## markers for Endothelial cells
  "ENDOTHELIAL" = c(
    "apnln", "hbegfb", "admb", "apln", "si:ch211-195b11.3", "sele", "itga2b",
    "thbs1b", "fgl2a", "FO681357.1", "pim2", "gp1bb", "cldn5b", "gcnt4a",
    "il11b", "si:ch211-105j21.9", "tln1", "si:ch211-214p13.9", "f5",
    "si:ch211-214p16.2", "si:ch211-153b23.5", "si:ch211-214p16.1", "rap1b",
    "blf", "Sep-12", "myh9a", "fam212ab", "plek", "lpp", "srgn", "pfn1",
    "hsp70l", "hacd3", "si:ch211-250g4.3", "fermt3b", "blvrb", "atp2b1a",
    "cbx7a", "hsp70.1", "hsp70.2", "crema", "zgc:171686", "mpl", "tpte", "lmo2",
    "arpc1b", "grasp", "hsp70.3", "myh11a", "si:ch211-222l21.1", "bmp16"),
  ## markers for Cardiomyocytes (dediff.)
  "CM_GENERAL" = c(
    "nppa", "nppb", "hsp90aa1.1", "ttn.2", "ttn.1", "myh6", "atp2a2a", "xirp1",
    "si:ch211-131k2.3", "hspb11", "her4.1", "mybpc3", "lmo7a",
    "si:ch211-270g19.5", "unc45b", "tcap", "ryr2b", "tnnc1b", "synpo2lb",
    "ldb3a", "pla1a", "sorbs1", "bves", "flnca", "ACTC1", "slc8a1a", "ndrg4",
    "trdn", "myh7l", "pdk2a", "pabpc4", "myom1b", "cox7c", "hspb8", "cited4b",
    "ATP5MD", "desma", "cox6a2", "zgc:101840", "mt-co1", "bag3", "ITGB1BP2",
    "kcnh6a", "tnni1b", "tnnt2a", "laptm4b", "smtnl1", "aldoab", "pygmb",
    "cox7a1"),
  ## markers for Cardiomyocytes (Atrium)
  "CM_ATRIUM" = c(
    "tnnc1b", "myh6", "si:ch211-270g19.5", "tcap", "tnni1b", "aldoab", "tnnt2a",
    "smtnl1", "cmlc1", "gapdh", "tnnc1a", "nppa", "ddx54", "ckma", "zgc:101840",
    "mybphb", "myl7", "tpm4a", "actc1a", "acta1b", "idh2", "cox5b2", "mpc1",
    "mdh2", "uqcrfs1", "mdh1aa", "ldhba", "ak1", "ptp4a3", "nmrk2", "atp5mc3b",
    "atp5mc1", "calm3a", "atp5pd", "atp5meb", "atp5f1b", "chchd10", "adprhl1",
    "fhl2a", "ckbb", "atp5mc3a", "lmod2b", "si:dkey-51e6.1", "cox5aa",
    "ndufb10", "mybpc3", "zgc:85722", "slc25a3b", "atp5if1b", "hbaa1"),
  ## markers for Cardiomyocytes (Ventricle)
  "CM_VENTRICLE" = c(
    "tnni4a", "CR926459.1", "fabp3", "myl7", "ak1", "cmlc1", "actc1a", "acta1b",
    "myh7l", "ndufa4", "gapdh", "ckma", "eno3", "tnnt2a", "atp5pd", "cox5b2",
    "cox7a1", "aldoaa", "tnnc1a", "tnni1b", "tpm4a", "cox4i1l", "mdh2",
    "ppdpfa", "mustn1b", "atp5mc1", "idh2", "nme2b.2", "rbpms2b", "pgam2",
    "ldb3b", "hspb11", "cox6b1", "atp5meb", "cox6a2", "ATP5MD", "mpc1", "gamt",
    "atp5po", "atp5mc3b", "zgc:193541", "tpi1b", "uqcr10", "slc25a5", "csrp3",
    "atp5if1b", "atp5f1d", "atp5f1e", "cox7b", "mdh1aa"),
  ## markers for Endocardium (Atrium)
  "ENDO_ATRIUM" = c(
    "zgc:158343", "spock3", "im:7152348", "ptgs2a", "ptgs2b", "id2b", "vcam1b",
    "aqp8a.1", "mycb", "her6", "krt18", "grb10a", "clic2", "si:ch211-248e11.2",
    "gpr182", "klf2a", "cxcl18b", "fosl1a", "tmem88b", "dusp5", "glulb",
    "si:ch73-335l21.4", "akap12b", "ctsla", "fosb", "ackr4b", "igf2b", "bzw1b",
    "si:ch211-153b23.5", "tmem88a", "zfp36l1b", "ramp2", "errfi1a", "ponzr1",
    "fosab", "si:ch211-145b13.6", "btg2", "junba", "cyp2ad2", "nppc", "lpl",
    "slc22a31", "sat1a.2", "jun", "junbb", "tmem2", "vwf", "kctd12.2", "spint2",
    "pim1"),
  ## markers for Endocardium (Ventricle)
  "ENDO_VENTRICLE" = c(
    "si:ch73-86n18.1", "aqp8a.1", "spock3", "mb", "fabp11a", "epas1b", "ramp2",
    "si:ch211-145b13.6", "f8", "her6", "cdh5", "fosab", "si:ch211-248e11.2",
    "eppk1", "socs3a", "tmem88b", "podxl", "cyp1b1", "grb10a",
    "si:ch73-335l21.4", "lpl", "jun", "si:ch211-160j14.3", "spon1b", "tmem88a",
    "akap12b", "klf2b", "id2b", "zgc:64106", "egfl7", "ecscr", "tspan4b",
    "fam174b", "kdrl", "junba", "ackr4b", "zgc:152791", "si:dkey-261h17.1",
    "si:dkeyp-97a10.2", "junbb", "nppc", "slc22a31", "klf6a", "clic2", "appa",
    "flt1", "tfpia", "ier2b", "im:7152348", "ccdc187"),
  ## markers for Epicardium (Atrium)
  "EPIC_ATRIUM" = c(
    "gstm.3", "s100a10a", "fn1b", "krt15", "mmp2", "si:ch211-105c13.3", "krt4",
    "frzb", "mmel1", "anxa2a", "ppdpfa", "krt94", "rbp4", "si:ch211-286o17.1",
    "anxa1a", "CR318588.3", "dkk3b", "tcf21", "postnb", "slc29a1b", "tmsb1",
    "phactr4a", "podxl", "c6.1", "myh10", "tuba1a", "f3b", "timp2a", "tbx18",
    "endouc", "cxcl14", "ckba", "aldh1a2", "col6a1", "c1qtnf5", "hsd11b2",
    "c7a", "si:ch211-156j16.1", "id2a", "rn7sk", "acvrl1", "spaca4l", "mxra8b",
    "mfap5", "s100a10b", "c4", "lxn", "flnb", "olfml3b", "cav1"),
  ## markers for Epicardium (Ventricle)
  "EPIC_VENTRICLE" = c(
    "si:ch211-106h4.12", "gstm.3", "krt15", "postnb", "s100a10a", "fn1b",
    "endouc", "mmp2", "tfa", "krt4", "krt94", "frzb", "si:ch211-105c13.3",
    "cldnc", "tmsb1", "si:ch211-198c19.3", "adh8a", "si:ch1073-459j12.1",
    "mdka", "c6.1", "anxa2a", "col1a2", "podxl", "si:ch211-137i24.10",
    "col1a1a", "spaca4l", "igfbp5b", "col1a1b", "rbp4", "eppk1", "sparc",
    "anxa1a", "dkk3b", "cygb1", "col6a3", "tgm2b", "si:dkey-8k3.2", "sostdc1a",
    "mgp", "aldh1a2", "c4", "crip1", "zgc:123068", "ppdpfa", "hsd11b2", "fn1a",
    "col5a1", "col6a1", "colec12", "si:ch211-286o17.1"),
  ## markers for Macrophages
  "MACROPHAGES" = c(
    "grn1", "grn2", "ccl35.1", "cd74a", "c1qb", "ctsd", "c1qc", "lygl1",
    "lgals3bpb", "si:ch211-147m6.2", "marco", "fcer1gl", "ctss2.2", "mfap4",
    "NPC2 (1 of many)", "si:ch211-147m6.1", "cd74b", "grna",
    "si:busm1-266f07.2", "ctss2.1", "atp6v0ca", "hmox1a", "zgc:174904",
    "si:dkey-5n18.1", "ctsba", "rasgef1ba", "lgmn", "c1qa", "psap", "sqstm1",
    "CR855311.1", "ctsz", "pfn1", "atp6v1g1", "si:ch211-212k18.7", "tspan36",
    "ctsc", "si:dkey-27i16.2", "mrc1b", "cd9b", "zgc:92066", "ccl34a.4", "tppp",
    "cotl1", "ncf1", "cfp", "laptm5", "coro1a", "si:ch211-194m7.4", "scpep1"),
  ## markers for Monocytes
  "MONOCYTES" = c(
    "ccl35.1", "epdl1", "si:ch211-214p16.1", "si:busm1-266f07.2",
    "si:ch211-214p16.2", "ccl35.2", "CU459094.3", "cd74a", "arpc1b", "cd74b",
    "tmsb1", "anxa3b", "sftpbb", "id2a", "mpeg1.1", "mhc2dab", "cxcr4b",
    "si:dkey-33i11.9", "zgc:92066", "pfn1", "hspbp1", "DNAJA4",
    "si:dkey-27i16.2", "stip1", "cxcl11.6", "aif1l", "hspa4a", "nr4a3",
    "syngr3b", "sqstm1", "zgc:152863", "rhbdd1", "hspe1", "BX004816.2",
    "pdcd4b", "zgc:103700", "chchd2", "actb2", "zdhhc18b", "psap", "rasgef1ba",
    "ctss2.2", "zgc:64051", "spi1b", "zgc:136870", "coro1a", "hspa8",
    "hsp90aa1.2", "wbp2", "ckbb"),
  ## markers for Neuronal cells
  "NEURONAL_CELLS" = c(
    "vipb", "elavl4", "stmn1b", "syt1a", "tuba1c", "gap43", "snap25a", "stmn2a",
    "sncga", "rtn1b", "mllt11", "crhbp", "slc5a7a", "anxa13l", "zgc:65894",
    "stmn2b", "ddc", "chata", "gng3", "FO907089.1", "phox2a", "pcsk1nl",
    "scg2b", "atp1b2a", "cpe", "nsg2", "vip", "vgf", "hpcal4", "atp1a3a",
    "cplx2", "phox2bb", "prph", "synpr", "map1b", "maptb", "rab6bb", "tuba2",
    "syt11a", "chrna3", "stxbp1a", "si:ch73-119p20.1", "ret", "vamp2", "cplx2l",
    "cntnap2a", "sypa", "eef1a1a", "stx1b", "elavl3"),
  ## markers for Neutrophils
  "NEUTROPHILS" = c(
    "lyz", "lect2l", "BX908782.2", "npsn", "si:ch211-117m20.5", "mmp13a.1",
    "si:ch211-9d9.1", "scpp8", "si:ch1073-67j19.1", "mmp9", "pnp5a", "mpx",
    "ncf1", "fcer1gl", "pfn1", "scinlb", "cfl1l", "cotl1", "plekhf1", "cpa5",
    "illr4", "myh9a", "adam8a", "si:dkey-102g19.3", "lta4h", "timp2b", "il6r",
    "si:ch73-54b5.2", "cfbl", "vsir", "alox5ap", "nccrp1", "arg2",
    "si:ch211-284o19.8", "CU499336.2", "si:ch1073-443f11.2", "s1pr4", "arpc1b",
    "DNAJA4", "antxr2a", "si:ch211-136m16.8", "tmsb1", "si:ch73-343l4.8",
    "srgn", "zgc:77112", "coro1a", "stip1", "cd7al", "si:ch73-248e21.5",
    "CT030188.1"),
  ## markers for Perivascular cells
  "PERIVASCULAR_CELLS" = c(
    "TCIM", "cxcl12b", "pdgfrb", "rasl12", "BX901920.1", "kcne4", "rgs5b",
    "rgs4", "agtr2", "TPM1", "si:ch211-198c19.3", "ifitm1", "fxyd6l",
    "slc20a1a", "calm2b", "notch3", "cd248a", "tagln", "si:ch73-193i22.1",
    "col18a1b", "mcamb", "atp1a1b", "stk17al", "rgs5a", "mustn1a", "myh11a",
    "pcdh18a", "plp1b", "fhl3b", "sh3d21", "ccdc3b", "prx", "snai1a", "sparc",
    "si:ch211-251b21.1", "si:ch211-270g19.5", "slc7a2", "myo1b", "jam3a",
    "arhgap36", "AL954322.2", "BX908750.1", "atp1b4", "en1b", "CABZ01068367.1",
    "oscp1a.1", "ppp1r14aa", "rftn2", "rassf4", "spdya"),
  ## markers for Proliferating cells
  "PROLIFERATING_CELLS" = c(
    "pcna", "stmn1a", "hmgb2b", "DUT", "sumo3b", "rrm2.1", "si:ch211-288g17.3",
    "rpa3", "banf1", "rpa2", "lig1", "fen1", "sfrp5", "dnajc9", "rbbp4", "tyms",
    "rrm1", "top2a", "nutf2l", "cks1b", "hells", "chaf1a", "selenoh", "mki67",
    "zgc:110540", "mcm5", "mcm6", "ubr7", "si:ch211-156b7.4", "smc2",
    "si:ch211-66k16.27", "mcm4", "mcm3", "rrm2", "cdk1", "dek", "rpa1", "h2afx",
    "mcm2", "nasp", "cks2", "gins2", "asf1ba", "dhfr", "prim2", "rfc5",
    "rnaseh2a", "chtf8", "si:dkey-6i22.5", "dctd"),
  ## markers for Smooth muscle cells
  "SMOOTH_MUSCLE_CELLS" = c(
    "rgs5a", "acta2", "TPM1", "C11orf96", "tagln", "myh11a", "krt91", "itih1",
    "myl6", "myl9a", "BX663503.3", "vim", "si:dkey-164f24.2", "ctgfa", "lmod1b",
    "flna", "cfd", "tpm2", "si:busm1-57f23.1", "zgc:77517", "mdkb", "tns1b",
    "si:ch211-1a19.3", "mylka", "gch2", "slmapb", "cyr61", "lpp", "dnajb4",
    "CR855337.1", "fhl1a", "fhl3b", "fbln5", "alcama", "tmsb4x", "c2cd4a",
    "htra1a", "phlda2", "anxa1a", "cav2", "cst3", "zgc:92162", "f3a", "sort1a",
    "cald1b", "kng1", "epas1a", "wfdc2", "tpm4b", "fhl2a"),
  ## markers for T-cells
  "T_CELLS" = c(
    "si:ch211-214p16.1", "ccl34b.4", "ccl36.1", "ccl38.6", "DNAJA4", "ccr9a",
    "CR936442.1", "cxcr4b", "hspbp1", "zgc:64051", "rgs13", "pfn1", "nkl.1",
    "pdcd4b", "hspa4a", "sla2", "cebpb", "cyb5a", "coro1a", "stip1", "b2m",
    "si:dkey-27i16.2", "si:dkey-109a10.2", "si:ch211-202a12.4", "chchd2",
    "hsp90aa1.2", "rhbdd1", "cxcr4a", "tnfb", "calm1b", "tnfrsf9b", "ptprc",
    "si:ch211-105j21.9", "dusp2", "arpc1b", "sh2d1ab", "FP236356.1", "wasb",
    "BX120005.1", "capgb", "si:dkey-262k9.2", "si:ch211-165d12.4", "BX005105.1",
    "laptm5", "oser1", "mcl1a", "wbp2", "fkbp4", "scinlb", "si:dkey-177p2.6")
)

annotate_clusters <- function(cluster_markers, cell_type_markers) {
  annotations <- sapply(unique(cluster_markers$cluster), function(cluster) {
    # Get the top 100 genes for the cluster
    cluster_genes <- cluster_markers %>% 
      filter(cluster == !!cluster) %>%
      pull(gene)
    
    # Score each cell type marker set by the number of top 100 genes that match
    scores <- sapply(cell_type_markers, function(markers) {
      sum(markers %in% cluster_genes)
    })
    
    # Return the cell type with the highest score
    names(which.max(scores))
  })
  
  # Return the annotations
  setNames(annotations, unique(cluster_markers$cluster))
}

# Get the cell type annotations
cluster_annotations <- annotate_clusters(markers, cell_type_markers)

# Assign cell type based on cluster IDs
cell_type_annotations <- cluster_annotations[
  as.character(
    integrated$seurat_clusters
  )
]

# Ensure the cell barcodes match the integrated data
names(cell_type_annotations) <- colnames(integrated)

# Add metadata to the Seurat object
integrated$cell_type <- cell_type_annotations

head(integrated@meta.data)

# Get the expression matrix (normalized data)
expression_matrix <- GetAssayData(
  integrated,
  layer = "data",
  assay = "RNA"
)

signature_matrix <- sapply(unique(integrated$cell_type), function(ct) {
  cells_of_type <- colnames(integrated)[integrated$cell_type == ct]
  rowMeans(expression_matrix[, cells_of_type, drop = FALSE])
})

# Convert to a data frame
signature_matrix <- as.data.frame(signature_matrix)

# Add Gene column
signature_matrix$Gene <- rownames(signature_matrix)

signature_matrix <- signature_matrix %>%
  select(Gene, everything())

# View the first few rows
head(signature_matrix)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Save the signature matrix
write.table(
  signature_matrix,
  file = "output/adult_signature_matrix.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ZFIN Biomarkers
##----------------------------------------------------------------------------

zfin_biomarkers <- list(
  "AVC" = c(
    "bmp4", "vcana", "notch1b", "tbx2b", "has2", "klf2a", "cdh5", "spp1",
    "snai1b", "nfatc1", "tbx2a", "klf2b", "myl7", "rhoua", "tbx3a", "alcama",
    "aldh1a2", "anxa5b", "apbb1ip", "crip2", "dchs1b", "egr1", "egr3", "fbln1",
    "fn1b", "foxn4", "gja3", "hapln1a", "hcn4", "hey2", "hnf1ba", "hspb1",
    "hspb12", "hspb7", "id4", "itga5", "itgb1b", "lama4", "lama5", "lamb1a",
    "lamb1b", "lamb2", "lamc1", "nppa", "pdcd4b", "pdgfaa", "pou6f1", "prkd2",
    "prss23", "si:ch73-105b23.6", "tbx1", "tbx20", "tbx5a", "thbs1a", "tln1",
    "twist1b", "vclb", "wnt2bb", "wnt9a"
  ),
  "ATRIAL_CM" = c(
    "myh6", "myl7", "nppa", "add3a", "bmp4", "cdh2", "cemip2", "gars1", "hspb12",
    "hspb7", "meis2b", "plcg1"
  ),
  "ATRIAL_ENDO" = c(
    "aldh1a2", "cdh5", "esama", "gata5", "nfatc1", "plcg1", "rxrab"
  ),
  "EPICARDIUM" = c(
    "tcf21", "wt1a", "tbx18", "adma", "cldn11a", "cxadr", "cxcl12a", "cxcr4b",
    "elnb", "jam2b", "ptprc", "wtip"
  ),
  "VENTRICAL_CM" = c(
    "myh7", "myl7", "nppa", "bmp4", "itgav", "itgb3b", "add3a", "bmp10l", "cdh2",
    "cemip2", "crip2", "hspa5", "hspb12", "hspb7", "nrp2b", "pkd2", "plcg1", "tbx2b"
  ),
  "TRABECULAE" = c(
    "col11a1a", "col11a2", "bglap", "col2a1a", "hey2", "mgp", "mstnb", "nppa"
  ),
  "VENTRICAL_ENDO" = c(
    "notch1b", "cdh5", "cemip2", "esama", "gata5", "nfatc1", "plcg1"
  ),
  "SINOATRIAL_RING" = c(
    "fzd9b", "isl1a", "shox2"
  )
)

zfin_annotations <- annotate_clusters(markers, zfin_biomarkers)

zfin_cell_type <- cluster_annotations[
  as.character(
    integrated$seurat_clusters
  )
]

names(zfin_cell_type) <- colnames(integrated)

integrated$zfin_cell_type <- zfin_cell_type

## ZFIN Signature Matrix
##------------------------------------------------------------------------------

signature_matrix <- sapply(unique(integrated$zfin_cell_type), function(ct) {
  cells_of_type <- colnames(integrated)[integrated$zfin_cell_type == ct]
  rowMeans(expression_matrix[, cells_of_type, drop = FALSE])
})

# Convert to a data frame
signature_matrix <- as.data.frame(signature_matrix)

# Add Gene column
signature_matrix$Gene <- rownames(signature_matrix)

signature_matrix <- signature_matrix %>%
  select(Gene, everything())

# View the first few rows
head(signature_matrix)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Save the signature matrix
write.table(
  signature_matrix,
  file = "output/adult_signature_matrix_zfin.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)