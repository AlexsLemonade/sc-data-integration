library(SingleCellExperiment)


citeseq_altexp_name <- "CITEseq"


sce <- readRDS("data/scpca/filtered_sce/SCPCP000007/SCPCS000221/SCPCL000295_filtered.rds")

if (citeseq_altexp_name %in% altExpNames(sce)) {
  
  cite_altexp <- altExp(sce, citeseq_altexp_name)

  # OSCA reference:
  #  http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#normalization
  
  # Note that although OSCA suggests to use the medianSizeFactors() function, it can perform poorly
  #  as described here under "Caveats": https://rdrr.io/bioc/scuttle/man/medianSizeFactors.html
  #  We therefore use librarySizeFactors() as described in these docs, which does not require a priori
  #  baseline calculation. The docs further explain that medianSizeFactors _is more appropriate_ for CITE-seq, 
  #  since those counts are "usually large enough to avoid zeroes." Alas, ours are not!
  
  # Calculate baseline:
  baseline <- DropletUtils::ambientProfileBimodal(cite_altexp)
  
  # Calculate size factors, and handle 0's
  size_factors <- scuttle::medianSizeFactors(cite_altexp, reference = baseline)
  
  # if there are 0s (there will never be negatives with median)
  if (sum(size_factors == 0) > 1) {
    
    # We can instead use geometricSizeFactors which should be more appropriate for antibody-derived tags;
    #  See Details section: https://rdrr.io/bioc/scuttle/man/geometricSizeFactors.html
    #size_factors <- scuttle::geometricSizeFactors(cite_altexp)
    
    # And/or, we can change all the 0's to ~0:
    #size_factors[size_factors == 0] <- 1e-100
    
  }
  
  
  scran::modelGeneVar(cite_altexp)
  
  # Perform normalization
  scater::logNormCounts(cite_altexp, size.factors = size_factors)
  
  

  
  
  sizeFactors(altExp(sce)) <- sf_amb
  
  # remove 0's?
  # https://github.com/edward130603/BayesSpace/issues/28#issuecomment-717547233
  

  
}

