# One-shot dependency install for the Pathview Shiny app.
cran <- c("shiny", "bslib", "DT", "shinycssloaders", "colourpicker", "pdftools", "BiocManager")
to_get <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")

bioc <- c("pathview", "clusterProfiler", "AnnotationDbi",
          "org.Hs.eg.db", "org.Mm.eg.db")  # add org.*.eg.db for more organisms
BiocManager::install(bioc, update = FALSE, ask = FALSE)
