# Pathview Shiny App
# Visualize KEGG pathways from a DEG (differential expression) list.
#
# Run locally:   shiny::runApp()
# Deploy to web: see README.md (shinyapps.io / Posit Connect)

library(shiny)
library(bslib)
library(DT)
library(shinycssloaders)
library(colourpicker)
library(pathview)
library(clusterProfiler)
library(AnnotationDbi)

# ---- Organism configuration -------------------------------------------------
# Maps a friendly name to the KEGG organism code and its Bioconductor OrgDb.
# Add more rows here (and install the matching org.*.eg.db) to support others.
ORGANISMS <- list(
  "Human (Homo sapiens)"      = list(kegg = "hsa", db = "org.Hs.eg.db"),
  "Mouse (Mus musculus)"      = list(kegg = "mmu", db = "org.Mm.eg.db"),
  "Rat (Rattus norvegicus)"   = list(kegg = "rno", db = "org.Rn.eg.db"),
  "Zebrafish (Danio rerio)"   = list(kegg = "dre", db = "org.Dr.eg.db"),
  "Fly (D. melanogaster)"     = list(kegg = "dme", db = "org.Dm.eg.db"),
  "Yeast (S. cerevisiae)"     = list(kegg = "sce", db = "org.Sc.sgd.db")
)

# Supported input gene ID types. pathview converts these to Entrez internally.
ID_TYPES <- c("SYMBOL", "ENTREZID", "ENSEMBL", "REFSEQ", "UNIGENE")

# ---- UI ---------------------------------------------------------------------
ui <- page_sidebar(
  title = "PathviewViz",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # Keep the colour-picker pop-up above neighbouring controls and stop the
  # sidebar's scroll container from clipping it.
  tags$head(tags$style(HTML("
    .colourpicker-panel { z-index: 1080 !important; }
    .bslib-sidebar-layout > .sidebar .sidebar-content { overflow: visible !important; }
  "))),

  sidebar = sidebar(
    width = 360,
    fileInput("deg_file", "DEG list (CSV/TSV)",
              accept = c(".csv", ".tsv", ".txt")),
    helpText("Needs a gene ID column and a numeric value column (e.g. log2FC)."),

    selectInput("organism", "Organism", choices = names(ORGANISMS)),
    selectInput("id_type", "Gene ID type", choices = ID_TYPES, selected = "SYMBOL"),

    selectInput("gene_col", "Gene ID column", choices = NULL),
    selectInput("value_col", "Value column (log2FC)", choices = NULL),

    textInput("pathway_id", "KEGG pathway ID(s)",
              placeholder = "e.g. 04110 or 04110,04151"),
    helpText("Digits only; organism prefix is added automatically."),

    radioButtons("render_mode", "Render style",
                 c("Native KEGG PNG" = "native", "Graphviz layout" = "graphviz"),
                 selected = "native"),
    sliderInput("limit", "Color scale limit (|value|)",
                min = 0.5, max = 10, value = 2, step = 0.5),

    tags$label("Node colors", class = "form-label fw-bold"),
    colourInput("col_low",  "Down-regulated", value = "#00B050",
                showColour = "both", width = "100%"),
    colourInput("col_mid",  "Midpoint",       value = "#DDDDDD",
                showColour = "both", width = "100%"),
    colourInput("col_high", "Up-regulated",   value = "#FF0000",
                showColour = "both", width = "100%"),
    helpText("Colors for down-regulated, midpoint, and up-regulated genes."),

    actionButton("run", "Render pathway", class = "btn-primary", icon = icon("play")),

    tags$hr(),
    tags$label("Suggest pathways (GSEA)", class = "form-label fw-bold"),
    numericInput("gsea_pcut", "Adj. p-value cutoff", value = 0.05,
                 min = 0, max = 1, step = 0.01),
    helpText("Ranks all mapped genes by the value column and runs KEGG GSEA.",
             "No threshold needed - the sign of the value sets the ranking direction."),
    actionButton("run_gsea", "Suggest top pathways", icon = icon("magnifying-glass"))
  ),

  navset_card_tab(
    nav_panel(
      "Pathway",
      uiOutput("status"),
      withSpinner(uiOutput("pathway_images")),
      uiOutput("download_ui")
    ),
    nav_panel(
      "Suggested pathways",
      uiOutput("gsea_status"),
      withSpinner(DTOutput("gsea_table"))
    ),
    nav_panel(
      "DEG table",
      DTOutput("deg_table")
    ),
    nav_panel(
      "Mapping summary",
      verbatimTextOutput("map_summary")
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {

  # Per-session working dir so concurrent web users don't clobber each other.
  work_dir <- file.path(tempdir(), paste0("pv_", session$token))
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  session$onSessionEnded(function() unlink(work_dir, recursive = TRUE))

  rendered <- reactiveVal(NULL)  # data frame of {pathway, png_path}

  # --- Read uploaded DEG file ---
  deg_data <- reactive({
    req(input$deg_file)
    path <- input$deg_file$datapath
    sep <- if (grepl("\\.tsv$|\\.txt$", input$deg_file$name, ignore.case = TRUE)) "\t" else ","
    df <- tryCatch(
      read.delim(path, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE),
      error = function(e) NULL
    )
    validate(need(!is.null(df) && ncol(df) >= 2,
                  "Could not parse the file, or it has fewer than 2 columns."))
    df
  })

  # Populate column selectors once a file is loaded.
  observeEvent(deg_data(), {
    cols <- names(deg_data())
    num_cols <- cols[vapply(deg_data(), is.numeric, logical(1))]
    updateSelectInput(session, "gene_col", choices = cols, selected = cols[1])
    updateSelectInput(session, "value_col", choices = cols,
                      selected = if (length(num_cols)) num_cols[1] else cols[min(2, length(cols))])
  })

  output$deg_table <- renderDT({
    datatable(deg_data(), options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE)
  })

  # --- Build the named numeric vector keyed by Entrez ID that pathview wants ---
  gene_vector <- reactive({
    df <- deg_data()
    req(input$gene_col, input$value_col)
    org <- ORGANISMS[[input$organism]]

    ids <- as.character(df[[input$gene_col]])
    vals <- suppressWarnings(as.numeric(df[[input$value_col]]))
    keep <- !is.na(ids) & ids != "" & !is.na(vals)
    ids <- ids[keep]; vals <- vals[keep]
    validate(need(length(ids) > 0, "No usable gene/value rows after cleaning."))

    if (!requireNamespace(org$db, quietly = TRUE)) {
      validate(need(FALSE, sprintf("OrgDb package '%s' is not installed.", org$db)))
    }
    db <- getFromNamespace(org$db, org$db)

    if (input$id_type == "ENTREZID") {
      entrez <- ids
    } else {
      mapped <- suppressMessages(AnnotationDbi::mapIds(
        db, keys = ids, column = "ENTREZID", keytype = input$id_type,
        multiVals = "first"))
      entrez <- unname(mapped[ids])
    }

    ok <- !is.na(entrez)
    validate(need(sum(ok) > 0,
                  "No genes could be mapped to Entrez IDs. Check ID type/organism."))

    v <- vals[ok]
    names(v) <- entrez[ok]
    # Collapse duplicate Entrez IDs by mean. tapply() yields a 1-D array;
    # pathview only accepts a plain named numeric vector, so coerce back.
    v <- tapply(v, names(v), mean)
    v <- stats::setNames(as.numeric(v), names(v))
    list(vec = v, n_in = length(ids), n_mapped = length(v))
  })

  output$map_summary <- renderPrint({
    gv <- gene_vector()
    cat("Organism:        ", input$organism, "\n")
    cat("Input ID type:   ", input$id_type, "\n")
    cat("Genes submitted: ", gv$n_in, "\n")
    cat("Mapped to Entrez:", gv$n_mapped, "\n")
    cat("Value range:     ", paste(round(range(gv$vec), 3), collapse = " to "), "\n")
  })

  # --- GSEA (suggest pathways) ---
  gsea_res <- reactiveVal(NULL)   # data frame of suggested pathways, or NULL

  observeEvent(input$run_gsea, {
    org <- ORGANISMS[[input$organism]]
    # GSEA needs every mapped gene ranked by the value, sorted decreasing.
    ranks <- sort(gene_vector()$vec, decreasing = TRUE)

    if (length(ranks) < 10) {
      gsea_res(list(error = sprintf(
        "Only %d genes mapped - GSEA needs a larger ranked list (~hundreds+).",
        length(ranks))))
      return()
    }

    res <- tryCatch(
      withProgress(message = "Running KEGG GSEA...", value = 0.5, {
        set.seed(1)  # fgsea is stochastic; fix for reproducibility
        gs <- clusterProfiler::gseKEGG(
          geneList = ranks, organism = org$kegg, keyType = "kegg",
          minGSSize = 10, maxGSSize = 500,
          pvalueCutoff = input$gsea_pcut, verbose = FALSE)
        if (is.null(gs) || nrow(as.data.frame(gs)) == 0) return(NULL)
        d <- as.data.frame(gs)
        d <- d[order(d$p.adjust, -abs(d$NES)), ]
        data.frame(
          Pathway     = d$Description,
          ID          = sub(paste0("^", org$kegg), "", d$ID),  # strip species prefix
          NES         = signif(d$NES, 3),
          `Adj. p`    = signif(d$p.adjust, 3),
          `Set size`  = d$setSize,
          Direction   = ifelse(d$NES > 0, "Up", "Down"),
          check.names = FALSE, stringsAsFactors = FALSE)
      }),
      error = function(e) list(error = conditionMessage(e)))

    if (is.null(res)) {
      gsea_res(list(error = "No gene sets passed the cutoff. Try a higher p-value cutoff."))
    } else {
      gsea_res(res)
    }
  })

  output$gsea_status <- renderUI({
    res <- gsea_res()
    if (is.null(res)) {
      return(div(class = "text-muted",
                 "Click 'Suggest top pathways' to rank KEGG pathways by GSEA of ",
                 "your value-ranked gene list. Positive NES = enriched in up-regulated ",
                 "genes, negative = down-regulated."))
    }
    if (!is.null(res$error)) {
      return(div(class = "alert alert-warning", res$error))
    }
    div(class = "text-muted",
        "Ranked by adjusted p-value. Click a row to load its ID into the pathway box.")
  })

  output$gsea_table <- renderDT({
    res <- gsea_res()
    req(is.data.frame(res))
    datatable(res, selection = "single", rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })

  # Clicking a suggested pathway fills the KEGG pathway ID field.
  observeEvent(input$gsea_table_rows_selected, {
    res <- gsea_res(); req(is.data.frame(res))
    pid <- res$ID[input$gsea_table_rows_selected]
    updateTextInput(session, "pathway_id", value = pid)
    showNotification(sprintf("Loaded pathway %s - click 'Render pathway'.", pid),
                     type = "message")
  })

  # --- Render on button click ---
  observeEvent(input$run, {
    org <- ORGANISMS[[input$organism]]
    ids_raw <- strsplit(gsub("\\s", "", input$pathway_id), ",")[[1]]
    ids_raw <- ids_raw[nzchar(ids_raw)]
    if (length(ids_raw) == 0) {
      showNotification("Enter at least one KEGG pathway ID (digits).", type = "error")
      return()
    }
    pathways <- gsub("[^0-9]", "", ids_raw)
    pathways <- pathways[nzchar(pathways)]

    gv <- gene_vector()$vec
    old_wd <- getwd(); setwd(work_dir); on.exit(setwd(old_wd), add = TRUE)

    native <- input$render_mode == "native"
    results <- lapply(pathways, function(pid) {
      tryCatch({
        suppressMessages(pathview(
          gene.data  = gv,
          pathway.id = pid,
          species    = org$kegg,
          gene.idtype = "ENTREZID",
          kegg.native = native,
          limit      = list(gene = input$limit, cpd = 1),
          low  = list(gene = input$col_low,  cpd = input$col_low),
          mid  = list(gene = input$col_mid,  cpd = input$col_mid),
          high = list(gene = input$col_high, cpd = input$col_high),
          out.suffix = "deg"
        ))
        # Native mode writes <species><pid>.deg.png; graphviz writes a PDF.
        f <- file.path(work_dir,
                       paste0(org$kegg, pid, if (native) ".deg.png" else ".deg.pdf"))
        if (!file.exists(f)) {
          return(list(pid = pid, display = NA_character_, download = NA_character_,
                      error = "No output produced (invalid pathway ID or no mapped genes?)."))
        }
        # Browsers can't show the graphviz PDF inline, so rasterize it for display.
        if (native) {
          display <- f
        } else {
          display <- file.path(work_dir, paste0(org$kegg, pid, ".deg.png"))
          pdftools::pdf_convert(f, format = "png", pages = 1, dpi = 150,
                                filenames = display, verbose = FALSE)
        }
        list(pid = pid, display = display, download = f, error = NA_character_)
      }, error = function(e)
        list(pid = pid, display = NA_character_, download = NA_character_,
             error = conditionMessage(e)))
    })
    rendered(results)
  })

  output$status <- renderUI({
    res <- rendered()
    if (is.null(res)) {
      return(div(class = "text-muted",
                 "Upload a DEG list, pick columns, enter a KEGG pathway ID, then Render."))
    }
    errs <- Filter(function(r) !is.na(r$error), res)
    if (length(errs)) {
      div(class = "alert alert-warning",
          lapply(errs, function(e) tags$div(sprintf("Pathway %s: %s", e$pid, e$error))))
    }
  })

  output$pathway_images <- renderUI({
    res <- rendered(); req(res)
    imgs <- Filter(function(r) !is.na(r$display), res)
    if (length(imgs) == 0) return(NULL)
    tagList(lapply(seq_along(imgs), function(i) {
      id <- paste0("pvimg_", i)
      output[[id]] <- renderImage({
        list(src = imgs[[i]]$display, contentType = "image/png",
             width = "100%", alt = imgs[[i]]$pid)
      }, deleteFile = FALSE)
      div(tags$h5(sprintf("Pathway %s", imgs[[i]]$pid)),
          imageOutput(id, height = "auto"), tags$hr())
    }))
  })

  output$download_ui <- renderUI({
    res <- rendered(); req(res)
    files <- Filter(function(r) !is.na(r$download), res)
    if (length(files) == 0) return(NULL)
    downloadButton("download_zip", "Download image(s)", class = "btn-success")
  })

  output$download_zip <- downloadHandler(
    filename = function() paste0("pathview_", Sys.Date(), ".zip"),
    content = function(file) {
      res <- rendered()
      files <- vapply(Filter(function(r) !is.na(r$download), res),
                      function(r) r$download, character(1))
      utils::zip(file, files, flags = "-j")
    }
  )
}

shinyApp(ui, server)
