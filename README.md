# PathviewViz

A Shiny app for visualizing KEGG pathways colored by your differential expression results. Upload a DEG table, pick a pathway, and see the genes painted by their fold change, with optional GSEA to help you decide which pathways are worth looking at.

## Demo

```{=html}
<!--
  The line below embeds demo.mp4 from this repo. GitHub plays MP4 inline in most
  browsers. If it ever shows up as a plain link instead of a player, the most
  reliable fix is to open the README in GitHub's web editor and drag demo.mp4
  into it — GitHub will replace the src with a guaranteed-playable attachment URL.
-->
```
<p align="center">

<video src="demo.mp4" controls width="100%">

Your browser can't play the video here, <a href="demo.mp4">download demo.mp4</a> to watch it.

</video>

</p>

## What it does

1.  Upload a DEG table (CSV or TSV) that has a gene-ID column and a numeric column such as `log2FC`.
2.  Choose the organism, the gene-ID type (SYMBOL, ENSEMBL, ENTREZID, and so on), and which two columns to use.
3.  Enter one or more KEGG pathway IDs, just the digits, for example `04110` for the cell cycle. Separate several with commas.
4.  Render. The KEGG map appears in the app colored by your values, and you can download it as a PNG.

You can also set the colors for up- and down-regulated genes, and run GSEA (gene set enrichment analysis) to get a ranked list of pathway names to explore. Behind the scenes, gene IDs are mapped to Entrez through the organism's `org.*.eg.db` package and then handed to pathview.

## Running it on your computer

The app runs locally in R; a browser tab opens automatically when it starts.

### 1. Install R

Download and install R from [cran.r-project.org](https://cran.r-project.org/). [RStudio Desktop](https://posit.co/download/rstudio-desktop/) is optional but makes launching the app a one-click affair.

### 2. Get the code

Clone the repository, or download it as a ZIP from GitHub and unzip it:

``` bash
git clone https://github.com/<your-username>/PathviewViz.git
cd PathviewViz
```

### 3. Install the dependencies (first time only)

From an R session started in the project folder, run:

``` r
source("install.R")
```

This installs the CRAN and Bioconductor packages the app needs (pathview, clusterProfiler, the annotation databases, and so on). It is a sizable download and can take several minutes the first time. The exact package list is also spelled out in [install.R](install.R) if you prefer to install things by hand.

### 4. Launch the app

Either open [app.R](app.R) in RStudio and click **Run App**, or start it from an R session:

``` r
shiny::runApp()      # run from inside the project folder
```

A browser window opens with the app. To stop it, press Esc in RStudio or Ctrl-C in the console.

### 5. Try it with the sample data

A small human example, `sample_deg.csv` (gene symbols with log2 fold changes), is included. Upload it, keep the defaults (Human, SYMBOL), and try pathway `04110` (cell cycle) or `04115` (p53 signaling).

## Notes

-   The first time you render a given pathway, the app downloads its data from KEGG and caches it for the rest of the session, so an internet connection is needed.
-   KEGG pathway IDs are the five-digit numbers; the organism prefix (hsa, mmu, and so on) is added for you.
-   Each session works in its own temporary folder, so several people can use a shared deployment without overwriting each other's images.

## 
