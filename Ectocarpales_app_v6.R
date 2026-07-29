# ============================
# Phylogenetic Tree Explorer
# ============================

library(shiny)
library(readr)
library(dplyr)
library(ape)
library(phangorn)
library(ggplot2)
library(sf)
library(leaflet)
library(RColorBrewer)
library(viridis)
library(familiar) # masks is.waive() dependency

# -------------------------
# Preload your sequence + metadata behind the app
# -------------------------
  fasta_all <- read.FASTA("data/Master_Ectocarpales_alignment_15i26.fasta")
meta_all <- read_csv("data/Master_data_cluster_biogeogr_5xii25.csv")
regions_shp <- st_read("data/regions.shp")

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  titlePanel("ᚠ Norwegian Ectocarpales Explororer ᚠ"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$div("The Fehu rune symbol ᚠ means wealth, prosperity, and abundance,
               a reflection of the dominance of Ectocarpales in terms of brown algal species diversity",
               style = "margin-bottom: 10px; color: #555; font-style: italic;"),
      h3("Filters"),
      
      selectInput("filter_datatype", "Data Type:", choices = NULL, multiple = TRUE),
      selectInput("filter_family", "Family:", choices = NULL, multiple = TRUE),
      selectInput("filter_genus", "Genus:", choices = NULL, multiple = TRUE),
      selectInput("filter_species", "Updated Species ID:", choices = NULL, multiple = TRUE),
      
      checkboxInput("enable_downsample", "Downsample sequences by species", FALSE),
      
      conditionalPanel(
        condition = "input.enable_downsample == true",
        
        sliderInput(
          "seq_per_species",
          "Sequences per species:",
          min = 1,
          max = 10,
          value = 3,
          step = 1
        ),
        
        checkboxInput("prefer_nti", "Prioritize Norwegian Taxonomic Initiative samples", TRUE)
      ),
      
      selectInput(
        "dist_model",
        "DNA Distance Model:",
        choices = c("JC69", "F84", "TN93"),
        selected = "TN93"
      ),
      
      actionButton("compute_tree", "Compute tree"),
      
      hr(),
      helpText("Filters update automatically based on previous selections."),
      hr(),
      strong("Filtering info:"),
      verbatimTextOutput("debug")
    ),
    
    mainPanel(
      width = 8,
      h3("Neighbour-joining phylogenetic tree"),
      tags$div("The neighbour joining tree uses coxI data and selected distance model, calculated using the dna.dist function of ape R package.
               An asterisk indicates a specimen sampled through the Norwegian Taxonomic Initiative (2022-2024).
               If plot looks messy (adjusting window, data filters), simply hit compute to reset with selected inputs.
               Note the coxI data are reliable for clustering species, but do not accurately reflect intergeneric relationships. The tree is not bootstrapped or rooted.
               Some genbank accessions species IDs were updated, see species metadata below for original IDs",
               style = "margin-bottom: 10px; color: #555; font-style: italic;"),
      plotOutput("tree_plot", height = 800),
      h4("Download Options"),
      numericInput("dl_width", "Width (px):", value = 1600, min = 400, max = 8000),
      numericInput("dl_height", "Height (px):", value = 1200, min = 400, max = 8000),
      selectInput("dl_format", "File format:",
                  choices = c("PNG", "PDF"),
                  selected = "PNG"),
      downloadButton("download_tree", "Download Tree"),
      h3("Biogeography of species"),
      tags$div("The map displays the number of species in the tree for broad oceanic basins",
               style = "margin-bottom: 10px; color: #555; font-style: italic;"),
          leafletOutput("map", height = 400),
      h3("Specimen metadata for sequences used in the tree"),
      tags$div("Genetic clusters were determined using a 2% threshold,
               determined using cophenetic.phylo [ape] and hclust[base R] R functions",
               style = "margin-bottom: 10px; color: #555; font-style: italic;"),
      dataTableOutput("nti_table")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # -------------------------
  # --- Initialize cascading filters ---
  # -------------------------
  
  # Helper to safely provide choices
  safe_choices <- function(x) if(length(x)) x else x
  
  # Initialize Data_type after UI is ready
  session$onFlushed(function() {
    datatypes <- sort(unique(meta_all$Data_type))
    updateSelectInput(session, "filter_datatype",
                      choices = safe_choices(datatypes),
                      )
  }, once = TRUE)
  
  # Families depend on Data_type
  observeEvent(input$filter_datatype, {
    req(input$filter_datatype)
    
    fams <- meta_all %>%
      filter(Data_type %in% input$filter_datatype) %>%
      pull(Family) %>% unique() %>% sort()
    
    updateSelectInput(session, "filter_family",
                      choices = safe_choices(fams),
                      )
  }, ignoreNULL = TRUE)
  
  # Genera depend on Family
  observeEvent(input$filter_family, {
    req(input$filter_family)
    
    genera <- meta_all %>%
      filter(Family %in% input$filter_family,
             if(!is.null(input$filter_datatype)) Data_type %in% input$filter_datatype else TRUE
             ) %>%
      pull(Genus) %>% unique() %>% sort()
    
    updateSelectInput(session, "filter_genus",
                      choices = safe_choices(genera),
                      )
  }, ignoreNULL = TRUE)
  
  # Species depend on Genus
  observeEvent(input$filter_genus, {
    req(input$filter_genus)
    
    spp <- meta_all %>%
      filter(Genus %in% input$filter_genus,
            if(!is.null(input$filter_family)) Family %in% input$filter_family else TRUE,
            if(!is.null(input$filter_datatype)) Data_type %in% input$filter_datatype else TRUE
            ) %>%
    pull(Updated_Species_ID) %>% unique() %>% sort()
    
    updateSelectInput(session, "filter_species",
                      choices = safe_choices(spp)
                      )
  }, ignoreNULL = TRUE)
  
  # -------------------------
  # --- Align metadata to sequences
  # -------------------------
  aligned_data <- reactive({
    seqs <- fasta_all
    meta <- meta_all
    meta_filtered <- meta %>% filter(seq_id %in% names(seqs))
    validate(need(nrow(meta_filtered) > 0, "No metadata matches sequence names."))
    meta_filtered <- meta_filtered[match(names(seqs), meta_filtered$seq_id), ]
    list(seqs = seqs, meta = meta_filtered)
  })
  
  # -------------------------
  # --- Filter sequences based on selections
  # -------------------------
  filtered_subset <- reactive({
    req(aligned_data())
    req(input$compute_tree)
    
    meta <- aligned_data()$meta
    seqs <- aligned_data()$seqs
    
    meta_f <- meta %>%
      filter(
        (is.null(input$filter_datatype) | Data_type %in% input$filter_datatype),
        (is.null(input$filter_family)   | Family %in% input$filter_family),
        (is.null(input$filter_genus)    | Genus %in% input$filter_genus),
        (is.null(input$filter_species)  | Updated_Species_ID %in% input$filter_species)
      )
    
    seqs_f <- seqs[names(seqs) %in% meta_f$seq_id]
    
    validate(need(length(seqs_f) > 0, "No sequences remain after filtering."))
    
    list(seqs = seqs_f, meta = meta_f)
  })
  
  # -------------------------
  # --- Downsample sequences based on selections
  # -------------------------
  downsample_sequences <- function(df, n_per_species = 3, prefer_nti = TRUE) {
    
    df %>%
      group_by(Updated_Species_ID) %>%
      group_modify(~{
        
        if (prefer_nti) {
          # Split into NTI vs non-NTI
          nti <- .x %>% filter(Data_type == "Norwegian Taxonomic Initiative")
          non_nti <- .x %>% filter(Data_type != "Norwegian Taxonomic Initiative")
          
          # 1. Take NTI first
          selected <- nti %>% slice_sample(n = min(nrow(nti), n_per_species))
          
          # 2. Fill remaining positions
          n_remaining <- n_per_species - nrow(selected)
          if (n_remaining > 0 && nrow(non_nti) > 0) {
            filler <- non_nti %>% slice_sample(n = min(n_remaining, nrow(non_nti)))
            selected <- bind_rows(selected, filler)
          }
          
        } else {
          # No priority — plain random sampling by species
          selected <- .x %>% slice_sample(n = min(nrow(.x), n_per_species))
        }
        
        selected
      }) %>%
      ungroup()
  }
  
  downsampled_seqs <- reactive({
    req(filtered_subset())
    
    seqs <- filtered_subset()$seqs
    meta <- filtered_subset()$meta
    
    # If the user does NOT downsample: return full set
    if (!input$enable_downsample) {
      return(list(seqs = seqs, meta = meta))
    }
    
    # Downsample metadata by species
    meta_ds <- downsample_sequences(
      df = meta,
      n_per_species = input$seq_per_species,
      prefer_nti = input$prefer_nti
    )
    
    # Reduce sequences to those selected in metadata
    seqs_ds <- seqs[names(seqs) %in% meta_ds$seq_id]
    
    list(seqs = seqs_ds, meta = meta_ds)
  })
  
  # -------------------------
  # --- Debug info
  # -------------------------
  output$debug <- renderPrint({
    nseq <- tryCatch(length(filtered_subset()$seqs), error = function(e) NA)
    list(
      filter_datatype = input$filter_datatype,
      filter_family   = input$filter_family,
      filter_genus    = input$filter_genus,
      filter_species  = input$filter_species,
      n_sequences_after_filtering = nseq
    )
  })
  
  # -------------------------
  # --- Species counts per region
  # -------------------------
  region_counts <- eventReactive(input$compute_tree, {
    req(regions_shp, filtered_subset())
    meta <- filtered_subset()$meta
    counts <- meta %>%
      group_by(Region) %>%
      summarise(species_count = n_distinct(Updated_Species_ID), .groups = "drop")
    
    sf_counts <- left_join(regions_shp, counts, by = "Region")
    sf_counts$species_count[is.na(sf_counts$species_count)] <- 0
    sf_counts
  })
  
  # -------------------------
  # --- Build NJ tree
  # -------------------------
  tree_data <- eventReactive(input$compute_tree, {
    req(downsampled_seqs())
    req(input$compute_tree)
    req(input$dist_model)
    
    seqs <- downsampled_seqs()$seqs
    meta <- downsampled_seqs()$meta
    
    if(length(seqs) > 500) {
      validate(need(FALSE, paste0("Too many sequences (", length(seqs),
                                  "). Narrow your filters or click compute.")))
    }
    
    # Ensure meta rows match the order of sequences
    meta <- meta[match(names(seqs), meta$seq_id), ]
    stopifnot(all(names(seqs) == meta$seq_id))  # sanity check
    
    # Compute distances and NJ tree
    dist_mat <- dist.dna(seqs, model = input$dist_model)
    tree <- nj(dist_mat)
    tree$edge.length[tree$edge.length < 1e-8] <- 1e-8
    tree <- phangorn::midpoint(tree)
    tree <- ape::ladderize(tree, right = TRUE)
    
    # Attach metadata for plotting
    list(tree = tree, meta = meta)
  })
  
  # -------------------------
  # --- Reusable render tree function
  # -------------------------
  plot_phylo_tree <- function(tree, meta, tip_cex = 1.2, symbol_cex = 1.5) {
    req(tree_data(), filtered_subset())
    
    tree <- tree_data()$tree
    meta <- tree_data()$meta
    
    # Align metadata to tree tips (ensures order matches)
    tip_meta <- meta[match(tree$tip.label, meta$seq_id), ]
    stopifnot(all(tree$tip.label == tip_meta$seq_id))  # sanity check
    
    # Replace tip labels with Tree_label
    tree$tip.label <- tip_meta$Tree_label
    
    # (Modify the condition if NTI is stored differently)
    tip_shapes <- ifelse(tip_meta$Data_type == "Norwegian Taxonomic Initiative", 22, 21)
    
    # --- Color mapping for symbols based on Region ---
    regions <- unique(tip_meta$Region)
    symbol_colors <- viridis(length(regions))
    names(symbol_colors) <- regions
    tip_symbols_colors <- symbol_colors[tip_meta$Region]
    
    # Plot tree
    par(mar = c(5, 2, 2, 8))
    plot(tree,
         tip.color = "black",
         cex = 1.2,
         label.offset = 0.001,
         use.edge.length = TRUE)
    
    # Add colored symbols at tips (now with shape & color)
    tiplabels(
      pch = tip_shapes,
      bg = tip_symbols_colors,
      col = "black",
      cex = 1.5
    )
    
    # Add legend: 1) Regions (color), 2) NTI vs GenBank (shape)
    legend("bottomleft",
           legend = regions,
           pch = 21,
           pt.bg = symbol_colors,
           pt.cex = 1.5,
           cex = 1.2,
           title = "Region",
           inset = c(0.01, 0))
    
    legend("topleft",
           legend = c("Norwegian Taxonomic Initiative", "GenBank"),
           pch = c(22, 21),     # square, circle
           pt.bg = "white",
           pt.cex = 1.5,
           cex = 1.2,
           inset = c(0.01, 0))
    
  }
  
  # -------------------------
  # --- Render Tree
  # -------------------------
  output$tree_plot <- renderPlot({
    req(tree_data())
    tree <- tree_data()$tree
    meta <- tree_data()$meta
    plot_phylo_tree(tree, meta)
  })
  
  # -------------------------
  # --- Download handler for trees
  # -------------------------
  output$download_tree <- downloadHandler(
    
    filename = function() {
      ext <- ifelse(input$dl_format == "PNG", "png", "pdf")
      paste0("Ectocarpales_NTI_phylogeny_", Sys.Date(), ".", ext)
    },
    
    content = function(file) {
      req(tree_data())
      
      width  <- input$dl_width
      height <- input$dl_height
      
      if (input$dl_format == "PNG") {
        png(file, width = width, height = height, res = 150)
        plot_phylo_tree(tree_data())
        dev.off()
        
      } else {
        pdf(file, width = width/100, height = height/100)  
        # convert px → inches assuming 100px per inch
        plot_phylo_tree(tree_data())
        dev.off()
      }
    }
  )
  
  # -------------------------
  # --- Render biogeography map + site locations
  # -------------------------
  output$map <- renderLeaflet({
    req(region_counts(), tree_data())

    df <- tree_data()$meta
    
    sf_counts <- region_counts()
    
    # ---- Filter NTI sample locations ----
    nti_points <- df %>%
      filter(Data_type == "Norwegian Taxonomic Initiative",
             !is.na(Latitude), !is.na(Longitude)) %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    pal <- colorNumeric(palette = rev(magma(256)), domain = sf_counts$species_count)
    
    leaflet(sf_counts) %>%
      addProviderTiles("Esri.OceanBasemap") %>%
      addPolygons(
        fillColor = ~pal(species_count),
        color = "black",
        weight = 1,
        opacity = 0.7,
        fillOpacity = 0.5,
        label = ~paste0(Region, ": ", species_count, " species")
      ) %>%
      addCircleMarkers(
        data = nti_points,
        radius = 4,
        color = "#57106e",
        fillOpacity = 0.8,
        popup = ~paste0("Norwegian Taxonomic Initiative Sample<br>",
                        "Accession:Species: ", Tree_label, "<br>",
                        "Latitude: ", round(st_coordinates(geometry)[2], 4),
                        "<br>Longitude: ", round(st_coordinates(geometry)[1], 4))
      ) %>%
      addLegend(
        pal = pal,
        values = ~species_count,
        title = "Species count",
        position = "bottomright"
      )
  })
  
  # -------------------------
  # --- Render specimen table
  # -------------------------
  specimen_table <- reactive({
    req(tree_data(), filtered_subset())
    
    tree <- tree_data()$tree
    meta <- tree_data()$meta
    
    present_tips <- tree$tip.label
    
    # Return ALL specimens used in the tree
    meta %>%
      filter(seq_id %in% present_tips)
  })
  
  output$nti_table <- renderDataTable({
    req(specimen_table())
    specimen_table()
  }, options = list(
    pageLength = 10,
    scrollX = TRUE
  ))
}

# -------------------------
# Run the app
# -------------------------
shinyApp(ui, server)
