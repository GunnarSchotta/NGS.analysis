library(shiny)
library(DESeq2)
library(ggplot2)
library(DT)
library(plotly)
library(shinyjs)
library(shinyFiles)

shinyApp(
  ui = fluidPage(
    useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        shinyFilesButton('rds', 'Choose RDS File', 'Select RDS File', multiple = FALSE),
        fluidRow(
          column(6,
            selectInput("select.study.a", label = ("Study A"),
                  choices = character(0)),
            selectInput("select.celltype.a", label = ("Cell type A"),
                  choices = character(0)),
            selectInput("select.target.a", label = ("Target A"),
                  choices = character(0)),
            selectInput("select.genotype.a", label = ("Genotype A"),
                  choices = character(0))
            ),
          column(6,
            selectInput("select.study.b", label = ("Study B"),
                   choices = character(0)),
            selectInput("select.celltype.b", label = ("Cell type B"),
                   choices = character(0)),
            selectInput("select.target.b", label = ("Target B"),
                   choices = character(0)),
            selectInput("select.genotype.b", label = ("Genotype B"),
                   choices = character(0))
            )
          ),
          actionButton("analyze", "Analyze"),
          actionButton("export.tables", "Export Tables"),
          plotOutput("pca.analysis")
        ),
      mainPanel(
        plotlyOutput("degenes"),
        DTOutput('detable'), style = "font-size:80%"
      )
    )
  ),
  server = function(input, output, session) {
    #max upload size
    options(shiny.maxRequestSize = 500 * 1024^2)
    shinyjs::disable("analyze")
    shinyjs::disable("export.tables")
    
    volumes <- getVolumes()()
    
    #RDS File Chooser
    shinyFileChoose(input = input, id = 'rds', session=session, 
                    roots=volumes, filetype = list(rds='rds'))
    
    #load rds file if requested
    se <- reactive (
      {
        if (!is.integer(input$rds))
        {
          rds.file <- parseFilePaths(roots = volumes, selection = input$rds)
          se <- readRDS(rds.file$datapath)
          return(list("se"=se, "sl"=colData(se)))
        }
      })
    
    #enable analyze when sample selection changes
    observeEvent(input$select.study.a, {
      req(input$select.study.a, input$select.study.a!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.a, ]
      updateSelectInput(session, "select.celltype.a", choices = unique(sl$cell_type))
      updateSelectInput(session, "select.target.a", choices = unique(sl$target))
      updateSelectInput(session, "select.genotype.a", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.study.b, {
      req(input$select.study.b, input$select.study.b!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.b, ]
      updateSelectInput(session, "select.celltype.b", choices = unique(sl$cell_type))
      updateSelectInput(session, "select.target.b", choices = unique(sl$target))
      updateSelectInput(session, "select.genotype.b", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.celltype.a, {
      req(input$select.celltype.a, input$select.celltype.a!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.a & se()$sl$cell_type == input$select.celltype.a, ]
      updateSelectInput(session, "select.target.a", choices = unique(sl$target))
      updateSelectInput(session, "select.genotype.a", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.celltype.b, {
      req(input$select.celltype.b, input$select.celltype.b!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.b & se()$sl$cell_type == input$select.celltype.b, ]
      updateSelectInput(session, "select.target.b", choices = unique(sl$target))
      updateSelectInput(session, "select.genotype.b", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.target.a, {
      req(input$select.target.a, input$select.target.a!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.a & 
                      se()$sl$cell_type == input$select.celltype.a &
                      se()$sl$target == input$select.target.a, ]
      updateSelectInput(session, "select.genotype.a", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.target.b, {
      req(input$select.target.b, input$select.target.b!="NA")
      sl <- se()$sl[se()$sl$study == input$select.study.b & 
                      se()$sl$cell_type == input$select.celltype.b &
                      se()$sl$target == input$select.target.b, ]
      updateSelectInput(session, "select.genotype.b", choices = unique(sl$genotype))
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.genotype.a, {
      req(input$select.genotype.a, input$select.genotype.a!="NA")
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })
    observeEvent(input$select.genotype.b, {
      req(input$select.genotype.b, input$select.genotype.b!="NA")
      if (samples()$a != samples()$b) {shinyjs::enable("analyze")} else {shinyjs::disable("analyze")}
    })

    #currently selected samples
    samples <- reactive (
      {
        a <- paste(input$select.study.a,input$select.celltype.a,input$select.target.a,input$select.genotype.a, sep="")
        b <- paste(input$select.study.b,input$select.celltype.b,input$select.target.b,input$select.genotype.b, sep="")
        return (list("a"=a, "b"=b))
      }
    )

    #observe if RDS file was selected - then load file
    observe(
      {
        l <- se()
        sl <- l$sl
        study <- as.list(unique(sl$study))
        names(study) <- study
        updateSelectInput (session, "select.study.a", choices = study)
        updateSelectInput (session, "select.study.b", choices = study)

        celltype <- as.list(unique(sl$cell_type))
        names(celltype) <- celltype
        updateSelectInput (session, "select.celltype.a", choices = celltype)
        updateSelectInput (session, "select.celltype.b", choices = celltype)
        
        target <- as.list(unique(sl$target))
        names(target) <- target
        updateSelectInput (session, "select.target.a", choices = target)
        updateSelectInput (session, "select.target.b", choices = target)
      }
    )

    #analysis upon click "Analyze" Button
    observeEvent(input$analyze, {
      #only perform analysis when data are loaded and proper selection took place
      req(input$select.study.a, input$select.study.a!="NA")
      req(input$select.study.b, input$select.study.b!="NA")
      req(input$select.celltype.a, input$select.celltype.a!="NA")
      req(input$select.celltype.b, input$select.celltype.b!="NA")
      req(input$select.target.a, input$select.target.a!="NA")
      req(input$select.target.b, input$select.target.b!="NA")
      req(input$select.genotype.a, input$select.genotype.a!="NA")
      req(input$select.genotype.b, input$select.genotype.b!="NA")

      #only submit if different samples are selected
      if (samples()$a != samples()$b)
      {
        shinyjs::disable("analyze")

        plots <- deseq.analysis ()

        output$pca.analysis <- renderPlot(plots$pca)
        output$degenes <- renderPlotly(plots$degenes)
        output$detable = renderDT(plots$detable)

      }
    })

    #export analysis tables upon clickin export tables
    observeEvent(input$export.tables, {
      plots <- deseq.analysis()
      data <- plots$detable
      contrast <- plots$contrast

      rds.file <- parseFilePaths(roots = volumes, selection = input$rds)
      export.path <- dirname(rds.file$datapath)
      export.file <- paste(export.path,"/",contrast,".lfcshrink.table.txt", sep="")
      write.table (data, file = export.file, sep="\t", col.names = NA)
    })
    
    #Deseq analysis, summarized experiment, study, target, comparison A vs B
    deseq.analysis <- reactive (
    {
      withProgress(message = 'prepare data', value = 0, {
      #assign input values to variables
      se <- se()$se
      study.a <- input$select.study.a
      celltype.a <- input$select.celltype.a
      target.a <- input$select.target.a
      gt.a <- input$select.genotype.a
      study.b <- input$select.study.b
      celltype.b <- input$select.celltype.b
      target.b <- input$select.target.b
      gt.b <- input$select.genotype.b

      contrast.a <- paste(study.a, celltype.a, target.a, gt.a, sep=".")
      contrast.b <- paste(study.b, celltype.b, target.b, gt.b, sep=".")
      
      se.a <- se[,se$study == study.a & se$cell_type == celltype.a & se$target == target.a & se$genotype == gt.a]
      se.b <- se[,se$study == study.b & se$cell_type == celltype.b & se$target == target.b & se$genotype == gt.b]
      se <- cbind(se.a, se.b)
      dt <- as.matrix(assays(se)$counts)
      sample.list <- as.data.frame(colData(se))
      sample.list$contrast <- paste(sample.list$study, sample.list$cell_type, sample.list$target, sample.list$genotype, sep=".")

      setProgress(message = "Deseq2 analysis", value = 0.1)
      
      dds <- DESeqDataSetFromMatrix(countData = dt, colData = sample.list, design = ~ contrast)
      dds <- dds[ rowMeans(counts(dds)) > 10, ]
      des <- DESeq(dds)

      setProgress(message = "vsd analysis", value = 0.6)
      
      vst_success <- TRUE
      
      # Attempt to apply vst, and if it fails, set the flag to FALSE
      tryCatch({
        vsd <- vst(dds, blind = FALSE)
      }, error = function(e) {
        warning("vst failed. Applying varianceStabilizingTransformation instead.")
        vst_success <<- FALSE
      })
      
      # If vst was unsuccessful, apply varianceStabilizingTransformation
      if (!vst_success) {
        vsd <- varianceStabilizingTransformation(dds)
      }
      
      mat <- assay(vsd)
      assay(vsd) <- mat

      p <- plotPCA(vsd, intgroup = "contrast",  returnData = TRUE)
      percentVar <- round(100 * attr(p, "percentVar"))

      pca <- ggplot(p, aes(x = PC1, y = PC2, color= contrast, label=name)) +
        geom_point(size =3) +
        geom_text(hjust = -0.2, nudge_y = 0.1, size=4) +
        scale_color_manual(values=c("darkorange", "grey30")) +
        theme(legend.position = "bottom") +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance"))

      setProgress(message = "calculate contrasts", value = 0.8)

      degenes <- results(des, contrast = c("contrast", contrast.b, contrast.a))
      degenes.sig <- subset(degenes, padj < 0.05)
      des.shrink <- as.data.frame(lfcShrink(type = "normal", contrast = c("contrast", contrast.b, contrast.a), des))
      up <- des.shrink[rownames(des.shrink) %in% rownames(degenes.sig) & des.shrink$log2FoldChange>0,]
      down <- des.shrink[rownames(des.shrink) %in% rownames(degenes.sig) & des.shrink$log2FoldChange< 0,]

      up <- rownames(up)
      down <- rownames(down)

      up.shrink <- des.shrink[rownames(des.shrink) %in% up,]
      up.shrink$geneID <- rownames(up.shrink)
      down.shrink <- des.shrink[rownames(des.shrink) %in% down,]
      down.shrink$geneID <- rownames(down.shrink)

      setProgress(message = "prepare plots", value = 0.9)

      degenes <- ggplot(data = des.shrink, aes(baseMean,log2FoldChange,label=rownames(des.shrink))) +
        scale_x_log10() + geom_point(color="gray") +
        geom_point (data = up.shrink, aes(baseMean,log2FoldChange,label=geneID), color = "coral1") +
        geom_point (data = down.shrink, aes(baseMean,log2FoldChange,label=geneID), color = "steelblue3") +
        labs (x="average normalized Expression",
              y=paste("log2fc ",contrast.b," vs. ",contrast.a, sep="")) +
        theme(axis.title.y = element_text(size = rel(0.8)))

      detable = as.data.frame(des.shrink)
      
      #enable table export of result files
      shinyjs::enable("export.tables")

      return (list("pca"=pca, "degenes"=degenes, "detable"=detable, 
                   "contrast"=paste(contrast.b,"_vs_",contrast.a,sep = "")))
      })
    })

    #Render plot according to table selection
    observeEvent(input$detable_rows_selected, ignoreNULL = FALSE, {
      if (!is.integer(input$rds))
      {
        plots <- deseq.analysis()
        data <- plots$detable
        selplot <- plots$degenes
        if (!is.null(input$detable_rows_selected))
          {
          selplot <- selplot +
            geom_point (data = data[input$detable_rows_selected,],
                      aes(baseMean,log2FoldChange,
                          label=rownames(data[input$detable_rows_selected,])),
                      color = "darkgreen", size = 2) +
            geom_text (data = data[input$detable_rows_selected,],
                      aes(baseMean,log2FoldChange,fontface = "bold",
                          label=rownames(data[input$detable_rows_selected,])),
                     color = "black", hjust = 0, nudge_x = 0.05)
          }
        selplotly <- ggplotly (selplot) %>% style(textposition = "right")
        output$degenes <- renderPlotly(selplotly)
      }
      else {return()}
    })
  }
)
