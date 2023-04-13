library(shiny)
library(shinyhelper)
library(bslib)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
if(!exists("sc1conf")) sc1conf = readRDS("sc1conf.rds")
if(!exists("sc1def")) sc1def  = readRDS("sc1def.rds")
if(!exists("sc1gene")) sc1gene = readRDS("sc1gene.rds")
if(!exists("sc1meta")) sc1meta = readRDS("sc1meta.rds")

### UI code
shinyUI(
fluidPage(style="margin:0;padding:0;",
tags$head(includeCSS("www/styles.css")),

theme = bslib::bs_theme(bootswatch = "flatly"),
navbarPage(
"IFNB",
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("sc1_civge_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_civge_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_civge_togL == true",
              selectInput("sc1_civge_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_civge_sub1.ui"),
              actionButton("sc1_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_civge_tog0 == true",
              sliderInput("sc1_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_civge_inp1", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civge_tog1 == true",
                  radioButtons("sc1_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc1_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("sc1_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("sc1_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.sc1_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("sc1_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("sc1_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civge_tog2 == true",
                  radioButtons("sc1_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc1_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("sc1_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge

,
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("sc1_civci_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_civci_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_civci_togL == true",
              selectInput("sc1_civci_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_civci_sub1.ui"),
              actionButton("sc1_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_civci_tog0 == true",
              sliderInput("sc1_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_civci_inp1", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civci_tog1 == true",
                  radioButtons("sc1_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("sc1_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("sc1_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_civci_inp2", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civci_tog2 == true",
                  radioButtons("sc1_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("sc1_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci

,
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("sc1_gevge_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_gevge_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc1_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_gevge_togL == true",
              selectInput("sc1_gevge_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_gevge_sub1.ui"),
              actionButton("sc1_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_gevge_tog0 == true",
              sliderInput("sc1_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_gevge_tog1 == true",
                  radioButtons("sc1_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc1_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("sc1_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("sc1_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_gevge_tog2 == true",
                  radioButtons("sc1_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc1_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("sc1_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("sc1_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("sc1_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("sc1_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge

,
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("sc1_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("sc1_gem_drX", "X-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[1]
            ),
            selectInput("sc1_gem_drY", "Y-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[2]
            ),
            checkboxInput("sc1_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_gem_togL == true",
              selectInput("sc1_gem_sub1", "Cell info to subset:",
                          choices = sc1conf[grp == TRUE]$UI,
                          selected = sc1def$grp1
              ),
              uiOutput("sc1_gem_sub1.ui"),
              actionButton("sc1_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_gem_tog0 == true",
              sliderInput("sc1_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("sc1_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_gem_txt", "Show axis text", value = FALSE),
              radioButtons("sc1_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc1_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("sc1_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("sc1_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc1_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc1_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("sc1_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("sc1_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("sc1_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem

,
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("sc1_gec_drX", "X-axis:",
                               choices = sc1conf[dimred == TRUE]$UI,
                               selected = sc1def$dimred[1]
                   ),
                   selectInput("sc1_gec_drY", "Y-axis:",
                               choices = sc1conf[dimred == TRUE]$UI,
                               selected = sc1def$dimred[2]
                   ),
                   selectInput("sc1_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("sc1_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   checkboxInput("sc1_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.sc1_gec_togL == true",
                    selectInput("sc1_gec_sub1", "Cell info to subset:", choices = sc1conf[grp == TRUE]$UI, selected = sc1def$grp1),
                     uiOutput("sc1_gec_sub1.ui"),
                     actionButton("sc1_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("sc1_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("sc1_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.sc1_gec_tog0 == true",
                     radioButtons("sc1_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("sc1_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("sc1_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("sc1_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("sc1_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("sc1_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("sc1_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("sc1_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("sc1_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("sc1_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("sc1_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("sc1_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("sc1_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec

,
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("sc1_vio_inp1", "Cell info (X-axis):",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("sc1_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("sc1_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("sc1_vio_pts", "Show data points", value = FALSE),
            checkboxInput("sc1_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_vio_togL == true",
              selectInput("sc1_vio_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_vio_sub1.ui"),
              actionButton("sc1_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_vio_tog == true",
              sliderInput("sc1_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.sc1_vio_typ == 'lineplot'",
              sliderInput("sc1_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc1_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc1_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("sc1_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc1_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio

,
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("sc1_pro_inp1", "Cell info to plot (X-axis):",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("sc1_pro_inp2", "Cell info to group / colour by:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("sc1_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("sc1_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("sc1_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_pro_togL == true",
              selectInput("sc1_pro_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_pro_sub1.ui"),
              actionButton("sc1_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_pro_tog == true",
              radioButtons("sc1_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc1_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc1_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("sc1_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc1_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro

,
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("sc1_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("sc1_hea_grp", "Group by:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("sc1_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("sc1_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("sc1_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("sc1_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("sc1_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_hea_togL == true",
              selectInput("sc1_hea_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_hea_sub1.ui"),
              actionButton("sc1_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_hea_tog == true",
              radioButtons("sc1_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc1_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("sc1_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("sc1_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("sc1_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("sc1_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("sc1_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("sc1_hea_oupTxt")),
          uiOutput("sc1_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea


,
# about ----
tabPanel(
  "About",
  fluidRow(
    class = "container page",
    column(
      12,
      includeMarkdown("about.md")
    )
  )
))
)
)
