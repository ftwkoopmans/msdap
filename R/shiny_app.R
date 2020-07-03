
#' open the MS-DAP graphical user interface (a R shiny app)
#'
#' @param ... optionally, pass some additional parameters to shiny::shinyApp()
#' @importFrom shiny shinyApp
#' @export
gui = function() {
  shiny::shinyApp(ui = msdap_shiny_ui(), server = msdap_shiny_server)
}





#' toggle shiny input elements
#'
#' https://stackoverflow.com/a/48852284
#'
#' @param input_list List of inputs, eg; input_list = reactiveValuesToList(input)
#' @param enable_components boolean flag; enable or disable?
#' @param exclude_names element names to exclude
#' @importFrom shinyjs enable disable
toggle_inputs <- function(input_list, enable_components = T, exclude_names = c("log", "btnScreenPrevious", "btnScreenNext")) {
  # print(sort(names(input_list)))
  for(x in setdiff(names(input_list), exclude_names)) {
    if(enable_components) {
      shinyjs::enable(x)
    } else {
      shinyjs::disable(x)
    }
  }
}


#' shiny function to show/hide elements
#'
#' @param idshow list of all ids (to hide)
#' @param idall single id to show after hiding all
#' @importFrom shinyjs show hide
shiny_show_screen = function(idshow, idall) {
  sapply(idall, shinyjs::hide)
  shinyjs::show(idshow)
}


validate_input_peptides = function(datastore_reactive) {
  return(length(datastore_reactive$dataset) > 0 && "peptides" %in% names(datastore_reactive$dataset) &&
           tryCatch({msdap::check_valid_tibble_peptides(datastore_reactive$dataset$peptides)}, error=function(err) {print(err); return(F)}) )
}

validate_input_samples = function(datastore_reactive) {
  return(length(datastore_reactive$dataset) > 0 &&
           "samples" %in% names(datastore_reactive$dataset) &&
           tryCatch({msdap::check_valid_tibble_samples(datastore_reactive$dataset$samples)}, error=function(err) {print(err); return(F)}) )
}

modalResetConfirm = function(id) {
  modalDialog(
    title = "confirm new dataset",
    "Changing these settings will unload current dataset. Continue?",
    easyClose = FALSE,
    footer = tagList(
      actionButton(paste0(id,"cancel"), "Cancel"),
      actionButton(paste0(id,"ok"), "OK")
    )
  )
}


#' MS-DAP shiny UI
#'
#' @importFrom shiny reactiveValuesToList a actionButton br checkboxGroupInput column conditionalPanel div fluidPage fluidRow h2 h3 h4 hr HTML htmlOutput numericInput observeEvent observe p radioButtons reactiveValues renderTable renderText selectInput span strong tableOutput tagAppendAttributes textOutput updateSelectInput showModal modalDialog shinyApp
#' @importFrom shinyFiles parseDirPath parseFilePaths shinyDirButton shinyDirChoose shinyFileChoose shinyFilesButton
#' @importFrom shinyjs click disable enable hide onclick show useShinyjs
#'
#' @export
msdap_shiny_ui = function() {

  # CSS left-handside text overflow cutoff  @  https://stackoverflow.com/a/9793669
  css_overflow_text = "white-space: nowrap; overflow: hidden; text-overflow: ellipsis; direction: rtl; text-align: left;"
  # css_overflow_text_multiline = "white-space: pre-wrap; overflow: hidden; text-overflow: ellipsis; direction: rtl; text-align: left; background-color: #DDDDDD;"
  # css_bootstrap_form_control = "border-color: #66afe9; outline: 0; -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.075),0 0 8px rgba(102,175,233,.6); box-shadow: inset 0 1px 1px rgba(0,0,0,.075),0 0 8px rgba(102,175,233,.6);"
  css_bootstrap_form_control = "border: 1px solid #ccc; border-radius: 4px; -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.075); box-shadow: inset 0 1px 1px rgba(0,0,0,.075);"

  help_text = parameter_documentation(func="analysis_quickstart", package="msdap")

  shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        body {
          max-width: 1024px;
        }
        .shiny-output-error-validation {
          color: orange;
        }
        .form-inline > input {
          width: 80px;
          max-width: 80px;
        }
      "))
    ),
    title = paste("MS-DAP", packageVersion("msdap")),
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shiny::column(6, shiny::h2("MS-DAP")),
      shiny::column(6, shiny::div(shiny::br(), shiny::span(paste("version:", packageVersion("msdap"))), shiny::HTML(" &nbsp; "), shiny::span(" additional resources:"), shiny::br(), shiny::span(shiny::a("https://github.com/ftwkoopmans/msdap")), class="text-right"))
      # shiny::column(6, shiny::div(shiny::br(), paste("version:", packageVersion("msdap"), " additional resources:"), shiny::br(), shiny::a("https://github.com/ftwkoopmans/msdap"), class="text-right")),
      # class="bg-info"
    ),

    shiny::fluidRow(
      id="screenData",
      shiny::h3("Step 1: import dataset"),
      shiny::hr(),
      shiny::selectInput("software", "upstream software:", choices = c("DIA-NN" = "diann", FragPipe = "fragpipe", MaxQuant = "maxquant", MetaMorpheus = "metamorpheus", OpenMS = "openms", OpenSWATH = "openswath", Peaks = "peaks", ProteomeDiscoverer = "pd", "Skyline - DDA" = "skyline_dda", "Skyline - DIA" = "skyline_dia", Spectronaut = "spectronaut", EncyclopeDIA = "encyclopedia"), multiple = F),

      shiny::p(shiny::strong("input data")),
      # dynamic UI: some input data requires a filepicker for directories while others just need the user to select a file
      shiny::conditionalPanel(condition = "input.software == 'maxquant' || input.software == 'metamorpheus' || input.software == 'encyclopedia'", # this is a javascript expression (...)
                              shiny::fluidRow(shiny::column(3, shinyFiles::shinyDirButton('dirDataset', label='browse files...', title='Please select a directory')),
                                              shiny::column(8, shiny::tagAppendAttributes(shiny::textOutput("dirDatasetResult"), style = paste(css_overflow_text, "background-color: #DDDDDD; color: #555;"), class="form-control"))
                              )),
      shiny::conditionalPanel(condition = "input.software != 'maxquant' && input.software != 'metamorpheus' && input.software != 'encyclopedia'", # this is a javascript expression (...)
                              shiny::fluidRow(shiny::column(3, shinyFiles::shinyFilesButton('fileDataset', label='browse files...', title='Please select a file', multiple=FALSE, viewtype = "detail")),
                                              shiny::column(8, shiny::tagAppendAttributes(shiny::textOutput("fileDatasetResult"), style = paste(css_overflow_text, "background-color: #DDDDDD; color: #555;"), class="form-control"))
                              )),

      shiny::br(),
      shiny::p(shiny::strong("fasta")),
      shiny::fluidRow(
        shiny::column(3, shinyFiles::shinyFilesButton('fileFasta', label='browse files...', title='hold the shift or ctrl key to select multiple', multiple=TRUE, viewtype = "detail")),
        shiny::column(8, shiny::tagAppendAttributes(shiny::tableOutput("fileFastaResult"), class="panel panel-default", style=paste(css_bootstrap_form_control, "background-color: #DDDDDD; color: #555;")))
      ),

      shiny::br(),
      shiny::p(shiny::strong("output directory")),
      shiny::fluidRow(
        shiny::column(3, shinyFiles::shinyDirButton('dirOutput', label='select dir...', title='Please select a directory')),
        shiny::column(8, shiny::tagAppendAttributes(shiny::textOutput("dirOutputResult"), style = paste(css_overflow_text, "background-color: #DDDDDD; color: #555;"), class="form-control"))
      )
    ),


    shiny::fluidRow(
      id="screenSamples",
      style = "display: none;",
      shiny::h3("Step 2: describe sample metadata"),
      shiny::hr(),
      shiny::fluidRow(
        shiny::column(width=3, shiny::actionButton("fileWriteSampleTemplate", label="1) create template file")),
        shiny::column(width=8, shiny::tags$p("open the template file in Excel and complete the sample metadata table")),
      ),
      shiny::br(),
      shiny::fluidRow(
        shiny::column(width=3, shinyFiles::shinyFilesButton('fileSamples', label='2) load sample metadata table', title='only files with .xlsx, .csv or .tsv extension are shown', multiple=FALSE, viewtype = "detail")),
        shiny::column(width=8, shiny::tagAppendAttributes(shiny::textOutput("fileSamplesResult"), style = paste(css_overflow_text, "background-color: #DDDDDD; color: #555;"), class="form-control"))
      ),
      shiny::hr(),
      shiny::br(),
      shiny::fluidRow(
        shiny::column(12, shiny::tagAppendAttributes(shiny::tableOutput('tblSamples'), style = "min-height:200px; max-height:600px; overflow-y: scroll;"))
      )
    ),

    shiny::fluidRow(
      id="screenContrasts",
      style = "display: none;",
      shiny::h3("Step 3: define contrasts"),
      shiny::hr(),
      module_contrast_ui()
    ),

    shiny::fluidRow(
      id="screenSettings",
      style = "display: none;",
      shiny::h3("Step 4: settings"),
      shiny::hr(),
      shiny::h4("Peptide filtering rules applied within each sample group; in how many samples should a peptide be confidently identified ?"),
      shiny::fluidRow(
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("mindetect", "identified sample count", value = 3, min = 0, max = NA, step = 1), class="form-inline"), create_help_icon(help_text["filter_min_detect"])) ),
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("fracdetect", "% of samples", value = 25, min = 0, max = 100, step = 1), class="form-inline"), create_help_icon(help_text["filter_fraction_detect"])) )
      ),
      shiny::fluidRow(
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("minquant", "quantified (MBR or identified)", value = 3, min = 0, max = NA, step = 1), class="form-inline"), create_help_icon(help_text["filter_min_quant"])) ),
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("fracquant", "% of samples", value = 75, min = 0, max = 100, step = 1), class="form-inline"), create_help_icon(help_text["filter_fraction_quant"])) )
      ),

      shiny::hr(),
      shiny::fluidRow(
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("filter_topn_peptides", "topN peptides per protein:", value = 0, min = 0, step = 1), class="form-inline"), create_help_icon(help_text["filter_topn_peptides"])) ),
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("filter_min_peptide_per_prot", "Peptides per protein", value = 1, min = 0, step = 1), class="form-inline"), create_help_icon(help_text["filter_min_peptide_per_prot"])) )
      ),

      shiny::hr(),
      shiny::fluidRow(
        shiny::column(12, insert_child_before(shiny::tagAppendAttributes(shiny::radioButtons("filter_by_contrast", "Apply filtering rules:", choices = c("to entire dataset"="all_group", "separately within each contrast"="by_contrast"), inline = T), class="form-inline"), create_help_icon(help_text["filter_by_contrast"])) )
      ),

      shiny::hr(),
      shiny::fluidRow(
        shiny::column(12, insert_child_before(shiny::tagAppendAttributes(shiny::radioButtons("norm_algorithm", "Normalization algorithm:", choices = c("none","median","vsn","loess","rlr","vwmb"), selected = "vwmb", inline = T), class="form-inline"), create_help_icon(help_text["norm_algorithm"])) )
      ),

      shiny::hr(),
      shiny::fluidRow(
        shiny::column(12, insert_child_before(shiny::tagAppendAttributes(shiny::checkboxGroupInput("dea_algorithm", "DEA algorithm", choices = c("ebayes","msqrob","msempire"), selected = "ebayes", inline = T), class="form-inline"), create_help_icon(help_text["dea_algorithm"])) )
      ),
      shiny::fluidRow(
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("dea_qvalue_threshold", "FDR threshold (%)", value = 1, min = 1, step = 1), class="form-inline"), create_help_icon(help_text["dea_qvalue_threshold"])) ),
        shiny::column(6, insert_child_before(shiny::tagAppendAttributes(shiny::numericInput("dea_log2foldchange_threshold", "log2-foldchange threshold", value = -1, min = -1), class="form-inline"), create_help_icon(help_text["dea_log2foldchange_threshold"])) )
      ),

      shiny::hr(),
      shiny::fluidRow(
        shiny::column(12, insert_child_before(shiny::tagAppendAttributes(shiny::checkboxGroupInput("output_files", "Output", choices = c("PDF report"="report", "data tables"="abundance_tables"), selected = "report", inline = T), class="form-inline"), create_help_icon(help_text["output_qc_report"])) )
       # shiny::column(12, insert_child_before(create_help_icon(paste(help_text["output_qc_report"], "\n\n", help_text["output_peptide_plots"], "\n\n", help_text["output_abundance_tables"])), shiny::tagAppendAttributes(shiny::checkboxGroupInput("output_files", "Output", choices = c("report","peptide_plots","abundance_tables"), selected = "report", inline = T), class="form-inline")) )
      ),

      shiny::hr(),
      shiny::br(),
      shiny::fluidRow(
        shiny::column(12, shiny::tagAppendAttributes(shiny::actionButton("btnAnalysisQuickstart", label="start analysis !", title="Running the pipeline may take a while, a popup message will notify you when ready !"), class="btn btn-lg btn-danger"), class="text-center" )
      )
    ),

    shiny::hr(),
    shiny::fluidRow(
      shiny::column(2, shiny::actionButton("btnScreenPrevious", "Previous", class="btn-info")),
      shiny::column(8, shiny::div()),
      shiny::column(2, shiny::actionButton("btnScreenNext", "Next", class="btn-primary"))
    ),
    shiny::hr(),
    shiny::p(shiny::strong('live logging your analysis (latest actions are on top)')),
    shiny::tagAppendAttributes(shiny::htmlOutput("log"), style="white-space:pre-wrap; background-color: #DDDDDD; height:400px; overflow-y: scroll;", class="form-control text-monospace")
  )
}



#' MS-DAP shiny server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @importFrom shiny renderText reactive
#' @importFrom testthat capture_output_lines
#' @export
msdap_shiny_server <- function(input, output, session) {

  resetPeptides = function() {
    # print("reset called")
    datastore_reactive$dataset = list()
    datastore_reactive$fileDataset=NULL
    datastore_reactive$dirDataset=NULL
    datastore_reactive$fileFasta=NULL
    datastore_reactive$output_dir = NULL
    output$fileDatasetResult = shiny::renderText("No file selected")
    output$dirDatasetResult = shiny::renderText("No file selected")
    output$dirOutputResult = shiny::renderText("No directory selected")
    output$fileFastaResult = shiny::renderTable(data.table(file="No file selected"), rownames = F, colnames = F)
    resetSamples()
    resetContrasts()
    # limit screen progress
    datastore_reactive$screens_validated_max = min(datastore_reactive$screens_validated_max, 1)

    shinyjs::disable("fileFasta")
    shinyjs::disable("btnScreenPrevious")
    shinyjs::disable("btnScreenNext")
  }

  resetSamples = function() {
    datastore_reactive$fileSamples=NULL
    datastore_reactive$sample_groups = NULL
    if(length(datastore_reactive$dataset) > 0) {
      datastore_reactive$dataset$samples = NULL
    }
    output$fileSamplesResult = shiny::renderText("No file selected")
    output$tblSamples = shiny::renderTable(data.frame(status="data has not been loaded yet..."))
    resetContrasts()
    # limit screen progress
    datastore_reactive$screens_validated_max = min(datastore_reactive$screens_validated_max, 2)
    shinyjs::disable("btnScreenNext")
  }

  resetContrasts = function() {
    # remove dynamic content
    shiny::removeUI(selector = "[id^='contrast_dyn_']", multiple = T)
    remove_shiny_inputs("^contrast_dyn_", input)
    # limit screen progress
    datastore_reactive$screens_validated_max = min(datastore_reactive$screens_validated_max, 3)
    shinyjs::disable("btnScreenNext")
    #
    contrast_reactive_values <<- shiny::reactiveValues() # use <<- to update in PARENT ENVIRONMENT, per https://adv-r.hadley.nz/environments.html#env-basics
  }


  update_shiny_file_pickers = function() {
    shinyFiles::shinyFileChoose(input, 'fileDataset', root = isolate(datastore_reactive$filepicker_roots), defaultRoot = "recent")
    shinyFiles::shinyFileChoose(input, 'fileFasta', root = isolate(datastore_reactive$filepicker_roots), defaultRoot = "recent", filetypes=c('fasta'))
    shinyFiles::shinyFileChoose(input, 'fileSamples', root = isolate(datastore_reactive$filepicker_roots), defaultRoot = "recent", filetypes=c('xlsx', 'tsv', 'csv'))
    shinyFiles::shinyDirChoose(input, 'dirDataset', root = isolate(datastore_reactive$filepicker_roots), defaultRoot = "recent")
    shinyFiles::shinyDirChoose(input, 'dirOutput', root = isolate(datastore_reactive$filepicker_roots), defaultRoot = "recent")
  }


  init_roots = function() {
    x = shinyFiles::getVolumes()()

    if(dir.exists("/exampledata")) {
      x["exampledata"] = "/exampledata"
    }

    if(dir.exists("/data")) {
      x = c(x, data = "/data", recent="/data")
    } else {
      x["recent"] = x[1]
    }

    return(x)
  }


  ######### reactive values
  datastore_reactive <- shiny::reactiveValues(screens_current_index = 1,
                                              screens_validated_max = 1,
                                              dataset = list(),
                                              filepicker_roots = init_roots(),
                                              log=NULL, history=NULL,
                                              fileDataset=NULL, dirDataset=NULL, fileFasta=NULL, fileSamples=NULL, output_dir = NULL,
                                              sample_groups = NULL) #reactiveVal(c("test1","test2","test3")))

  # init file pickers
  update_shiny_file_pickers()


  ############ custom code for jumping between screens, because shinyglide has a javascript bug / incompatability with the sortablejs code for drag&drop
  screens_all = c("screenData", "screenSamples", "screenContrasts", "screenSettings")

  # datastore_reactive$sample_groups = letters[1:3]
  # shiny::observe({
  #   Sys.sleep(2)
  #   shinyjs::click("btnAddContrast")
  # })
  # shiny_show_screen(screens_all[3], screens_all) # debug; while testing, jump directly to some screen


  shiny::observeEvent(input$btnScreenNext, {
    if(datastore_reactive$screens_current_index < length(screens_all)) {
      datastore_reactive$screens_current_index = datastore_reactive$screens_current_index + 1
      shiny_show_screen(screens_all[datastore_reactive$screens_current_index], screens_all)

      # just jumped to last screen, disable next button
      if(datastore_reactive$screens_current_index == length(screens_all)) {
        shinyjs::disable("btnScreenNext")
      }

      # we haven't validated input beyond current, so disable next
      if(datastore_reactive$screens_validated_max <= datastore_reactive$screens_current_index) {
        shinyjs::disable("btnScreenNext")
      }

      # just jumped to second screen, enable back button
      if(datastore_reactive$screens_current_index == 2) {
        shinyjs::enable("btnScreenPrevious")
      }
    }
  }, ignoreInit = TRUE)


  shiny::observeEvent(input$btnScreenPrevious, {
    # cat(sprintf("screens_validated_max:%s screens_current_index%s\n", datastore_reactive$screens_validated_max, datastore_reactive$screens_current_index))
    if(datastore_reactive$screens_current_index > 1) {
      datastore_reactive$screens_current_index = datastore_reactive$screens_current_index - 1
      shiny_show_screen(screens_all[datastore_reactive$screens_current_index], screens_all)

      # just jumped to first screen, disable previous button
      if(datastore_reactive$screens_current_index == 1) {
        shinyjs::disable("btnScreenPrevious")
      }

      # going to previous screen, we can always jump forward again
      shinyjs::enable("btnScreenNext")
    }
  }, ignoreInit = TRUE)


  ######### LOG
  output$log <- shiny::renderText({
    datastore_reactive$log
  }, sep = "")

  ### register a logger for the msdap pipeline so we can redirect log to the GUI
  # based on msdap::logger.default()
  log_ <<- list() # reset log on each run
  logger_msdap <<- function(x, type = "info") {
    if (!exists("log_")) {
      log_ <<- list()
    }
    log_[[length(log_)+1]] = c(x, type, Sys.time())
    log_ <<- log_

    # if(is.character(type) && type == "progress") { return() } # example; disable some log type
    clr = log_type_to_color(type)
    if(type == "gui") {
      datastore_reactive$log = c(paste0("<font color=\"", clr, "\">", x, "</font>\n"), datastore_reactive$log)
    } else {
      datastore_reactive$log = c(paste0("<font color=\"", clr, "\">[", format(Sys.time(), "%H:%M:%S"), "]  ", x, "</font>\n"), datastore_reactive$log)
    }
  }


  ################################################################################################
  ###### input software -->> update UI
  # input_software = c("MaxQuant"="maxquant", "MetaMorpheus"="metamorpheus", "FragPipe"="fragpipe", "Skyline - DDA"="skyline_dda", "ProteomeDiscoverer"="pd", "Peaks"="peaks", "OpenMS"="openms", "Spectronaut"="spectronaut", "DIA-NN"="diann", "Skyline - DIA"="skyline_dia", "OpenSWATH"="openswath")
  # input_software = input_software[order(names(input_software))]
  #
  # # populate UI, once
  # shiny::observeEvent(input$software, {
  #   shiny::updateSelectInput(session, "software", choices = input_software)
  # }, once = TRUE)

  # what to do after the user selects software? reset whatever data was previously loaded
  shiny::observeEvent(input$software, {
    if(length(datastore_reactive$dataset) > 0) {
      append_log("\nunloaded current dataset \n\n", type = "warning")
      # resetPeptides()
    } else {
      # datastore_reactive$fileDataset=NULL
      # datastore_reactive$dirDataset=NULL
      # output$fileDatasetResult = shiny::renderText("No file selected")
      # output$fileFastaResult = shiny::renderTable(data.table(file="No file selected"), rownames = F, colnames = F)
      # # output$fileFastaResult = shiny::renderTable(data.table(file=c("No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1No file selected 1","No file selected 2","No file selected 3")), rownames = F, colnames = F, width = "100%", spacing = "xs") # striped = F, hover = F, bordered = F,
    }

    resetPeptides()
  })


  shiny::observeEvent(input$fileWriteSampleTemplate, {
    tryCatch({
      msdap::write_template_for_sample_metadata(datastore_reactive$dataset, sprintf("%s/samples_%s.xlsx", datastore_reactive$output_dir, format(Sys.time(), "%Y-%m-%d_%H;%M;%S")))
      msdap::append_log(paste("sample metadata template file was created at:", datastore_reactive$output_dir), type = "progress")
    }, error = function(err) { print(err); warning(err) })

  }, ignoreInit = TRUE)



  shiny::observeEvent(input$fileDataset, {
    is_reload = length(datastore_reactive$dataset) > 0
    filename_input_file <- shinyFiles::parseFilePaths(roots = datastore_reactive$filepicker_roots, selection = input$fileDataset)

    if(!is.null(filename_input_file) && nrow(filename_input_file) == 1 && "datapath" %in% names(filename_input_file)) {
      filename_input_file$datapath = normalizePath(filename_input_file$datapath, winslash = "/")

      if(!is_reload || datastore_reactive$fileDataset != filename_input_file$datapath) {

        input_list <- shiny::reactiveValuesToList(input)
        tryCatch({
          # reset
          datastore_reactive$fileDataset = NULL
          # datastore_reactive$output_dir = NULL
          datastore_reactive$dataset = list()
          # store recent file dir
          datastore_reactive$filepicker_roots["recent"] = dirname(filename_input_file$datapath)
          update_shiny_file_pickers()

          # block UI and parse data
          toggle_inputs(input_list, enable_components = FALSE)
          if(length(datastore_reactive$log) > 0) {
            append_log("\n\n --------------------------------------------------------------- \n\n", type = "gui")
          }

          datastore_reactive$dataset = msdap::import_dataset__generic(filename_input_file$datapath, type = input$software)
          datastore_reactive$history = c(datastore_reactive$history, sprintf('dataset = msdap::import_dataset__generic("%s", type="%s")', filename_input_file$datapath, input$software))

          # ##################### LOAD MOCK DATA
          # cat("loading...")
          # f <- system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")
          # datastore_reactive$dataset = import_dataset_skyline(f, confidence_threshold = 0.01, return_decoys = F, acquisition_mode = "dia")
          # cat("done!\n")
          # ##################### LOAD MOCK DATA

          # no errors thus successfull, store values and update UI
          datastore_reactive$fileDataset = filename_input_file$datapath
          output$fileDatasetResult <- shiny::renderText(filename_input_file$datapath)

          # default output dir is the parent of selected directory
          if(length(datastore_reactive$output_dir) == 0) {
            datastore_reactive$output_dir = dirname(filename_input_file$datapath)
            output$dirOutputResult <- shiny::renderText(datastore_reactive$output_dir)
          }
        }, error = function(err) { print(err); warning(err) })

        toggle_inputs(input_list, enable_components = TRUE)

        if(length(datastore_reactive$output_dir) == 1 && validate_input_peptides(datastore_reactive)) {
          shinyjs::enable("fileFasta")
          shinyjs::enable("btnScreenNext")
          datastore_reactive$screens_validated_max = 2
        } else {
          shinyjs::disable("fileFasta")
          shinyjs::disable("btnScreenNext")
          datastore_reactive$screens_validated_max = 1
        }
      }
    }
  }, ignoreInit = TRUE)



  shiny::observeEvent(input$dirDataset, {
    is_reload = length(datastore_reactive$dataset) > 0
    input_dataset_path = shinyFiles::parseDirPath(roots = datastore_reactive$filepicker_roots, selection = input$dirDataset)

    if(!is.null(input_dataset_path) && length(input_dataset_path) == 1) {
      input_dataset_path = normalizePath(paste0(input_dataset_path, "/"), winslash = "/")

      if(!is_reload || datastore_reactive$dirDataset != input_dataset_path) {

        input_list <- shiny::reactiveValuesToList(input)
        tryCatch({
          # reset
          datastore_reactive$dirDataset = NULL
          # datastore_reactive$output_dir = NULL
          datastore_reactive$dataset = list()
          # store recent file dir
          datastore_reactive$filepicker_roots["recent"] = input_dataset_path
          update_shiny_file_pickers()

          # block UI and parse data
          toggle_inputs(input_list, enable_components = FALSE)
          if(length(datastore_reactive$log) > 0) {
            append_log("\n\n --------------------------------------------------------------- \n\n", type = "gui")
          }

          datastore_reactive$dataset = msdap::import_dataset__generic(input_dataset_path, type = input$software)
          datastore_reactive$history = c(datastore_reactive$history, sprintf('dataset = msdap::import_dataset__generic("%s", type="%s")', input_dataset_path, input$software))

          # no errors thus successfull, store values and update UI
          datastore_reactive$dirDataset = input_dataset_path
          output$dirDatasetResult <- shiny::renderText(input_dataset_path)

          # default output dir is the parent of selected directory
          if(length(datastore_reactive$output_dir) == 0) {
              datastore_reactive$output_dir = dirname(datastore_reactive$dirDataset)
              output$dirOutputResult <- shiny::renderText(datastore_reactive$output_dir)
          }
          # if(nchar(input_dataset_path) <= nchar(root_path) || substr(input_dataset_path, 1, nchar(root_path)) != root_path) { # test; root_path = "C:/data/"; input_dataset_path = "C:/";;;;;root_path = "C:/data/"; input_dataset_path = "C:/data/a/"
          #   datastore_reactive$output_dir = root_path
          # } else {
          #   datastore_reactive$output_dir = dirname(datastore_reactive$dirDataset)
          # }
          # msdap::append_log(paste("output_dir:", datastore_reactive$output_dir), type = "info")
        }, error = function(err) { print(err); warning(err) })

        toggle_inputs(input_list, enable_components = TRUE)

        if(length(datastore_reactive$output_dir) == 1 && validate_input_peptides(datastore_reactive)) {
          shinyjs::enable("fileFasta")
          shinyjs::enable("btnScreenNext")
          datastore_reactive$screens_validated_max = 2
        } else {
          shinyjs::disable("fileFasta")
          shinyjs::disable("btnScreenNext")
          datastore_reactive$screens_validated_max = 1
        }
      }
    }
  }, ignoreInit = TRUE)



  shiny::observeEvent(input$dirOutput, {
    output_dataset_path = shinyFiles::parseDirPath(roots = datastore_reactive$filepicker_roots, selection = input$dirOutput)
    if(!is.null(output_dataset_path) && length(output_dataset_path) == 1) {
      output_dataset_path = normalizePath(paste0(output_dataset_path, "/"), winslash = "/")

      datastore_reactive$output_dir = output_dataset_path
      output$dirOutputResult <- shiny::renderText(output_dataset_path)

      if(length(datastore_reactive$output_dir) == 1 && validate_input_peptides(datastore_reactive)) {
        shinyjs::enable("fileFasta")
        shinyjs::enable("btnScreenNext")
        datastore_reactive$screens_validated_max = 2
      } else {
        shinyjs::disable("fileFasta")
        shinyjs::disable("btnScreenNext")
        datastore_reactive$screens_validated_max = 1
      }
    }
  }, ignoreInit = TRUE)



  shiny::observeEvent(input$fileFasta, {
    fasta_input_file <- shinyFiles::parseFilePaths(roots = datastore_reactive$filepicker_roots, selection = input$fileFasta)

    if(!is.null(fasta_input_file) && nrow(fasta_input_file) > 0 && "datapath" %in% names(fasta_input_file)) {
      fasta_input_file$datapath = normalizePath(fasta_input_file$datapath, winslash = "/")

      input_list <- shiny::reactiveValuesToList(input)
      tryCatch({
        datastore_reactive$fileFasta = NULL

        # block UI and parse data
        toggle_inputs(input_list, enable_components = FALSE)

        datastore_reactive$dataset = msdap::import_fasta(datastore_reactive$dataset, files = fasta_input_file$datapath)
        # log the R command
        datastore_reactive$history = c(datastore_reactive$history, sprintf('dataset = msdap::import_fasta(dataset, files=%s)',
                                                                           r_variable_to_string(as.character(fasta_input_file$datapath))))

        # no errors thus successfull, store values and update UI
        datastore_reactive$fileFasta = fasta_input_file$datapath

        output$fileFastaResult = shiny::renderTable(data.table(file=basename(fasta_input_file$datapath)), rownames = F, colnames = F, width = "100%", spacing = "xs") # striped = F, hover = F, bordered = F,
      }, error = function(err) { print(err); warning(err) })

      toggle_inputs(input_list, enable_components = TRUE)
    }
  }, ignoreInit = TRUE)







  ####################### CONTRASTS SCREEN

  contrast_reactive_values = shiny::reactiveValues()

  # dynamically add/remove contrasts
  # use shinyjs to listen for event in order to be compatible with; observe({shinyjs::click("btnAddContrast")})
  shinyjs::onclick("btnAddContrast", {
    if(length(datastore_reactive$sample_groups) == 0) {
      append_log("cannot add a contrast, there are no groups in the dataset", type = "warning")
      return()
    }
    contr_id = sprintf('contrast_dyn_%04d', input$btnAddContrast)
    shiny::insertUI(
      selector = '#rowContrastValidation',
      where = "beforeBegin",
      ui = contrast_creator_ui(contr_id, item_list=datastore_reactive$sample_groups)
    )
    contrast_reactive_values[[contr_id]] = shiny::callModule(contrast_creator_server, id = contr_id)
    shiny::observeEvent(input[[paste0(contr_id, '-btnRemove')]], {
      shiny::removeUI(selector = sprintf("[id^='%s']", contr_id), multiple = T)
      remove_shiny_inputs(contr_id, input)
      contrast_reactive_values[[contr_id]] = NULL
    })
  })

  ### input validation
  shiny::observe({
    if(datastore_reactive$screens_validated_max >= 3 && datastore_reactive$screens_current_index == 3) {
      # from list of reactives to a list of contrast strings
      l = lapply(shiny::reactiveValuesToList(contrast_reactive_values), function(contr_mod_result) {
        lapply(contr_mod_result, function(x) x())
      })

      # remove empty contrasts
      l = l[lengths(lapply(l,unlist,use.names=F,recursive=F)) > 0]

      # only validate if there are any values in UI
      if(length(l) > 0) {
        # returns an array of error strings
        errors = contrast_definitions_validate(l)
        if(length(errors) > 0) {
          output$contrastErrors = shiny::renderText(paste(errors, collapse="\n"))
          shinyjs::disable("btnScreenNext")
          datastore_reactive$screens_validated_max = 3
        } else {
          output$contrastErrors = shiny::renderText("")
          shinyjs::enable("btnScreenNext")
          datastore_reactive$screens_validated_max = 4
        }
      } else {
        output$contrastErrors = shiny::renderText("")
        shinyjs::enable("btnScreenNext")
        datastore_reactive$screens_validated_max = 4
      }
    }
  }, label = "observeVALIDATEscreen")



  ########### SAMPLES

  shiny::observeEvent(input$fileSamples, {
    metadata_input_file <- shinyFiles::parseFilePaths(roots = datastore_reactive$filepicker_roots, selection = input$fileSamples)


    if(!is.null(metadata_input_file) && nrow(metadata_input_file) == 1 && "datapath" %in% names(metadata_input_file)) {
      metadata_input_file$datapath = normalizePath(metadata_input_file$datapath, winslash = "/")

      input_list <- shiny::reactiveValuesToList(input)
      tryCatch({
        # reset
        datastore_reactive$dataset$samples = NULL
        datastore_reactive$fileSamples = NULL
        datastore_reactive$sample_groups = NULL

        # block UI and parse data
        toggle_inputs(input_list, enable_components = FALSE)

        datastore_reactive$dataset = msdap::import_sample_metadata(datastore_reactive$dataset, filename = metadata_input_file$datapath)
        datastore_reactive$history = c(datastore_reactive$history, sprintf('dataset = msdap::import_sample_metadata(dataset, filename="%s")', metadata_input_file$datapath))

        datastore_reactive$fileSamples = metadata_input_file$datapath

        if(validate_input_samples(datastore_reactive)) {
          output$fileSamplesResult <- shiny::renderText(metadata_input_file$datapath)
          # remove extra information added by msdap in this view
          print(datastore_reactive$dataset$samples)
          cols_remove = c("sample_index", grep("^(detected|all)_(peptides|proteins)$", colnames(datastore_reactive$dataset$samples), value=T))
          tib_samples_plot = datastore_reactive$dataset$samples %>% select(!!setdiff(colnames(datastore_reactive$dataset$samples), cols_remove))
          output$tblSamples = shiny::renderTable(tib_samples_plot)
          datastore_reactive$sample_groups = datastore_reactive$dataset$samples %>% distinct(group) %>% pull()

          shinyjs::enable("btnScreenNext")
          datastore_reactive$screens_validated_max = 3

          # create contrast by simulating a click (works even though that part of the UI is hidden until user clicks 'next')
          shiny::observe({
            shinyjs::click("btnAddContrast")
          }, label = "OBSERVEbtnAddContrast")

        } else {
          shinyjs::disable("btnScreenNext")
          datastore_reactive$screens_validated_max = 2
          output$fileSamplesResult = shiny::renderText("No file selected")
          output$tblSamples = shiny::renderTable(data.frame(status="data has not been loaded yet..."))
        }

        # ##################### LOAD MOCK DATA
        # cat("loading...")
        # dataset = sample_metadata_custom(dataset, group_regex_array = c(A = "007|009|011", B = "008|010|012") )
        # datastore_reactive$dataset = setup_contrasts(dataset, contrast_list = list(c("A", "B")))
        # print(dataset$samples %>% select(sample_id, group))
        # cat("done!\n")
        # append_log("HARDCODED CONTRAST", type = "warning")
        # datastore_reactive$dataset = setup_contrasts(dataset, contrast_list = list(c("wt", "ko")))
        # append_log("HARDCODED CONTRAST", type = "warning")
        # ##################### LOAD MOCK DATA
      }, error = function(err) { print(err); warning(err) })

      toggle_inputs(input_list, enable_components = TRUE)
    }
  }, ignoreInit = TRUE)



  ####################### SETTINGS SCREEN

  shiny::observeEvent(input$btnAnalysisQuickstart, {
    print("input$btnAnalysisQuickstart")
    input_list <- shiny::reactiveValuesToList(input)
    tryCatch({
      # step 1) input validation
      if(is.null(datastore_reactive$dataset) || length(datastore_reactive$dataset) == 0) {
        append_log("ERROR: no dataset", type="error")
        return()
      }

      ## disable UI controls
      toggle_inputs(input_list, enable_components = FALSE)


      ## validate dataset integrity (should be fine as we have checks at intermediate step already)
      msdap::check_dataset_integrity(datastore_reactive$dataset)


      ## extract contrasts from UI and inject into dataset (don't do this runtime while user is fiddling with UI, waste of resources)
      # from list of reactives to a list of contrast strings
      l = lapply(shiny::reactiveValuesToList(contrast_reactive_values), function(contr_mod_result) {
        lapply(contr_mod_result, function(x) x())
      })
      # remove empty contrasts
      l = l[lengths(lapply(l,unlist,use.names=F,recursive=F)) > 0]

      if(length(l) > 0) {
        # actual contrast setup, includes most stringent contrast validation
        datastore_reactive$dataset = msdap::setup_contrasts(datastore_reactive$dataset, l)
        # from a list to the R command required to generate this list
        datastore_reactive$history = c(datastore_reactive$history, sprintf('dataset = msdap::setup_contrasts(dataset, %s)', r_variable_to_string(l)) )
      }



      ## analysis_quickstart()


      cmd = sprintf('dataset = msdap::analysis_quickstart(
        dataset,
        filter_min_detect = %d,
        filter_fraction_detect = %s,
        filter_min_quant = %d,
        filter_fraction_quant = %s,
        filter_by_contrast = %s,
        filter_topn_peptides = %d,
        filter_min_peptide_per_prot = %d,
        norm_algorithm = %s,
        dea_algorithm = %s,
        dea_qvalue_threshold = %s,
        dea_log2foldchange_threshold = %s,
        output_qc_report = %s,
        output_peptide_plots = "%s",
        output_abundance_tables = %s,
        output_dir = "%s",
        output_within_timestamped_subdirectory = TRUE
      )',
                    as.integer(input$mindetect),
                    as.numeric(input$fracdetect) / 100,
                    as.integer(input$minquant),
                    as.numeric(input$fracquant) / 100,
                    input$filter_by_contrast %in% "by_contrast",
                    as.integer(input$filter_topn_peptides),
                    as.integer(input$filter_min_peptide_per_prot),
                    r_variable_to_string(input$norm_algorithm),
                    r_variable_to_string(input$dea_algorithm),
                    as.numeric(input$dea_qvalue_threshold) / 100,
                    ifelse(as.numeric(input$dea_log2foldchange_threshold)<0, NA, as.numeric(input$dea_log2foldchange_threshold)),
                    "report" %in% input$output_files,
                    ifelse("peptide_plots" %in% input$output_files, "signif", "none"),
                    "abundance_tables" %in% input$output_files,
                    datastore_reactive$output_dir
      )
      datastore_reactive$history = c(datastore_reactive$history, cmd)
      history_ <<- isolate(datastore_reactive$history)

      datastore_reactive$dataset = msdap::analysis_quickstart(
        datastore_reactive$dataset,
        filter_min_detect = as.integer(input$mindetect),
        filter_fraction_detect = as.numeric(input$fracdetect) / 100,
        filter_min_quant = as.integer(input$minquant),
        filter_fraction_quant = as.numeric(input$fracquant) / 100,
        filter_by_contrast = input$filter_by_contrast %in% "by_contrast",
        filter_topn_peptides = as.integer(input$filter_topn_peptides),
        filter_min_peptide_per_prot = as.integer(input$filter_min_peptide_per_prot),
        norm_algorithm = input$norm_algorithm,
        dea_algorithm = input$dea_algorithm,
        dea_qvalue_threshold = as.numeric(input$dea_qvalue_threshold) / 100,
        dea_log2foldchange_threshold = ifelse(as.numeric(input$dea_log2foldchange_threshold)<0, NA, as.numeric(input$dea_log2foldchange_threshold)),
        output_qc_report = "report" %in% input$output_files,
        output_peptide_plots = ifelse("peptide_plots" %in% input$output_files, "signif", "none"),
        output_abundance_tables = "abundance_tables" %in% input$output_files,
        output_dir = datastore_reactive$output_dir,
        # dump_all_data = T,
        output_within_timestamped_subdirectory = TRUE
      )


      ### log summary to console
      x = testthat::capture_output_lines(print_dataset_summary(datastore_reactive$dataset), width = 1000)
      x = grep("^\\s*[<#].*", x, value=T, invert=T) # remove comments
      x = paste0(x, collapse="\n")
      datastore_reactive$log = c("\n---------------------------------------------------------------------------------------------------\n\n",
                                 paste0("<font color=\"#018023\">[", format(Sys.time(), "%H:%M:%S"), "]</font>\n"),
                                 paste0("<font color=\"#018023\">", x, "</font>\n"),
                                 "\n---------------------------------------------------------------------------------------------------\n\n",
                                 datastore_reactive$log)


      # ##################### LOAD MOCK DATA
      # cat("loading...")
      # ds = dataset
      # ds$proteins$classification = regex_classification(ds$proteins$fasta_headers, regex=c(human="_HUMA", yeast="_YEAS", discard="_ECOL"))
      # print(table(ds$proteins$classification))
      # ds = filter_dataset(ds,
      #                     filter_min_detect = 3,
      #                     norm_algorithm = "vwmb",
      #                     by_group = F, all_group = T, by_contrast = F)
      # ds = dea(ds, algo_de = "ebayes", qval_signif = 0.05, fc_signif = NA)
      # tib_plot = left_join(ds$de_proteins, ds$proteins, by="protein_id") %>% filter(classification %in% c("human", "yeast"))
      # print(tib_plot %>%
      #         group_by(classification) %>%
      #         summarise(`5% FDR` = sum(qvalue <= 0.05),
      #                   `5% FDR AND foldchange threshold` = sum(signif)))
      # datastore_reactive$dataset = ds
      # cat("done!\n")
      ##################### LOAD MOCK DATA
    }, error = function(err) { print(err); warning(err) })

    toggle_inputs(input_list, enable_components = TRUE)

    shiny::showModal(shiny::modalDialog(
      title = "analysis done",
      "check the output log at the bottom for a status report"
    ))

  }, ignoreInit = TRUE)


  onStop(msdap_shiny_server_stop)
}


# cleanup global variables we created
msdap_shiny_server_stop = function() {
  # importantly, get rid of the reactive context @ logger
  suppressWarnings(try(rm(logger_msdap, envir = .GlobalEnv)))
  suppressWarnings(try(rm(history_, envir = .GlobalEnv)))
  suppressWarnings(try(rm(log_, envir = .GlobalEnv)))
  suppressWarnings(try(rm(contrast_reactive_values, envir = .GlobalEnv)))
}



### for local testing
## setwd("C:\\PHD\\code\\R\\msdap")
## setwd("C:\\PHD\\code\\R\\shiny_docker\\dc") # setwd('/data')
# library(shiny); library(shinyFiles); library(shinyjs); library(sortable)
# devtools::load_all(); shiny::shinyApp(msdap_shiny_ui(), msdap_shiny_server)

# app = shinyApp(msdap_shiny_ui, msdap_shiny_server)
# runApp(app, host ="0.0.0.0", port = 3841, launch.browser = T)
