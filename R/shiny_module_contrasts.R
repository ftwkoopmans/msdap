
#'
#' @importFrom shiny tagList fluidRow column NS tags div actionButton br HTML textOutput
#' @importFrom sortable sortable_js sortable_options sortable_js_capture_input
#' @importFrom htmlwidgets JS
contrast_creator_ui = function(id, item_list = c("ctrl1","ctrl2","cond1","cond2"), random_param = "bla") {
  ns = shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      id = ns("fluidrow"), # make sure the entire row has an ID, so we can cleanly delete from UI
      class = "panel panel-primary",
      shiny::div(class="panel-body",
          shiny::column(
            width = 3,
            shiny::tags$div(
              class = "panel panel-default",
              shiny::tags$div(class = "panel-heading", "unused"),
              shiny::tags$ul(class = "list-group", id = ns("sort1"), items_to_html_list(item_list))
            )
          ),
          shiny::column(
            width = 3,
            shiny::tags$div(
              class = "panel panel-default",
              shiny::tags$div(class = "panel-heading", "A"),
              shiny::tags$ul(class = "list-group", id = ns("sort2")),
              shiny::HTML("&nbsp;") # some extra whitespace at bottom to make dragging easier
            )
          ),
          shiny::column(
            width = 3,
            shiny::tags$div(
              class = "panel panel-default",
              shiny::tags$div(class = "panel-heading","B"),
              shiny::tags$ul(class = "list-group",id = ns("sort3")),
              shiny::HTML("&nbsp;") # some extra whitespace at bottom to make dragging easier
            )
          ),
          shiny::column(
            width = 3,
            shiny::tags$div(
              shiny::actionButton(ns("btnRemove"), "remove contrast", class="btn btn-danger"),
              shiny::tags$p(shiny::br()),
              shiny::tags$p(shiny::textOutput(ns("txtOut")))
            )
          )
      )
    ),
    sortable::sortable_js(
      css_id = ns("sort1"),
      options = sortable::sortable_options(
        group = list(
          name = ns("sortGroup"),
          put = TRUE
        ),
        swap = TRUE,
        swapClass = "sortable-swap-highlight",
        sort = FALSE#,
        # onSort = sortable_js_capture_input("sort_vars")
      )
    ),
    sortable::sortable_js(
      css_id = ns("sort2"),
      options = sortable::sortable_options(
        group = list(
          name = ns("sortGroup"),
          put = htmlwidgets::JS("function (to) { return true; }"),
          pull = TRUE
        ),
        swap = TRUE,
        swapClass = "sortable-swap-highlight",
        onSort = sortable::sortable_js_capture_input(ns("sort2_values"))
      )
    ),
    sortable::sortable_js(
      css_id = ns("sort3"),
      options = sortable::sortable_options(
        group = list(
          name = ns("sortGroup"),
          put = htmlwidgets::JS("function (to) { return true; }"),
          pull = TRUE
        ),
        swap = TRUE,
        swapClass = "sortable-swap-highlight",
        onSort = sortable::sortable_js_capture_input(ns("sort3_values"))
      )
    )
  )
}



#'
#' @importFrom shiny reactive renderText validate need
contrast_creator_server = function(input, output, session) {
  a = shiny::reactive({
    input$sort2_values %>% trimws()
  })

  b = shiny::reactive({
    input$sort3_values %>% trimws()
  })

  output$txtOut = shiny::renderText({
    # only throw error warning if either is missing, not for both
    hasA = length(a()) > 0
    hasB = length(b()) > 0

    shiny::validate( # note; need() takes FALSE as the error state. so we invert the conditional for error state(eg; hasB and not hasA)
      shiny::need(!(hasB && !hasA), "A is empty"),
      shiny::need(!(hasA && !hasB), "B is empty")
    )

    if(hasA && hasB) {
      paste(paste(a(), collapse=","), " vs ", paste(b(), collapse=","))
    }
  })

  list(A=a, B=b)
}



#'
#' @importFrom shiny tags tagList fluidRow column tagAppendAttributes actionButton textOutput
module_contrast_ui = function() {
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(width=8, shiny::tags$div("If you want to perform A/B testing, drag&drop from the list of available sample groups (initially positioned in the 'unused' column) into 'A' and 'B' accordingly", style="display: inline-block; vertical-align: middle; float: none;")),
      shiny::column(width=3, shiny::actionButton("btnAddContrast", label = "add contrast", class = "btn btn-success"))
    ),
    shiny::fluidRow(
      id="rowContrastValidation",
      shiny::tagAppendAttributes(shiny::textOutput("contrastErrors"), style='color: red; font-weight: bold;')
    )
  )
}
