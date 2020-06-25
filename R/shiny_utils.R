#'
#' @importFrom shiny tags
items_to_html_list <- function(x) {
  sapply(x, function(n) {
    shiny::tags$li(
      class = "list-group-item",
      tags$span(class = "glyphicon glyphicon-move"),
      tags$b(n)
    )
  }, simplify = F)
}



# https://www.r-bloggers.com/shiny-add-removing-modules-dynamically/
remove_shiny_inputs <- function(id, .input) {
  invisible(
    lapply(grep(id, names(.input), value = TRUE), function(i) {
      .subset2(.input, "impl")$.values$remove(i)
    })
  )
}


#'
#' @importFrom shiny span tagAppendAttributes icon
create_help_icon = function(x) {
  if(length(x) == 1 && !is.na(x) && x != "") {
    return(shiny::span(shiny::tagAppendAttributes(shiny::icon("question-circle"), title=x)))
  }
}


#'
#' @importFrom htmltools tagSetChildren
insert_child_before = function(tag, new_element) {
  htmltools::tagSetChildren(tag, new_element, tag$children)
}
