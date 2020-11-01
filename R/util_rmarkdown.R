
#' split a ggplot object into 2 ggplots; original plot minus legend and only the legend
#'
#' @param p ggplot object
#' @importFrom ggpubr get_legend as_ggplot
ggplot_split_legend <- function(p) {
  return(list(plot=p + theme(legend.position = "none"),
              legend=ggpubr::as_ggplot(ggpubr::get_legend(p))))
}

# #' if there are more than N legend entries, ggplot objects will be split
# #'
# #' @param plotlist list of ggplot objects
# #' @param max_items if there are more than N legend entries, the plot and legend will be split
# #' @importFrom ggpubr get_legend as_ggplot
# ggplot_list_split_legend <- function(plotlist, n = 10) {
#   result = list()
#   for(p in plotlist) {
#     leg = ggpubr::get_legend(p)
#     n_legend =  0  # TODO: how can we extract the number of legend items from a ggplot object?
#     if(n_legend > n) {
#       # disable legend + store plot
#       p = p + theme(legend.position = "none")
#       result[[length(result) + 1]] = p
#       # store legend as a new plot
#       result[[length(result) + 1]] = ggpubr::as_ggplot(leg) # TODO: guide label size
#     } else {
#       # don't split, store plot as-is
#       result[[length(result) + 1]] = p
#     }
#   }
#
#   return(result)
# }


#' placeholder title
#' https://github.com/yihui/knitr/issues/1494
#' http://michaeljw.com/blog/post/subchunkify/
#' @param g todo
#' @param unique_chunk_id todo
#' @param fig_height todo
#' @param fig_width todo
subchunkify <- function(g, unique_chunk_id, fig_height=7, fig_width=5) {
  unique_chunk_id = gsub("[^a-zA-Z0-9_]+", "", unique_chunk_id)

  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')

  sub_chunk <- paste0("
  `","``{r sub_chunk_", unique_chunk_id, ", echo=F, message=F, warning=F, fig.height=", fig_height, ", fig.width=", fig_width, "}",
                      "\n(",
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")

  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}



#' placeholder title
#' https://stackoverflow.com/a/16580525
#' @param x matrix/data.frame
#' @param caption todo
#' @param align todo
#' @param scalebox todo
rmarkdown_xtable_custom = function(x, caption = NULL, align=NULL, scalebox = 0.9) {
  # stripes = NULL
  # if(nrow(x) > 1) {
  #   xtbl_rws <- seq(1, (nrow(x)-1), by = 2)
  #   xtbl_col <- rep("\\rowcolor[gray]{0.95}", length(xtbl_rws))
  #   stripes = list(pos = as.list(xtbl_rws), command = xtbl_col)
  # }

  ## handle long tables, assuming latex packages are not available
  # eg;  tabular.environment = "longtable"  @   https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
  if(nrow(x) > 55) {
    i = 1
    while(i < nrow(x)) {
      if(i != 1) {
        cat('\n\\newpage \n')
      }

      j = min(nrow(x), i + 49)
      print(xtable::xtable(x[i:j,,drop=F], caption = caption, align=align, booktabs = TRUE), scalebox = scalebox, floating = FALSE, latex.environments = "center", include.rownames = FALSE, comment=FALSE) #, add.to.row = stripes)
      i = i + 50
    }
  } else {
    print(xtable::xtable(x, caption = caption, align=align, booktabs = TRUE), scalebox = scalebox, floating = FALSE, latex.environments = "center", include.rownames = FALSE, comment=FALSE) #, add.to.row = stripes)
  }
}
