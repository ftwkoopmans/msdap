### original code from this commit: https://github.com/statOmics/MSqRobSum/commit/2875dcde1f89578685ed0a3316d6f66fa510b732
### here adapted to get optimal multiprocessing for vastly reduced computation times  &  compatability with latest dplyr

### FRANK: removed all documentation/examples/"addition functions we don't need" and kept minimal set of imports
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import lme4
#' @import limma
### FRANK: add 2 imports below
#' @importFrom foreach %dopar%
#' @importFrom parallel clusterExport
#
### FRANK: below code is adapted only on lines/sections tagged with my name
msqrobsum <- function(
  data, formulas, group_vars = 'protein', contrasts = NULL
  , mode = c('msqrobsum','msqrob', 'sum')
  , robust_lmer_iter = 'auto'
  , squeeze_variance = TRUE
  , p_adjust_method = c("BH",p.adjust.methods)
  , keep_model = FALSE
  , rlm_args = list(maxit = 20L)
  , lmer_args = list(control = lmerControl(calc.derivs = FALSE))
  , parallel_args = list(strategy = 'multisession')
  ## EXPERIMENTAL OPTIONS
  , type_df = 'traceHat'
  ## , type_df = c('traceHat','conservative')
  , squeeze_covariate = FALSE
  , fit_fun = do_mm
){
  # rlang::exec(plan, !!!parallel_args) ### FRANK: remove this line, updated multiprocessing below
  mode = match.arg(mode)

  if (robust_lmer_iter == 'auto')
    robust_lmer_iter <- ifelse(mode == 'msqrob', 20L, 1L)
  type_df = match.arg(type_df)
  p_adjust_method = match.arg(p_adjust_method)
  if (missing(formulas) & mode == 'sum') formulas = ''
  if ((mode != 'sum') & length(unlist(map(formulas,findbars))) == 0)
    stop('msqrobsum only works with formulas with at least 1 random effect specified')
  if (is(data,'MSnSet')) data <- MSnSet2df(data)
  if (is(contrasts, 'character')){
    if (!all(c(formulas) %>% map_lgl(~{contrasts %in% map_chr(findbars(.),all.vars)})))
      stop('Contrasts can only be constructed from variable name when the variable is specified as random effect')
    contrasts <- make_simple_contrast(data, contrasts)}
  group_vars <- map_chr(group_vars, ~match.arg(.,names(data)))
  model_vars <- c(formulas) %>% map(all.vars) %>% unlist %>% unique


  ### FRANK: begin multiprocessing fix
  ### begin original code
  # df <- data %>% select(group_vars,model_vars, 'expression', 'feature','sample') %>%
  #   group_by_at(group_vars) %>% nest %>%
  #   mutate(mm = future_lapply(data, fit_fun, robust_lmer_iter = robust_lmer_iter
  #                             , rlm_args = rlm_args, lmer_args =lmer_args
  #                             , contrasts = contrasts, keep_model = keep_model
  #                             , formulas = formulas, mode = mode)) %>%
  #   unnest(mm) %>% ungroup
  ### end original code
  fit_fun = match.fun(fit_fun)
  # prep data
  df <- data %>% select(group_vars,model_vars, 'expression', 'feature','sample') %>% group_by_at(group_vars) %>% nest
  ## actual code change, from future package to foreach. Internal version number for updates to this codebase; v2.2
  if(!exists("cl")) stop("first, setup a cluster for multiprocessing")
  ## functions from global/parent environment
  parallel::clusterExport(cl, varlist = c("do_mm", "do_lmerfit", "summarizeRobust", "getVcovBetaBUnscaled", "getBetaB", "calculate_df", "getDf", "get_contrasts"), envir = parent.env(environment())) # .GlobalEnv)
  # arguments from local environment
  parallel::clusterExport(cl, varlist = c("robust_lmer_iter", "rlm_args", "lmer_args", "contrasts", "keep_model", "formulas", "mode"), envir = environment()) # environment(fun = NULL)
  df$mm = foreach::foreach(d = df$data, .packages = c("lme4", "purrr", "dplyr", "tibble", "MASS")) %dopar% {
    fit_fun(d, robust_lmer_iter = robust_lmer_iter
            , rlm_args = rlm_args, lmer_args =lmer_args
            , contrasts = contrasts, keep_model = keep_model
            , formulas = formulas, mode = mode)
  }
  ### FRANK: end multithreading fix

  df = df %>% unnest(cols = mm) %>% ungroup() ### FRANK; add ungroup

  if (mode == 'sum') return(df)
  ## Return also failed ones afterward
  rows_fail = apply(df, 1, function(x) any(lengths(x) == 0 | is.na(x)) ) ### FRANK; fix fail group filter @ dplyr 1.0 compatability
  df_failed <- df[rows_fail,] ### FRANK
  df <- df[!rows_fail,] ### FRANK
  # df_failed <- filter(df, is.na(df))### FRANK
  # df <- filter(df, !is.na(df))### FRANK
  if(!nrow(df)) {warning("No models could be fitted"); return(df_prot_failed)}
  ## Squeeze variance
  df <- mutate(df, sigma_post = sigma, df_prior = 0, sigma_prior = 0)
  if(squeeze_variance) {
    ## no shrinkage if df < 1 because greatly influences emp. bay.
    ## TODO check only works when at least 3 protein
    id = df$df >= 1L
    sq = squeezeVar(df$sigma[id]^2, df$df[id])
    ## experimental: use intercept of model as covariate to
    ## squeeze variance to correct for  mean variance
    if(squeeze_covariate) sq = squeezeVar(df$sigma[id]^2, df$df[id], covariate = df$intercept[id])
    df[id,] = mutate(df[id,],
                     sigma_prior = sqrt(sq$var.prior), sigma_post = sqrt(sq$var.post), df_prior = sq$df.prior)} ### FRANK
  if(type_df == "conservative"){
    ## Calculate df on protein level, assumption is that there is only one protein value/run,
    ## TODO remove conservative option? Now broken, fix vars = colnames(...)
    df <- mutate(df, df_prior = 0
                 , df = map2_dbl(data, model,~calculate_df(.x,.y, vars = colnames(pData(msnset)))))}
  df <- mutate(df, df_post = df + df_prior)

  ## Calculate q and p values for contrasts (if contrasts are found)
  df <- df %>% select(group_vars ,df_post,sigma_post, contrasts) %>% unnest(contrasts) %>% ### FRANK: specify the unnest
    # df <- df %>% select(group_vars ,df_post,sigma_post, contrasts) %>% unnest %>% ### FRANK: orig line
    ## skip if there are no contrasts
    {if(nrow(.)) {
      mutate(., se = sigma_contrast * sigma_post
             , t = logFC/se
             , pvalue = pt(-abs(t), df_post) * 2) %>%
        select(-sigma_post, - df_post) %>%
        group_by(contrast) %>%
        mutate(qvalue = p.adjust(pvalue, method = p_adjust_method)) %>%
        group_by_at(group_vars) %>% group_nest(.key="contrasts") ### FRANK: adjust to dplyr updates
      # group_by_at(group_vars) %>% nest(contrasts = - group_vars) ### FRANK: orig code
    } else .} %>%
    select(group_vars, contrasts) %>%
    left_join(select(df, -contrasts),., by = group_vars)

  ## mutate(df_failed, contrasts = list(tibble())) %>% bind_rows(df,.)
  ## give empty tibble instead of NULL when no contrast (nicer to work with posthoc)
  # bind_rows(df,df_failed) %>% ### FRANK
  bind_rows(df, df_failed %>% select(-contrasts)) %>% ### FRANK
      mutate(contrasts = map(contrasts, ~{if (is(.x,'tbl')) .x else tibble()}))
}

MSnSet2df <- function(msnset){
  ## Converts Msnset to a tidy dataframe
  ## Always creates feature and vector column so these shouldn't be defined by user.
  ## convenient for downstream analysis steps.
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(any(c("sample", "feature", "expression") %in% c(colnames(fData(msnset)),colnames(pData(msnset))))){
    stop("Column names in the \"fData\" or \"pData\" slot of the \"msnset\" object cannot be named
         \"sample\", \"feature\" or \"expression\". Please rename.")
  }


  dt <- as_tibble(Biobase::exprs(msnset), rownames="feature") %>% gather(sample, expression, - feature, na.rm=TRUE) ### FRANK: update deprecated dplyr code, use as_tibble to convert matrix while maintaining rownames
  # dt <- as.data.frame(Biobase::exprs(msnset)) %>% mutate(feature = rownames(.)) %>% gather(sample, expression, - feature, na.rm=TRUE) ### FRANK: orig line
  dt <- fData(msnset) %>% mutate(feature = (rownames(.))) %>% left_join(dt,. , by = 'feature')
  dt <- pData(msnset) %>% mutate(sample = rownames(.)) %>% left_join(dt,. , by = 'sample')
  return(dt) ### FRANK: as_data_frame is a deprecated function. downstream code uses only tibbles, we return dt for now. if data.frame is desired, use as.data.frame(dt)
  # as_data_frame(dt) ### FRANK: orig line
}

#' @import lme4
#' @import tibble
do_mm <- function(d, formulas, contrasts, mode ='msqrobsum', robust_lmer_iter = 1L
                  , rlm_args = list(), lmer_args = list(), keep_model = TRUE){
  out <- list()
  if (mode != 'msqrob') {
    d <- mutate(d,expression = rlang::exec(summarizeRobust, expression, feature, sample, !!!rlm_args)) %>%
      select(-feature) %>% distinct
    out$data_summarized <- list(d)
  }
  if (mode == 'sum') return(new_tibble(out ,nrow = 1L))
  for (form in c(formulas)){
    ## TODO do lm fit if there are no random effects (all output inside loop, do while loop)
    model <- try(do_lmerfit(d, form,  robust_lmer_iter,lmer_args), silent = TRUE)
    if (is(model,"lmerMod")) break #else message(paste0('Cannot fit ', format(form)))
  }
  if (keep_model) out$model <- list(model)
  if (!is(model,"lmerMod")) return(new_tibble(out ,nrow = 1L))
  out$formula <- as.character(list(attributes(model@frame)$formula))
  out$df <- getDf(model)
  out$sigma <- sigma(model)
  out$contrasts <- list(get_contrasts(model,contrasts))
  out$intercept <- model@beta[1]
  new_tibble(out , nrow = 1L)
}

#' @importFrom MASS psi.huber
#' @import lme4
do_lmerfit <- function(df, form, robust_lmer_iter, args_lmer = list()){
  tol <- 1e-6
  fit <- rlang::exec(lmer,form, data = df, !!!args_lmer)
  ##Initialize SSE
  sseOld <- fit@devcomp$cmp['pwrss']
  while (robust_lmer_iter > 0){
    robust_lmer_iter= robust_lmer_iter-1
    res <- resid(fit)
    fit@frame$`(weights)` <- psi.huber(res/(mad(res,0)))
    fit <- refit(fit)
    sse <- fit@devcomp$cmp['pwrss']
    if(abs(sseOld-sse)/sseOld <= tol) break
    sseOld <- sse
  }
  fit
}

#' @importFrom MASS rlm
summarizeRobust <- function(expression, feature, sample, ...) {
  ## Assumes that intensities mx are already log-transformed
  ## characters are faster to construct and work with then factors
  feature <- as.character(feature)
  ##If there is only one 1 peptide for all samples return expression of that peptide
  if (length(unique(feature)) == 1L) return(expression)
  sample <- as.character(sample)
  ## modelmatrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(sample)) == 1L) sample <- rep(1,length(sample))

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  X <- model.matrix(~ -1 + sample + feature,contrasts.arg = list(feature = 'contr.sum'))
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced modelmatrix untill no singular values are present
  repeat {
    fit <- .lm.fit(X,expression)
    id <- fit$coefficients != 0
    X <- X[ , id, drop =FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm
  fit <- rlm(X, expression, ...)
  ## Only return the estimated effects effects as summarised values
  ## sampleid = seq_along(unique(sample))
  ## return(X[,sampleid,drop = FALSE] %*% fit$coefficients[sampleid])
  fit$coefficients[paste0('sample',sample)]
}

#' @import tibble
get_contrasts <- function(model, contrasts){
  ## TODO allow for contrasts from fixed effects
  ## tricky because reference level can change
  betaB <- getBetaB(model)
  vcov <- getVcovBetaBUnscaled(model)
  coefficients <- names(betaB)
  id <- coefficients %in% rownames(contrasts)
  if (!any(id)) return(new_tibble(list(), nrow = 0))
  coefficients <- coefficients[id]
  vcov <- vcov[id,id]
  betaB <- betaB[id]

  ## check for which contrasts I have data
  missing_levels <- !(rownames(contrasts) %in% coefficients)
  id <- !apply(contrasts,2,function(x){any(x[missing_levels] != 0)})
  contrasts <- contrasts[coefficients, id, drop = FALSE]
  ## If no contrasts could be found, terminate
  if (is.null(colnames(contrasts))) return(new_tibble(list(), nrow = 0))
  new_tibble(list(contrast = colnames(contrasts)
                  , logFC = logFC <- (t(contrasts)%*%betaB)[,1]
                  , sigma_contrast = sqrt(diag(t(contrasts)%*%vcov%*%contrasts))
  ),nrow = sum(id))
}

#' @import purrr
getBetaB <- function(model) {
  betaB <- c(as.vector(getME(model,"beta")),as.vector(getME(model,"b")))
  names(betaB) <- c(colnames(model@pp$X), unlist(imap(model@flist,~{paste0(.y,levels(.x))})))
  betaB
}

#' @import lme4
getVcovBetaBUnscaled <- function(model){
  X <- getME(model,"X")
  Z <- getME(model,"Z")
  vcovInv <- Matrix::crossprod(cbind2(X,Z))
  Ginv <- Matrix::solve(Matrix::tcrossprod(getME(model,"Lambda")) +
                          Matrix::Diagonal(ncol(Z),1e-18))
  i <- -seq_len(ncol(X))
  vcovInv[i,i] <- vcovInv[i,i]+Ginv
  as.matrix(Matrix::solve(vcovInv))
}

calculate_df <- function(df, model, vars){
  ## Get all the variables in the formula that are not defined in vars
  form <- attributes(model@frame)$formula
  vars_formula <- all.vars(form)
  vars_drop <- vars_formula[!vars_formula %in% vars]
  ## Sum of number of columns -1 of Zt mtrix of each random effect that does not involve a variable in vars_drop
  mq <- getME(model,'q_i')
  id <- !map_lgl(names(mq),~{any(stringr::str_detect(.x,vars_drop))})
  p <- sum(mq[id]) - sum(id)
  ## Sum of fixed effect parameters that do not involve a variable in vars_drop
  mx <- getME(model,'X')
  id <- !map_lgl(colnames(mx),~{any(stringr::str_detect(.x,vars_drop))})
  p <- p + sum(id)

  ## n is number of sample because 1 protein defined per sample
  n <- n_distinct(df$sample)
  n-p
}

#' @import lme4
getDf <- function(object){
  w <- object@frame$"(weights)"
  if (is.null(w)) w <- 1
  sigma <- sigma(object)
  sum((resid(object)* sqrt(w))^2)/sigma^2
}

make_simple_contrast <- function(data, contrast_var){
  c <- pull(data,contrast_var) %>% unique %>% paste0(contrast_var, .) %>% sort %>% as.factor
  comp <- combn(c,2,simplify = FALSE)
  condIds <- map(comp, ~which(c %in% .x))
  L <- rep(0,nlevels(c))
  L <- sapply(comp,function(x){L[x]=c(-1,1);L})
  rownames(L) <- levels(c)
  colnames(L) <- map_chr(comp, ~paste(rev(.x),collapse = '-'))
  L
}
