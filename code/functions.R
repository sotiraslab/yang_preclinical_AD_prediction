filter_progressor <- function(.data, amypos_col, max_time_to_progression = 3, recovered = FALSE) {

    if (recovered) {
        progressor_class <- c("progressor", "progressor_recovered")
    } else {
        progressor_class <- c("progressor")
    }

    .data %>% filter(
        .progressor %in% progressor_class,
        .is_valid_baseline,
        .data[[amypos_col]],
        .time_to_progression <= max_time_to_progression,
    )

}

filter_stable <- function(.data, amypos_col, min_time_stable = 3) {

    .data %>% filter(
        .progressor == "stable",
        .is_valid_baseline,
        .data[[amypos_col]],
        .time_stable >= min_time_stable
    )

}

compute.pacc <- function(df, pacc.columns,
                         cn.mask, min.required = 2,
                         higher.better = NULL) {
  cn.data <- df[cn.mask, ]
  n = nrow(df)
  k = length(pacc.columns)
  
  normed.scores <- matrix(data=NA, nrow=n, ncol=k)
  normed.scores <- as.data.frame(normed.scores)
  colnames(normed.scores) <- pacc.columns
  
  if (is.null(higher.better)) {
    higher.better <- rep(TRUE, k)
  }
  
  for (i in 1:k) {
    col <- pacc.columns[i]
    mu <- mean(cn.data[[col]], na.rm = T)
    s <- sd(cn.data[[col]], na.rm = T)
    z <- (df[[col]] - mu) / s
    
    if (! higher.better[i]) {
      z <- (-1 * z)
    }
    
    normed.scores[, i] <- z
  }
  
  pacc.score <- rowMeans(normed.scores, na.rm = T)
  count.present <- rowSums(! is.na(normed.scores))
  pacc.score <- ifelse(count.present >= min.required, pacc.score, NA)
  
  return(pacc.score)
}
