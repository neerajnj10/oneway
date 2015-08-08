#' Perform a oneway analysis of variance
#'
#' @param z A list of responses grouped by factor level
#' @param ... Addition arguments used by other S3 methods
#'
#' @seealso \code{\link{oneway.factor}}, \code{\link{oneway.formula}},
#' \code{\link{summary.oneway}}, \code{\link{plot.oneway}}, \code{\link{lsmeans}}
#' @examples
#' data(coagulation)
#' attach(coagulation)
#' oneway(split(coag, diet))
#' oneway(diet, coag)
#' oneway(coag~diet)
#' summary.oneway(oneway(coag~diet))
#'
#' lsmeans(oneway(coag~diet))
#' plot(oneway(coag~diet))
#'
#' @export
oneway <- function(z, ...) UseMethod("oneway")
#' @export
oneway.default <- function(z, ...) {


  l <- length(z)
  N <- length(unlist(z))
  sums <- sum(unlist(z)^2)
  g <- sum(sapply(z, function(x) sum(x)^2/length(x)))
  h <- sum(unlist(z))^2/N

  sum.wit <- g - h ## within
  sum.bet <- sums - g ## between
  df.wit <- l - 1
  df.bet <- N - l
  m <- unlist(lapply(z, mean))
  n <- unlist(lapply(z, length))

  if(is.null(names(z))){
    groups <- as.character(1:length(z))
  }
  else{
    groups <- names(z)
  }
  names(m) <- groups
  ## create object, add class, & return
  res <- list(df = c(df.wit, df.bet), SS = c(sum.wit, sum.bet),
              groups = groups, call = match.call(), data = z,
              means = m, n = n, N = N, l = l)
  class(res) <- "oneway"
  return(res)
}


#' oneway.factor
#'
#' S3 method for \code{\link{oneway}} using a response vector and factor
#'
#' @param z A factor of levels for each observation
#' @param y A vector of responses
#' @param ... Addition arguments used by other S3 methods
#'
#' @seealso \code{link{oneway}}
#'
#' @export
oneway.factor <- function(z, y, ...) {
  responses <- oneway.default(split(y, z))
  responses$call <- match.call()
  responses
}
#' oneway.formula
#'
#' An S3 method for \code{\link{oneway}} for formulas.
#'
#' @param formula A formula of the form \code{response~factor}
#' @param data An (optional) data frame used by the formula.
#' @param ... Addition arguments used by other S3 methods
#'
#' @seealso \code{\link{oneway}}
#'
#' @export
oneway.formula <- function(formula, data=list(), ...) {
  frames <- model.frame(formula, data)
  responses <- oneway.factor(frames[,2], frames[,1])
  responses$call <- match.call()
  responses
}

#' @export
print.oneway <- function(x, ...) {
  print(x$call)
  cat("\nWithin SS:", x$SS[1], "on", x$df[1],
      "degrees of freedom.\n")
  cat("Between SS:", x$SS[2], "on", x$df[2],
      "degrees of freedom.\n")
}


#' summary.oneway
#'
#' Creates an Analysis of Variance table for a \code{oneway} object
#'
#' @param object An object of class \code{oneway}
#' @param ... Addition arguments used by other S3 methods
#'
#' @export
summary.oneway <- function(object, ...) {
  attach(object)

  SS.tot <- SS[1] + SS[2]
  DF.tot <- df[1] + df[2]

  MS.wit <- SS[1]/df[1]
  MS.bet <- SS[2]/df[2]

  F <- MS.wit/MS.bet
  p <- pf(F, df[1], df[2], lower.tail = FALSE)

  aov.table <- with(object, cbind(DF = c(df, DF.tot),
                                  SS = c(SS, SS.tot),
                                  MS = c(MS.wit, MS.bet, NA),
                                  F = c(F, NA, NA),
                                  "Pr(>F)" = c(p, NA, NA)))
  rownames(aov.tab) <- c("Among Group", "Within Group",
                         "Total")
  res <- list(call=call, aov.table=aov.table, groups=groups, means=means,
              P = p, MS = c(MS.wit,MS.bet), n = n, N = N, l = l)
  class(res) <- "summary.oneway"
  detach(object)
  return(res)
}

#' @export
print.summary.oneway <- function(x, ...) {

  cat("print: \n\t")
  print(x$call)
  cat("\nMeans:\n")
  print(x$means)
  cat("\n")


  # AOV Table

  printCoefmat(x$aov.table, P.values=TRUE, has.Pvalue=TRUE, signif.stars=TRUE, na.print="")
}


#' Perform Fisher's LSD
#'
#' Testing  pairwise differences using Fisher's LSD proceduce on a
#' \code{oneway} object
#'
#' @param object A \code{oneway} object
#'
#' @export
lsmeans <- function(object) UseMethod("lsmeans")


#' @export
lsmeans.oneway <- function(object, ...) {
  object <- summary(object)

  comparison <- function(i, j){
    d <- object$means[i] - object$means[j]
    SE <- sqrt(object$MS[2]*(1/object$n[i] + 1/object$n[j]))
    d.se <- d/SE
    round(2*pt(abs(d.se), object$N-object$l, lower.tail=FALSE),4)
  }
  p.se <- pairwise.table(compare.levels=compare,
                         level.names=object$groups,
                         p.adjust.method="none")
  result <- list(p.value=p.se, call=match.call())
  class(result) <- "lsmeans"
  result
}

#' Creating boxplot of groups in a \code{oneway} object
#'
#' @param x A \code{oneway} object
#' @param xlab X label
#' @param ylab Y label
#' @param main Main plot title
#' @param ... Optional graphing arguments to be passed to
#' \code{boxplot}
#'@param names Names of factor levels
#' @seealso \code{link{boxplot}}
#'
#' @export
plot.oneway <- function(x, names=x$groups, xlab="Grouping variable",ylab="Response variable", main=capture.output(x$call)){
  boxplot(x=x$data, names=names, xlab=xlab,  main=main)
}


#' @export
print.lsmeans <- function(x, ...){
  cat("Call:\n\t")
  print(x$call)
  cat("\nFisher's LSD Table\n")
  cat("\nP-Values:\n")
  print.table(x$p.value, na.print="-")
}
