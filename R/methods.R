#' Methods for FMR fits
#'
#' \code{print} and \code{summary} show the estimates of \eqn{\beta} (FMR with the
#' bootstrap-\eqn{t} and the \eqn{\tau}-Wald interval, and the comparators), the first stage
#' and the test of \eqn{H_0:\alpha\equiv0}; \code{coef} returns \eqn{\hat\beta};
#' \code{confint} the FMR intervals; \code{plot} draws \eqn{\hat\alpha(t)} with its pointwise
#' band.
#' @param x,object An \code{"fmr"} object.
#' @param digits Significant digits.
#' @param parm,level Ignored (the level is that of the fit).
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Further arguments passed to \code{plot}.
#' @name fmr-methods
NULL

#' @rdname fmr-methods
#' @export
print.fmr <- function(x, digits = 3, ...) {
  cat("Functional Mendelian randomization\n")
  cat(sprintf("n = %d, positions p = %d, FPCA scores k = %d, lambda (GCV) = %s\n",
              x$n, x$p, x$k, format(x$lambda, digits = 3)))
  cat(sprintf("First stage: F = %.1f, partial R^2 = %.3f\n\n", x$first_stage[["F"]],
              x$first_stage[["partial_R2"]]))
  tab <- x$coef
  tab[, -1] <- lapply(tab[, -1], function(v) signif(v, digits))
  names(tab)[4:5] <- sprintf("%g%% %s", 100 * x$level, c("lower", "upper"))
  print(tab, row.names = FALSE)
  if (is.finite(x$talpha[["T_alpha"]]))
    cat(sprintf("\nTest of no direct effect, H0: alpha(t) = 0: T_alpha = %s, p = %s\n",
                format(x$talpha[["T_alpha"]], digits = digits), format(x$talpha[["p"]], digits = digits)))
  invisible(x)
}

#' @rdname fmr-methods
#' @export
summary.fmr <- function(object, ...) {
  print(object, ...)
  cat(sprintf("\nBootstrap-t: %d of %d draws usable, critical value c-hat = %.3f (V-hat SE %.4g)\n",
              object$n_boot, object$B, object$c_hat, object$se_vhat))
  cat(sprintf("Wald interval: beta-hat(lambda) +- z SE(%g lambda)\n", object$tau))
  cat(sprintf("Second-stage effective degrees of freedom: %.2f; FPCA retained %d components\n",
              object$df, object$K_retained))
  invisible(object)
}

#' @rdname fmr-methods
#' @export
coef.fmr <- function(object, ...) c(beta = object$beta)

#' @rdname fmr-methods
#' @export
confint.fmr <- function(object, parm, level, ...) {
  tab <- object$coef[grepl("^FMR", object$coef$method), c("method", "lower", "upper")]
  rownames(tab) <- NULL
  tab
}

#' @rdname fmr-methods
#' @export
plot.fmr <- function(x, xlab = "position", ylab = expression(hat(alpha)(t)),
                     main = "Direct (pleiotropic) effect", ...) {
  a <- x$alpha
  graphics::plot(a$pos, a$alpha, type = "l", lwd = 2, col = "firebrick", ylim = range(a$lower, a$upper, 0),
                 xlab = xlab, ylab = ylab, main = main, ...)
  graphics::polygon(c(a$pos, rev(a$pos)), c(a$lower, rev(a$upper)), col = grDevices::adjustcolor("firebrick", 0.15),
                    border = NA)
  graphics::abline(h = 0, lty = 3)
  graphics::legend("topleft", c(expression(hat(alpha)(t)), sprintf("pointwise %g%% band", 100 * x$level)),
                   col = c("firebrick", grDevices::adjustcolor("firebrick", 0.3)), lwd = c(2, 8), bty = "n")
  invisible(x)
}
