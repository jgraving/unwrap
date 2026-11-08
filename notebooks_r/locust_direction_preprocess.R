# Locust direction mixture: preprocessing helpers.
#
# Sourced by notebooks_r/locust_direction_mixture_model.ipynb (no
# standalone output). R twin of the Python data-prep step
# (unwrap.preprocess.rediscretize_directional); the R stack has no
# package equivalent, so the shared circular helpers and the
# rediscretization are local.

wrap <- function(x) (x + pi) %% (2 * pi) - pi

circ_mean <- function(x) atan2(mean(sin(x)), mean(cos(x)))

resultant <- function(x) Mod(mean(exp(1i * x)))

rediscretize_trajectory <- function(x, y, step_length = 6 * 0.04, tolerance = 0.0001) {
  if (length(x) < 2L) return(data.frame(X = x, Y = y))
  counts <- pmax(as.integer(sqrt(diff(x)^2 + diff(y)^2) / tolerance), 2L)
  n_dense <- sum(counts) - length(counts) + 1L
  dense_x <- dense_y <- numeric(n_dense)
  cursor <- 1L
  for (i in seq_along(counts)) {
    sx <- seq(x[i], x[i + 1L], length.out = counts[i])
    sy <- seq(y[i], y[i + 1L], length.out = counts[i])
    if (i > 1L) {
      sx <- sx[-1L]
      sy <- sy[-1L]
    }
    take <- cursor + seq_along(sx) - 1L
    dense_x[take] <- sx
    dense_y[take] <- sy
    cursor <- cursor + length(sx)
  }
  new_x <- x[1]
  new_y <- y[1]
  while (length(dense_x) > 0L) {
    distances <- sqrt((dense_x - tail(new_x, 1))^2 + (dense_y - tail(new_y, 1))^2)
    j <- which.min(abs(distances - step_length))
    if (distances[j] < step_length - tolerance) break
    new_x <- c(new_x, dense_x[j])
    new_y <- c(new_y, dense_y[j])
    keep <- seq.int(j, length(dense_x))
    dense_x <- dense_x[keep]
    dense_y <- dense_y[keep]
    if (j == 1L) break
  }
  data.frame(X = new_x, Y = new_y)
}

