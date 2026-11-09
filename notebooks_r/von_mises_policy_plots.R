# von Mises navigation policy: figure construction.
#
# Sourced by notebooks_r/von_mises_policy.ipynb (no standalone
# output). R twin of notebooks/von_mises_policy_plots.py: the
# two-row navigation figure (exact density as a radial modulation of
# the policy plane on top, the Monte-Carlo trajectories below). The
# notebook defines the policy and runs the simulation; only the figure
# construction lives here.

navigation_panels <- function(x_start, target, kappas, theta_grid, radials,
  trajectories, T_steps, num_traj, mu0, R_policy, amp) {
  policy_panels <- trajectory_panels <- vector("list", length(kappas))
  # Policy panels follow the Python twin (von_mises_policy_plots.py,
  # plot_navigation): the density curve is embedded in the policy plane as
  # x_start + radius * cos/sin(theta) - the same plane as the trajectory
  # panels below - with the policy circle (radius R_policy around the
  # start), a dot at the start and the goal-direction marker drawn
  # directly, instead of a separate polar panel.
  Rmax_policy <- R_policy + amp
  ref_ring <- data.frame(
    x = x_start[1] + R_policy * cos(seq(0, 2 * pi, length.out = 361L)),
    y = x_start[2] + R_policy * sin(seq(0, 2 * pi, length.out = 361L)))
  for (i in seq_along(kappas)) {
    polar <- data.frame(x = x_start[1] + radials[[i]] * cos(theta_grid),
      y = x_start[2] + radials[[i]] * sin(theta_grid))
    policy_panels[[i]] <- ggplot(polar, aes(x, y)) +
      geom_path(aes(group = 1), colour = "black", linewidth = 0.6) +
      geom_path(data = ref_ring, inherit.aes = FALSE, aes(x, y),
        colour = "black", alpha = 0.95, linewidth = 0.3) +
      annotate("point", x = x_start[1], y = x_start[2], colour = "black", size = 2) +
      annotate("point",
        x = x_start[1] + (R_policy + amp + 0.22) * cos(mu0),
        y = x_start[2] + (R_policy + amp + 0.22) * sin(mu0),
        colour = "black", shape = 17, size = 2.5) +
      scale_x_continuous(limits = x_start[1] + c(-1, 1) * (Rmax_policy + 0.3),
        expand = c(0, 0)) +
      scale_y_continuous(limits = x_start[2] + c(-1, 1) * (Rmax_policy + 0.3),
        expand = c(0, 0)) +
      coord_equal() +
      labs(title = paste("kappa =", kappas[i]), x = "", y = "Policy density") +
      theme_void(base_size = 12)
    pos <- trajectories[[i]]
    seg <- do.call(rbind, lapply(seq_len(num_traj), function(j) {
      data.frame(iteration = seq_len(T_steps), trajectory = factor(j),
        x = pos[-(T_steps + 1L), j, 1], y = pos[-(T_steps + 1L), j, 2],
        xend = pos[-1L, j, 1], yend = pos[-1L, j, 2])
    }))
    # Segments encode consecutive time points explicitly: no sorting of an
    # x/y path can change a trajectory's temporal order.
    trajectory_panels[[i]] <- ggplot(seg, aes(x, y, colour = trajectory)) +
      geom_segment(aes(xend = xend, yend = yend), linewidth = 0.3, alpha = 0.9) +
      annotate("point", x = x_start[1], y = x_start[2], colour = "black", size = 2) +
      annotate("point", x = target[1], y = target[2], colour = "black", shape = 17, size = 2.5) +
      scale_colour_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")) +
      coord_equal(xlim = c(-7.5, 7.5), ylim = c(-3.5, 10.8)) +
      theme_void(base_size = 12) + theme(legend.position = "none")
  }
  list(policy = policy_panels, trajectory = trajectory_panels)
}
