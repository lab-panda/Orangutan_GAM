library(gamm4)
library(ggplot2)

# Time varying associations, as per Fig 3,5 ,8
plot_tv_gam <- function(outcome, data = use) {
  # Fit the GAMM with the chosen outcome
  form <- as.formula(paste0(outcome, ' ~ s(year, bs = "cr", k = 6)'))
  mod  <- gamm4(form, data = data, random = ~ (1 | Animal_ID))
  
  # Extract the smooth, as in the original code
  fit_se <- plot(mod$gam)
  fplot <- data.frame(cbind(unlist(fit_se[[1]][1]),
                            unlist(fit_se[[1]][10]),
                            unlist(fit_se[[1]][3])))
  names(fplot) <- c("TSB", "Beta", "SE")
  
  # yields pointwise alpha ~ 0.0001, a conservative interval given the number of models fitted.
  fplot$upper <- fplot$Beta + (1.96 * fplot$SE)
  fplot$lower <- fplot$Beta - (1.96 * fplot$SE)
  
  # Title: outcome name with trailing isotope mass removed (e.g. "Ba138" -> "Ba")
  title <- sub("[0-9]+$", "", outcome)
  
  ggplot(data = fplot, aes(x = TSB, y = Beta)) +
    geom_line(size = 1.5) +
    theme_bw() +
    geom_hline(yintercept = 0, size = 0.6, lty = 2) +
    xlab("Time Since Birth (Years)") +
    ylab("Time-Varying Concentration (SD)") +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 2, alpha = 0.35) +
    theme(text = element_text(size = 14)) +
    ggtitle(title)
}


# GAM to estimate seasonal effects in urine
fit_seasonal_gam <- function(outcome, data = urine) {
  
  formula <- as.formula(paste0(outcome, " ~ s(doy, bs = 'cc', k = 4) + s(Day, bs = 'cr')"))
  
  m <- gamm(formula,
            data = data,
            correlation = corCAR1(form = ~Day))
  
  # --- Plot 1: Raw data + smooth full model fit ---
  pred_full <- data.frame(
    Day = seq(min(data$Day), max(data$Day), length.out = 200),
    doy = seq(min(data$doy), max(data$doy), length.out = 200)
  )
  pred_f <- predict(m$gam, newdata = pred_full, se.fit = TRUE)
  pred_full$fit <- pred_f$fit
  pred_full$lwr <- pred_f$fit - 1.96 * pred_f$se.fit
  pred_full$upr <- pred_f$fit + 1.96 * pred_f$se.fit
  pred_full$Date_Collected <- seq(min(data$Date_Collected), max(data$Date_Collected), length.out = 200)
  
  p_raw <- ggplot() +
    geom_ribbon(data = pred_full, aes(x = Date_Collected, ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(data = pred_full, aes(x = Date_Collected, y = fit)) +
    geom_point(data = data, aes(x = Date_Collected, y = .data[[outcome]]), colour = "steelblue") +
    labs(x = "Date collected", y = outcome, title = "Observed and fitted values") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  # --- Plot 2: Seasonal partial effect ---
  pred_doy <- data.frame(doy = 1:365, Day = mean(data$Day))
  pred <- predict(m$gam, newdata = pred_doy, type = "terms", se.fit = TRUE)
  pred_doy$fit <- pred$fit[, "s(doy)"]
  pred_doy$lwr <- pred_doy$fit - 1.96 * pred$se.fit[, "s(doy)"]
  pred_doy$upr <- pred_doy$fit + 1.96 * pred$se.fit[, "s(doy)"]
  
  p_doy <- ggplot(pred_doy, aes(x = doy)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = fit)) +
    labs(x = "Day of year", y = "Partial effect", title = "Seasonal effect") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  # --- Plot 3: Annual partial effect ---
  pred_day <- data.frame(Day = seq(min(data$Day), max(data$Day), length.out = 200),
                         doy = mean(data$doy))
  pred2 <- predict(m$gam, newdata = pred_day, type = "terms", se.fit = TRUE)
  pred_day$fit <- pred2$fit[, "s(Day)"]
  pred_day$lwr <- pred_day$fit - 1.96 * pred2$se.fit[, "s(Day)"]
  pred_day$upr <- pred_day$fit + 1.96 * pred2$se.fit[, "s(Day)"]
  
  p_day <- ggplot(pred_day, aes(x = Day)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = fit)) +
    labs(x = "Sequential day", y = "Partial effect", title = "Annual trend") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  print(summary(m$gam))
  print(p_raw / (p_doy + p_day) + plot_annotation(title = outcome))
  
  invisible(m)
}

# GAM to estimate time-varying associations in urine, adjusting for seasonal and annual trends
fit_varying_gam <- function(outcome, covariate, data = urine) {
  
  formula <- as.formula(paste0(outcome, " ~ s(doy, bs = 'cc', k = 4) + s(doy, by = ", covariate, ", bs = 'cc', k = 4) + s(Day, bs = 'cr')"))
  
  m <- gamm(formula,
            data = data,
            correlation = corCAR1(form = ~Day))
  
  # --- Plot 1: Time-varying association between covariate and outcome ---
  pred_vc <- data.frame(doy = 1:365, Day = mean(data$Day))
  pred_vc[[covariate]] <- 1
  
  pred_v <- predict(m$gam, newdata = pred_vc, type = "terms", se.fit = TRUE)
  term_name <- paste0("s(doy):", covariate)
  pred_vc$fit <- pred_v$fit[, term_name]
  pred_vc$lwr <- pred_vc$fit - 1.96 * pred_v$se.fit[, term_name]
  pred_vc$upr <- pred_vc$fit + 1.96 * pred_v$se.fit[, term_name]
  
  p_vc <- ggplot(pred_vc, aes(x = doy)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = fit)) +
    labs(x = "Day of year", y = paste0("Effect of ", covariate, " on ", outcome),
         title = paste0("Time-varying association: ", covariate, " → ", outcome)) +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  # --- Plot 2: Seasonal partial effect (main effect of doy) ---
  pred_doy <- data.frame(doy = 1:365, Day = mean(data$Day))
  pred_doy[[covariate]] <- mean(data[[covariate]], na.rm = TRUE)
  
  pred <- predict(m$gam, newdata = pred_doy, type = "terms", se.fit = TRUE)
  pred_doy$fit <- pred$fit[, "s(doy)"]
  pred_doy$lwr <- pred_doy$fit - 1.96 * pred$se.fit[, "s(doy)"]
  pred_doy$upr <- pred_doy$fit + 1.96 * pred$se.fit[, "s(doy)"]
  
  p_doy <- ggplot(pred_doy, aes(x = doy)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = fit)) +
    labs(x = "Day of year", y = "Partial effect", title = "Seasonal effect") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  # --- Plot 3: Annual partial effect ---
  pred_day <- data.frame(Day = seq(min(data$Day), max(data$Day), length.out = 200),
                         doy = mean(data$doy))
  pred_day[[covariate]] <- mean(data[[covariate]], na.rm = TRUE)
  
  pred2 <- predict(m$gam, newdata = pred_day, type = "terms", se.fit = TRUE)
  pred_day$fit <- pred2$fit[, "s(Day)"]
  pred_day$lwr <- pred_day$fit - 1.96 * pred2$se.fit[, "s(Day)"]
  pred_day$upr <- pred_day$fit + 1.96 * pred2$se.fit[, "s(Day)"]
  
  p_day <- ggplot(pred_day, aes(x = Day)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(aes(y = fit)) +
    labs(x = "Sequential day", y = "Partial effect", title = "Annual trend") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  print(summary(m$gam))
  print(p_vc / (p_doy + p_day) + plot_annotation(title = paste0(outcome, " ~ ", covariate)))
  
  invisible(m)
}