


mad <- function(x) {
    median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
}



glmnetReplicate <- function(model, reps = 1000, type = c("lambda.min",
    "lambda.1se")) {
    
# Checks the class of the object
    if (!class(model) == "cv.glmnet") {
        stop("'model' must be an output object of cv.glmnet!\n")
    }
    
# Checks the type
    if (! type[[1]] %in% c("lambda.min", "lambda.1se")) {
        stop("'type' must either be lambda.min or lambda.1se!\n")
    }

# Extracts the model call value
    mod_call <- model$call
    
# Replicates the model n times
    mod_reps <- replicate(reps, {
        mod_x <- eval(mod_call)
        list("coefs" = coef(mod_x, s = mod_x$lambda.min),
            "per_dev" = mod_x$glmnet.fit$dev.ratio[which(mod_x$glmnet.fit$lambda == mod_x$lambda.min)])        
        
    }, simplify = FALSE)

# Extracts the percent deviance and coefficient values
    per_devs <- sapply(mod_reps, function(x) x$per_dev)
    coefs <- sapply(mod_reps, function(x) x$coefs[, 1])
    
    
# Output
    mod_out <- list(per_dev = per_devs, coefs = coefs, call = mod_call,
        reps = reps, type = type[[1]])
    class(mod_out) <- "cv.glmrep"
    mod_out

}



print.cv.glmrep <- function(x) {
    
cat("\n\nReplicated glmnet() call: ")
cat("\nNumber of replications: ")
cat(x$reps)
cat("\nModel selection type: ")
cat(x$type)
cat("\nModel median % deviance (+- 1MAD): ")
cat(paste0(round(median(x$per_dev * 100), 2), " (+-",
    round(mad(x$per_dev * 100), 8), ")"))
cat("\n\nMedian parameter estimates (95% CI interval): \n")

mod_95CIs <- round(t(apply(x$coefs, 1, quantile, c(0.025, 0.5, 0.975))), 5)[, c(2, 1, 3)]
colnames(mod_95CIs) <- c("Estimate", "2.5%", "97.5%")

print(mod_95CIs[order(mod_95CIs[, 1]), ])


}






plot.cv.glmrep <- function(x, ordered = TRUE) {
    
# Stores median and 95% CI values
    mod_95CIs <- round(t(apply(x$coefs, 1, quantile, c(0.025, 0.5, 0.975))), 5)[, c(2, 1, 3)]
    colnames(mod_95CIs) <- c("Estimate", "2.5%", "97.5%")
    
    if (ordered == TRUE) {
        mod_95CIs <- mod_95CIs[order(mod_95CIs[, 1]), ]
    }

# Stores a good range for xlim values
    x_max_mag <- max(abs(mod_95CIs))
    x_range <- c(-x_max_mag * 1.1, x_max_mag * 1.1)
    
# Stores the ylim values
    y_vals <- nrow(mod_95CIs):1
    y_range <- c(nrow(mod_95CIs), 1)
    
# Takes a guess at the best Y axis margin size
    margin_2 <- max(nchar(rownames(x$coefs))) / 1.5
    
# Sets the framing parameters
    par(mar = c(4.5, margin_2, 1, 1))
    
# Opens an empty plot
    plot(NULL, xlim = x_range, ylim = y_range,
        xlab = "Effect Size (Median + 95% CI)", ylab = "", yaxt = "n")
    
# Adds the Y axis
    axis(side = 2, at = y_vals, labels = rownames(mod_95CIs), las = 2,
        cex = 0.8)
    mtext("Variable", side = 2, line = margin_2 * 0.85)
    
# Adds gridlines
    abline(h = y_vals, lwd = 0.75, col = "grey95")
    abline(v = seq(-1, 1, by = 0.05), lwd = 0.75, col = "grey95")
    
# Adds a zero line
    abline(v = 0, lty = 3)
    
# Adds the 95% CI intervals
    segments(mod_95CIs[, 2], y_vals, mod_95CIs[, 3], y_vals, lwd = 5,
        col = "grey80", lend = 1)
    
# Adds the median points
    points(mod_95CIs[, 1], y_vals, pch = 18)
 
# Bounding box
    box()
       
}





