compute_power_type1error <- function(pvalues, alpha = 0.05) {
    return((sum(pvalues <= 0.05))/length(pvalues))
}

extract_statistics <- function(res_simul, method, imputation_methods, statistics) {

    stats_data <- list()
    for(i in 1:length(imputation_methods)) {
        stats_data[[i]] <- lapply(
            res_simul,
            function(x) x[[2]][[method]][[imputation_methods[i]]]$pars$trt_6[statistics]
        )
        stats_data[[i]] <- as.data.frame(t(sapply(
            stats_data[[i]],
            function(x) unlist(x)
        )))
    }
    names(stats_data) <- imputation_methods

    return(stats_data)

}

summary_stats <- function(H0, N_sim = NULL, read_eachsim = FALSE) {

    if(read_eachsim) {
        res_simul <- list()
        for(i in 1:N_sim) {
            res_simul[[i]] <- readRDS(paste0("Results/Res_each_sim/res_",i,"_",H0,"_.rds"))
        }
    } else {
        res_simul <- readRDS(paste0("Results/results_", H0, ".rds"))
    }

    methods <- c("bootstrap", "jackknife", "bayesian")
    imputation_methods <- c("res_MAR", "res_JR", "res_CR", "res_CIR")
    statistics <- c("est", "se", "pvalue")

    stats_data <- sapply(
        methods,
        function(x) extract_statistics(res_simul, x, imputation_methods, statistics),
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    summary_stats <- sapply(
        stats_data,
        function(x) sapply(
            x,
            function(y) {

                ret_obj <- list("mean_est" = mean(y$est),
                                "sd_est" = sd(y$est),
                                "mean_se" = mean(y$se),
                                "type1error" = compute_power_type1error(y$pvalue)
                )

                # adjust name of type1error/power according to H0
                if(!H0) names(ret_obj)[which(names(ret_obj) == "type1error")] = "power"

                return(ret_obj)
            },
            simplify = FALSE,
            USE.NAMES = TRUE),
        simplify = FALSE,
        USE.NAMES = TRUE)

}

H0 <- c(TRUE, FALSE)
N_sim <- 10000

res <- sapply(H0,
              function(x) summary_stats(x),
              simplify = FALSE,
              USE.NAMES = TRUE)

saveRDS(res, file = "Results/table_res.rds")
