library(readr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(DT)
library(tidyr)
library(gtools)
library(knitr)
rm(list = ls())

hms_span <- function(start, end) {
    dsec <- as.numeric(difftime(end, start, unit = "secs"))
    hours <- floor(dsec / 3600)
    minutes <- floor((dsec - 3600 * hours) / 60)
    seconds <- dsec - 3600 * hours - 60 * minutes
    paste0(sapply(c(hours, minutes, seconds), function(x) {
        formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

# parameters

beta.1 <- log(1.2)
beta.2 <- log(2)
beta.3 <- 0
MAF <- 0.2
f2 <- 0.1
para_prev <- 0.01
N_cases_expected <- 10000
n_simu <- 10000
output_folder <- "./"

# calculate some parameters
N_all <- round(N_cases_expected / para_prev)
f1 <- MAF * MAF + 2 * MAF * (1 - MAF)

my_fun <- function(x) { (1-f1)*(1-f2)*exp(x)/(1+exp(x)) + f1*(1-f2)*(exp(x+beta.1)/(1+exp(x+beta.1))) + (1-f1)*f2*(exp(x+beta.2)/(1+exp(x+beta.2))) + f1*f2*(exp(x+beta.1+beta.2+beta.3)/(1+exp(x+beta.1+beta.2+beta.3))) - para_prev}
estimate_beta0 <- uniroot(f = my_fun, interval = c(-10, 0))$root
group_00 = (1-f1)*(1-f2)*exp(estimate_beta0)/(1+exp(estimate_beta0))
group_10 = f1*(1-f2)*(exp(estimate_beta0+beta.1)/(1+exp(estimate_beta0+beta.1))) 
group_01 = (1-f1)*f2*(exp(estimate_beta0+beta.2)/(1+exp(estimate_beta0+beta.2)))
group_11 = f1*f2*(exp(estimate_beta0+beta.1+beta.2+beta.3)/(1+exp(estimate_beta0+beta.1+beta.2+beta.3)))
estimated_prevalence = group_00+group_10+group_01+group_11


output_path <- output_folder
simulation_result_prefix <- paste0("simu_b1_", format(beta.1, nsmall = 3, digits = 3), "_b2_", format(beta.2, nsmall = 3, digits = 3), "_b3_", format(beta.3, nsmall = 3, digits = 3), "_MAF_", format(MAF, nsmall = 3, digits = 3), "_freq.E_", format(f2, nsmall = 3, digits = 3), "_prev_", format(para_prev, nsmall = 3, digits = 3),"_ncases_", as.character(N_cases_expected))


######################## SIMULATIONS ################################
time_1 <- Sys.time()
print(paste("Now we start at", time_1))

print("This is case-control study")
print("[====parameter setting====]")
print(paste0("b1=", as.character(beta.1)))
print(paste0("b2=", as.character(beta.2)))
print(paste0("b3=", as.character(beta.3)))
print(paste0("MAF=", as.character(MAF)))
print(paste0("freq.E=", as.character(f2)))
print(paste0("disease prevalence=", as.character(para_prev)))
print(paste0("estimated prevalence (based on estimate beta0)=", as.character(estimated_prevalence)))
print(paste0("estimated beta=", as.character(exp(estimate_beta0)/(1+exp(estimate_beta0)))))
print(paste0("baseline risk=", as.character(exp(estimate_beta0)/(1+exp(estimate_beta0)))))
print(paste0("N_cases_expected=", as.character(N_cases_expected)))
print(paste0("n_simu=", as.character(n_simu)))

output_file_simulation_table <- paste0(output_path, "/", simulation_result_prefix, ".csv")

write("true_beta_1, true_beta_2, true_beta_3, para_MAF, para_freq.x2, para_prevalence, simulated_prevalence, estimate_beta0, N_case_expected, N_cases_simulated, N_all_simulated, n_x1_0_x2_0_case, n_x1_1_x2_0_case, n_x1_0_x2_1_case, n_x1_1_x2_1_case, n_x1_0_x2_0_ctrl, n_x1_1_x2_0_ctrl, n_x1_0_x2_1_ctrl, n_x1_1_x2_1_ctrl, cohort_beta_0, cohort_beta_1, cohort_beta_2, cohort_beta_3, cohort_p_beta_3, case_alpha_0, case_alpha_1, case_p_alpha_1,  asct_beta_0_ratio_1, asct_beta_1_ratio_1, asct_beta_2_ratio_1, asct_beta_3_ratio_1, asct_p_gamma_ratio_1,  asct_beta_0_ratio_2, asct_beta_1_ratio_2, asct_beta_2_ratio_2, asct_beta_3_ratio_2, asct_p_gamma_ratio_2, asct_beta_0_ratio_3, asct_beta_1_ratio_3, asct_beta_2_ratio_3, asct_beta_3_ratio_3, asct_p_gamma_ratio_3,  asct_alpha_0_ratio_1, asct_alpha_1_ratio_1, asct_p_alpha_ratio_1, asct_alpha_0_ratio_2, asct_alpha_1_ratio_2, asct_p_alpha_ratio_2, asct_alpha_0_ratio_3, asct_alpha_1_ratio_3, asct_p_alpha_ratio_3", output_file_simulation_table, append = F)

i <- 0
while (i < n_simu) {
    # generation
    x1_a <- rbinom(N_all, size = 1, prob = MAF)
    x1_b <- rbinom(N_all, size = 1, prob = MAF)
    x1 <- as.numeric(x1_a + x1_b >= 1)
    x2 <- rbinom(N_all, size = 1, prob = f2)
    logit.p.case <- estimate_beta0 + beta.1 * x1 + beta.2 * x2 + beta.3 * x1 * x2
    p.case <- exp(logit.p.case) / (1 + exp(logit.p.case))
    y <- rbinom(N_all, 1, p.case)
    simulate_data <- data.frame(cbind(y, x1, x2))
    simulated_prevalence <- as.numeric(table(y)[2] / N_all)
    simulate_data <- simulate_data %>%
        mutate(
            x1 = factor(x1, levels = c("0", "1")),
            x2 = factor(x2, levels = c("0", "1"))
        )
    # select cases and controls
    case_data <- simulate_data %>%
        filter(y == 1)

    case_stat <- table(case_data$x1, case_data$x2)
    n_x1_0_x2_0_case <- case_stat[1, 1]
    n_x1_1_x2_0_case <- case_stat[2, 1]
    n_x1_0_x2_1_case <- case_stat[1, 2]
    n_x1_1_x2_1_case <- case_stat[2, 2]
    ctrl_data <- simulate_data %>%
        filter(y == 0)
    ctrl_stat <- table(ctrl_data$x1, ctrl_data$x2)
    n_x1_0_x2_0_ctrl <- ctrl_stat[1, 1]
    n_x1_1_x2_0_ctrl <- ctrl_stat[2, 1]
    n_x1_0_x2_1_ctrl <- ctrl_stat[1, 2]
    n_x1_1_x2_1_ctrl <- ctrl_stat[2, 2]
    if (nrow(ctrl_data) >= nrow(case_data) * 3 & length(unique(case_data$x1)) != 1 & length(unique(case_data$x2)) != 1) {
        if (table(simulate_data$x1, simulate_data$x2)[2, 2] != 0) {
            # cohort data
            cohort.regress.y.x1.x2.x1x2 <- glm(y ~ x1 + x2 + x1 * x2, data = simulate_data, family = binomial)
            cohort_beta_0 <- cohort.regress.y.x1.x2.x1x2$coefficients[1]
            cohort_beta_1 <- cohort.regress.y.x1.x2.x1x2$coefficients[2]
            cohort_beta_2 <- cohort.regress.y.x1.x2.x1x2$coefficients[3]
            cohort_beta_3 <- cohort.regress.y.x1.x2.x1x2$coefficients[4]
            cohort_p_beta_3 <- coef(summary(cohort.regress.y.x1.x2.x1x2))[4, 4]
        } else {
            cohort_beta_0 <- NA
            cohort_beta_1 <- NA
            cohort_beta_2 <- NA
            cohort_beta_3 <- NA
            cohort_p_beta_3 <- NA
        }
        # case-only data
        case.regress.x1.x2 <- glm(x1 ~ x2, data = case_data, family = binomial)
        case_alpha_0 <- case.regress.x1.x2$coefficients[1]
        case_alpha_1 <- case.regress.x1.x2$coefficients[2]
        case_p_alpha_1 <- coef(summary(case.regress.x1.x2))[2, 4]
        # case-control data
        asct_beta_0 <- asct_beta_1 <- asct_beta_2 <- asct_beta_3 <- asct_p_gamma <- c()
        asct_alpha_0 <- asct_alpha_1 <- asct_p_alpha <- c()
        for (ratio in 1:3) {
            sig <- 1
            while (sig == 1) {
                asct_case_data <- case_data
                asct_ctrl_data <- ctrl_data[sample(nrow(ctrl_data), nrow(asct_case_data) * ratio), ]
                asct_data <- rbind(asct_case_data, asct_ctrl_data)
                if (length(unique(asct_case_data$x1)) != 1 & length(unique(asct_case_data$x2)) != 1) {
                    sig <- 0
                    # regression using case-control data
                    if (table(asct_data$x1, asct_data$x2)[2, 2] != 0) {
                        asct.regress.y.x1.x2.x1x2 <- glm(y ~ x1 + x2 + x1 * x2, data = asct_data, family = binomial)
                        asct_beta_0 <- c(asct_beta_0, asct.regress.y.x1.x2.x1x2$coefficients[1])
                        asct_beta_1 <- c(asct_beta_1, asct.regress.y.x1.x2.x1x2$coefficients[2])
                        asct_beta_2 <- c(asct_beta_2, asct.regress.y.x1.x2.x1x2$coefficients[3])
                        asct_beta_3 <- c(asct_beta_3, asct.regress.y.x1.x2.x1x2$coefficients[4])
                        asct_p_gamma <- c(asct_p_gamma, coef(summary(asct.regress.y.x1.x2.x1x2))[4, 4])
                    } else {
                        asct_beta_0 <- c(asct_beta_0, NA)
                        asct_beta_1 <- c(asct_beta_1, NA)
                        asct_beta_2 <- c(asct_beta_2, NA)
                        asct_beta_3 <- c(asct_beta_3, NA)
                        asct_p_gamma <- c(asct_p_gamma, 1)
                    }
                    # regression using case-only data
                    asct.regress.x1.x2 <- glm(x1 ~ x2, data = asct_case_data, family = binomial)
                    asct_alpha_0 <- c(asct_alpha_0, asct.regress.x1.x2$coefficients[1])
                    asct_alpha_1 <- c(asct_alpha_1, asct.regress.x1.x2$coefficients[2])
                    asct_p_alpha <- c(asct_p_alpha, coef(summary(asct.regress.x1.x2))[2, 4])
                }
            }
        }
        # write to output
        output_line <- paste(
            beta.1, beta.2, beta.3,
            MAF, f2, para_prev, simulated_prevalence, estimate_beta0,
            N_cases_expected, nrow(case_data), N_all,
            n_x1_0_x2_0_case, n_x1_1_x2_0_case, n_x1_0_x2_1_case, n_x1_1_x2_1_case,
            n_x1_0_x2_0_ctrl, n_x1_1_x2_0_ctrl, n_x1_0_x2_1_ctrl, n_x1_1_x2_1_ctrl,
            cohort_beta_0, cohort_beta_1, cohort_beta_2, cohort_beta_3, cohort_p_beta_3,
            case_alpha_0, case_alpha_1, case_p_alpha_1,
            asct_beta_0[1], asct_beta_1[1], asct_beta_2[1], asct_beta_3[1], asct_p_gamma[1],
            asct_beta_0[2], asct_beta_1[2], asct_beta_2[2], asct_beta_3[2], asct_p_gamma[2],
            asct_beta_0[3], asct_beta_1[3], asct_beta_2[3], asct_beta_3[3], asct_p_gamma[3],
            asct_alpha_0[1], asct_alpha_1[1], asct_p_alpha[1],
            asct_alpha_0[2], asct_alpha_1[2], asct_p_alpha[2],
            asct_alpha_0[3], asct_alpha_1[3], asct_p_alpha[3],
            sep = ","
        )
        write(output_line, output_file_simulation_table, append = TRUE)
        i <- i + 1
    }
}

time_2 <- Sys.time()
print(paste("Now we end at", time_2, " (the whole process takes", hms_span(time_1, time_2), ")."))


######################## EVALUATION ################################

############### functions required for evaluation ###############
import_data <- function(data_path, name_pattern = "simu_*.csv", test_code = FALSE, print = TRUE) {
    # import simulated data
    simu_res_filelist <- list.files(path = data_path, pattern = glob2rx(name_pattern), full.names = TRUE)
    if (print) {
        print(paste0("There are all together ", as.character(length(simu_res_filelist)), " files under ", data_path))
    }
    time_1 <- Sys.time()
    # significance_level <- 0.05
    simu_res_df <- data.frame()
    if (test_code) {
        for (simulation_result_file in simu_res_filelist) {
            simulation_result_tmp <- read.csv(file = simulation_result_file, nrows = 10)
            simu_res_df <- rbind(simu_res_df, simulation_result_tmp)
        }
    } else {
        for (simulation_result_file in simu_res_filelist) {
            simulation_result_tmp <- read.csv(file = simulation_result_file)
            simu_res_df <- rbind(simu_res_df, simulation_result_tmp)
        }
    }
    time_2 <- Sys.time()
    if (print){
        print(paste("Importing", as.character(length(simu_res_filelist)), "csv files takes", hms_span(time_1, time_2), "."))
    }
    simu_res_df <- simu_res_df %>%
        mutate(para_beta = paste0("(", format(true_beta_1, digits = 3, nsmall = 3), ",", format(true_beta_2, digits = 3, nsmall = 3), ",", format(true_beta_3, digits = 3, nsmall = 3), ")"))
    if (print){
        print("============ Number of replicates ============")
        print(table(simu_res_df$para_beta, simu_res_df$para_prev))
    }
    return(simu_res_df)
}

format_data <- function(df) {
    simu_res_df <- df %>%
        mutate(
            true_beta_1 = ifelse(abs(true_beta_1 - 0.182321556793955) <= 0.01, log(1.2), true_beta_1), # this is for accuracy in numerical results
            true_beta_2 = ifelse(abs(true_beta_2 - 0.693147180559945) <= 0.01, log(2), true_beta_2), # this is for accuracy in numerical results
            true_beta_3 = ifelse(abs(true_beta_3 - 0.182321556793955) <= 0.01, log(1.2), true_beta_3) # this is for accuracy in numerical results
        ) %>%
        select(
            true_beta_1, true_beta_2, true_beta_3, para_prevalence,
            case_alpha_1, asct_beta_3_ratio_1, asct_beta_3_ratio_2, asct_beta_3_ratio_3, # this is the beta
            case_p_alpha_1, asct_p_gamma_ratio_1, asct_p_gamma_ratio_2, asct_p_gamma_ratio_3 # this is the p-value
        ) %>%
        rename(
            est_beta_case_only = case_alpha_1,
            est_beta_case_control_R1 = asct_beta_3_ratio_1,
            est_beta_case_control_R2 = asct_beta_3_ratio_2,
            est_beta_case_control_R3 = asct_beta_3_ratio_3,
            pvalue_case_only = case_p_alpha_1,
            pvalue_case_control_R1 = asct_p_gamma_ratio_1,
            pvalue_case_control_R2 = asct_p_gamma_ratio_2,
            pvalue_case_control_R3 = asct_p_gamma_ratio_3
        )
    simu_res_df_reshaped <- simu_res_df %>%
        pivot_longer(cols = starts_with("est_beta_"), names_to = "design", values_to = "beta") %>%
        pivot_longer(cols = starts_with("pvalue_"), names_to = "pvalue_design", values_to = "pvalue") %>%
        mutate(
            design = gsub("est_beta_", "", design),
            pvalue_design = gsub("pvalue_", "", pvalue_design)
        ) %>%
        filter(design == pvalue_design) %>%
        select(true_beta_1, true_beta_2, true_beta_3, para_prevalence, design, beta, pvalue)
    simu_res_df_reshaped <- simu_res_df_reshaped %>%
        mutate(Design = case_when(
            design == "case_only" ~ "\ncase-only\nN = 10 000 cases\n",
            design == "case_control_R1" ~ "\ncase-control\nN = 10 000 cases    \n + 10 000 controls\n",
            design == "case_control_R2" ~ "\ncase-control\nN = 10 000 cases    \n + 20 000 controls\n",
            design == "case_control_R3" ~ "\ncase-control\nN = 10 000 cases    \n + 30 000 controls\n"
        )) %>%
        mutate(
            OR_G = exp(true_beta_1),
            OR_E = exp(true_beta_2),
            OR_GE = exp(true_beta_3),
            Scenario = paste0(
                "OR[G]==", format(OR_G, nsmall = 3, digits = 3),
                "~ and~ OR[E]==", format(OR_E, nsmall = 3, digits = 3),
                "~ and~ OR[GE]==", format(OR_GE, nsmall = 3, digits = 3)
            )
        )
    simu_res_df_reshaped$Design <- factor(simu_res_df_reshaped$Design,
        levels = c(
            "\ncase-only\nN = 10 000 cases\n",
            "\ncase-control\nN = 10 000 cases    \n + 10 000 controls\n",
            "\ncase-control\nN = 10 000 cases    \n + 20 000 controls\n",
            "\ncase-control\nN = 10 000 cases    \n + 30 000 controls\n"
        )
    )
    return(simu_res_df_reshaped)
}



calculate_type_I_err <- function(df, significance_level_threshold = 0.05) {
    type_i_err_res <- df %>%
        mutate(
            significance_level = significance_level_threshold,
            reject_sign = ifelse(pvalue < significance_level, 1, 0)
        ) %>%
        group_by(true_beta_1, true_beta_2, true_beta_3, para_prevalence, significance_level, design) %>%
        summarise(
            count = n(),
            type_i_err = mean(reject_sign),
            sem = sqrt(type_i_err * (1 - type_i_err) / count)
        )
    return(type_i_err_res)
}

calculate_stat_power <- function(df, significance_level_threshold = alpha) {
    stat_power_res <- df %>%
        mutate(
            significance_level = significance_level_threshold,
            reject_sign = ifelse(pvalue < significance_level, 1, 0)
        ) %>%
        group_by(true_beta_1, true_beta_2, true_beta_3, para_prevalence, significance_level, Design) %>%
        summarise(
            count = n(),
            stat_power = mean(reject_sign),
            sem = sqrt(stat_power * (1 - stat_power) / count)
        )
    return(stat_power_res)
}

draw_stat_power_figure <- function(df, output_path = "./", output_name) {
    options(repr.plot.width = 30, repr.plot.height = 20, repr.plot.res = 100)
    p <- ggplot(df, aes(x =  factor(paste0(format(para_prevalence * 100, digits = 3, nsmall = 1), "%")), y = stat_power, group = Design, fill = Design)) +
        # geom_point(size = 10, alpha = 0.8, aes(color = estimator)) +
        geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
        geom_errorbar(aes(ymin = stat_power - sem, ymax = stat_power + sem), width = 0.3, position = position_dodge(0.7)) +
        geom_hline(aes(yintercept = 0.5), color = "red", linewidth = 1, linetype = "dashed") +
        geom_text(aes(label = format(stat_power, digits = 2, nsmall = 2)), position = position_dodge(width = 0.70), size = 5, vjust = -0.5) +
        # facet_wrap(~scenario, nrow = 2, ncol = 2, labeller = label_parsed, scales = "fixed") +
        labs(
            x = "Disease prevalence",
            y = "Statistical power",
            title = "Statistical power for the case-only and case-control designs"
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1.05)) +
        theme_bw() +
        theme(
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank(),
            plot.subtitle = element_text(size = 28),
            axis.title.x = element_text(size = 30, face = "bold"),
            axis.text.x = element_text(size = 20, face = "bold"),
            axis.title.y = element_text(size = 30, face = "bold"),
            axis.text.y = element_text(size = 20, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_blank(),
            strip.text = element_text(size = 25, face = "bold"),
            legend.position = "none"
            )
    print(p)
    ggsave(paste0(output_path, output_name, ".pdf"), width = 30, height = 20, dpi = 300)
    ggsave(paste0(output_path, output_name, ".png"), width = 30, height = 20, dpi = 300)
}


gg_qqplot_input <- function(ps, ci = 0.95) {
    n <- length(ps)
    df <- data.frame(
        observed = -log10(sort(ps)),
        expected = -log10(ppoints(n)),
        clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
}


draw_qq_plot <- function(df, output_path = "./", output_name){
    case_only_df <- df %>% filter(design == "case_only")
    case_control_R1_df <- df %>% filter(design == "case_control_R1")
    qq_input_1 <- gg_qqplot_input(case_only_df$pvalue)
    qq_input_2 <- gg_qqplot_input(case_control_R1_df$pvalue)
    qq_input_1$Design <- "\ncase-only\nN = 10,000 cases\n"
    qq_input_2$Design <- "\ncase-control\nN = 10 000 cases    \n + 10 000 controls\n"
    qq_input <- rbind(qq_input_1, qq_input_2)
    qq_input$Design <- factor(qq_input$Design, levels = c("\ncase-only\nN = 10,000 cases\n", "\ncase-control\nN = 10 000 cases    \n + 10 000 controls\n"))
    p <- ggplot() +
            geom_ribbon(
                data = qq_input_1,
                mapping = aes(x = expected, ymin = clower, ymax = cupper),
                alpha = 0.1
            ) +
            geom_point(data = qq_input, aes(expected, observed, color = Design), shape = 16, size = 2) +
            geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
            labs(
                x = expression(paste("Expected -log"[10], plain(P))),
                y = expression(paste("Observed -log"[10], plain(P)))
            ) +
            theme_bw() +
            theme(
                panel.background = element_blank(),
                panel.border = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 25, face = "bold"),
                plot.subtitle = element_text(size = 28),
                axis.title.x = element_text(size = 25, face = "bold"),
                axis.text.x = element_text(size = 25, face = "bold"),
                axis.title.y = element_text(size = 25, face = "bold"),
                axis.text.y = element_text(size = 25, face = "bold"),
                legend.title = element_blank(),
                legend.text = element_blank(),
                strip.text = element_text(size = 30, face = "bold"),
                legend.position = "none"
            ) +
            guides(color = guide_legend(override.aes = list(size = 8)), fill = guide_legend(byrow = TRUE))
        print(p)
    ggsave(paste0(output_path, output_name, ".pdf"), width = 30, height = 20, dpi = 300)
    ggsave(paste0(output_path, output_name, ".png"), width = 30, height = 20, dpi = 300)
    }

draw_violin_plot <- function(df, true_beta_GE, output_path = "./", output_name){
    cols <- c("#F8766D", "#4A536B")
    p <- ggplot(df, aes(x = factor(paste0(format(para_prevalence * 100, digits = 3, nsmall = 1), "%")), y = beta - true_beta_GE, fill = Design)) +
        geom_boxplot(width = 0.7, color = "black", alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = 2, color = "red") +
        labs(
        x = "Disease prevalence",
        y = "Bias"
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(-0.6, 0.6)) +
        theme_bw() +
        theme(
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 35, face = "bold"),
        plot.subtitle = element_text(size = 28),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.text.y = element_text(size = 25, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text = element_text(size = 30, face = "bold"),
        legend.position = "none"
        ) +
        guides(color = guide_legend(override.aes = list(size = 5)), fill = guide_legend(byrow = TRUE)) + 
        scale_fill_manual(values = cols)
    print(p)
    ggsave(paste0(output_path, output_name, ".pdf"), width = 30, height = 20, dpi = 300)
    ggsave(paste0(output_path, output_name, ".png"), width = 30, height = 20, dpi = 300)
    }


calculate_analytical_bias <- function(OR_G, OR_E, OR_GE, MAF, freq.E, output_path = "./", prefix = "Table"){
    beta.1 <- log(OR_G)
    beta.2 <- log(OR_E)
    beta.3 <- log(OR_GE)
    f1 <- MAF * MAF + 2 * MAF * (1 - MAF)
    f2 <- freq.E
    result <- data.frame()
    for (disease_prevalence in c(seq(from = 0.01, to = 0.05, by = 0.005), 0.10, 0.15, 0.20)) {
        my_fun <- function(x) {
            (1 - f1) * (1 - f2) * exp(x) / (1 + exp(x)) + f1 * (1 - f2) * (exp(x + beta.1) / (1 + exp(x + beta.1))) + (1 - f1) * f2 * (exp(x + beta.2) / (1 + exp(x + beta.2))) + f1 * f2 * (exp(x + beta.1 + beta.2 + beta.3) / (1 + exp(x + beta.1 + beta.2 + beta.3))) - disease_prevalence
        }
        estimated.beta.0 <- uniroot(f = my_fun, interval = c(-10, 0))$root
        baseline_risk <- 1 / (1 + exp(-estimated.beta.0))
        group_00 <- (1 - f1) * (1 - f2) * exp(estimated.beta.0) / (1 + exp(estimated.beta.0))
        group_10 <- f1 * (1 - f2) * (exp(estimated.beta.0 + beta.1) / (1 + exp(estimated.beta.0 + beta.1)))
        group_01 <- (1 - f1) * f2 * (exp(estimated.beta.0 + beta.2) / (1 + exp(estimated.beta.0 + beta.2)))
        group_11 <- f1 * f2 * (exp(estimated.beta.0 + beta.1 + beta.2 + beta.3) / (1 + exp(estimated.beta.0 + beta.1 + beta.2 + beta.3)))
        estimated_disease_prevalence <- group_00 + group_10 + group_01 + group_11
        numerator <- exp(2 * estimated.beta.0 + beta.1 + beta.2) + exp(estimated.beta.0 + beta.1) + exp(estimated.beta.0 + beta.2) + 1
        denominator <- exp(2 * estimated.beta.0 + beta.1 + beta.2 + beta.3) + exp(estimated.beta.0 + beta.1 + beta.2 + beta.3) + exp(estimated.beta.0) + 1
        bias <- log(numerator / denominator)
        result <- rbind(result, c(OR_G, beta.1, OR_E, beta.2, OR_GE, beta.3, MAF, f1, f2, disease_prevalence, estimated_disease_prevalence, estimated.beta.0, baseline_risk, bias))
    }
    colnames(result) <- c("OR_G", "beta.1", "OR_E", "beta.2", "OR_GE", "beta.3", "MAF", "f1", "freq.E", "expected_disease_prevalence", "estimated_disease_prevalence", "estimated.beta.0", "baseline_risk", "bias")

    write.table(result,
        file = paste0(
            output_path,
            prefix, "_bias_b1_", format(beta.1, nsmall = 3, digits = 3),
            "_b2_", format(beta.2, nsmall = 3, digits = 3),
            "_b3_", format(beta.3, nsmall = 3, digits = 3),
            "_MAF_", format(MAF, nsmall = 3, digits = 3),
            "_freq.E_", format(f2, nsmall = 3, digits = 3), ".txt"
        ),
        append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
    )
    display_table <- result %>% select(expected_disease_prevalence, baseline_risk, bias)
    custom_column_names <- c("Disease prevalence", as.character(c(seq(from = 0.01, to = 0.05, by = 0.005), 0.10, 0.15, 0.20)))
    knitr::kable(t(display_table), col.names = custom_column_names)
}


############### Evaluate Type I error ###############

for (alpha in c(0.001, 0.01, 0.05)){
    header = c("true_beta_1, true_beta_2, true_beta_3, para_prevalence, significance_level, design, count, type_i_err, sem")
    write("true_beta_1, true_beta_2, true_beta_3, para_prevalence, significance_level, design, count, type_i_err, sem",
          paste0(output_path, "Table_type_i_err_alpha_",format(alpha, ndigits=3, nsmall=3), ".csv"), 
          append = F)
    }

for (disease_prevalence in c(seq(from = 0.01, to = 0.05, by = 0.005), 0.10, 0.15, 0.20)){
    simulated_data = import_data(data_path = "./", 
                                 name_pattern = paste0("simu_","*_prev_", format(disease_prevalence, ndigits=3, nsmall=3), "*.csv"),
                                 test = FALSE, print = FALSE)
    simulated_data_formatted <- format_data(simulated_data)
    for (alpha in c(0.001, 0.01, 0.05)){
        type_i_res <- calculate_type_I_err(simulated_data_formatted, significance_level_threshold = alpha)
        write.table(type_i_res, file = paste0(output_path, "Table_type_i_err_alpha_",format(alpha, ndigits=3, nsmall=3), ".csv"), append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

for (alpha in c(0.001, 0.05)){ # select the alpha levels of interest
    res = fread(paste0(output_path, "Table_type_i_err_alpha_",format(alpha, ndigits=3, nsmall=3), ".csv"))
    type_i_res_case_only <- res %>% 
        filter(design %in% c("case_only", "case_control_R1")) %>% arrange(design, para_prevalence)
    print(paste0("When significance level is ", format(alpha, ndigits = 3, nsmall = 3), ", the type I error of case-only design (based on ", unique(res$count), " replicates) is "))
    print(knitr::kable(type_i_res_case_only, digits = 10, align = "c"))
}


############### Evaluating power ###############

# import simulated data
simulated_data <- import_data(data_path = "./", name_pattern = "simu_*.csv", test = FALSE)

simulated_data_formatted <- format_data(simulated_data)

significance_level_threshold = 0.05
stat_power_res <- calculate_stat_power(simulated_data_formatted, significance_level_threshold)
draw_stat_power_figure(stat_power_res, output_path = output_path, output_name = "Fig_stat_power")


############### Calculating bias ###############
OR_G <- 1.2
OR_E <- 2
OR_GE <- 1.2
MAF=0.2
freq.E=0.1
print(paste0("beta1 = ", format(log(OR_G), digits=3, nsmall=3)))
print(paste0("beta2 = ", format(log(OR_E), digits=3, nsmall=3)))
print(paste0("beta3 = ", format(log(OR_GE), digits=3, nsmall=3)))
calculate_analytical_bias(OR_G, OR_E, OR_GE, MAF, freq.E, output_path = output_path, prefix = "Table_analytical_bias")
