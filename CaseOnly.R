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

# calculate
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
