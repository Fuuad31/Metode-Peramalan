library("tseries")
library("forecast")
library("TTR")
library("TSA")
library("graphics")
library("astsa")
library("car")
library("portes")
library("ggplot2")

dt = ts(Data_Responsi_1_[2])

# Create a time series plot
plot(dt, main = "Proses Industri",
     xlab = "Year-Month", ylab = "Proses Industri")

# Plot ACF
ggAcf(dt,lag.max = 15) + ggtitle("ACF")

seasonplot(dt,6,main="Proses Produksi", ylab="Produksi",year.labels = TRUE, col=rainbow(18))

#Diff 1 with Log-Trans
dtrans1 = diff(log(dt), differences=1)
adf.test(dtrans1)

#Diff 2 with Log-Trans
dtrans2 = diff(log(dt), differences=2)
adf.test(dtrans2)

# Log-Transformasi d = 1 dan D = 1
dtrans11 = diff(diff(log(dt), lag=6), differences=1)
adf.test(dtrans11)
autoplot(dtrans11, main = "Time series plot of logtrans d=1, D=6")

# Log-Transformasi d = 2 dan D = 1
dtrans21 = diff(diff(log(dt), lag=6), differences=2)
adf.test(dtrans21)
autoplot(dtrans21, main = "Time series plot of logtrans d=2, D=6")

# Log-Transformasi d = 1 dan D = 2
dtrans12 = diff(diff(log(dt), lag=12), differences=1)
adf.test(dtrans12)
autoplot(dtrans12, main = "Time series plot of logtrans d=1, D=12")

# Log-Transformasi d = 2 dan D = 2
dtrans22 = diff(diff(log(dt), lag=12), differences=2)
adf.test(dtrans22)
autoplot(dtrans22, main = "Time series plot of logtrans d=2, D=12")

# Plot ACF
ggAcf(dtrans21,lag.max = 15) + ggtitle("ACF dtrans21")

# Plot PACF
ggPacf(dtrans21,lag.max = 15) + ggtitle("PACF dtrans21")

create_model_name <- function(p, d, q, P, D, Q, constant, multiplicative) {
  paste0(ifelse(constant, "c", "nc"), p, d, q, P, D, Q, ifelse(multiplicative, "m", "a"))
}

results <- data.frame(Model = character(), Coef = character(), P_Value = numeric(), Conclusion = character(), stringsAsFactors = FALSE)
model_list <- list()
max_p <- 2
max_d <- 2
max_q <- 1
max_P <- 2
max_D <- 1
max_Q <- 2
seasonal_period <- 6
multiplicative <- TRUE  # Set multiplicative to TRUE as specified

for (p in 0:max_p) {
  for (d in max_d:max_d) {
    for (q in 0:max_q) {
      if (p == 0 && q == 0) next
      for (P in 0:max_P) {
        for (D in max_D:max_D) {
          for (Q in 0:max_Q) {
            if (!multiplicative) {
              if (p == P && d == D && q == Q) next
            } else {
              if ((P == 0 && Q == 0) || (P != 0 && Q != 0)) next
            }
            for (constant in c(FALSE)) {
              model_name <- create_model_name(p, d, q, P, D, Q, constant, multiplicative)
              
              tryCatch({
                model <- Arima(dt, order = c(p, d, q), 
                               seasonal = list(order = c(P, D, Q), period = seasonal_period), 
                               include.constant = constant, lambda = 0)
                test_results <- coeftest(model)
                if (any(is.nan(test_results[ , "Pr(>|z|)"]))) next
                
                significant <- TRUE
                for (i in 1:nrow(test_results)) {
                  coef_name <- rownames(test_results)[i]
                  p_value <- round(test_results[i, "Pr(>|z|)"], 3)
                  conclusion <- ifelse(p_value < 0.05, "Significant", "Not Significant")
                  results <- rbind(results, data.frame(Model = model_name, Coef = coef_name, P_Value = p_value, Conclusion = conclusion, stringsAsFactors = FALSE))
                  if (conclusion == "Not Significant") {
                    significant <- FALSE
                  }
                }
                if (significant) {
                  model_list[[model_name]] <- model
                }
              }, error = function(e) {
                cat("Error fitting model:", model_name, "\n")
              })
            }
          }
        }
      }
    }
  }
}
print(results)

### Detect all significant models
all_significant_models <- data.frame(Model = names(model_list), stringsAsFactors = FALSE)
print(all_significant_models)

# Diagnostic Checking
run_diagnostic_tests <- function(model, model_name) {
  results <- data.frame(Model = character(), Diagnostic = character(), P_Value = numeric(), Result = character(), stringsAsFactors = FALSE)
  
  ## No-Autocorrelation Residuals Test
  test <- Box.test(model$residuals, type = "Ljung")
  p_value <- round(test$p.value, 3)
  result <- ifelse(p_value > 0.05, "No Autocorrelation", "Autocorrelated")
  results <- rbind(results, data.frame(Model = model_name, Diagnostic = "Autocorrelation", P_Value = p_value, Result = result))
  
  ## Homoscedasticity Test
  test <- Box.test((model$residuals)^2, type = "Ljung")
  p_value <- round(test$p.value, 3)
  result <- ifelse(p_value > 0.05, "Homoscedastic", "Heteroscedastic")
  results <- rbind(results, data.frame(Model = model_name, Diagnostic = "Homoscedasticity", P_Value = p_value, Result = result))
  
  ## Normality Test
  test <- jarque.bera.test(model$residuals)
  p_value <- round(test$p.value, 3)
  result <- ifelse(p_value > 0.05, "Normally Distributed", "Not Normally Distributed")
  results <- rbind(results, data.frame(Model = model_name, Diagnostic = "Normality Residuals", P_Value = p_value, Result = result))
  
  return(results)
}

all_diagnostics <- data.frame(Model = character(), Diagnostic = character(), P_Value = numeric(), Result = character(), stringsAsFactors = FALSE)
for (model_name in all_significant_models$Model) {
  model <- model_list[[model_name]]
  diagnostics <- run_diagnostic_tests(model, model_name)
  all_diagnostics <- rbind(all_diagnostics, diagnostics)
}
print(all_diagnostics)

# Model Selection by Log-Likelihood, AIC, and BIC
model_selection <- data.frame(Model = character(), LogLik = numeric(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)
for (model_name in all_significant_models$Model) {
  model <- model_list[[model_name]]
  model_info <- data.frame(Model = model_name, 
                           LogLik = logLik(model), 
                           AIC = AIC(model), 
                           BIC = BIC(model))
  model_selection <- rbind(model_selection, model_info)
}
print(model_selection)

forecast(nc221210m,1)
