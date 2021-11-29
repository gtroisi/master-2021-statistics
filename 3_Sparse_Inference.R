#############################################
# Packages
load("Screening_pt2.RData")

install.packages("e1071")
install.packages("glmnet")
install.packages("caret")
library("tibble")
library("dplyr")
library("glmnet")
library("caret")

#############################################
## 2 - SPARSE INFERENCE

# Estrazione (del nome) delle variabili individuate durante la fase di screening
var.red <- summary.pval %>% 
    filter(Benjamin.Hochberg == "reject H0") %>%
    pull(Var) %>% as.character

# Estrazione degli identificatori di colonna (boolean) della matrice ridotta
id.red <- is.element(colnames(X), var.red)     # Verifica la presenza del gene all'interno di X.
                                               # Se esiste, id = TRUE, viceversa id = FALSE.

# Design Matrix
X.red <- as.matrix(X[, id.red])
dim(X.red)

###########################################
## 2.1 Regularized Regression Models

# Lasso Regression
out.lasso <- glmnet(x = X.red, y = y, family = "binomial", alpha = 1)    # alpha = 1 -> per Lasso

## 2.2 Selezione del parametro di penalizzazione ottimale tramite Information Criteria

# Parametri
p <- ncol(X.red)                # Numero di predittori totali
npar <- out.lasso$df + 1        # Numero di parametri stimati (df: numero di coefficienti non nulli)
lambda <- out.lasso$lambda      # Estrazione dei parametri di tuning
nlambda <- length(lambda)       # Numero dei parametri di tuning

# Costruzione della matrice contenente le misure delle bontà di adattamento
IC.lasso <- matrix(0, nrow = nlambda, ncol = 5)
colnames(IC.lasso) <- c("AIC", "AICc", "GIC", "BIC", "BICg")

IC.lasso[, "AIC"]  <- deviance(out.lasso) + 2 * npar
IC.lasso[, "AICc"] <- deviance(out.lasso) + 2 * npar * n / (n - npar - 1)
IC.lasso[, "GIC"]  <- deviance(out.lasso) + npar * log(log(n)) * log(p)
IC.lasso[, "BIC"]  <- deviance(out.lasso) + log(n) * npar
g <- 0.5
IC.lasso[, "BICg"] <-  deviance(out.lasso) + log(n) * npar + 2 * npar * g * log(p)

# Identificazione dei valori ottimi del parametro di tuning
lambda.opt.IC <- out.lasso$lambda[apply(IC.lasso, 2L, which.min)]   # Si selezionano i lambda minimi per i 5 metodi
names(lambda.opt.IC) <- colnames(IC.lasso)

# Stima dell'insieme A_hat
A.hat.lasso.IC <- predict(out.lasso, type = "nonzero", s = lambda.opt.IC)
names(A.hat.lasso.IC) <- colnames(IC.lasso)

# Variabili selezionate
nms.IC <- lapply(A.hat.lasso.IC, function(id) names(out.lasso$beta[, 1L])[id])

# Numero di variabili selezionate
df.IC <- sapply(A.hat.lasso.IC, length)    # Calcola la lunghezza (delle liste) di A_hat (5 liste per ogni metodo)

###########################################
## 2.3 Risultati

# Plot del percorso della soluzione
plot(out.lasso,
     xvar = "lambda",
     xlab = expression(log(lambda)),
     lwd = 2,
     cex.lab = 1.6,
     col = rainbow(20))

title("Modello di Regressione Logistica col metodo di Penalizzazione di Lasso", line = 2.5)

# Grafico
matplot(log(out.lasso$lambda),
        IC.lasso,
        type = "l",
        lwd = 3,
        lty = 1,
        main = "Parametri di tuning valutati con metodi ICs",
        xlab = expression(log(lambda)),
        ylab = "Information Criteria",
        col = rainbow(5))

# Retta che mostra i minimi delle curve
abline(v = log(lambda.opt.IC),
       lty = 2,
       lwd = 1,
       col = rainbow(5))

legend("topright", legend = colnames(IC.lasso), lty = 1, lwd = 3, col = rainbow(5), bty = "n")

# Visualizzazione delle variabili e del numero totale di variabili per i metodi IC
nms.IC
df.IC

# Stima della variabile di risposta
y.hat.aic <- predict(out.lasso, newx = X.red, type = "response", s = lambda.opt.IC["AIC"])

y.hat.bic <- predict(out.lasso, newx = X.red, type = "response", s = lambda.opt.IC["BIC"])

# Stima della capacità predittiva
prediction <- function(y.hat) {
    result <- confusionMatrix(factor(ifelse(y.hat >= .5, 1, 0), levels = c(0,1)),    # Si approssimano i valori ottenuti a:
                                  factor(y, levels = c(0,1)))                        # 0 (se y < 0.5); 1 (se y >= 0.5)
    
    out <- list(method = colnames(y.hat), summary = result)
    class(out) <- "prediction"                                   # Attribuzione alla classe "prediction"
    return(out)
}

# Generazione dei risultati mediante funzione
result.aic = prediction(y.hat.aic)
result.bic = prediction(y.hat.bic)

# Metodi della classe "prediction"
print.prediction <- function(x, ...) {
    method.name = x$method[1L]
    matrix = x$summary$table
    accuracy = round(x$summary$overall[1L],3)
    kappa = round(x$summary$overall[2L],3)
    cat(c("Method: ", method.name, "\n \n"))
    print(matrix)
    cat(c("\n", "Accuracy:", accuracy*100, "%", "\n", "--------------------", "\n \n"))
}

# Visualizzazione dell'accuratezza dei modelli
result.aic
result.bic

#############################################
# Salvataggio dei risultati
#save.image(file = "SparseInference_results.RData")
