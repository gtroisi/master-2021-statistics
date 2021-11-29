#############################################
# Packages
install.packages("pbmcapply")    # Pacchetto per la computazione in parallelo
library("tibble")
library("dplyr")
library("pbmcapply")

#############################################
### 1 SCREENING

## 1.1 Acquisizione e pulizia del Dataset

# Lettura del dataset
cancer.df <- read.table("Dataset_3.txt", header = TRUE)

# Analisi della struttura del dataset
str(cancer.df)    # n = 102 -> numero di osservazioni; p = 6033 -> numero di variabili esplicative

# Estrazione delle variabili esplicative e di risposta
X <- as.matrix(cancer.df[ , -1L])
y <- cancer.df$Status

# La variabile di risposta è convertita in una variabile numerica (0 = "control" - 1 = "cancer")
y <- ifelse(y == "control","0","1")
y <- as.numeric(y)

# Definizione dei parametri (dimensioni) del modello
n <- nrow(X)         # numero di osservazioni
p <- ncol(X)         # numero di variabili

#############################################
## 1.2 Calcolo del P-Value

# Funzione per il calcolo del P-Value

marginPval <- function(X, y, permutation = FALSE, B = 500L, n.threads = 1L) {

  # Calcolo del P-Value col metodo della permutazione
  if(permutation) {
    
    # Fix seed per la generazione di numeri casuali
    set.seed(1)
    
    # Estrazione dei coefficienti del summary per ogni variabile
    temp <- t(sapply(seq.int(p), function(i) summary(lm(y ~ X[, i]))$coef[2L, ]))
    
    # Creazione della variabile temporanea "id" da usare come indice per la permutazione
    id <- split(c(replicate(B, sample(n))), rep(seq.int(B), each = n))
  
    # Applicazione della permutazione (tramite computazione in parallelo)
    t.perm <- pbmclapply(seq.int(B), function(i) fun.permutation(id = id[[i]], X = X, y = y), mc.cores = n.threads)
    t.perm <- matrix(unlist(t.perm), ncol = B)
    
    # Calcolo del P-Value permutato
    Pval.cmptd <- rowMeans(abs(t.perm) >= abs(temp[, 3L]))
  }
  
  # Calcolo del P-Value con approccio teorico
  else {
    Pval.cmptd <- apply(X, 2L, fun.theoretical, y = y)
  }
  
  # Creazione della matrice di Output
  temp.Pval <- data.frame(Var = colnames(X), p.val = Pval.cmptd)    # Combinazione del nome delle osservazioni col nuovo P-Value calcolato
  Pval <- as_tibble(temp.Pval)                                      # Conversione in un tibble
  out <- list(call = match.call(), summary = Pval)                  # Creazione variabile output
  
  ## Definizione della classe per la variabile in output
  class(out) <- "marginPval"
  
  return(out)
}

# Funzione per il metodo teorico
fun.theoretical <- function(x, y){
  out.summary <- summary(lm(y ~ x))$coefficients
  out <- out.summary[2L, "Pr(>|t|)"]               # Estrazione del P-Value
  return(out)
}

# Funzione per il metodo di permutazione
fun.permutation <- function(id, X, y){
  y.perm <- y[id]
  sapply(seq.int(p), function(i) summary(lm(y.perm ~ X[, i]))$coef[2L, 3L])   # Estrazione del T-value (calcolato con la variabile di risposta permutata)
}

#############################################
# Approccio iterativo per il calcolo del P-Value
permutation = TRUE                                 # FALSE: Theoretical Approach,
                                                   # TRUE:  Permutation Approach
out.Pval <- marginPval(X = X, y = y, 
                       permutation = permutation,
                       B = 1000L,
                       n.threads = 8L)

#############################################
# Metodi della classe "marginPval"
hist.marginPval <- function(x, ...) {
  if(permutation){
    main = "P-Value distribution - permutation null distribution"
  }
  else{
    main = "P-Value distribution - theoretical null distribution"
  }
  xlab = "P-Values"
  col = "lightskyblue3"
  p.val <- x$summary %>% pull(p.val)
  col = "lightskyblue3"
  hist(p.val, xlab = xlab, col = col, main = main, nclass = 15, ylim = c(0,600))
}

print.marginPval <- function(x, ...) {
  out <- out.Pval$summary
  print(out, n = 20)
}

hist(out.Pval)
print(out.Pval)

#############################################
# Salvataggio dei risultati
#save.image("Screening_pt1.RData")
