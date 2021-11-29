#############################################
# Packages

load("Screening_pt1.RData")

library("tibble")
library("dplyr")

############################################
## 1.3 Riduzione delle variabili tramite Multiple Test Procedures

summary.pval <- out.Pval$summary

MTP <- function(out, method = c("Bonferroni",
                                       "Sidak",
                                       "Bonferroni.Holm",
                                       "Sidak.Holm",
                                       "Hochberg",
                                       "Benjamin.Hochberg",
                                       "Benjamin.Yekutieli",
                                       "Lehmann.Romano.1",
                                       "Lehmann.Romano.2"),
                            alpha = 0.05,                     # Livello di significatività: 5%
                            c = 0.1) {
  this.call <- match.call()
  method <- match.arg(method)                                 # Selezione del metodo
  
  out <- out %>%
    arrange(p.val)                                            # Ordinamento del p-value (in ordine crescente)
                                                              # Questo si applica già la prima volta che la funzione
                                                              # viene chiamata. Non influenza i metodi One-Step
  
  if (method == "Bonferroni") {
    out[method] <- ifelse(out$p.val <= alpha / p,
                          "reject H0", "do not reject H0")
  }
  
  else if (method == "Sidak") {
    out[method] <- ifelse(out$p.val <= 1 - (1 - alpha)^(1/p),  
                          "reject H0", "do not reject H0")
  }
  
  else if (method == "Bonferroni.Holm") {
    out[method] <- rep("do not reject H0", times = p)

    for(j in seq_len(p)) {
      temp.p.val <- out$p.val[j]
      if (temp.p.val <= alpha / (p - j + 1))
        out[j, method] <- "reject H0"
      else break
    }
  }
  
  else if (method == "Sidak.Holm") {
    out[method] <- rep("do not reject H0", times = p)
    
    for(j in seq_len(p)) {
      temp.p.val <- out$p.val[j]
      if (temp.p.val <= 1 - (1 - alpha)^(1 / (p - j + 1)) )
        out[j, method] <- "reject H0"
      else break
    }
  }
  
  else if (method == "Hochberg") {
    out[method] <- rep("reject H0", times = p)

    for(j in p:1) {
      temp.p.val <- out$p.val[j]
      if (temp.p.val > alpha / (p - j + 1))
        out[j, method] <- "do not reject H0"
      else break
    }
  }
  
  else if (method == "Benjamin.Hochberg") {
    out[method] <- rep("reject H0", times = p)
      
    for(j in p:1) {
      temp.p.val <- out$p.val[j]
      if (temp.p.val > j * alpha / p)
        out[j, method] <- "do not reject H0"
      else break
    }
  }

  else if (method == "Benjamin.Yekutieli") {
    out[method] <- rep("reject H0", times = p)
    
    cost <- sum(1 / p:1)
    
    for(j in p:1) {
      temp.p.val <- out$p.val[j]
      if (temp.p.val > j * alpha / (p * cost))
        out[j, method] <- "do not reject H0"
      else break
    }
  }
  
  else if (method == "Lehmann.Romano.1") {
    out[method] <- rep("do not reject H0", times = p)
    
    for(j in seq_len(p)) {
      temp.p.val <- out$p.val[j]
      alpha_j <- alpha * (floor(c * j) + 1) / (p - j + floor(c * j) + 1)
      if (temp.p.val <= alpha_j)
        out[j, method] <- "reject H0"
      else break
    }
  }
  
  else if (method == "Lehmann.Romano.2") {
    out[method] <- rep("do not reject H0", times = p)
    
    cost <- sum(1 / seq_len(floor(c * p) + 1))
    
    for(j in seq_len(p)) {
      temp.p.val <- out$p.val[j]
      alpha_j <- alpha * (floor(c * j) + 1) / (cost * (p - j + floor(c * j) + 1))
      if (temp.p.val <= alpha_j)
        out[j, method] <- "reject H0"
      else break
    }
  }
  return(out)
}

# Ogni chiamata alla funzione MTP, aggiunge una nuova colonna "method" al summary
summary.pval <- MTP(summary.pval, "Bonferroni")
summary.pval <- MTP(summary.pval, "Sidak")
summary.pval <- MTP(summary.pval, "Bonferroni.Holm")
summary.pval <- MTP(summary.pval, "Sidak.Holm")
summary.pval <- MTP(summary.pval, "Hochberg")
summary.pval <- MTP(summary.pval, "Benjamin.Hochberg")
summary.pval <- MTP(summary.pval, "Benjamin.Yekutieli")
summary.pval <- MTP(summary.pval, "Lehmann.Romano.1")
summary.pval <- MTP(summary.pval, "Lehmann.Romano.2")

#############################################
# Riepilogo
MTP.summary <- summary.pval[c(-1,-2)] %>%                 # Si conta il numero di occorrenze delle due categorie (reject H0/do not reject H0)
  apply(2, table)                                         # e si dispongono in una tabella

MTP.summary = t(MTP.summary)                              # Si traspone la tabella per facilitare la lettura

MTP.summary

#############################################
# Salvataggio dei risultati
# save.image("Screening_pt2.RData")