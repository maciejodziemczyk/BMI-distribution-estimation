#dane: https://www.kaggle.com/yersever/500-person-gender-height-weight-bodymassindex
library("maxLik") #Biblioteka do metody NR

# Instrukcja, po ustawieniu œcie¿ki dostêpu, mo¿na uruchomiæ ca³y skrpyt
# wyœwietlone zostan¹ wyniki i wykresy, które opisane s¹ w pracy
# do regulacji wykresów kwantylowych s³u¿y zmienna wiersz, mog¹ca przyjmowaæ parametry
# od 1 do 3 (1-Mê¿czyŸni, 2-Kobiety, 3-Ca³a próbka)

getwd()
setwd('F:/Studia/II stopieñ/WNE/Programowanie narzêdzi analitycznych/Projekt')


#********************************************************************
#                       Dane: zmienna - BMI
#********************************************************************
dane = read.csv(file = "BMI.csv", header = TRUE, sep = ";" , dec = ",")

indeksy = which(dane$Gender=="Male")

# Mê¿czyŸni
x1 = dane$BMI[indeksy]
# Kobiety
x2 = dane$BMI[-indeksy]
# Razem
x3 = dane$BMI
# lista próbek
X = list(x1, x2, x3)

# ***Histogramy***
par(mfrow=c(1,3))
hist(x1, main = "BMI - mê¿czyŸni", xlab = "BMI", ylab = "czêstotliwoœæ")
hist(x2, main = "BMI - kobiety", xlab = "BMI", ylab = "czêstotliwoœæ")
hist(x3, main = "BMI - razem", xlab = "BMI", ylab = "czêstotliwoœæ")

#***Tabele do zapisywania wyników (macierze)***
# macierz statystyk populacji
statystyki = matrix(nrow = 3, ncol = 6)
row.names(statystyki) = c("Mê¿czyŸni", "Kobiety", "Razem")
colnames(statystyki) = c("min", "max", "œrednia", "mediana", "dominanta", "odchylenie std.")

# macierz wyników estymacji
estymacje = matrix(nrow = 3, ncol = 8)
row.names(estymacje) = c("Mê¿czyŸni", "Kobiety", "Razem")
colnames(estymacje) = c("mu_MNW", "sigma_MNW", "mu_GMM", "sigma_GMM", "k_MNW", "lambda_MNW", "k_GMM", "lambda_GMM")

# macierz wynikóW hipotez z³o¿onych
zlozone = matrix(nrow = 2, ncol = 3)
row.names(zlozone) = c("LR_test ","test Walda")
colnames(zlozone) = c("Mê¿czyŸni", "Kobiety", "Razem")

# macierz wynikóW hipotez prostych
proste = matrix(nrow = 2, ncol = 3)
row.names(proste) = c("z test MNW","z test UMM")
colnames(proste) = c("Mê¿czyŸni", "Kobiety", "Razem")
proste

#***Obliczanie statystyk opisowych***
Dominanta = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# zapis do tabeli wyników
i = 1
for(x in X){
  statystyki[i,1] = min(x)            # min
  statystyki[i,2] = max(x)            # max
  statystyki[i,3] = mean(x)           # œrednia
  statystyki[i,4] = median(x)         # mediana
  statystyki[i,5] = Dominanta(x)      # dominanta
  statystyki[i,6] = sd(x)             # odchylenie standardowe
  i = i+1
}
#********************************************************************
#                         Estymacje
#********************************************************************
#----------------------ROZK£AD NORMALNY-----------------------------
#---------------------------MNW-------------------------------------
# funkcja log-wiarygodnoœci
lnl_norm = function(parms){
  mu = parms[1]
  s = parms[2]
  ll = -N/2*log(2*pi)-N*log(s)-1/2*sum((x-mu)^2/(s^2))
  return(ll)
}
# gradient
gr_norm = function(params){
  mu = params[1]
  s = params[2]
  g = rep(0, times=2)
  g[1] = sum((x-mu)/(s^2))
  g[2] = -N/s+sum((x-mu)^2/(s^3))
  return(g)
}
# hessian
he_norm = function(params){
  mu = params[1]
  s = params[2]
  h = matrix(0, nrow=2, ncol=2)
  h[1,1] = -N/(s^2)
  h[2,2] = N/(s^2)-3*sum((x-mu)^2/(s^4))
  h[1,2] = -2*sum((x-mu)/(s^3))
  h[2,1] = h[1,2]
  return(h)
}
# zapisanie do tabeli wyników
i = 1
for(x in X){
  N = length(x)
  wynik = maxNR(fn=lnl_norm, grad=gr_norm, hess=he_norm, start=c(20, 1))
  estymacje[i,1] = wynik$estimate[1]
  estymacje[i,2] = wynik$estimate[2]
  i = i+1
}

#---------------------------GMM-------------------------------------
Qmin_norm = function(parms){
  mu = parms[1]
  s = parms[2]
  m = rbind(x-mu, (x-mu)^2-s^2, (x-mu)^4-3*s^4)
  W = m%*%t(m)/N
  M = rowMeans(m)
  val = -t(M)%*%solve(W)%*%M
  return(val)
}
# zapisanie do tabeli wyników
i = 1
for(x in X){
  N = length(x)
  wynik = maxNR(fn = Qmin_norm,  start = c(mean(x), sd(x)))
  estymacje[i,3] = wynik$estimate[1]
  estymacje[i,4] = wynik$estimate[2]
  i = i+1
}
#----------------------ROZK£AD WEIBULLA-----------------------------
#---------------------------MNW-------------------------------------
# funkcja log-wiarygodnoœci
lnl_wbl = function(parms){
  k = parms[1]
  lambda = parms[2]
  ll = N*log(k)-N*k*log(lambda)+(k-1)*sum(log(x))-sum((x/lambda)^k)
  return(ll)
}
# gradient
gr_wbl = function(parms){
  k = parms[1]
  lambda = parms[2]
  g = rep(0, times = 2)
  g[1] = N/k-N*log(lambda)+sum(log(x))-sum((x/lambda)^k*log(x/lambda))
  g[2] = -N*k/lambda+sum(k*(x^k)/(lambda^(k+1)))
  return(g)
}
# hessian
he_wbl = function(parms){
  k = parms[1]
  lambda = parms[2]
  h = matrix(0, nrow = 2, ncol = 2)
  h[1,1] = -N/(k^2)-sum((x/lambda)^k*log(x/lambda)*log(x/lambda))
  h[2,2] = N*k/(lambda^2)+sum(k*lambda^(-k-2)*x^k*(-k-1))
  h[1,2] = -N/lambda-sum(-k*(x^k)/(lambda^(k+1))*log(x/lambda)+((x/lambda)^k)*(-1/lambda))
  h[2,1] = h[1,2]
  return(h)
}
# zapisanie do tabeli wyników
i = 1
for(x in X){
  N = length(x)
  wynik = maxNR(fn=lnl_wbl, grad=gr_wbl, hess=he_wbl, start=c(0.1, 0.1))
  estymacje[i,5] = wynik$estimate[1]
  estymacje[i,6] = wynik$estimate[2]
  i = i+1
}

#---------------------------GMM-------------------------------------
Qmin_wbl = function(parms){
  k = parms[1]
  lambda = parms[2]
  m = rbind(x-(lambda*gamma(1+1/k)),
            (x-(lambda*gamma(1+1/k)))^2-(lambda^2*(gamma(1+2/k)-(gamma(1+1/k))^2)),
            (x-(lambda*gamma(1+1/k)))^4-(lambda*gamma(1+3/k))^3
  )
  W = m%*%t(m)/N
  M = rowMeans(m)
  val = -t(M)%*%solve(W)%*%M
  return(val)
}
# zapisanie do tabeli wyników
i = 1
for(x in X){
  N = length(x)
  wynik = maxNR(fn = Qmin_wbl,  start = c(3, 40))
  estymacje[i,7] = wynik$estimate[1]
  estymacje[i,8] = wynik$estimate[2]
  i = i+1
}

#********************************************************************
#                       Wykresy kwantylowe
#********************************************************************
# kwantyle empiryczne
emp_quantiles = quantile(x = x1, probs = seq(0.1, 0.99, 0.01))
# kwanyle teoretyczne
wiersz = 3 # 1 - wykresy dla mê¿czyzn, 2 - dla kobiet, 3 - ca³a próba
norm_quantiles_mnw = qnorm(p=seq(0.1, 0.99, 0.01), mean = estymacje[wiersz,1], sd = estymacje[wiersz,2])
weibull_quantiles_mnw = qweibull(p=seq(0.1, 0.99, 0.01), shape = estymacje[wiersz,5], scale = estymacje[wiersz,6])
norm_quantiles_gmm = qnorm(p=seq(0.1, 0.99, 0.01), mean = estymacje[wiersz,3], sd = estymacje[wiersz,4])
weibull_quantiles_gmm = qweibull(p=seq(0.1, 0.99, 0.01), shape = estymacje[wiersz,7], scale = estymacje[wiersz,8])

# rysowanie wykresów
par(mfrow=c(2,2))
plot(norm_quantiles_mnw, emp_quantiles, main = "norm_mnw")
abline(0,1, col='red')
plot(norm_quantiles_gmm, emp_quantiles, main = "norm_gmm")
abline(0,1, col='red')
plot(weibull_quantiles_mnw, emp_quantiles, main = "weibull_mnw")
abline(0,1, col='red')
plot(weibull_quantiles_gmm, emp_quantiles, main = "weibull_gmm")
abline(0,1, col='red')

#********************************************************************
#                         Hipotezy 
#********************************************************************
# z³o¿one
i = 1 
for(x in X){
  # LR_test
  wynik = maxNR(fn=lnl_norm, grad=gr_norm, hess=he_norm, start=c(20, 1))
  summary(wynik)
  lnl_U = wynik$maximum             # model bez ograniczeñ
  lnl_R = lnl_norm(c(37.5, 2.5))    # model z ograniczeniami
  LR_test = 2*(lnl_U-lnl_R)         # nale¿y do X^2 df=2
  p = 1-pchisq(q = LR_test, df = 2) # p_value
  zlozone[1,i] = p
  # test walda
  wynik = maxNR(fn = Qmin_norm,  start = c(mean(x), sd(x)))
  R = diag(2)
  theta = wynik$estimate
  q = rbind(37.5, 2.5)
  S = R%*%theta-q
  vcov = -solve(wynik$hessian)/N
  w_test = t(S)%*%solve(R%*%vcov%*%t(R))%*%S
  p = 1-pchisq(q = w_test, df = 2)     
  zlozone[2,i] = p
  i = i+1
}
# proste
i = 1
for(x in X){
  # z_test MNW
  k_0 = 3
  wynik = maxNR(fn=lnl_wbl, grad=gr_wbl, hess=he_wbl, start=c(0.1, 0.1))
  vcov = -solve(wynik$hessian)/N                    # macierz wariancji kowariancji
  std.err.k = sqrt(vcov[1,1])                       # b³¹d standardowy lambda
  z_test = (wynik$estimate[1]-k_0)/std.err.k        # z_test ~N(0,1)
  p = 2*(1-pnorm(q=abs(z_test), mean = 0, sd = 1))  # p_value
  proste[1,i] = p
  # z_test GMM
  wynik = maxNR(fn = Qmin_wbl,  start = c(3, 40))
  vcov = -solve(wynik$hessian)/N                    # macierz wariancji kowariancji
  std.err.k = sqrt(vcov[1,1])                       # b³¹d standardowy lambda
  z_test = (wynik$estimate[1]-k_0)/std.err.k        # z_test ~N(0,1)
  p = 2*(1-pnorm(q=abs(z_test), mean = 0, sd = 1))  # p_value
  proste[2,i] = p
  i = i+1
}

#********************************************************************
#                         Wyniki 
#********************************************************************
statystyki
estymacje
zlozone
proste











