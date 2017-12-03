#'Funcao para gerar valores da Weibull Discreta
#'
#'Esta funcao gera valores aleatorios com distribuicao Weibull Discreta.
#'
#' @param n tamanho da amostra
#' @param q parametro do modelo
#' @param b parametro do modelo
#' 
#' @details 
#' \code{q} e um numero entre 0 e 1 
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' rweidc(100,0.5,0.9)
#' 
#' @export
#' 
rweidc<-function(n,q,b) {
  mt<-matrix(ncol=n,nrow=1,0)
  m<-weidc(0:20000,q,b)
  for (i in 2:20001) {
    m[i]<-sum(m[i],m[i-1])
  }
  for (j in 1:n) {
    k<-runif(1,0,1)
    mt[j]<-which(m>k)[1]-1
  }
  as.vector(mt)
}


#' Distribuicao de probabilidades da Weibull Discreta
#' 
#' Esta funcao nos fornece os valores das probabilidades da Distribuicao Weibull Discreta
#'
#' @param x amostra
#' @param q parametro do modelo
#' @param b parametro do modelo
#' 
#' @details 
#' \code{q} e um numero entre 0 e 1 
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' x = rweidc(100,0.5,0.9)
#' weidc(x,0.5,0.9)
#' 
#' @export
#' 
weidc<-function(x,q,b){ ((q)^(x^b)) - (q)^((x+1)^b) }


#'Funcao de sobrevivencia da Weibull Discreta
#'
#'Esta funcao nos fornece os valores das probabilidades de sobrevivencia da Distribuicao Weibull Discreta
#' @param x amostra
#' @param q parametro do modelo
#' @param b parametro do modelo
#' 
#' @details 
#' \code{q} e um numero entre 0 e 1
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' x = rweidc(100,0.5,0.9)
#' sob.weidc(x,0.5,0.9)
#' 
#' @export
#' 
sob.weidc<-function(x,q,b){ (q)^((x+1)^b) }


#' Funcao de risco ou taxa de falha da Weibull Discreta
#' 
#' Esta funcao nos fornece os valores as taxas de falha da Distribuicao Weibull Discreta
#' @param t amostra (tempos)
#' @param q parametro do modelo
#' @param b parametro do modelo
#' 
#' @details 
#' \code{q} e um numero entre 0 e 1
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' t = rweidc(100,0.5,0.9)
#' HWD(t,0.5,0.9)
#' 
#' @export
#' 
HWD<-function(t,q,b){ 1-(q)^((t+1)^b-(t)^b) }



#' Funcao de verossimilhanca da Weibull Discreta
#' 
#' Esta funcao nos fornece o valor Verossimilhanca da Distribuicao Weibull Discreta
#' 
#' @param t amostra (tempos)
#' @param par parametro do modelo (q,b)
#' @param censura censuras do banco
#' 
#' @details 
#' \code{q} e um numero entre 0 e 1
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' t = rweidc(100,0.5,0.9)
#' par = c(0.5,0.9)
#' censura = c(rep(10,0))
#' LWD(t,par,censura)
#' 
#' @export
#' 
LWD <-function(par,t,censura){
  q<-par[1]
  b<-par[2]
  
  if ( (b>0)&&(q>0)&&(q<1) )
    return (-1*sum( censura*log(weidc(t,q,b)) + (1-censura)*log(sob.weidc(t,q,b)) ) )
  else return (-Inf)
} 


#' Log-verossiminhanca da Weibull Discreta
#' 
#' Esta funcao nos fornece o valor da Log-Verossimilhanca da Distribuicao Weibull Discreta
#' 
#' @param ti amostra (tempos)
#' @param theta parametro do modelo (q,b)
#' @param c censuras do banco
#'  
#' @details 
#' \code{q} e um numero entre 0 e 1
#' 
#' \code{b} e um numero real
#' 
#' @examples 
#' ti = rweidc(100,0.5,0.9)
#' theta = c(0.5,0.9)
#' c = c(rep(10,0))
#' logweidc(ti,theta,c)
#' 
#' @export
#' 
#'  
logweidc <- function(c,ti,theta){
  l = sum(c * log((theta[1]^(ti^theta[2])) - (theta[1]^((ti+1)^theta[2])))) +
    sum((1-c) * ((ti+1)^theta[2]) * log(theta[1]) )
  return(l)
}

