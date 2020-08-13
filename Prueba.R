# Ejercicio 1 -------------------------------------------------------------
#Primero necesitamos cargar los datos 
mis.colores <- colorRampPalette(c("yellow", "red", "yellow"))
muestra <- read.csv("datos/muestra.csv")[1]$muestra
(mu_hat<-mean(muestra))
#Para ver que forma tienen los datos, veamos un histograma de la muestra.
hist(muestra,freq = F,col = mis.colores(14),xlab = "Muestra observada")
lines(density(muestra),lwd=3)
#Pongamos la densidad de la normal(0,1),para tener una referencia 
curve(dnorm(x),col="green",add = T,lwd = 3)
#Propongamos una N(mu_hat,sd(muestra)) 
curve(dnorm(x,mu_hat,sd(muestra)),col="red",add = T,lwd = 3)
#Vemos que tiene similitud la esta normal,vemos que pasa si modicamos sd
#Sd=1
curve(dnorm(x,mu_hat,1),col="blue",add = T,lwd = 3)
#Pongamos otra normal con promedio de sd de sd(muestra) y 1 
curve(dnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd = 3)
#graficas_densidades<-function(){
#hist(muestra,freq = F,col = mis.colores(14),xlab = "Muestra observada")
#lines(density(muestra),lwd=3)
#curve(dnorm(x),col="green",add = T,lwd = 3)
#curve(dnorm(x,mu_hat,sd(muestra)),col="red",add = T,lwd = 3)
#curve(dnorm(x,mu_hat,1),col="blue",add = T,lwd = 3)
#curve(dnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd = 3)
#legend(max(muestra)-3.5,0.35,c("N(0,1)","N(mu_had,1)","N(mu_had,sd(muestra))",
#"N(mu_hat,sum(sd(muestra),1)/2)"),
#col =c("green","red","blue","purple"),lty=1,lwd=2,cex=0.5)
#}
#Esta funcion se utilizara  Para mostrar las comparaciones en el proyecto 
#Ahora veamos como se comportan las distribuciones 
f_gorro<-ecdf(muestra)
curve(f_gorro(x),col="yellow",-5,5,lwd=3)
curve(pnorm(x),col="green",add = T,lwd=2.5)
curve(pnorm(x,mu_hat,sd(muestra)),col="red",add = T,lwd=2.5)
curve(pnorm(x,mu_hat,1),col="blue",add = T,lwd=2.5)
curve(pnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd=2.5)
#graficas_distribuciones<-function(){
#       curve(f_gorro(x),col="yellow",-5,5,lwd=3)
#        curve(pnorm(x),col="green",add = T,lwd=2.5)
#        curve(pnorm(x,mu_hat,sd(muestra)),col="red",add = T,lwd=2.5)
#        curve(pnorm(x,mu_hat,1),col="blue",add = T,lwd=2.5)
#        curve(pnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd=2.5)
#        legend(min(muestra)-6,1,c("Muestra","N(0,1)","N(mu_had,1)",
#                                  "N(mu_had,sd(muestra))",
#                                  "N(mu_hat,sum(sd(muestra),1)/2)"),
#               col =c("yellow","green","red","blue","purple"),lty=1,lwd=2,cex=.6)
#}
#Observese que las tres normales centradas en Mu son muy similares a la muestra
#Asi que aplicaremos la prueba de Anderson-Darling ya que es la muestra mas 
#estricta de las pruebas vistas en clase.La significancia la tomamos de 0.05
#Pruba para N(mu,sd(muestra))
#Esta funcion se utilizara  Para mostrar las comparaciones en el proyecto 
prueba<-function(){
        N_sdmustra<-goftest::ad.test(muestra, pnorm, mean =mu_hat, sd =sd(muestra))
        N_1<-goftest::ad.test(muestra, pnorm, mean =mu_hat, sd =1)
        N_1.088<-goftest::ad.test(muestra, pnorm, mean =mu_hat, sd =1.088)
        Resultados<-list(Normal_Sdmuestra=N_sdmustra$p.value,Normal_1=N_1$p.value,
                         Normal_1.088=N_1.088$p.value)
        return(Resultados)
}
prueba()
#Vemos que en ninguna prueba rechazamos, asi que escojamos
#el que tenga mayor p-value.
#Vemos que la ultima pureba es la que tiene mayor p-value 
#Entonces podemos decir que Nuestros datos se distribuye como N(mu,1.088)
#La aceptacionse puede deber que a que las desviaciones en el redondeo son 1, 
#pero el p-value mas grande se dio con 1.088
#Asi Concluimos que la muestra de distribuya como una N(mu,sd=1.088)
hist(muestra,freq = F,col = mis.colores(14),xlab = "Muestra observada")
lines(density(muestra),lwd=3)
curve(dnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd = 3)
legend(max(muestra)-3.5,0.35,c("Distribucion de la muestra","Distribucion Propuesta"),
       col =c("Black","Purple"),lty=1,lwd=2)
curve(f_gorro(x),col="yellow",-5,5,lwd=3)
curve(pnorm(x,mu_hat,sum(sd(muestra),1)/2),col="purple",add = T,lwd=2.5)
legend(min(muestra)-6,1,c("Distribucion de la muestra","N(mu,sd=1.088)"),
       col =c("yellow","purple"),lty=1,lwd=2)

# Ejercicio 2 -------------------------------------------------------------
#boot parametrico
#Hacemos una funcion teniendo en cuenta que los datos se distribuyen de una 
#N(mu,1.088)
Bootparametrico<-function(datos,replicas=1000,confianza=.90){
        set.seed(432)
        mu_bp<-vector()
        for (i in 1:replicas) {
                mu_bp[i]<-mean(rnorm(length(datos),mu_hat,(sd(muestra)+1)/2))
                
        }
        mu_bp_hat<-mean(mu_bp)
        desviacion<-sd(mu_bp)
        alpha<-1-confianza
        vecpivotal<-c(mu_bp-mu_bp_hat)
        nor<-mean(datos)+c(-qnorm((alpha)/2,lower.tail = F)*desviacion,
                           qnorm((alpha)/2,lower.tail = F)*desviacion)
        ampN<-nor[2]-nor[1]
        resultados<-list(Estadistico = mu_hat,Estimacion = mu_bp_hat,
                         Sd_hat_bp = desviacion ,
                         min_max = range(mu_bp),
                         Int.Normal = nor,
                         Amplitud.Normal = ampN)
        hist(mu_bp,freq = F,col = mis.colores(14),main = "Bootstrap Parametrico",
             xlab = "mu_hat", border = mis.colores(14),lty="blank")
        rug(mean(datos), lwd = 2)
        rug(mu_bp_hat, lwd = 2, col = "blue")
        abline(v = 2.5, lty = 2, col = "black", lwd = 2)
        abline(v = nor[1],lty = 2, col = "Tomato", lwd = 2)
        abline(v = nor[2],lty = 2, col = "Tomato", lwd = 2)
        lines(density(mu_bp),lwd=3,col="blue")
        legend(min(mu_bp)-.025,3.5,col = c("Tomato",
                                           "Black"),legend = c("Int.Normal","Mu real"),
               lty =2,lwd = 2, cex=0.7)
        
        return(resultados)
}
(Parametrico<-Bootparametrico(muestra))
#Ahora haremos una funcion para calcular el bootstap no parametrico

Bootnoparametrico<-function(datos,replicas=1000,confianza=.90){
        set.seed(432)
        mu_bnp<-vector()
        for (i in 1:replicas) {
                mu_bnp[i]<-mean(sample(datos,size=length(datos),replace = T))
                
        }
        mu_bnp_hat<-mean(mu_bnp)
        desviacion<-sd(mu_bnp)
        alpha<-1-confianza
        vecpivotal<-c(mu_bnp-mu_bnp_hat)
        # El intervalo pivotal está incorrecto. 
        #  Los cuantiles son 1-alpha/2 y alpha/2
        piv<-c(mu_hat-quantile(vecpivotal,1-alpha),
               mu_hat-quantile(vecpivotal,alpha))
        per<-quantile(mu_bnp,probs = c(alpha/2,1-alpha/2))
        nor<-mean(datos)+c(-qnorm((alpha)/2,lower.tail = F)*desviacion,
                           qnorm((alpha)/2,lower.tail = F)*desviacion)
        ampN<-nor[2]-nor[1]
        ampPiv<-piv[2]-piv[1]
        ampPer<-per[2]-piv[1]
        resultados<-list(Estadistico = mean(datos),Estimacion = mu_bnp_hat,
                         Sd_hat_bnp = desviacion ,
                         min_max = range(mu_bnp),
                         Int.Normal = nor,
                         Amplitud.Normal =ampN,
                         Int.Pivotal=piv,
                         Amplitud.Pivotal =ampPiv,
                         Int.Percentil = per,
                         Amplitud.Percentil =ampPer)
        hist(mu_bnp,freq = F,col = mis.colores(14),main = "Bootstrap No Parametrico",
             xlab = "mu_hat", border = mis.colores(14),lty="blank")
        rug(mean(datos), lwd = 2)
        rug(mu_bnp_hat, lwd = 2, col = "blue")
        abline(v = 2.5, lty = 2, col = "black", lwd = 2)
        abline(v = per[1],lty = 2, col = "Purple", lwd = 2)
        abline(v = per[2],lty = 2, col = "Purple", lwd = 2)
        abline(v = nor[1],lty = 2, col = "Tomato", lwd = 2)
        abline(v = nor[2],lty = 2, col = "Tomato", lwd = 2)
        abline(v = piv[1],lty = 2, col = "dodgerblue3", lwd = 2)
        abline(v = piv[2],lty = 2, col = "dodgerblue3", lwd = 2)
        lines(density(mu_bnp),lwd=3,col="blue")
        legend(min(mu_bnp)-.025,3,col = c("Tomato","dodgerblue3","Purple",
                                          "Black"),legend = c("Int.Normal","Int.Pivotal",
                                                              "Int.Percelil","Mu real"),lty =2,lwd = 2, cex=0.7)
        return(resultados)
}
(BNoparametrico<-Bootnoparametrico(muestra))
#Ahora haremos una fucion para Jackknife
Jackknife<-function(datos,confianza=.90){
        set.seed(432)
        n<-length(datos)
        mu_j<-vector()
        for(i in 1:n){
                mu_j[i]<-mean(datos[-i])
        }
        mu_j_hat<-mean(mu_j)
        alpha<-1-confianza
        desviacion<-sqrt(((n-1)/n)*sum((mu_j-mu_j_hat)^2))
        nor<-mean(datos)+c(-qnorm((alpha)/2,lower.tail = F)*desviacion,
                           qnorm((alpha)/2,lower.tail = F)*desviacion)
        ampN<-nor[2]-nor[1]
        resultados<-list(Estadistico = mean(datos),Estimacion =mu_j_hat,
                         Sd_hat_jackk = desviacion ,
                         min_max = range(mu_j) ,
                         Int.Normal = nor,
                         Amplitud.Normal = ampN)
        hist(mu_j,freq = F,col = mis.colores(14),main = "Jackknife",
             xlab = "mu_j_hat", border = mis.colores(14),lty="blank")
        rug(mean(datos), lwd = 2)
        rug(mu_j_hat, lwd = 2, col = "blue")
        abline(v = 2.5, lty = 2, col = "black", lwd = 2)
        abline(v = nor[1],lty = 2, col = "Tomato", lwd = 2)
        abline(v = nor[2],lty = 2, col = "Tomato", lwd = 2)
        lines(density(mu_j),lwd=3,col="blue")
        legend(2.39,30,col = c("Tomato","Black"),legend = c("Int.Normal",
                                                            "Mu real"),lty =2,lwd = 2, cex=0.7)
        return(resultados)
}
(Jackk<-Jackknife(muestra))
# Ejercicio 3 -------------------------------------------------------------
#Todoslos intervalos ya fueron programados anteriomente en las funciones pasadas
#Aqui solamente mostraremos los resultados 
#Intervalos Normales----
#BootsParametrico
Parametrico$Int.Normal
#Boots no parametrico
BNoparametrico$Int.Normal
#Jackknife
Jackk$Int.Normal
#Intervalos Pivotal----
BNoparametrico$Int.Pivotal
#Intervalo Percentil----
BNoparametrico$Int.Percentil
#Amplitudes de intervalos -----
#Intervalos Normales----
#BootsParametrico
Parametrico$Amplitud.Normal
#Boots no parametrico
BNoparametrico$Amplitud.Normal
#Jackknife
Jackk$Amplitud.Normal
#Intervalos Pivotal----
BNoparametrico$Amplitud.Pivotal
#Intervalo Percentil----
BNoparametrico$Amplitud.Percentil
