main_folder <- 'C:/Users/39346/OneDrive/Desktop/THESIS'
setwd(main_folder)

###AUTOCORRELATION FUNCTION###
##daily##
data <- read.csv("daily_278.csv", header = TRUE)
stations <- colnames(data)
m <- dim(data)[2]
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/ACF_plot_278_d.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:m){
  col_data <- na.omit(data[,i])
  emp_acf <- acf(col_data,lag.max = 365,ci = 0.99,plot = FALSE)
  plot(emp_acf$lag,emp_acf$acf,type = "h",ylim = c(-0.1,1), main = banks[i], lwd = 3, las = 1,
       ylab = "", xlab = "", cex.main = 1.6, cex.axis = 1.3)
  title(xlab="Lag", line=1, cex.lab=1.6)
  title(ylab="ACF", line=2, cex.lab=1.6)
  abline(h=0, lwd = 3)
  abline(h=qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
  abline(h=-qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
}
dev.off()

##biweekly##
data <- read.csv("biweekly_278.csv", header = TRUE)
stations <- colnames(data)
m <- dim(data)[2]
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/ACF_plot_278_bw.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:m){
  col_data <- na.omit(data[,i])
  emp_acf <- acf(col_data,lag.max = 52 ,ci = 0.99,plot = FALSE)
  plot(emp_acf$lag,emp_acf$acf,type = "h",ylim = c(-0.1,1), main = stations[i], lwd = 3, las = 1,
       ylab = "", xlab = "", cex.main = 1.6, cex.axis = 1.3)
  title(xlab="Lag", line=1, cex.lab=1.6)
  title(ylab="ACF", line=2, cex.lab=1.6)
  abline(h=0, lwd = 3)
  abline(h=qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
  abline(h=-qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
}
dev.off()

##monthly##
data <- read.csv("monthly_278.csv", header = TRUE)
stations <- colnames(data)
m <- dim(data)[2]
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/ACF_plot_278_m.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:m){
  col_data <- na.omit(data[,i])
  emp_acf <- acf(col_data,lag.max = 24,ci = 0.99,plot = FALSE)
  plot(emp_acf$lag,emp_acf$acf,type = "h",ylim = c(-0.1,1), main = stations[i], lwd = 3, las = 1,
       ylab = "", xlab = "", cex.main = 1.6, cex.axis = 1.3)
  title(xlab="Lag", line=1, cex.lab=1.6)
  title(ylab="ACF", line=2, cex.lab=1.6)
  abline(h=0, lwd = 3)
  abline(h=qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
  abline(h=-qnorm((1 + 0.99)/2)/sqrt(emp_acf$n.used), col = "blue", lty = 2, lwd = 3)
}
dev.off()

###HEAVY TAIL DISTRIBUTION###
library(ExtremeRisks)
library(ReIns)
library(laeken)
data <- read.csv("monthly_278.csv", header = TRUE)
data <- sapply(data[,], as.numeric)

##exponential QQ plot##
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/EXP_QQ_278.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:278){
	vector <- as.vector(data[,i])
	ExpQQ(vector, plot = TRUE, main = colnames(data)[i])
}
dev.off()
	
##mean excess plot##
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/Mean_Excess_plot_278.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:278){
	vector <- as.vector(data[,i])
	meanExcessPlot(vector, main = colnames(data)[i])
}
dev.off()

###PARETO-TYPE DISTRIBUTION### 
##pareto qq##
data <- read.csv("monthly_aggregate_282_381_160.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/PARETO_QQ_160.pdf")
par(mai=c(.60,.65,.35,.1), mgp=c(2.1,.6,0))
for(i in 1:160){
	vector <- as.vector(data[,i])
	ParetoQQ(vector, plot = TRUE, main = colnames(data)[i])
}
dev.off()


###FIRST STEP OF MES ESTIMATION - EV INDEX ESTIMATION###
##estimation through hill estimator##
data <- read.csv("monthly_aggregate_282_381_160.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/HILL_tailindex_160.pdf")
for(j in 1:160){
	vector <- as.vector(data[, j])
	kk <- seq(5,45, by=1)   #largest order statistics to use 
	nk <- length(kk)
	SPTindexHill <- array(0, c(nk,3))

	for(i in 1:nk){
  	#estimates the tail index by the Hill estimator
	  temp <- HTailIndex(vector, kk[i], TRUE, varType = "asym-Ind")
	  SPTindexHill[i, ] <- c(temp$CIgamHat[1], temp$gammaHat, temp$CIgamHat[2])
	} 
	plot(kk, SPTindexHill[,2], type="l", ylim=c(-0.6,1), xlim=c(0,45), lwd=3, ylab="EV index", 
     xlab="k", cex.lab=1.3, cex.axis=1.3, main=colnames(data)[j], las=1)
	minor_tick_positions <- seq(0, 45, by = 1)
	axis(1, at = minor_tick_positions, labels = FALSE, tck = -0.015)
	abline(h=0, col=2, lty=4, lwd=2)
	abline(h=0.5, col=2, lty=4, lwd=2)
	#abline(v = list_ranges[j,1], col = "red", lty = 2)  #to check after the stability of the range
	#abline(v = list_ranges[j,2], col = "red", lty = 2)
	}
dev.off()

##regression approach - linear ##
library(ExtremeRisks) 
vect <- c(c(22,27),c(17,22),c(25,30), c(30,35), c(21,26), c(26,31), c(22,27), c(20,25), c(21,26), c(21,26),
c(18,23), c(24,29), c(22,27), c(18,23), c(21,26), c(30,35), c(30,35), c(22,27), c(27,32), c(29,34),
c(27,32), c(30,35), c(25,30), c(24,30), c(22,27), c(17,22), c(21,26), c(30,35), c(30,35), c(30,35),
c(30,35), c(25,30), c(25,30), c(19,24), c(28,33), c(20,25), c(35,40), c(25,30), c(26,31), c(25,30),
c(20,25), c(30,35), c(23,28), c(21,26), c(35,40), c(35,40), c(22,27), c(18,23), c(31,36), c(28,33),
c(20,25), c(20,25), c(35,40), c(30,35), c(21,26), c(20,25), c(30,35), c(35,40), c(33,38), c(34,39),
c(33,38), c(28,33), c(25,30), c(30,35), c(30,35), c(25,30), c(30,35), c(30,35), c(30,35), c(23,28),
c(30,35), c(26,31), c(15,20), c(18,23), c(25,30), c(17,22), c(28,33), c(30,35), c(25,30), c(32,37),
c(23,28), c(28,33), c(28,33), c(20,25), c(18,23), c(20,25), c(26,31), c(20,25), c(28,33), c(17,22),
c(19,24), c(20,25), c(22,27), c(28,33), c(23,28), c(23,28), c(28,33), c(20,25), c(22,27), c(28,33),
c(30,35), c(33,38), c(27,32), c(30,35), c(30,35), c(16,21), c(27,32), c(22,27), c(35,40), c(25,30),
c(25,30), c(35,40), c(30,35), c(20,25), c(20,25), c(33,38), c(17,22), c(20,25), c(30,35), c(30,35),
c(32,37), c(18,23), c(20,25), c(15,20), c(23,28), c(32,37), c(30,35), c(32,37), c(33,38), c(24,29),
c(26,31), c(25,30), c(22,27), c(30,35), c(23,28), c(23,28), c(28,33), c(35,40), c(28,33), c(30,35),
c(27,32), c(22,27), c(22,27), c(21,26), c(20,25), c(30,35), c(31,36), c(25,30), c(20,25), c(20,25),
c(17,22), c(31,36), c(25,30), c(30,35), c(15,20), c(30,35), c(23,28), c(30,35), c(23,28), c(22,27))
list_ranges <-  matrix(vect, ncol = 2, byrow = TRUE)  #intervals of lenght 5 chosen for the MES estimation

#seasons
data <- read.csv("monthly_aggregate_282_381_160_seasons_dummy_eautunno.csv", header = TRUE)
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/EV_index_seasons.pdf", width=20, height=6)
for(j in 1:160){
	vector <- as.vector(data[, j])
	kk <- seq(8, 41, by=1) 
	nk <- length(kk)
	SPTindexHill <- array(0, c(nk,1))
	tk <-  array(0, c(nk,1))
	sorted_vector <- sort(vector, decreasing = TRUE)

	for(i in 1:nk){
  	# estimates the tail index by the Hill estimator to compare with regression EV index estimates
	  #tk[i] <- sorted_vector[kk[i]]
	  temp <- HTailIndex(vector, kk[i], TRUE, varType = "asym-Ind")
	  SPTindexHill[i, ] <-  temp$gammaHat
	}

	matrix_data <- data[, j]
	matrix_data <- as.data.frame(matrix_data)
	
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML_0 <- array(0, c(nt,1))  #estimates the tail index by the ML estimator without covariate (also used as a comparison) 

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(vector,  matrix_data , threshold = t[i], type="GP", time.units = "months"),
			error = function(e) {
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML_0[i] <- c(NA)
			next
		} else {
		indexML_0[i] <- tryCatch(
			temp$results$par[2],
			error = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}

	
	matrix_data <- cbind(data[, j], data[, "autunno_prim"])  #estimates the tail index by the ML estimator with autumn_spring
	matrix_data <- as.data.frame(matrix_data)
	
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML <- array(0, c(nt,2))  

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML[i, ] <- c(NA, NA)
			next
		} else {
		indexML[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"])    #estimates the tail index by the ML estimator with autumn_winter
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML1 <- array(0, c(nt,2))

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML1[i, ] <- c(NA, NA)
			next
		} else {
		indexML1[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno_est"])    #estimates the tail index by the ML estimator with autumn_summer
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML2 <- array(0, c(nt,2))
	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML2[i, ] <- c(NA, NA)
			next
		} else {
		indexML2[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}	

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno"])  #estimates the tail index by the ML estimator with autumn
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2) 
	nt <- length(t)
	indexML3 <- array(0, c(nt,2)) 

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML3[i, ] <- c(NA, NA)
			next
		} else {
		indexML3[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}	
		
	
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	plot(k1, indexML[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Spring"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML1[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Winter"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML2[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Summer"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML3[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	mtext(colnames(data)[j],                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
} #for j 
dev.off()

#temperature, wind speed, relative humidity, solar radiation
data <- read.csv("monthly_aggregate_282_381_160_seasons_dummy_eautunno.csv", header = TRUE) 
data_t<- read.csv("temperatura_monthly_aggregate_282_381_132.csv", header = TRUE)
data_v <- read.csv("vento_monthly_aggregate_282_381_72.csv", header = TRUE) 
data_r<- read.csv("radiazione_monthly_aggregate_282_381_51.csv", header = TRUE)               
data_u<- read.csv("umidità_monthly_aggregate_282_381_91.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/EV_index_atmospheric.pdf", width=20, height=6)

for(j in 1:160){
	vector <- as.vector(data[, j])
	kk <- seq(8, 41, by=1) 
	nk <- length(kk)
	SPTindexHill <- array(0, c(nk,1))
	tk <-  array(0, c(nk,1))
	sorted_vector <- sort(vector, decreasing = TRUE)

	for(i in 1:nk){
  	# estimates the tail index by the Hill estimator
	  #tk[i] <- sorted_vector[kk[i]]
	  temp <- HTailIndex(vector, kk[i], TRUE, varType = "asym-Ind")
	  SPTindexHill[i, ] <-  temp$gammaHat
	}

	matrix_data <- data[, j]
	matrix_data <- as.data.frame(matrix_data)

	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(8, 42, by = 2)
	nt <- length(t)
	indexML_0 <- array(0, c(nt,1))

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(vector,  matrix_data , threshold = t[i], type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML_0[i] <- c(NA)
			next
		} else {
		indexML_0[i] <- tryCatch(
			temp$results$par[2],
			error = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}


		indexML1_t <- array(0, c(nt,1))   #estimates for temperature at different levels (1.5,2.5, 5, 7.5, 10, 12.5)
		indexML2_t <- array(0, c(nt,1))	
		indexML3_t <- array(0, c(nt,1))
		indexML4_t <- array(0, c(nt,1))
		indexML5_t <- array(0, c(nt,1))	
		if (colnames(data)[j] %in% colnames(data_t)) {
			matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_t[y ] <- c(NA)
			indexML2_t[y ] <- c(NA)
			indexML3_t[y ] <- c(NA)
			indexML4_t[y ] <- c(NA)
			indexML5_t[y ] <- c(NA)
			indexML6_t[y ] <- c(NA)
			next
		} else {
		indexML1_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*1.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*2.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML3_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML4_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*7.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML5_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*10,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML6_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*12.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		}}} else {
		indexML1_t <- array(NA, c(nt,1))
		indexML2_t <- array(NA, c(nt,1))
		indexML3_t <- array(NA, c(nt,1))
		indexML4_t <- array(NA, c(nt,1))
		indexML5_t <- array(NA, c(nt,1))
		indexML6_t <- array(NA, c(nt,1))
		}

		indexML1_u <- array(0, c(nt,1))  #estimates for humidity at different levels 75%, 80%, 85%, 90%
		indexML2_u <- array(0, c(nt,1))	
		indexML3_u <- array(0, c(nt,1))
		indexML4_u <- array(0, c(nt,1))
		if (colnames(data)[j] %in% colnames(data_u)) {
			matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_u[y ] <- c(NA)
			indexML2_u[y ] <- c(NA)
			next
		} else {
		indexML1_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*75,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*80,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML3_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*85,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML4_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*90,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)

		}}} else {
		indexML1_u <- array(NA, c(nt,1))
		indexML2_u <- array(NA, c(nt,1))
		indexML3_u <- array(NA, c(nt,1))
		indexML4_u <- array(NA, c(nt,1))
		}


		indexML1_v <- array(0, c(nt,1))
		indexML2_v <- array(0, c(nt,1))
		indexML3_v <- array(0, c(nt,1))
		indexML4_v <- array(0, c(nt,1))
	
		if (colnames(data)[j] %in% colnames(data_v)) {
			matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_v[y ] <- c(NA)  #estimates for the wind speed at different levels 10m/s, 15m/s, 25m/s, 30m/s
			indexML2_v[y ] <- c(NA)
			indexML3_v[y ] <- c(NA)
			indexML4_v[y ] <- c(NA)
			next
		} else {
		indexML1_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*10,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*15,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*25,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*30,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)

		}}} else {
		indexML1_v <- array(NA, c(nt,1))
		indexML2_v <- array(NA, c(nt,1))
		indexML3_v <- array(NA, c(nt,1))
		indexML4_v <- array(NA, c(nt,1))
		}

		indexML1_r <- array(0, c(nt,1))   #estimates for solar radiation at different levels 25 W/m2, 50 W/m2
		indexML2_r <- array(0, c(nt,1))		
		if (colnames(data)[j] %in% colnames(data_r)) {
			matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)
		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_r[y ] <- c(NA)
			indexML2_r[y ] <- c(NA)
			next
		} else {
		indexML1_r[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*25,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_r[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*50,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}} else {
		indexML1_r <- array(NA, c(nt,1))
		indexML2_r <- array(NA, c(nt,1))
		}

	
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	if (colnames(data)[j] %in% colnames(data_t)) {
	plot(k1, indexML1_t, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="EV index",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="TEMPERATURE", las=1)
	lines(k1, indexML2_t, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_t, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_t, lty=1, lwd=2, col="red")
	lines(k1, indexML5_t, lty=1, lwd=2, col="grey")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "grey", "green", "magenta"), lwd=c(2,2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"temper_1.5c"),
       expression(hat(gamma)[n]~"temper_2.5c"), expression(hat(gamma)[n]~"temper_5c"),
	 expression(hat(gamma)[n]~"temper_10c"), expression(hat(gamma)[n]~"temper_12.5c"),
	expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill") ))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for temperature")
	}


	if (colnames(data)[j] %in% colnames(data_v)) {
	plot(k1, indexML1_v, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="WIND SPEED", las=1)
	lines(k1, indexML2_v, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_v, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_v, lty=1, lwd=2, col="red")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "green", "magenta"), lwd=c(2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"wind_10m/s"),
       expression(hat(gamma)[n]~"wind_15m/s"), expression(hat(gamma)[n]~"wind_25m/s"), 
	expression(hat(gamma)[n]~"wind_30m/s"), expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")
	))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for wind")
	}


	if (colnames(data)[j] %in% colnames(data_u)) {
	plot(k1, indexML1_u, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="REL HUMIDITY", las=1)
	lines(k1, indexML2_u, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_u, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_u, lty=1, lwd=2, col="red")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
legend("top", lty=c(1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "green", "magenta"), lwd=c(2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"humid_75%"),
       expression(hat(gamma)[n]~"humid_80%"), expression(hat(gamma)[n]~"humid_85%"), expression(hat(gamma)[n]~"humid_90%"),
	expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")
	))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for humidity")
	}

	if (colnames(data)[j] %in% colnames(data_r)) {
	plot(k1, indexML1_r, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="SOLAR RADIATION", las=1)
	lines(k1, indexML2_r, lty=1, lwd=2, col="blue")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
legend("top", lty=c(1,1,2,1), col=c("black", "blue", "green", "magenta"), lwd=c(2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"rad_25W/m2"),
       expression(hat(gamma)[n]~"rad_50W/m2"), expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for solar radiat")
	}

	mtext(colnames(data)[j],                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
} #for j 
dev.off()



##regression approach - exponential ##
#seasons
data <- read.csv("monthly_aggregate_282_381_160_seasons_dummy_eautunno.csv", header = TRUE)
pdf("C:/Users/39346/OneDrive/Desktop/THESIS/EV_index_seasons.pdf", width=20, height=6)
for(j in 1:160){
	vector <- as.vector(data[, j])
	kk <- seq(8, 41, by=1) 
	nk <- length(kk)
	SPTindexHill <- array(0, c(nk,1))
	tk <-  array(0, c(nk,1))
	sorted_vector <- sort(vector, decreasing = TRUE)

	for(i in 1:nk){
  	# estimates the tail index by the Hill estimator
	  #tk[i] <- sorted_vector[kk[i]]
	  temp <- HTailIndex(vector, kk[i], TRUE, varType = "asym-Ind")
	  SPTindexHill[i, ] <-  temp$gammaHat
	}

	matrix_data <- data[, j]
	matrix_data <- as.data.frame(matrix_data)
	
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML_0 <- array(0, c(nt,1))  #estimates the tail index by the ML estimator without covariate

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(vector,  matrix_data , threshold = t[i], type="GP", time.units = "months"),
			error = function(e) {
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML_0[i] <- c(NA)
			next
		} else {
		indexML_0[i] <- tryCatch(
			temp$results$par[2],
			error = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}

	
	matrix_data <- cbind(data[, j], data[, "autunno_prim"])  
	matrix_data <- as.data.frame(matrix_data)
	
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML <- array(0, c(nt,2)) #estimates the tail index by the ML estimator with autumn_spring 

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML[i, ] <- c(NA, NA)
			next
		} else {
		indexML[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"])  
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML1 <- array(0, c(nt,2)) #estimates the tail index by the ML estimator with autumn_winter

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML1[i, ] <- c(NA, NA)
			next
		} else {
		indexML1[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno_est"])  
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2)
	nt <- length(t)
	indexML2 <- array(0, c(nt,2)) #estimates the tail index by the ML estimator with autumn_summer
	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML2[i, ] <- c(NA, NA)
			next
		} else {
		indexML2[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}	

	vector <- as.vector(data[, j])
	matrix_data <- cbind(data[, j], data[, "autunno"])  
	matrix_data <- as.data.frame(matrix_data)
	sorted_vector <- sort(vector, decreasing = TRUE)
	#max_value <- sorted_vector[10]
	#min_value <- sorted_vector[40]
	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(10, 42, by = 2) 
	nt <- length(t)
	indexML3 <- array(0, c(nt,2)) #estimates the tail index by the ML estimator with autumn

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[i], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML3[i, ] <- c(NA, NA)
			next
		} else {
		indexML3[i,] <- tryCatch(
			c(temp$results$par[2], temp$results$par[2]+temp$results$par[3]),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA, NA))  # Return a default value or 'NA'
      		}
    		)}}	
		
	
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	plot(k1, indexML[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Spring"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML1[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Winter"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML2[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn_Summer"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	plot(k1, indexML3[,2], type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="autumn_spring" , las=1)
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,2,1), col=c("black", "green", "magenta"), lwd=c(2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"Autumn"),
       expression(hat(gamma)[n]~"ML_no_cov",  expression(hat(gamma)[n]~"Hill")))

	mtext(colnames(data)[j],                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
} #for j 
dev.off()

#temperature, wind speed, relative humidity, solar radiation
data <- read.csv("monthly_aggregate_282_381_160_seasons_dummy_eautunno.csv", header = TRUE) 
data_t<- read.csv("temperatura_monthly_aggregate_282_381_132.csv", header = TRUE)
data_v <- read.csv("vento_monthly_aggregate_282_381_72.csv", header = TRUE) 
data_r<- read.csv("radiazione_monthly_aggregate_282_381_51.csv", header = TRUE)               
data_u<- read.csv("umidità_monthly_aggregate_282_381_91.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/EV_index_atmospheric.pdf", width=20, height=6)

for(j in 1:160){
	vector <- as.vector(data[, j])
	kk <- seq(8, 41, by=1) 
	nk <- length(kk)
	SPTindexHill <- array(0, c(nk,1))
	tk <-  array(0, c(nk,1))
	sorted_vector <- sort(vector, decreasing = TRUE)

	for(i in 1:nk){
  	# estimates the tail index by the Hill estimator
	  #tk[i] <- sorted_vector[kk[i]]
	  temp <- HTailIndex(vector, kk[i], TRUE, varType = "asym-Ind")
	  SPTindexHill[i, ] <-  temp$gammaHat
	}

	matrix_data <- data[, j]
	matrix_data <- as.data.frame(matrix_data)

	t <- sorted_vector[seq(10, 42, by = 2)] 
	k1 <- seq(8, 42, by = 2)
	nt <- length(t)
	indexML_0 <- array(0, c(nt,1))

	for(i in 1:nt){
  		temp <- tryCatch(
			fevd(vector,  matrix_data , threshold = t[i], type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)	
	
		if (is.null(temp)) {
			indexML_0[i] <- c(NA)
			next
		} else {
		indexML_0[i] <- tryCatch(
			temp$results$par[2],
			error = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}


		indexML1_t <- array(0, c(nt,1))
		indexML2_t <- array(0, c(nt,1))	
		indexML3_t <- array(0, c(nt,1))
		indexML4_t <- array(0, c(nt,1))
		indexML5_t <- array(0, c(nt,1))
		indexML6_t <- array(0, c(nt,1))	
		if (colnames(data)[j] %in% colnames(data_t)) {
			matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_t[y ] <- c(NA)
			indexML2_t[y ] <- c(NA)
			indexML3_t[y ] <- c(NA)
			indexML4_t[y ] <- c(NA)
			indexML5_t[y ] <- c(NA)
			indexML6_t[y ] <- c(NA)
			next
		} else {
		indexML1_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*1.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*2.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML3_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML4_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*7.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML5_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*10,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML6_t[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*12.5,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		}}} else {
		indexML1_t <- array(NA, c(nt,1))
		indexML2_t <- array(NA, c(nt,1))
		indexML3_t <- array(NA, c(nt,1))
		indexML4_t <- array(NA, c(nt,1))
		indexML5_t <- array(NA, c(nt,1))
		indexML6_t <- array(NA, c(nt,1))
		}

		indexML1_u <- array(0, c(nt,1))
		indexML2_u <- array(0, c(nt,1))	
		indexML3_u <- array(0, c(nt,1))
		indexML4_u <- array(0, c(nt,1))
		if (colnames(data)[j] %in% colnames(data_u)) {
			matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_u[y ] <- c(NA)
			indexML2_u[y ] <- c(NA)
			next
		} else {
		indexML1_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*75,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*80,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML3_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*85,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML4_u[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*90,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)

		}}} else {
		indexML1_u <- array(NA, c(nt,1))
		indexML2_u <- array(NA, c(nt,1))
		indexML3_u <- array(NA, c(nt,1))
		indexML4_u <- array(NA, c(nt,1))
		}


		indexML1_v <- array(0, c(nt,1))
		indexML2_v <- array(0, c(nt,1))
		indexML3_v <- array(0, c(nt,1))
		indexML4_v <- array(0, c(nt,1))
	
		if (colnames(data)[j] %in% colnames(data_v)) {
			matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)

		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_v[y ] <- c(NA)
			indexML2_v[y ] <- c(NA)
			indexML3_v[y ] <- c(NA)
			indexML4_v[y ] <- c(NA)
			next
		} else {
		indexML1_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*10,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*15,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*25,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_v[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*30,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)

		}}} else {
		indexML1_v <- array(NA, c(nt,1))
		indexML2_v <- array(NA, c(nt,1))
		indexML3_v <- array(NA, c(nt,1))
		indexML4_v <- array(NA, c(nt,1))
		}

		indexML1_r <- array(0, c(nt,1))
		indexML2_r <- array(0, c(nt,1))		
		if (colnames(data)[j] %in% colnames(data_r)) {
			matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
			matrix_data <- as.data.frame(matrix_data)
		
		for(y in 1:nt){
  		temp <- tryCatch(
			fevd(V1,  matrix_data , threshold = t[y], shape.fun = ~exp(V2), type="GP", time.units = "months"),
			error = function(e) {
				return(NULL)  # Return a default value or 'NA'
      		}
    		)
		if (is.null(temp)) {
			indexML1_r[y ] <- c(NA)
			indexML2_r[y ] <- c(NA)
			next
		} else {
		indexML1_r[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*25,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)
		indexML2_r[y] <- tryCatch(
			temp$results$par[2]+temp$results$par[3]*50,
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		},
			warning = function(e) {
				 #ML_title <- "ML_extRemes_PROBLEM"
        			return(c(NA))  # Return a default value or 'NA'
      		}
    		)}}} else {
		indexML1_r <- array(NA, c(nt,1))
		indexML2_r <- array(NA, c(nt,1))
		}

	
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	if (colnames(data)[j] %in% colnames(data_t)) {
	plot(k1, indexML1_t, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="EV index",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="TEMPERATURE", las=1)
	lines(k1, indexML2_t, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_t, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_t, lty=1, lwd=2, col="red")
	lines(k1, indexML5_t, lty=1, lwd=2, col="grey")
	lines(k1, indexML6_t, lty=1, lwd=2, col="aquamarine")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,1,1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "grey", "aquamarine", "green", "magenta"), lwd=c(2,2,2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"temper_1.5c"),
       expression(hat(gamma)[n]~"temper_2.5c"), expression(hat(gamma)[n]~"temper_5c"),
	expression(hat(gamma)[n]~"temper_7.5c"),expression(hat(gamma)[n]~"temper_10c"), 
	expression(hat(gamma)[n]~"temper_12.5c"),
	expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill") ))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for temperature")
	}


	if (colnames(data)[j] %in% colnames(data_v)) {
	plot(k1, indexML1_v, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="WIND SPEED", las=1)
	lines(k1, indexML2_v, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_v, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_v, lty=1, lwd=2, col="red")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
	legend("top", lty=c(1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "green", "magenta"), lwd=c(2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"wind_10m/s"),
       expression(hat(gamma)[n]~"wind_15m/s"), expression(hat(gamma)[n]~"wind_25m/s"), 
	expression(hat(gamma)[n]~"wind_30m/s"), expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")
	))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for wind")
	}


	if (colnames(data)[j] %in% colnames(data_u)) {
	plot(k1, indexML1_u, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="REL HUMIDITY", las=1)
	lines(k1, indexML2_u, lty=1, lwd=2, col="blue")
	lines(k1, indexML3_u, lty=1, lwd=2, col="yellow")
	lines(k1, indexML4_u, lty=1, lwd=2, col="red")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
legend("top", lty=c(1,1,1,1,2,1), col=c("black", "blue", "yellow", "red", "green", "magenta"), lwd=c(2,2,2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"humid_75%"),
       expression(hat(gamma)[n]~"humid_80%"), expression(hat(gamma)[n]~"humid_85%"), expression(hat(gamma)[n]~"humid_90%"),
	expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")
	))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for humidity")
	}

	if (colnames(data)[j] %in% colnames(data_r)) {
	plot(k1, indexML1_r, type="l", ylim=c(-0.6,1), xlim=c(8,43), lwd=2, ylab="",
     xlab="k", cex.lab=1.3, cex.axis=1.3, main="SOLAR RADIATION", las=1)
	lines(k1, indexML2_r, lty=1, lwd=2, col="blue")
	lines(k1, indexML_0, lty=2, lwd=1, col="green")
	lines(kk, SPTindexHill, lty=1, lwd=2, col="magenta")
	abline(h=0, col=2, lty=4, lwd=2)
	abline(v = list_ranges[j,1], col = "red", lty = 2)
	abline(v = list_ranges[j,2], col = "red", lty = 2)
	axis(2, at=seq(-0.6, 1, by=0.1), labels=rep("", length(seq(-0.6, 1, by=0.1))))
legend("top", lty=c(1,1,2,1), col=c("black", "blue", "green", "magenta"), lwd=c(2,2,1,2),
       bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"rad_25W/m2"),
       expression(hat(gamma)[n]~"rad_50W/m2"), expression(hat(gamma)[n]~"ML_no_cov"), expression(hat(gamma)[n]~"Hill")))
	} else {
	  plot(1, type = "n", ylim=c(-0.6,1), xlab = "", ylab = "", main = "No data available for solar radiat")
	}

	mtext(colnames(data)[j],                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
} #for j 
dev.off()


###SECOND STEP OF MES ESTIMATION - UNWEIGHTED CASE ###
##MES with Hill vs with Coov - unweighted##
coov <- read.csv("coovariate choice.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/MES_HillvsCoov_unweighted.pdf", width=20, height=6)
kk <- seq(8,41, by=1)  #by 1 may be too large 
nk <- length(kk)
MES_1 <- array(0, c(nk,159,160)) #hill MES
MES_1_avg <- array(0, c(nk,160))
MES_2 <- array(0, c(nk,160))
MES_3 <- array(0, c(nk,159,160)) #coov MES
MES_3_avg <- array(0, c(nk,160))
MES_4 <- array(0, c(nk,160))

for(j in 1:160){
#print(j)

#first step

	vector_x <- as.vector(data[,j])
	sorted_vector_x <- sort(vector_x, decreasing = TRUE)
	v_ML <- list()
	for (k1 in list_ranges[j,1]:list_ranges[j,2]){
		temp <- HTailIndex(vector_x, k1, TRUE, varType = "asym-Ind")
		v_ML <- append(v_ML, temp$gammaHat)
	}
	v_ML <- unlist(v_ML)
	mean_tail_ML <- mean(v_ML)

	#1 autumn_spring
	if (!is.na(coov[j, 2]) && coov[j, 2] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_prim"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#2 autumn_winter
	if (!is.na(coov[j, 3]) && coov[j, 3] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#3 autumn_spring
	if (!is.na(coov[j, 4]) && coov[j, 4] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_est"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#4 autumn
	if (!is.na(coov[j, 5]) && coov[j, 5] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	
	} 

	#5 temperature 1.5
	if (!is.na(coov[j, 6]) && coov[j, 6] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*1.5)
		
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#6 temperature 2.5
	if (!is.na(coov[j, 7]) && coov[j, 7] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*2.5)
		#print(temp$results$par[2]+temp$results$par[3]*2.5)
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#7 temperature 5
	if (!is.na(coov[j, 8]) && coov[j, 8] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#8 temperature 7.5
	if (!is.na(coov[j, 9]) && coov[j, 9] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*7.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#9 temperature 10
	if (!is.na(coov[j, 10]) && coov[j, 10] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]])
	matrix_data <- as.data.frame(matrix_data) 
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#10 temperature 12.5
	if (!is.na(coov[j, 11]) && coov[j, 11] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*12.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#11 wind speed 10
	if (!is.na(coov[j, 12]) && coov[j, 12] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#12 wind speed 15
	if (!is.na(coov[j, 13]) && coov[j, 13] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*15)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#13 wind speed 25
	if (!is.na(coov[j, 14]) && coov[j, 14] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#14 wind speed 30
	if (!is.na(coov[j, 15]) && coov[j, 15] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*30)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#15 solar radiation 25
	if (!is.na(coov[j, 16]) && coov[j, 16] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#16 solar radiation 50
	if (!is.na(coov[j, 17]) && coov[j, 17] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*50)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#17 humidity 75
	if (!is.na(coov[j, 18]) && coov[j, 18] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*75)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	
	#18 humidity 80
	if (!is.na(coov[j, 19]) && coov[j, 19] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*80)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#19 humidity 85
	if (!is.na(coov[j, 20]) && coov[j, 20] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*85)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#20 humidity 90
	if (!is.na(coov[j, 21]) && coov[j, 21] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*90)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	# case of no covariate
	if (all(is.na(coov[j, 2:21]))) {
	mean_tail_ML_c <- NA
	}

#}

#second part	
	Y <- data[,-j]

	#paiwise HILL&ML
	for(z in 1:159){  
		vector <- as.vector(Y[, z])		
		sorted_vector <- sort(vector, decreasing = TRUE)
		for(i in 1:nk){
			thresh <- sorted_vector[kk[i]]
			indices_vector <- which(vector > thresh)
			selected_values <- vector_x[indices_vector]
			sum_s <- sum(selected_values)
			MES_1[i, z, j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*(sum_s) 			#if p=1/n
			if (!is.na(mean_tail_ML_c)) {
			MES_3[i, z, j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*(sum_s)} else {MES_3[i, z, j] <- NA}
		}} #}
	MES_1_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = mean)
	if (!is.na(mean_tail_ML_c)) {
	MES_3_avg[,j] <- apply(MES_3[,,j], MARGIN = 1, FUN = mean)} else {MES_3_avg[,j] <- NA}
	
	#average 
	Y_mean <- rowMeans(Y)  
	sorted_vector <- sort(Y_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_2[i,j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*sum_s 			#if p=1/n
		if (!is.na(mean_tail_ML_c)) {
		MES_4[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s}  else {MES_4[i,j] <- NA}
}

	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall(mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="PAIRWISE MES HILL without cov", las=1)
	lines(kk, MES_1_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES HILL without cov"),
       expression("average of pairwise MES HILL without cov")))

	
	plot(kk, MES_2[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",    
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES with AVERAGE as cond event", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      legend=expression("MES HILL with AVERAGE as cond event without cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML , col = "red") 
	text(44, 27, "highest obs:" , col = "red")
	text(44, 26, round(sorted_vector_x[1],1) , col = "red")
	text(44, 25, round(sorted_vector_x[2],1) , col = "red")
	text(44, 24, "tail thresh:" , col = "red")
	text(44, 23, round(sorted_vector_x[k1],1) , col = "red")

	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_3[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML with cov", las=1)
	lines(kk, MES_3_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML with cov"),
       expression("average of pairwise MES ML with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_4[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event without cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event with cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML_c , col = "red")

	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}

	mtext(colnames(data)[j], side = 3, line = - 2, outer = TRUE)
}
dev.off()

#get the difference between prediction and real data to check if overfitting or underfitting or both
Hill1 <- array(0, c(160))
Hill2 <- array(0, c(160))
Coov1 <- array(0, c(160))
Coov2 <- array(0, c(160))

for (i in 1:160) {
	Hill1[j] <- (max(MES_1_avg[,j])-max(data[,j]))
	Hill2[j] <- (max(MES_2[,j])-max(data[,j]))
	Coov1[j] <- (max(MES_3_avg[,j])-max(data[,j]))
	Coov2[j] <- (max(MES_4[,j])-max(data[,j]))
}
#get average absolute error
for (i in 1:160) {
	Hill1[j] <- abs(max(MES_1_avg[,j])-max(data[,j]))
	Hill2[j] <- abs(max(MES_2[,j])-max(data[,j]))
	Coov1[j] <- abs(max(MES_3_avg[,j])-max(data[,j]))
	Coov2[j] <- abs(max(MES_4[,j])-max(data[,j]))
}
mH1 <- mean(Hill1)
mH2 <- mean(Hill1)
mC1 <- mean(Hill1)
mC2 <- mean(Hill1)
#get mean squared error
for (i in 1:160) {
	Hill1[j] <- (abs(max(MES_1_avg[,j])-max(data[,j])))^2
	Hill2[j] <- (abs(max(MES_2[,j])-max(data[,j])))^2
	Coov1[j] <- (abs(max(MES_3_avg[,j])-max(data[,j])))^2
	Coov2[j] <- (abs(max(MES_4[,j])-max(data[,j])))^2
}
mH1 <- mean(Hill1)
mH2 <- mean(Hill1)
mC1 <- mean(Hill1)
mC2 <- mean(Hill1)

	

### DERIVATION OF WEIGHT - spatial correlation ###
##spatial model fitting##
library(RandomFields)
library(CompRandFld)
main_folder <- 'C:/Users/39346/OneDrive/Desktop/THESIS'
setwd(main_folder)
data <- read.csv("scaled_282_382_160.csv", header = TRUE)  #use scaled data to get gaussian distribution 
data <- sapply(data[,], as.numeric)
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
x <- as.vector(coo[,7])
y <- as.vector(coo[,6])
coords<-cbind(x,y)

corrmodel <- "matern"
mean <- 0
sill <- 1
nugget <- 0
scale <- 5
smooth <- 8
# Fixed parameters
fixed<-list(mean=mean)
 
# Starting value for the estimated parameters
start<-list(scale=scale, smooth=smooth, sill=sill, nugget=nugget)

fit1 <- FitComposite(data, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, replicates=100, corrmodel=corrmodel, distance = "Geod", likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit1)

vario1 <- EVariogram(data, coords, distance="Geod", replicates=100)
print(vario1)
par(mfrow=c(1,2))
Covariogram(fit1, show.cov=TRUE, show.range=TRUE, show.vario=TRUE, vario=vario1,pch=20)
cov1 <- Covariogram(fit1,answer.cov=TRUE,answer.vario=TRUE, vario=vario1,pch=20)
print(cov1)
par(mfrow=c(1,1))
xyg <- seq(0,1000, length=250)
plot(xyg, Covariogram(fit1, xyg, answer.cov=TRUE)$covariance, type="l", ylim=c(0,1), xlim=c(0,200), main="Estimated covariance function", xlab="lag", ylab="Cov")

par(mfrow=c(1,1)
xyg <- seq(0,1000, length=200)
covariance <- Covariogram(fit1, xyg, answer.cov=TRUE)$covariance + 0.46456 
plot(xyg, covariance , type="l", ylim=c(0,1), xlim=c(0,200), main="Estimated covariance function", xlab="lag", ylab="Cov")

##derivation of weights ##
#for hill#
dist <- read.csv("dist_159_R.csv", header = TRUE) 
weights_n <- array(0, c(159,160))
for (i in 1:160) {
	w <- Covariogram(fit1, dist[,i], answer.cov = TRUE)$covariance
	weights_n[,i] <- w/(sum(w))	
} 
write.csv(weights_n, file = "norm_weights_160.csv", row.names = FALSE)

#for coov#
dist <- read.csv("dist_124_R.csv", header = TRUE) 
weights_n <- array(0, c(159,160))
for (i in 1:160) {
	w <- Covariogram(fit1, dist[,i], answer.cov = TRUE)$covariance
	weights_n[,i] <- w/(sum(w))	
} 
write.csv(weights_n, file = "norm_weights_125.csv", row.names = FALSE)

###MES ESTIMATION - weighted ###
##MES with Hill vs with Coov - weighted##
coov <- read.csv("coovariate choice.csv", header = TRUE)
weights_160 <- read.csv("norm_weights_160.csv", header = TRUE)
weights_125 <- read.csv("norm_weights_125.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/MES_HillvsCoov_weighted.pdf", width=20, height=6)
kk <- seq(8,41, by=1)  #by 1 may be too large 
nk <- length(kk)
MES_1 <- array(0, c(nk,159,160)) #hill MES
MES_1_avg <- array(0, c(nk,160))
MES_2 <- array(0, c(nk,160))
MES_3 <- array(0, c(nk,159,160)) #coov MES
MES_3_avg <- array(0, c(nk,160))
MES_4 <- array(0, c(nk,160))

for(j in 1:160){

#first step

	vector_x <- as.vector(data[,j])
	sorted_vector_x <- sort(vector_x, decreasing = TRUE)
	v_ML <- list()
	for (k1 in list_ranges[j,1]:list_ranges[j,2]){
		temp <- HTailIndex(vector_x, k1, TRUE, varType = "asym-Ind")
		v_ML <- append(v_ML, temp$gammaHat) #called v_ML but it is the Hill EV index 
	}
	v_ML <- unlist(v_ML)
	mean_tail_ML <- mean(v_ML)

	#1 autumn_spring
	if (!is.na(coov[j, 2]) && coov[j, 2] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_prim"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#2 autumn_winter
	if (!is.na(coov[j, 3]) && coov[j, 3] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#3 autumn_spring
	if (!is.na(coov[j, 4]) && coov[j, 4] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_est"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#4 autumn
	if (!is.na(coov[j, 5]) && coov[j, 5] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	
	} 

	#5 temperature 1.5
	if (!is.na(coov[j, 6]) && coov[j, 6] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*1.5)
		
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#6 temperature 2.5
	if (!is.na(coov[j, 7]) && coov[j, 7] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*2.5)
		#print(temp$results$par[2]+temp$results$par[3]*2.5)
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#7 temperature 5
	if (!is.na(coov[j, 8]) && coov[j, 8] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#8 temperature 7.5
	if (!is.na(coov[j, 9]) && coov[j, 9] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*7.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#9 temperature 10
	if (!is.na(coov[j, 10]) && coov[j, 10] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]])
	matrix_data <- as.data.frame(matrix_data) 
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#10 temperature 12.5
	if (!is.na(coov[j, 11]) && coov[j, 11] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*12.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#11 wind speed 10
	if (!is.na(coov[j, 12]) && coov[j, 12] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#12 wind speed 15
	if (!is.na(coov[j, 13]) && coov[j, 13] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*15)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#13 wind speed 25
	if (!is.na(coov[j, 14]) && coov[j, 14] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#14 wind speed 30
	if (!is.na(coov[j, 15]) && coov[j, 15] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*30)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#15 solar radiation 25
	if (!is.na(coov[j, 16]) && coov[j, 16] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#16 solar radiation 50
	if (!is.na(coov[j, 17]) && coov[j, 17] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*50)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#17 humidity 75
	if (!is.na(coov[j, 18]) && coov[j, 18] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*75)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	
	#18 humidity 80
	if (!is.na(coov[j, 19]) && coov[j, 19] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*80)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#19 humidity 85
	if (!is.na(coov[j, 20]) && coov[j, 20] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*85)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#20 humidity 90
	if (!is.na(coov[j, 21]) && coov[j, 21] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*90)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	# case of no covariate
	if (all(is.na(coov[j, 2:21]))) {
	mean_tail_ML_c <- NA
	}

#}

#second part	
	Y <- data[,-j]

	#paiwise HILL&ML
	weights_160 <- weights_160[,j]
	weights_125 <- weights_125[,j]
	for(z in 1:159){  
		vector <- as.vector(Y[, z])		
		sorted_vector <- sort(vector, decreasing = TRUE)
		for(i in 1:nk){
			thresh <- sorted_vector[kk[i]]
			indices_vector <- which(vector > thresh)
			selected_values <- vector_x[indices_vector]
			sum_s <- sum(selected_values)
			MES_1[i, z, j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*(sum_s) 			#if p=1/n
			if (!is.na(mean_tail_ML_c)) {
			MES_3[i, z, j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*(sum_s)} else {MES_3[i, z, j] <- NA}
		}} 
	MES_1_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = function(x) weighted.mean(x, weights_160))
	if (!is.na(mean_tail_ML_c)) {
	MES_3_avg[,j] <- apply(MES_3[,,j], MARGIN = 1, FUN = function(x) weighted.mean(x, weights_125))} else {MES_3_avg[,j] <- NA}
	
	#average weighted 
 	Y_w <- sweep(Y,2, weights_160, `*`)
	Y_w_mean <- rowMeans(Y_w)
	sorted_vector <- sort(Y_w_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_w_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_2[i,j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*sum_s 			#if p=1/n

	if (!is.na(mean_tail_ML_c)) {
	Y_w <- sweep(Y,2, weights_125, `*`)
	Y_w_mean <- rowMeans(Y_w)
	sorted_vector <- sort(Y_w_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_w_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_4[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s}  else {MES_4[i,j] <- NA}}

	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall(mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="PAIRWISE MES HILL - weighted", las=1)
	lines(kk, MES_1_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES HILL without cov"),
       expression("average of pairwise MES HILL without cov")))
	
	plot(kk, MES_2[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",    
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES with AVERAGE as cond event", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      legend=expression("MES HILL with AVERAGE as cond event without cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML , col = "red") 
	text(44, 27, "highest obs:" , col = "red")
	text(44, 26, round(sorted_vector_x[1],1) , col = "red")
	text(44, 25, round(sorted_vector_x[2],1) , col = "red")
	text(44, 24, "tail thresh:" , col = "red")
	text(44, 23, round(sorted_vector_x[k1],1) , col = "red")

	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_3[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML with cov - weighted", las=1)
	lines(kk, MES_3_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML with cov"),
       expression("average of pairwise MES ML with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_4[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event without cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event with cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML_c , col = "red")

	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}

	mtext(colnames(data)[j], side = 3, line = - 2, outer = TRUE)
}
dev.off()

write.csv(MES_1_avg[37,], file = "MES_hill1_weight.csv", row.names = FALSE)
write.csv(MES_2[37,], file = "MES_hill2_weight.csv", row.names = FALSE)
write.csv(MES_3_avg[37,], file = "MES_cov1_weight.csv", row.names = FALSE)
write.csv(MES_4[37,], file = "MES_cov2_weight.csv", row.names = FALSE)

##check the goodness of fit##
Hill1 <- array(0, c(160))
Hill2 <- array(0, c(160))
Coov1 <- array(0, c(160))
Coov2 <- array(0, c(160))
#get the difference between prediction and real data to check if overfitting or underfitting or both
for (i in 1:160) {
	Hill1[j] <- (max(MES_1_avg[,j])-max(data[,j]))
	Hill2[j] <- (max(MES_2[,j])-max(data[,j]))
	Coov1[j] <- (max(MES_3_avg[,j])-max(data[,j]))
	Coov2[j] <- (max(MES_4[,j])-max(data[,j]))
}
#get average absolute error
for (i in 1:160) {
	Hill1[j] <- abs(max(MES_1_avg[,j])-max(data[,j]))
	Hill2[j] <- abs(max(MES_2[,j])-max(data[,j]))
	Coov1[j] <- abs(max(MES_3_avg[,j])-max(data[,j]))
	Coov2[j] <- abs(max(MES_4[,j])-max(data[,j]))
}
mH1 <- mean(Hill1)
mH2 <- mean(Hill1)
mC1 <- mean(Hill1)
mC2 <- mean(Hill1)
#get mean squared error
for (i in 1:160) {
	Hill1[j] <- (abs(max(MES_1_avg[,j])-max(data[,j])))^2
	Hill2[j] <- (abs(max(MES_2[,j])-max(data[,j])))^2
	Coov1[j] <- (abs(max(MES_3_avg[,j])-max(data[,j])))^2
	Coov2[j] <- (abs(max(MES_4[,j])-max(data[,j])))^2
}
mH1 <- mean(Hill1)
mH2 <- mean(Hill1)
mC1 <- mean(Hill1)
mC2 <- mean(Hill1)

##unweighted MES coov vs weighted MES coov##

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/MES_coov_unweighted_vs_weighted.pdf", width=20, height=6)
kk <- seq(5,41, by=1)  #by 1 may be too large 
nk <- length(kk)
MES_1 <- array(0, c(nk,159,160))
MES_1_avg <- array(0, c(nk,160))
MES_2 <- array(0, c(nk,160))
MES_3 <- array(0, c(nk,159,160))
MES_3_avg <- array(0, c(nk,160))
MES_4 <- array(0, c(nk,160))

for(j in 1:160){
#print(j)

#prima parte

	vector_x <- as.vector(data[,j])
	sorted_vector_x <- sort(vector_x, decreasing = TRUE)

	#1
	if (!is.na(coov[j, 2]) && coov[j, 2] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_prim"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#2
	if (!is.na(coov[j, 3]) && coov[j, 3] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#3
	if (!is.na(coov[j, 4]) && coov[j, 4] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_est"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#4
	if (!is.na(coov[j, 5]) && coov[j, 5] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	
	} 

	#5
	if (!is.na(coov[j, 6]) && coov[j, 6] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*1.5)
		
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#6
	if (!is.na(coov[j, 7]) && coov[j, 7] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*2.5)
		#print(temp$results$par[2]+temp$results$par[3]*2.5)
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#7
	if (!is.na(coov[j, 8]) && coov[j, 8] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#8
	if (!is.na(coov[j, 9]) && coov[j, 9] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*7)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#9 
	if (!is.na(coov[j, 10]) && coov[j, 10] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]])
	matrix_data <- as.data.frame(matrix_data) 
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#10
	if (!is.na(coov[j, 11]) && coov[j, 11] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*12.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#11
	if (!is.na(coov[j, 12]) && coov[j, 12] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#12
	if (!is.na(coov[j, 13]) && coov[j, 13] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*15)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#13
	if (!is.na(coov[j, 14]) && coov[j, 14] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#14
	if (!is.na(coov[j, 15]) && coov[j, 15] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*30)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}




	#15
	if (!is.na(coov[j, 16]) && coov[j, 16] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#16
	if (!is.na(coov[j, 17]) && coov[j, 17] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*50)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#17
	if (!is.na(coov[j, 18]) && coov[j, 18] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*75)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	
	#18
	if (!is.na(coov[j, 19]) && coov[j, 19] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*80)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#19
	if (!is.na(coov[j, 20]) && coov[j, 20] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*85)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#20
	if (!is.na(coov[j, 21]) && coov[j, 21] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*90)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

 
	#15  no coviates
	if (all(is.na(coov[j, 2:21]))) {
	mean_tail_ML_c <- NA
	}

#}

#seconda parte	
	Y <- data[,-j]

	#paiwise HILL&ML
	weights_125 <- weights_125[,j]
	for(z in 1:159){  
		vector <- as.vector(Y[, z])		
		sorted_vector <- sort(vector, decreasing = TRUE)
		for(i in 1:nk){
			thresh <- sorted_vector[kk[i]]
			indices_vector <- which(vector > thresh)
			selected_values <- vector_x[indices_vector]
			sum_s <- sum(selected_values)
			if (!is.na(mean_tail_ML_c)) {
			MES_1[i, z, j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*(sum_s)} else {MES_1[i, z, j] <- NA}			
	if (!is.na(mean_tail_ML_c)) {
	MES_1_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = mean)
	MES_3_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = function(x) weighted.mean(x, weights_125))
	} 
	else {
	MES_1_avg[,j] <- NA
	MES_3_avg[,j] <- NA }
	

	#average 
	Y_mean <- rowMeans(Y)  
	sorted_vector <- sort(Y_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values) 	#if p=1/n
		if (!is.na(mean_tail_ML_c)) {
		MES_2[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s}  else {MES_2[i,j] <- NA}
		}

	#average weighted 
 	Y_w <- sweep(Y,2, weights_125, `*`)
	Y_w_mean <- rowMeans(Y_w)
	sorted_vector <- sort(Y_w_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_w_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_4[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s 			
		} #if p=1/n

	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall (mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML UNWEIGHTED with cov", las=1)
	lines(kk, MES_1_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML unweight with cov"),
       expression("average of pairwise MES ML unweight with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_2[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event UNWEIGHTED with cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event unweight with cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML_c , col = "red") 
	text(44, 27, "highest obs:" , col = "red")
	text(44, 26, round(sorted_vector_x[1],1) , col = "red")
	text(44, 25, round(sorted_vector_x[2],1) , col = "red")
	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}
	
	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall (mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML WEIGHTED with cov", las=1)
	lines(kk, MES_3_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML weight with cov"),
       expression("average of pairwise MES ML weight with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_4[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event UNWEIGHTED with cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event unweight with cov"))
	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}
	mtext(colnames(data)[j], side = 3, line = - 2, outer = TRUE)
}
dev.off()


##check the goodness of fit##
Coov1 <- array(0, c(160))
Coov2 <- array(0, c(160))
#get the difference between prediction and real data to check if overfitting or underfitting or both
for (i in 1:160) {
	Coov1[j] <- (max(MES_3_avg[,j]-max(MES_1_avg[,j])
	Coov2[j] <- (max(MES_4[,j])-max(MES_2[,j]))
}

###SPATIAL KRIGING###
##fitting with 130 station - Hill 1 weighted##
library(RandomFields)
library(CompRandFld)
data <- read.csv("MES_hill1_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)

#create training and test sample
vect_test <- c(8, 10, 12,  17, 19, 20, 25, 26, 28, 32, 33, 34, 49, 50, 58, 64,
 75, 79, 84, 93, 94, 95, 96, 102, 109, 112,116, 133, 149, 160)

train <- data[-vect_test]
test <- data[vect_test]

#estimate the covariance function
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo_train <- coo[-vect_test,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

CorrelationParam("gauss")
corrmodel <- "gauss"
mean <- mean(train)
sill <- 1
nugget <- 0
scale <- 0.5
# Fixed parameters
fixed<-list(mean=mean)
 
# Starting value for the estimated parameters
start<-list(scale=scale, sill=sill, nugget=nugget)

fit <- FitComposite(train, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

#fit the kriging model and predict for test set 
coo_test <- coo[vect_test,]
x <- as.vector(coo_test[,7])
y <- as.vector(coo_test[,6])
coords_pred<-cbind(x,y)

result_matrix1 <- matrix(nrow = 1, ncol = 30)
for (i in 1:30) {

pr<-Kri(data=train,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix2 <- matrix(nrow = 1, ncol = 30)
for (i in 1:30) {
	result_matrix2[1,i] <- abs(result_matrix1[1,i]-test[i])
}
mean(result_matrix2[1,]) #mean absolute error

##fitting with 130 station - Hill 2 weighted##
library(RandomFields)
library(CompRandFld)
data <- read.csv("MES_hill2_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)

#create training and test sample
vect_test <- c(8, 10, 12,  17, 19, 20, 25, 26, 28, 32, 33, 34, 49, 50, 58, 64,
 75, 79, 84, 93, 94, 95, 96, 102, 109, 112,116, 133, 149, 160)

train <- data[-vect_test]
test <- data[vect_test]

#estimate the covariance function
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo_train <- coo[-vect_test,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

CorrelationParam("gauss")
corrmodel <- "gauss"
mean <- mean(train)
sill <- 1
nugget <- 0
scale <- 0.5
# Fixed parameters
fixed<-list(mean=mean)
 
# Starting value for the estimated parameters
start<-list(scale=scale, sill=sill, nugget=nugget)

fit <- FitComposite(train, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

#fit the kriging model and predict for test set 
coo_test <- coo[vect_test,]
x <- as.vector(coo_test[,7])
y <- as.vector(coo_test[,6])
coords_pred<-cbind(x,y)

result_matrix1 <- matrix(nrow = 1, ncol = 30)
for (i in 1:30) {

pr<-Kri(data=train,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix2 <- matrix(nrow = 1, ncol = 30)
for (i in 1:30) {
	result_matrix2[1,i] <- abs(result_matrix1[1,i]-test[i])
}
mean(result_matrix2[1,]) #mean absolute error


##fitting with 160 station - 29 new MES - Hill ##
#estimate the MES for the new 29 stations#
data <- read.csv("new_stations_29.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
					
vect <- c(c(31,36), c(27,32), c(11,16), c(14,19), c(29,34), c(30,35),
c(22,27), c(18,23), c(21,26), c(20,25), c(35,40), c(31,36), 
c(31,36), c(35,40), c(16,21), c(36,41), c(31,36), c(17,23),
c(24,29), c(25,30), c(35,40), c(23,28), c(35,40), c(33,38),
c(25,30), c(36,41), c(26,31), c(20,25), c(30,35)) 
						
list_ranges <-  matrix(vect, ncol = 2, byrow = TRUE)
weights <- read.csv("norm_weights_29.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/MES_hill_pred_29.pdf", width=20, height=6)
kk <- seq(5,41, by=1)  #by 1 may be too large 
nk <- length(kk)
MES_1 <- array(0, c(nk,28,29))
MES_1_avg <- array(0, c(nk,29))
MES_2 <- array(0, c(nk,29))
MES_3_avg <- array(0, c(nk,29))
MES_4 <- array(0, c(nk,29))

for(j in 1:29){
#print(j)

#prima parte

	vector_x <- as.vector(data[,j])
	sorted_vector_x <- sort(vector_x, decreasing = TRUE)
	v_ML <- list()
	for (k1 in list_ranges[j,1]:list_ranges[j,2]){
		temp <- HTailIndex(vector_x, k1, TRUE, varType = "asym-Ind")
		v_ML <- append(v_ML, temp$gammaHat)
	}
	v_ML <- unlist(v_ML)
	mean_tail_ML <- mean(v_ML)


#seconda parte	
	Y <- data[,-j]

	#paiwise HILL&ML
	weights <- weights[,j]
	for(z in 1:28){  
		vector <- as.vector(Y[, z])		
		sorted_vector <- sort(vector, decreasing = TRUE)
		for(i in 1:nk){
			thresh <- sorted_vector[kk[i]]
			indices_vector <- which(vector > thresh)
			selected_values <- vector_x[indices_vector]
			sum_s <- sum(selected_values)
			MES_1[i, z, j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*(sum_s) 			#if p=1/n
		}} #}
	MES_1_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = mean)
	
	#paiwise ML
	MES_3_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = function(x) weighted.mean(x, weights))

	#average 
	Y_mean <- rowMeans(Y)  
	sorted_vector <- sort(Y_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_2[i,j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*sum_s 			#if p=1/n
}

	#average weighted 
 	Y_w <- sweep(Y,2, weights, `*`)
	Y_w_mean <- rowMeans(Y_w)
	sorted_vector <- sort(Y_w_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_w_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_4[i,j] <- ((kk[i]^mean_tail_ML)*(1/kk[i]))*sum_s 			
		}  
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	#par(mfrow=c(1, 2), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall(mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="UNWEIGHTED pairwise MES hill without cov", las=1)
	lines(kk, MES_1_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES HILL without cov"),
       expression("average of pairwise MES HILL without cov")))

	
	plot(kk, MES_2[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",    
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="UNWEIGHTED MES with AVERAGE as cond event", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      legend=expression("MES HILL with AVERAGE as cond event without cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML , col = "red") 
	text(44, 27, "highest obs:" , col = "red")
	text(44, 26, round(sorted_vector_x[1],1) , col = "red")
	text(44, 25, round(sorted_vector_x[2],1) , col = "red")
	text(44, 24, "tail thresh:" , col = "red")
	text(44, 23, round(sorted_vector_x[k1],1) , col = "red")

	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall(mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="WEIGHTED pairwise MES hill without cov", las=1)
	lines(kk, MES_3_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES HILL without cov"),
       expression("average of pairwise MES HILL without cov")))
	
	plot(kk, MES_4[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",    
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="WEIGHTED MES with AVERAGE as cond event", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      legend=expression("MES HILL with AVERAGE as cond event without cov"))

	mtext(colnames(data)[j], side = 3, line = - 2, outer = TRUE)
}
dev.off()

#fitting kriging model and prediction for 29 - hill 1 weighted#
data <- read.csv("MES_hill1_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
x <- as.vector(coo[,7])
y <- as.vector(coo[,6])
coords<-cbind(x,y)

CorrelationParam("gauss")
corrmodel <- "gauss"
mean <- mean(data)
sill <- 1
nugget <- 0
scale <- 0.5
# Fixed parameters
fixed<-list(mean=mean)
 
# Starting value for the estimated parameters
start<-list(scale=scale, sill=sill, nugget=nugget)

fit <- FitComposite(data, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

loc_p <- read.csv("loc_new_stations_29.csv", header = TRUE)
x <- as.vector(loc_p[,3])
y <- as.vector(loc_p[,2])
coords_pred<-cbind(x,y)
result_matrix1 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
pr<-Kri(data=data,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix3 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
	result_matrix3[1,i] <- abs(result_matrix1[1,i]-MES_3_avg[37,i])
}
mean(result_matrix3[1,]) 


#fitting kriging model and prediction for 29 - hill 1 weighted#
data <- read.csv("MES_hill2_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
x <- as.vector(coo[,7])
y <- as.vector(coo[,6])
coords<-cbind(x,y)

CorrelationParam("gauss")
corrmodel <- "gauss"
mean <- mean(data)
sill <- 1
nugget <- 0
scale <- 0.5
# Fixed parameters
fixed<-list(mean=mean)
 
# Starting value for the estimated parameters
start<-list(scale=scale, sill=sill, nugget=nugget)

fit <- FitComposite(data, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

loc_p <- read.csv("loc_new_stations_29.csv", header = TRUE)
x <- as.vector(loc_p[,3])
y <- as.vector(loc_p[,2])
coords_pred<-cbind(x,y)
result_matrix1 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
pr<-Kri(data=data,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix3 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
	result_matrix3[1,i] <- abs(result_matrix1[1,i]-MES_4[37,i])
}
mean(result_matrix3[1,]) 

##training with 97 stations - Coov 1 weighted ##
data <- read.csv("MES_cov1_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
max_values <- apply(data, 2, max)
names(max_values) <- colnames(data)
na_indices <- which(is.na(max_values))
max_values <- (max_values[-na_indices])
train <- (max_values[vect_train])
test <- (max_values[vect_test])
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
coo_train <- coo[vect_train,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

corrmodel <- "stable"  
mean <- mean(data)
sill <- 1
nugget <- 0
scale <- 0.5
power <-2
# Fixed parameters
fixed<-list(mean=mean)
# Starting value for the estimated parameters
start<-list(scale=scale, power=power, sill=sill, nugget=nugget)

fit <- FitComposite(train, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

coo_test <- coo[vect_test,]
x <- as.vector(coo_test[,7])
y <- as.vector(coo_test[,6])
coords_pred<-cbind(x,y)

result_matrix1 <- matrix(nrow = 1, ncol = 28)
for (i in 1:28) {
pr<-Kri(data=train,coordx=coords ,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix2 <- matrix(nrow = 1, ncol = 28)
for (i in 1:28) {
	result_matrix2[1,i] <- abs(result_matrix1[1,i]-test[i])
	}
mean(result_matrix2[1,]) #absolute mean error

##training with 97 stations - Coov 2 weighted ##
data <- read.csv("MES_cov2_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
max_values <- apply(data, 2, max)
names(max_values) <- colnames(data)
na_indices <- which(is.na(max_values))
max_values <- (max_values[-na_indices])
train <- (max_values[vect_train])
test <- (max_values[vect_test])
coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
coo_train <- coo[vect_train,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

corrmodel <- "stable"  
mean <- mean(data)
sill <- 1
nugget <- 0
scale <- 0.5
power <-2
# Fixed parameters
fixed<-list(mean=mean)
# Starting value for the estimated parameters
start<-list(scale=scale, power=power, sill=sill, nugget=nugget)

fit <- FitComposite(train, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

coo_test <- coo[vect_test,]
x <- as.vector(coo_test[,7])
y <- as.vector(coo_test[,6])
coords_pred<-cbind(x,y)

result_matrix1 <- matrix(nrow = 1, ncol = 28)
for (i in 1:28) {
pr<-Kri(data=train,coordx=coords ,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix2 <- matrix(nrow = 1, ncol = 28)
for (i in 1:28) {
	result_matrix2[1,i] <- abs(result_matrix1[1,i]-test[i])
	}
mean(result_matrix2[1,]) #absolute mean error

##training 125 - 15 new - coov 1 weighted##
#estimating MES for 15 new stations#
library(ExtremeRisks)
library(extRemes)
main_folder <- 'C:/Users/39346/OneDrive/Desktop/THESIS'  #versione corretta
setwd(main_folder)
data <- read.csv("new_stations_29_seasons.csv", header = TRUE) 
data_t<- read.csv("new_stations_29_temp.csv", header = TRUE)
data_v <- read.csv("new_stations_29_vento.csv", header = TRUE) 
data_r<- read.csv("new_stations_29_rad.csv", header = TRUE)               
data_u<- read.csv("new_stations_29_umid.csv", header = TRUE)

#modificati					
vect <- c(
c(31,36), c(27,32), c(25,30), c(14,19), c(16,21), c(30,35),
c(22,27), c(18,23), c(21,26), c(20,25), c(35,40), c(31,36), 
c(31,36), c(35,40), c(16,21), c(36,41), c(29,34), c(17,23),
c(24,29), c(25,30), c(35,40), c(23,28), c(35,40), c(33,38),
c(24,29), c(36,41), c(26,31), c(17,23), c(30,35)) 
						
list_ranges <-  matrix(vect, ncol = 2, byrow = TRUE)


weights <- read.csv("norm_weights_29.csv", header = TRUE)
coov <- read.csv("COOV29.csv", header = TRUE)

pdf("C:/Users/39346/OneDrive/Desktop/THESIS/rainfall_analysis_unscaled_MES_WEIGHT_coov29.pdf", width=20, height=6)
kk <- seq(5,41, by=1)  #by 1 may be too large 
nk <- length(kk)
MES_1 <- array(0, c(nk,28,29))
MES_1_avg <- array(0, c(nk,29))
MES_2 <- array(0, c(nk,29))
MES_3 <- array(0, c(nk,28,29))
MES_3_avg <- array(0, c(nk,29))
MES_4 <- array(0, c(nk,29))

for(j in 1:29){
#print(j)

#prima parte

	vector_x <- as.vector(data[,j])
	sorted_vector_x <- sort(vector_x, decreasing = TRUE)

	#1
	if (!is.na(coov[j, 2]) && coov[j, 2] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_prim"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#2
	if (!is.na(coov[j, 3]) && coov[j, 3] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_inverno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#3
	if (!is.na(coov[j, 4]) && coov[j, 4] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno_est"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#4
	if (!is.na(coov[j, 5]) && coov[j, 5] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data[, "autunno"]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3])}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	
	} 

#4.5
	if (!is.na(coov[j, 6]) && coov[j, 6] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*1.5)
		
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#5
	if (!is.na(coov[j, 7]) && coov[j, 7] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
		#print(sorted_vector_x[k1])
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1] , shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)
		#print(fevd)
		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*2.5)
		#print(temp$results$par[2]+temp$results$par[3]*2.5)
		}
		
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#6
	if (!is.na(coov[j, 8]) && coov[j, 8] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#6
	if (!is.na(coov[j, 9]) && coov[j, 9] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*7)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#7
	if (!is.na(coov[j, 10]) && coov[j, 10] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]])
	matrix_data <- as.data.frame(matrix_data) 
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#8
	if (!is.na(coov[j, 11]) && coov[j, 11] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_t[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*12.5)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#9
	if (!is.na(coov[j, 12]) && coov[j, 12] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*10)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#9
	if (!is.na(coov[j, 13]) && coov[j, 13] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*15)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#10
	if (!is.na(coov[j, 14]) && coov[j, 14] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#10
	if (!is.na(coov[j, 15]) && coov[j, 15] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_v[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1],shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*30)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}




	#11
	if (!is.na(coov[j, 16]) && coov[j, 16] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*25)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#12
	if (!is.na(coov[j, 17]) && coov[j, 17] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_r[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*50)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#12
	if (!is.na(coov[j, 18]) && coov[j, 18] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*75)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	
	#13
	if (!is.na(coov[j, 19]) && coov[j, 19] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*80)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}


	#13
	if (!is.na(coov[j, 20]) && coov[j, 20] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*85)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}



	#14
	if (!is.na(coov[j, 21]) && coov[j, 21] == 1) {
	v_ML_c <- list()
	matrix_data <- cbind(data[, j], data_u[, colnames(data)[j]]) 
	matrix_data <- as.data.frame(matrix_data)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], shape.fun = ~V2, type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2]+temp$results$par[3]*90)}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	}

	#14
	if (!is.na(coov[j, 22]) && coov[j, 22] == 1) {
	v_ML_c <- list()
	matrix_data <- as.data.frame(vector_x)
	if ( list_ranges[j,1] %% 2 == 0) {
  	p <- seq( list_ranges[j,1],  list_ranges[j,2], by = 2)
	} else {
 	 p <- seq( list_ranges[j,1]+1,  list_ranges[j,2], by = 2)
	}
	for (k1 in p){
  		temp <- tryCatch(
			fevd(vector_x,  matrix_data, threshold = sorted_vector_x[k1], type="GP", time.units = "months"),
			error = function(e) {
				#ML_title <- "ML_extRemes_PROBLEM"
        			return(NULL)  # Return a default value or 'NA'
      		}
    		)

		if (is.null(temp)) {
			next
		} else {	
		v_ML_c <- append(v_ML_c, temp$results$par[2])}
	}
	v_ML_c <- unlist(v_ML_c)
	mean_tail_ML_c <- mean(v_ML_c)
	} 

	#15 
	if (all(is.na(coov[j, 2:22]))) {
	mean_tail_ML_c <- NA
	}

#}

#seconda parte	
	Y <- data[,-j]

	#paiwise HILL&ML
	weights <- weights[,j]
	for(z in 1:28){  
		vector <- as.vector(Y[, z])		
		sorted_vector <- sort(vector, decreasing = TRUE)
		for(i in 1:nk){
			thresh <- sorted_vector[kk[i]]
			indices_vector <- which(vector > thresh)
			selected_values <- vector_x[indices_vector]
			sum_s <- sum(selected_values)
			if (!is.na(mean_tail_ML_c)) {
			MES_1[i, z, j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*(sum_s)} else {MES_1[i, z, j] <- NA}			
			#if (!is.na(mean_tail_ML_c)) {
			#MES_3[i, z, j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*(sum_s)} else {MES_3[i, z, j] <- NA}
		}} #} #if p=1/n
	if (!is.na(mean_tail_ML_c)) {
	MES_1_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = mean)
	MES_3_avg[,j] <- apply(MES_1[,,j], MARGIN = 1, FUN = function(x) weighted.mean(x, weights))
	} 
	else {
	MES_1_avg[,j] <- NA
	MES_3_avg[,j] <- NA }
	

	#average 
	Y_mean <- rowMeans(Y)  
	sorted_vector <- sort(Y_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values) 	#if p=1/n
		if (!is.na(mean_tail_ML_c)) {
		MES_2[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s}  else {MES_2[i,j] <- NA}
		}

	#average weighted 
 	Y_w <- sweep(Y,2, weights, `*`)
	Y_w_mean <- rowMeans(Y_w)
	sorted_vector <- sort(Y_w_mean, decreasing = TRUE)
	for(i in 1:nk){
		thresh <- sorted_vector[kk[i]]
		indices_vector <- which(Y_w_mean > thresh)
		selected_values <- vector_x[indices_vector]
		sum_s <- sum(selected_values)
		MES_4[i,j] <- ((kk[i]^mean_tail_ML_c)*(1/kk[i]))*sum_s 			
		} #if p=1/n

	#}
	par(mfrow=c(1, 4), mai=c(.4,.3,.4,.3), mgp=c(1.6,.6,0))
	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall (mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML UNWEIGHTED with cov", las=1)
	lines(kk, MES_1_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML unweight with cov"),
       expression("average of pairwise MES ML unweight with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_2[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event UNWEIGHTED with cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event unweight with cov"))
	text(44, 30, "tail index:" , col = "red")
	text(44, 29, mean_tail_ML_c , col = "red") 
	text(44, 27, "highest obs:" , col = "red")
	text(44, 26, round(sorted_vector_x[1],1) , col = "red")
	text(44, 25, round(sorted_vector_x[2],1) , col = "red")
	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}
	
	if (!is.na(mean_tail_ML_c)) {
	matplot(kk, MES_1[,,j], type="l", ylim=c(0,30), xlim=c(0,45), lty=2, lwd=0.5, ylab="rainfall (mm)",
	xlab="k", cex.lab=1.3, cex.axis=1.3, main="pairwise MES ML WEIGHTED with cov", las=1)
	lines(kk, MES_3_avg[,j], lty=1, lwd=3, col="black")
	legend("bottomright", lty=c(2,1), col=c("magenta", "black"), lwd=c(0.5,3 ),
       bty="n", cex=1.3, legend=c(expression("pairwise MES ML weight with cov"),
       expression("average of pairwise MES ML weight with cov")))} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")}

	if (!is.na(mean_tail_ML_c)) {
	plot(kk, MES_4[,j], type="l", ylim=c(0,30), xlim=c(0,45), lwd=3, ylab="",     
      	xlab="k", cex.lab=1.3, cex.axis=1.3, main="MES ML with AVERAGE as cond event UNWEIGHTED with cov", las=1)
	legend("bottomright", col=c(1), lwd=c(3), bty="n", cex=1.3, lty=c(1),
      	legend=expression("MES ML with AVERAGE as cond event unweight with cov"))
	} else {
	 plot(1, type = "n", ylim=c(0,30), xlab = "", ylab = "", main = "No coovariate good enough")
	}
	mtext(colnames(data)[j], side = 3, line = - 2, outer = TRUE)
}
dev.off()


#fitting kriging model and predicting - coov 1 weighted #
data <- read.csv("MES_cov1_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
max_values <- apply(data, 2, max)
names(max_values) <- colnames(data)
na_indices <- which(is.na(max_values))
max_values <- (max_values[-na_indices])

coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
coo_train <- coo[-na_indices,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

corrmodel <- "stable"  
mean <- mean(max_values)
sill <- 1
nugget <- 0
scale <- 0.5
power <-2
# Fixed parameters
fixed<-list(mean=mean)
# Starting value for the estimated parameters
start<-list(scale=scale, power=power, sill=sill, nugget=nugget)
fit <- FitComposite(max_values, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

loc_p <- read.csv("loc_new_stations_29.csv", header = TRUE)
x <- as.vector(loc_p[,3])
y <- as.vector(loc_p[,2])
coords_pred<-cbind(x,y)
result_matrix1 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
pr<-Kri(data=max_values,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix3 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
result_matrix3[1,i] <- abs(result_matrix1[1,i]-MES_3_avg[37,i])}
mean(result_matrix3[1,], na.rm = TRUE)

#fitting kriging model and predicting - coov 2 weighted #
data <- read.csv("MES_cov1_weight.csv", header = TRUE)
data <- sapply(data[,], as.numeric)
max_values <- apply(data, 2, max)
names(max_values) <- colnames(data)
na_indices <- which(is.na(max_values))
max_values <- (max_values[-na_indices])

coo <- read.csv("location_282_382_160.csv", header = TRUE)
coo[,6:7] <- sapply(coo[,6:7], as.numeric)
coo_train <- coo[-na_indices,]
x <- as.vector(coo_train[,7])
y <- as.vector(coo_train[,6])
coords<-cbind(x,y)

corrmodel <- "stable"  
mean <- mean(max_values)
sill <- 1
nugget <- 0
scale <- 0.5
power <-2
# Fixed parameters
fixed<-list(mean=mean)
# Starting value for the estimated parameters
start<-list(scale=scale, power=power, sill=sill, nugget=nugget)
fit <- FitComposite(max_values, coordx=coords, coordy=NULL, coordt=NULL, grid=FALSE, corrmodel=corrmodel,likelihood="Full",type="Standard",varest=TRUE,start=start,fixed=fixed)
print(fit)

loc_p <- read.csv("loc_new_stations_29.csv", header = TRUE)
x <- as.vector(loc_p[,3])
y <- as.vector(loc_p[,2])
coords_pred<-cbind(x,y)
result_matrix1 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
pr<-Kri(data=max_values,coordx=coords,corrmodel=corrmodel, loc=coords_pred[i,], param= as.list(c(fit$param,fit$fixed)))
result_matrix1[1,i] <- pr$pred}
result_matrix1[1,]

result_matrix3 <- matrix(nrow = 1, ncol = 29)
for (i in 1:29) {
result_matrix3[1,i] <- abs(result_matrix1[1,i]-MES_4[37,i])}
mean(result_matrix3[1,], na.rm = TRUE)



