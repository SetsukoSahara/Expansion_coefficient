# script_common_Tc.R
# 
# Original script was written by Charles F Stevens
# Modified by Takashi Kodama
#
# (3/21/2019)
#
# This script implements a function "expansion.coefficient", which computes
# expansion coefficients of NE and RG based on common cell-cycle length.
#
# [How to use]
# 1. Make sure two csv files (cell_fraction.csv, Tc.csv) and this script file are 
#	all in the current working directory,
# 2. In console, run "source("script_common_Tc.R")"
# 3. In console, run "expansion.coefficients()"
#

expansion.coefficients = function(cell_fraction_csv_name = "cell_fraction.csv",
							Tc_csv_name = "Tc.csv",
							noise_amplitude = 0.01, # noise amplitude used in robustness test
							nTest = 10000, # number of test in robustness test
							sd.threshold = 0.1 # threshold value for robustness test
							) 
{
	#................................................................................................
	# SECTION 1: read in cell fraction data from csv files
	#................................................................................................
	
	# Matrix of fractions (%) of differnt cell types at different developmental stages
	C_raw = read.csv(cell_fraction_csv_name) 
	C_raw = as.matrix(C_raw)
	
	
	# Cell fraction (excluding GABAergic interneuron and glial progenitor)
	fraction_sum = apply(C_raw[,c("f_neuron", "f_NE", "f_RG_SNP", "f_bIP", "f_gliaP", "f_GABA")],1,sum)
	fraction_sum.target = apply(C_raw[,c("f_neuron", "f_NE", "f_RG_SNP", "f_bIP")],1,sum)	#w/o GABA and gliaP 
	fraction.adj.factor = fraction_sum/fraction_sum.target
	
	C = C_raw[, c("stage", "f_neuron", "f_NE", "f_RG_SNP", "f_bIP")] 
	C[, -1] = fraction.adj.factor * C[, -1]/100
	
	
	# SEM (excluding GABAergic interneuron and glial progenitor)
	C_sem = C_raw[, c("sem_neuron", "sem_NE", "sem_RG_SNP", "sem_bIP")] 
	C_sem = fraction.adj.factor * C_sem/100
	
	
	# Cell cycle duration [hour] at different developmental stages
	Tc = read.csv(Tc_csv_name) 
	Tc = as.matrix(Tc)


	nCellType = 4 # number of cell type consiered (NE, RG, bIP, neuron)


	stage_in_day = C[,"stage"] # E10.5 ~ 18.5
	

	#................................................................................................
	# SECTION 2: Cell-cycle lengths (Tc), modeled by sigmoid function
	#................................................................................................

	x = Tc[, "stage"] * 24 # [hour]
	y = Tc[, "Tc_all_dissociated"] # [hour]
	Tc.data = data.frame(x, y)
	Tc.coef = nls(y ~ a + b/(1 + exp(-(x - c)/d)), Tc.data, start = list(a = 10, b = 8, c = 330, d = 7)) # sigmoid fitting


	# The function "Tc.func" estimates Tc at arbitrary embryonic stage in hour (measuring from the time of firtilization), 
	# based on the regression model (sigmoid)
	Tc.func = function(t_in_hour) {
		a = coef(Tc.coef)["a"]
		b = coef(Tc.coef)["b"]
		c = coef(Tc.coef)["c"]
		d = coef(Tc.coef)["d"]

		return(a + b/(1 + exp(-(t_in_hour - c)/d)))
	}		


	#................................................................................................
	# SECTION 3: Functions to convert time embryonic stage in hour <-> cell division number
	#................................................................................................

	# Given that Tc.func(t) is cell cycle time at time t [hour],
	# 1/Tc.func(t) is "velocity" of cell cycle at time t [cell cycle/hour]
	# Integral of 1/Tc.func(t) therefore corresponds to cumulative cell cycle at time t [hour]

	# cell cycle "velocity" at time t in hour
	inv_Tc.func = function(t_in_hour) {
		return(1/Tc.func(t_in_hour))
	}

	# developmental stage t in hour -> cell division number (counting from E10.5)
	Ct.func = function(t_in_hour) { # cumulative cell cycle at time t_in_hour
		count_from = 10.5 * 24 #count from 10.5 days (228 hours) after the firtilization	
		return(integrate(inv_Tc.func, count_from, t_in_hour)$value)
	}

	# cell division number t (counting from E10.5) -> developmental stage in hour (inverse function of Ct.func)
	inv_Ct.func = function(t_in_cellDiv) {
		temp_func = function(x) {
			return(Ct.func(x) - t_in_cellDiv)
		}

		result = uniroot(temp_func, c(0, 20 * 24),extendInt = "yes") # developmental stage in hour ranges from 0 to 20*24 hours
		return(result$root)
	}


	#................................................................................................
	# SECTION 4: modeling cell fractions 
	#................................................................................................
	
	
	# modeled cell fraction as a function of developmental stage in hour, stored in C_in_hour	
	range_hour = 252:444 # [hour], corresponding to E10.5~18.5
	C_in_hour = matrix(NA, nrow = length(range_hour), ncol = nCellType)
	colnames(C_in_hour) = c("f_NE", "f_RG_SNP", "f_bIP", "f_neuron")
	

	# modeled cell fraction as a function of cell division number, stored in C_in_cellDiv
	range_cellDiv = seq(from = 0, to = floor(Ct.func(444)), by = 1) #  counting from E10.5 to E18.5(444 hour)
	C_in_cellDiv = matrix(NA, nrow = length(range_cellDiv), ncol = nCellType)
	colnames(C_in_cellDiv) = c("f_NE", "f_RG_SNP", "f_bIP", "f_neuron")


	# developmental stage in hour corresponding to each cell division number
	range_cellDiv_in_hour = sapply(range_cellDiv, inv_Ct.func)

	
	#
	# modeling by 3rd polynominal function
	#
	
	x = stage_in_day * 24 #[hour]	

	# NE
	y = C[,"f_NE"]
	t_fixed.NE = 18.5*24 #[hour]
	f_fixed.NE = 0 # cell fraction
	NE.data = data.frame(x, y)
	NE.coef = nls(y ~ f_fixed.NE + a*(x-t_fixed.NE) + b*(x-t_fixed.NE)^2 + c*(x-t_fixed.NE)^3, 
				NE.data, start = list(a = -0.010019, b = 2.2719e-5, c = 1.2723e-8))
	NE.func = function(t_in_hour){

		a = coef(NE.coef)["a"]
		b = coef(NE.coef)["b"]
		c = coef(NE.coef)["c"]
	
		return(f_fixed.NE + a*(t_in_hour-t_fixed.NE) + b*(t_in_hour-t_fixed.NE)^2 + c*(t_in_hour-t_fixed.NE)^3)
	}


	# RG
	y = C[, "f_RG_SNP"]
	t_fixed.RG = 10 * 24 #[hour]
	f_fixed.RG = 0 # cell fraction
	RG.data = data.frame(x, y)
	RG.coef = nls(y ~ f_fixed.RG + a * (x - t_fixed.RG) + b * (x - t_fixed.RG)^2 + c * (x - t_fixed.RG)^3, 
				RG.data, start = list(a = 0.010379, b = -9.3445e-05, c = 2.1323e-07))
	RG.func = function(t_in_hour) {

		a = coef(RG.coef)["a"]
		b = coef(RG.coef)["b"]
		c = coef(RG.coef)["c"]

		if (t_in_hour < t_fixed.RG) {
			return(0)
		} else {
			return(f_fixed.RG + a * (t_in_hour - t_fixed.RG) + b * (t_in_hour - t_fixed.RG)^2 + c * (t_in_hour - t_fixed.RG)^3)
		}
	}


	# bIP
	y = C[, "f_bIP"]
	t_fixed.IP = 12.5 * 24 #[hour]
	f_fixed.IP = 0 # cell fraction
	IP.data = data.frame(x, y)
	IP.coef = nls(y ~ f_fixed.IP + a * (x - t_fixed.IP) + b * (x - t_fixed.IP)^2 + c * (x - t_fixed.IP)^3,
				IP.data, start = list(a = 0.0007721, b = 1.0156e-05, c = -8.7322e-08))
	IP.func = function(t_in_hour) {

		a = coef(IP.coef)["a"]
		b = coef(IP.coef)["b"]
		c = coef(IP.coef)["c"]

		if (t_in_hour < t_fixed.IP) {
			return(0)
		} else {
			return(f_fixed.IP + a * (t_in_hour - t_fixed.IP) + b * (t_in_hour - t_fixed.IP)^2 + c * (t_in_hour - t_fixed.IP)^3)
		}
	}


	# Neuron (EN)
	y = C[, "f_neuron"]
	t_fixed.EN = 11 * 24 #[hour]
	f_fixed.EN = 0 # cell fraction
	EN.data = data.frame(x, y)
	EN.coef = nls(y ~ f_fixed.EN + a * (x - t_fixed.EN) + b * (x - t_fixed.EN)^2 + c * (x - t_fixed.EN)^3,
				EN.data, start = list(a = 3e-04, b = 1e-05, c = -1e-07))
	EN.func = function(t_in_hour) {

		a = coef(EN.coef)["a"]
		b = coef(EN.coef)["b"]
		c = coef(EN.coef)["c"]

		if (t_in_hour < t_fixed.EN) {
			return(0)
		} else {
			return(f_fixed.EN + a * (t_in_hour - t_fixed.EN) + b * (t_in_hour - t_fixed.EN)^2 + c * (t_in_hour - t_fixed.EN)^3)
		}
	}


	# modeled cell fraction as a function of developmental stage in hour (C_in_hour)
	C_in_hour[, "f_NE"] = sapply(range_hour, NE.func)
	C_in_hour[, "f_RG_SNP"] = sapply(range_hour, RG.func)
	C_in_hour[, "f_bIP"] = sapply(range_hour, IP.func)
	C_in_hour[, "f_neuron"] = sapply(range_hour, EN.func)
	

	# modelded cell fraction as a function of cell division number (C_in_cellDiv)
	C_in_cellDiv[, "f_NE"] = sapply(range_cellDiv_in_hour, NE.func)
	C_in_cellDiv[, "f_RG_SNP"] = sapply(range_cellDiv_in_hour, RG.func)
	C_in_cellDiv[, "f_bIP"] = sapply(range_cellDiv_in_hour, IP.func)
	C_in_cellDiv[, "f_neuron"] = sapply(range_cellDiv_in_hour, EN.func)

	
	
	# Check and report deviation of fraction sum from 1
	C_in_hour_fractionSum = apply(C_in_hour,1,sum)
	printf("Fraction sum deviation from 1: Avg: %f (%f ~ %f)\n", 
			mean(C_in_hour_fractionSum)-1, min(C_in_hour_fractionSum)-1, max(C_in_hour_fractionSum)-1)
	

	#................................................................................................
	# SECTION 5: Expansion coefficient (alpha) for NE (neuroepithelial cells)
	#................................................................................................

	alpha_NE = function(C_in_cellDiv) {

		# Making  sum of cell fraction one
		C_in_cellDiv.rowsum = apply(C_in_cellDiv, 1, sum)
		C_in_cellDiv.temp = C_in_cellDiv/C_in_cellDiv.rowsum
		
		
		as = rep(NA, times = dim(C_in_cellDiv)[1])
		s = C_in_cellDiv.temp[, "f_NE"]
		n = C_in_cellDiv.temp[, "f_neuron"]

		for (k in 1:dim(C_in_cellDiv)[1]) {
			as[k] = (2 - n[k]) * s[k + 1]/s[k]
		}

		return(as)
	}

	as = alpha_NE(C_in_cellDiv)

	#................................................................................................
	# SECTION 6: Expansion coeficient (alpha) for RG (radial glial cells), method 1
	#................................................................................................

	alpha_RG1 = function(C_in_cellDiv) {

		# Making  sum of cell fraction one
		C_in_cellDiv.rowsum = apply(C_in_cellDiv, 1, sum)
		C_in_cellDiv.temp = C_in_cellDiv/C_in_cellDiv.rowsum


		ar = rep(NA,times = dim(C_in_cellDiv)[1])
		s = C_in_cellDiv.temp[, "f_NE"]
		r = C_in_cellDiv.temp[, "f_RG_SNP"]
		t = C_in_cellDiv.temp[, "f_NE"] + C_in_cellDiv.temp[, "f_RG_SNP"]
		n = C_in_cellDiv.temp[, "f_neuron"]
		
		for (k in 1:dim(C_in_cellDiv)[1]) {
			ar[k] = ((2 - n[k]) * t[k + 1] - 2 * s[k])/r[k]
		}

		return(ar)
	}

	ar1 = alpha_RG1(C_in_cellDiv)

	#................................................................................................
	# SECTION 7: Expansion coeficient (alpha) for RG, method 2
	#................................................................................................

	alpha_RG2 = function(C_in_cellDiv) {

		# Making  sum of cell fraction one
		C_in_cellDiv.rowsum = apply(C_in_cellDiv, 1, sum)
		C_in_cellDiv.temp = C_in_cellDiv/C_in_cellDiv.rowsum

		br = rep(NA, times = dim(C_in_cellDiv)[1])
		r = C_in_cellDiv.temp[, "f_RG_SNP"]
		p = C_in_cellDiv.temp[, "f_bIP"]
		t = C_in_cellDiv.temp[, "f_neuron"] + C_in_cellDiv.temp[, "f_bIP"]
		n = C_in_cellDiv.temp[, "f_neuron"]

		for (k in 1:dim(C_in_cellDiv)[1]) {
			br[k] = ((2 - n[k]) * t[k + 1] - n[k] - 2 * p[k])/r[k]
		}

		return(2 - br)
	}

	ar2 = alpha_RG2(C_in_cellDiv)

	#................................................................................................
	# SECTION 8: roubustness of expansion coefficinet estimation
	#................................................................................................

	measure.SD.byNoise = function(C_in_cellDiv, alpha.func, noise_amplitude, nTest) {
		a.sum = rep(0, dim(C_in_cellDiv)[1])
		a.sqsum = rep(0, dim(C_in_cellDiv)[1])

		for (i in 1:nTest) {
			C_in_cellDiv_noise = C_in_cellDiv + runif(length(C_in_cellDiv), min = -noise_amplitude, max = noise_amplitude)
			a.noise = alpha.func(C_in_cellDiv_noise)
			a.sum = a.sum + a.noise
			a.sqsum = a.sqsum + a.noise^2
		}

		a.SD = sqrt(a.sqsum/nTest - (a.sum/nTest)^2)

		return(a.SD)
	}

	as_sd = measure.SD.byNoise(C_in_cellDiv, alpha_NE, noise_amplitude, nTest)
	ar1_sd = measure.SD.byNoise(C_in_cellDiv, alpha_RG1, noise_amplitude, nTest)
	ar2_sd = measure.SD.byNoise(C_in_cellDiv, alpha_RG2, noise_amplitude, nTest)

	
	#................................................................................................
	# SECTION 9: plot results
	#................................................................................................
	sizeGrWindow(8, 8)
	par(mfrow = c(2, 2))


	#
	# (Top left)
	# plots of cell fraction modeling results vs embryonic stage
	#
	plot(NA, xlab = "embryonic stage in day", ylab = "cell fraction",
		xlim = c(min(stage_in_day), max(stage_in_day)), 
		ylim = c(0, 1), main = "Cell fraction vs Embryonic stage")
	
	# raw data
	points(stage_in_day, C[, "f_NE"], col = "red",pch=19, cex = 2)
	points(stage_in_day, C[, "f_RG_SNP"], col = "green",pch=19, cex = 2)
	points(stage_in_day, C[, "f_bIP"], col = "blue",pch=19, cex = 2)
	points(stage_in_day, C[, "f_neuron"], col = "purple",pch=19, cex = 2)	
	arrows(stage_in_day, C[, "f_NE"]-C_sem[, "sem_NE"], stage_in_day, C[, "f_NE"]+C_sem[, "sem_NE"], length=0.1, angle=90, code=3,col = "red") # error bar
	arrows(stage_in_day, C[, "f_RG_SNP"]-C_sem[, "sem_RG_SNP"], stage_in_day, C[, "f_RG_SNP"]+C_sem[, "sem_RG_SNP"], length=0.1, angle=90, code=3,col = "green") # error bar
	arrows(stage_in_day, C[, "f_bIP"]-C_sem[, "sem_bIP"], stage_in_day, C[, "f_bIP"]+C_sem[, "sem_bIP"], length=0.1, angle=90, code=3,col = "blue") # error bar
	arrows(stage_in_day, C[, "f_neuron"]-C_sem[, "sem_neuron"], stage_in_day, C[, "f_neuron"]+C_sem[, "sem_neuron"], length=0.1, angle=90, code=3,col = "purple") # error bar
		
	# modeling results
	points(range_hour/24, C_in_hour[, "f_NE"], type = "l", col = "red")
	points(range_hour/24, C_in_hour[, "f_RG_SNP"], type = "l", col = "green")
	points(range_hour/24, C_in_hour[, "f_bIP"], type = "l", col = "blue")
	points(range_hour/24, C_in_hour[, "f_neuron"], type = "l", col = "purple")


	#
	# (Top right)
	# plot modeled cell cycle length
	#
	cellCycleTime.x = range_hour/24
	cellCycleTime.y = sapply(range_hour,Tc.func)
	
	plot(NA, xlab = "embryonic stage in day", ylab = "cell cycle length [hr]",
		xlim = c(min(cellCycleTime.x), max(cellCycleTime.x)), 
		ylim = c(11,22.5), main = "Cell-cycle length vs Embryonic stage")	
	points(cellCycleTime.x, cellCycleTime.y, type = "l", col = "red")
	points(Tc[, "stage"], Tc[, "Tc_all_dissociated"], col = "black",pch=19, cex = 2)
	arrows(Tc[, "stage"], Tc[, "Tc_all_dissociated"]-Tc[, "sem_Tc_all_dissociated"], Tc[, "stage"], Tc[, "Tc_all_dissociated"]+Tc[, "sem_Tc_all_dissociated"], length=0.1, angle=90, code=3,col = "black") # error bar


	#
	# (Bottom left)
	# plots of cell fraction modeling results vs cell division
	#
	stage_in_cellDiv = sapply(stage_in_day*24,Ct.func)
	
	plot(NA, xlab = "cell division number", ylab = "cell fraction",
		xlim = c(min(stage_in_cellDiv), max(stage_in_cellDiv)), 
		ylim = c(0, 1), main = "Cell fraction vs Division number")
	
	# raw data
	points(stage_in_cellDiv, C[, "f_NE"], col = "red",pch=19, cex = 2)
	points(stage_in_cellDiv, C[, "f_RG_SNP"], col = "green",pch=19, cex = 2)
	points(stage_in_cellDiv, C[, "f_bIP"], col = "blue",pch=19, cex = 2)
	points(stage_in_cellDiv, C[, "f_neuron"], col = "purple",pch=19, cex = 2)
	arrows(stage_in_cellDiv, C[, "f_NE"]-C_sem[, "sem_NE"], stage_in_cellDiv, C[, "f_NE"]+C_sem[, "sem_NE"], length=0.1, angle=90, code=3,col = "red") # error bar
	arrows(stage_in_cellDiv, C[, "f_RG_SNP"]-C_sem[, "sem_RG_SNP"], stage_in_cellDiv, C[, "f_RG_SNP"]+C_sem[, "sem_RG_SNP"], length=0.1, angle=90, code=3,col = "green") # error bar
	arrows(stage_in_cellDiv, C[, "f_bIP"]-C_sem[, "sem_bIP"], stage_in_cellDiv, C[, "f_bIP"]+C_sem[, "sem_bIP"], length=0.1, angle=90, code=3,col = "blue") # error bar
	arrows(stage_in_cellDiv, C[, "f_neuron"]-C_sem[, "sem_neuron"], stage_in_cellDiv, C[, "f_neuron"]+C_sem[, "sem_neuron"], length=0.1, angle=90, code=3,col = "purple") # error bar
		
	# modeling results
	range_hour_in_cellDiv = sapply(range_hour,Ct.func)
	points(range_hour_in_cellDiv, C_in_hour[, "f_NE"], type = "l", col = "red")
	points(range_hour_in_cellDiv, C_in_hour[, "f_RG_SNP"], type = "l", col = "green")
	points(range_hour_in_cellDiv, C_in_hour[, "f_bIP"], type = "l", col = "blue")
	points(range_hour_in_cellDiv, C_in_hour[, "f_neuron"], type = "l", col = "purple")


	#
	# (Bottom right)
	# Expansion coefficients
	# * Shown are robust estimations only *
	#
	plot(NA, xlab = "embryonic stage in day", ylab = "expansion coefficient", xlim = c(9, 18), ylim = c(0, 
		2), main = "Red:NE, Green:RG (NE+RG), Blue:RG (RG+IP+EN)", cex.main = 0.8)
	points(range_cellDiv_in_hour[ar1_sd < sd.threshold]/24, ar1[ar1_sd < sd.threshold], pch = 20, col = "green", cex = 2)
	points(range_cellDiv_in_hour[as_sd < sd.threshold]/24, as[as_sd < sd.threshold], pch = 20, col = "red", cex = 2)
	points(range_cellDiv_in_hour[ar2_sd < sd.threshold]/24, ar2[ar2_sd < sd.threshold], pch = 4, lwd=2,col = "blue", cex = 1)


}



#
# Misc functions
#
sizeGrWindow = function(width, height) {
	din = par("din")
	if ((din[1] != width) | (din[2] != height)) {
		dev.off()
		dev.new(width = width, height = height)
	}
}


printf <- function(...) cat(sprintf(...))
