library(mgcv)
# Example data setup (from example in code)
set.seed(123); n <- 200
dat <- gamSim(1, n = n, scale = 2)
dat$year <- factor(sample(2001:2010, n, replace = TRUE))
dat$fac <- factor(sample(letters[1:4], n, replace = TRUE))
dat$log_catch <- rnorm(n, mean = (dat$y / 5) + as.numeric(dat$year)*0.05, sd = 0.5)

# Fit model
m1 <- gam(log_catch ~ year + ns(x0, 4) + s(x1) + s(x2) + fac, data = dat)



data$year <- data$f.myear

m1 <- gam(lcpue ~ year * Stock + timestep + f.month + vessel * f.NNets +
            s(BottomDepth, bs = "ts", k = 6) + s(log(FishingDuration), bs = "ts", k = 6) + 
            s(time.mid.hr, bs = "cc", k = 6) + Bycatch +
            f.NNets + moon.phase + 
            ti(lat, k = 10) + ti(long, k = 10) + ti(long, lat, k = c(6, 6)), 
            data = data, gamma = 1.4, select = TRUE) #family = Gamma(link = "log"), select = TRUE) #, method = "REML")

influ_obj <- influence_gam(m1, focus = "year", data = data)

influ_obj <- calculate(influ_obj, islog = T)
# Check progress with print:
print(influ_obj)

summary(influ_obj)

# Standardization plot
plot(influ_obj, type = "stan")
# or plot_stan(influ_obj)

# Step plot (faceted)
plot(influ_obj, type = "step")
# or plot_step(influ_obj)

# Influence plot (faceted)
plot(influ_obj, type = "influ")
# or plot_influ(influ_obj)

# Combined Step and Influence plot
plot(influ_obj, type = "step_influ")
# or plot_step_influ(influ_obj)

# CDI plot for a specific term (e.g., smooth of x0)
plot(influ_obj, type = "cdi", term = "s(x0)")
# or plot_cdi(influ_obj, term = "s(x0)")

# CDI plot for a factor term
plot(influ_obj, type = "cdi", term = "fac")
# or plot_cdi(influ_obj, term = "fac")

