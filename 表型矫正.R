library(data.table)
dat_dis <- read.csv("../phe.csv", header = TRUE, na.strings = "NA")

fit <- lm(survival~Tank,data = dat_dis)
survivalmean <- mean(dat_dis$survival,na=TRUE)
survival <- dat_dis[,c(1,2,6)]
survival$adjust <- survival$survival - predict(fit) + survivalmean
write.csv(survival,file="../survivaladjust.csv",row.names = FALSE)

fit2 <- lm(bw~Tank*Sex,data = data)
bwmean <- mean(data$bw,na=TRUE)
bw <- data[,c(1,2,3,6)]
bw$adjust <- bw$bw - predict(fit2) + bwmean
write.csv(bw,file="bwadjust.csv",row.names = FALSE)
