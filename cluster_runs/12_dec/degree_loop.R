
degrees <- sample(1:9, 10000, prob = rev(c(seq(0.1,0.9,0.1))), replace = T)

breaks <- c()
i <- 1

while (i <= max(degrees)){
        k <- i
        if(length(x[x %in% i:k]) <= 1000){
                while (length(x[x %in% i:k]) <= 1000){
                        k <- k + 1
                }
                breaks[length(breaks)+1] <- k              
                i <- i + k
        }else{
                breaks[length(breaks)+1] <- i
                i <- i + 1
        }
}
