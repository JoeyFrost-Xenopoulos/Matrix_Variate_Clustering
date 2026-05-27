sim_function<-function(Ea,Eb,Ec,x,basis,fun){
  a=rnorm(30,Ea,0.4^2)
  b=rnorm(30,Eb,0.4^2)
  c=rnorm(30,Ec,0.4^2)
  y=a*fun(x-b)+c
  return(data.frame(basis,y=y))}

plotImage<-function(m){
  df <- as.data.frame(as.table(m))  
  df$Var1 <- as.numeric(df$Var1)
  df$Var2 <- as.numeric(df$Var2)
  
  ggplot(df, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile() +
    scale_y_reverse() +       
    coord_fixed() +          
    scale_fill_gradient(low = "white", high = "black") +  
    theme_void()   +theme(legend.position = "none")
}

make_shape<-function(shape){
  if(shape=="heart"){
    m <- matrix(0, 12, 12)
    m[3, 4:5] <- 0.8
    m[3, 8:9] <- 0.8
    m[4, 3:6] <- 0.9
    m[4, 7:10] <- 0.9
    m[5, 2:11] <- 1.0
    m[6, 2:11] <- 1.0
    m[7, 3:10] <- 0.9
    m[8, 4:9] <- 0.8
    m[9, 5:8] <- 0.6
    m[10, 6:7] <- 0.4
  }else{
    m <- matrix(0, 12, 12)
    for (i in 4:9) {
      for (j in 4:9) {
        m[i, j] <- runif(1, 0.7, 1)  # random grayscale inside square
      }
    }}
  noise <- matrix(rnorm(length(m), mean = 0, sd = 0.1), nrow = nrow(m))
  noisy <- pmin(pmax(m + noise, 0), 1)  # clamp between 0 and 1
  return(noisy)
}
