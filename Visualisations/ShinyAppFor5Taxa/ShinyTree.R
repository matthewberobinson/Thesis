#1. Copy the following into the console without the double quotation marks and press enter "install.packages('shiny')"
install.packages("shiny")
install.packages("ggplot2")
install.packages("plyr")
install.packages("ape")
library(shiny)
library(ggplot2)
library(plyr)
library(ape)
options(shiny.error = NULL)

TrueTree = read.tree(text="((t1:0.1,t2:0.1):0.1,((t3:0.1,t4:0.1):0.1,t5:0.1):0.4);")
e <- read.table(file=paste('theHMCtestMattR.txt'))
load("petersen.Rdata")
if(TRUE){
e$V19 <- NULL
e$V18 <- NULL
e$V9 <- NULL
e$V10 <- NULL
e$V11 <- NULL
e$V12 <- NULL
e$V13 <- NULL
e$V14 <- NULL
e$V15 <- NULL
e$V1 <- NULL
e$V2 <- NULL
e$V3 <- NULL
e$V4 <- NULL
e$V5 <- NULL
colnames(e)[1] <- "V1"
colnames(e)[2] <- "V2"
colnames(e)[3] <- "V3"
colnames(e)[4] <- "V7"
colnames(e)[5] <- "V6"
}

#colnames(e)[1] <- "V1"
#colnames(e)[2] <- "V2"
#colnames(e)[3] <- "V3"
#colnames(e)[4] <- "V4"
#colnames(e)[5] <- "V5"
#colnames(e)[6] <- "V7"
#colnames(e)[7] <- "V6"

petersenplot <- function(){
  plot(0,0,col=0,ylim=c(-0.1,1.5),xlim=c(-0.4,1.5),asp=1,main="Petersen graph")
  for(i in 1:4){
    lines(c(petersen[[i]][1],petersen[[i+1]][1]),c(petersen[[i]][2],petersen[[i+1]][2]))
  }
  for(i in 1:5){
    lines(c(petersen[[i]][1],petersen[[i+5]][1]),c(petersen[[i]][2],petersen[[i+5]][2]))
  }
  lines(c(petersen[[5]][1],petersen[[1]][1]),c(petersen[[5]][2],petersen[[1]][2]))
  lines(c(petersen[[6]][1],petersen[[8]][1]),c(petersen[[6]][2],petersen[[8]][2]))
  lines(c(petersen[[6]][1],petersen[[9]][1]),c(petersen[[6]][2],petersen[[9]][2]))
  lines(c(petersen[[7]][1],petersen[[9]][1]),c(petersen[[7]][2],petersen[[9]][2]))
  lines(c(petersen[[7]][1],petersen[[10]][1]),c(petersen[[7]][2],petersen[[10]][2]))
  lines(c(petersen[[8]][1],petersen[[10]][1]),c(petersen[[8]][2],petersen[[10]][2]))
}

ui <- fluidPage(
  
  titlePanel("Tree Plots"),
  sliderInput("Iteration", label="Iteration", 1, length(e$V1), value=1, step = 1),
  checkboxInput("animation",label="Animate",T),
  fluidRow(
    splitLayout(cellWidths = c("50%", "50%"), plotOutput("internalsplot2"), plotOutput("internalsplot",width = "95%"))
  ),
  fluidRow(
    splitLayout(cellWidths = c("50%", "50%"), plotOutput("treeplot"), plotOutput("petersen"))
  ),
  fluidRow(
    splitLayout(cellWidths = c("50%", "50%"), plotOutput("acceptancerate"), plotOutput("likelihood"))
  ),
  plotOutput("densityplot")
)

server <- function(input, output) {
  
  output$internalsplot <- renderPlot({
    colourtopology <- c(rep(dist.topo(TrueTree,read.tree(text=toString(e$V1[1])))+1,length(e$V1)))
    for(i in 1:(length(e$V1)-1)){
      colourtopology[i+1] <- colourtopology[i] + dist.topo(read.tree(text=toString(e$V1[i])),read.tree(text=toString(e$V1[i+1])))
    }
    plot(e$V2,e$V3,type='l',xlim=c(0,1),ylim=c(0,1),xlab="first internal edge",ylab="second internal edge",cex.main=0.75,main="Plot of the Internal Edges: colour changes when you cross into another orthant",col=0)
    for(i in 1:(length(e$V1)-1)){
      if(input$animation){
        if(i<input$Iteration){
          lines(c(e$V2[i],e$V2[i+1]),c(e$V3[i],e$V3[i+1]),col=colourtopology[i])
        }
      }
      else{ lines(c(e$V2[i],e$V2[i+1]),c(e$V3[i],e$V3[i+1]),col=colourtopology[i])}
    }
    points(e$V2[1],e$V3[1])
    points(e$V2[input$Iteration],e$V3[input$Iteration],col=colourtopology[input$Iteration],pch=21,  bg=colourtopology[input$Iteration])
  })
  
  output$internalsplot2 <- renderPlot({
    colourtopology <- c(rep(dist.topo(TrueTree,read.tree(text=toString(e$V1[1])))+1,length(e$V1)))
    for(i in 1:(length(e$V1)-1)){
      colourtopology[i+1] <- dist.topo(TrueTree,read.tree(text=toString(e$V1[i+1])))+1
    }
    plot(e$V2,e$V3,type='l',xlim=c(0,1),ylim=c(0,1),xlab="first internal edge",ylab="second internal edge",cex.main=0.75,main="Plot of the Internal Edges: colour determines topological distance",col=0)
    for(i in 1:(length(e$V1)-1)){
      if(input$animation){
        if(i<input$Iteration){
          lines(c(e$V2[i],e$V2[i+1]),c(e$V3[i],e$V3[i+1]),col=colourtopology[i])
        }
      }
      else{ lines(c(e$V2[i],e$V2[i+1]),c(e$V3[i],e$V3[i+1]),col=colourtopology[i])}
     
    }
    points(e$V2[1],e$V3[1])
    points(e$V2[input$Iteration],e$V3[input$Iteration],col=colourtopology[input$Iteration],pch=21,  bg=colourtopology[input$Iteration])
  })
  
  output$treeplot <- renderPlot({
    plot(read.tree(text=toString(e$V1[input$Iteration])))
  })
  
  output$petersen <- renderPlot({
    petersenplot()
    if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                 read.tree(text=toString(
                   "(t2:0.1000000,(t3:0.1000000,(t5:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[5]][1],petersen[[1]][1]),c(petersen[[5]][2],petersen[[1]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t2:0.1000000,(t1:0.1000000,(t5:0.1000000,t4:0.1000000):0.1):0.1,t3:0.1000000);")))==0){
      lines(c(petersen[[1]][1],petersen[[6]][1]),c(petersen[[1]][2],petersen[[6]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t3:0.1000000,(t2:0.1000000,(t5:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[1]][1],petersen[[2]][1]),c(petersen[[1]][2],petersen[[2]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t3:0.1000000,(t5:0.1000000,(t2:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[3]][1],petersen[[2]][1]),c(petersen[[3]][2],petersen[[2]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t3:0.1000000,(t1:0.1000000,(t2:0.1000000,t4:0.1000000):0.1):0.1,t5:0.1000000);")))==0){
      lines(c(petersen[[3]][1],petersen[[4]][1]),c(petersen[[3]][2],petersen[[4]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t3:0.1000000,(t4:0.1000000,(t2:0.1000000,t1:0.1000000):0.1):0.1,t5:0.1000000);")))==0){
      lines(c(petersen[[4]][1],petersen[[5]][1]),c(petersen[[4]][2],petersen[[5]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t3:0.1000000,(t4:0.1000000,(t5:0.1000000,t2:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[2]][1],petersen[[7]][1]),c(petersen[[2]][2],petersen[[7]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t5:0.1000000,(t3:0.1000000,(t2:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[3]][1],petersen[[8]][1]),c(petersen[[3]][2],petersen[[8]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t4:0.1000000,(t2:0.1000000,(t5:0.1000000,t3:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[4]][1],petersen[[9]][1]),c(petersen[[4]][2],petersen[[9]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t2:0.1000000,(t5:0.1000000,(t3:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[5]][1],petersen[[10]][1]),c(petersen[[5]][2],petersen[[10]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t4:0.1000000,(t3:0.1000000,(t5:0.1000000,t2:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[7]][1],petersen[[9]][1]),c(petersen[[7]][2],petersen[[9]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t2:0.1000000,(t1:0.1000000,(t3:0.1000000,t4:0.1000000):0.1):0.1,t5:0.1000000);")))==0){
      lines(c(petersen[[7]][1],petersen[[10]][1]),c(petersen[[7]][2],petersen[[10]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t5:0.1000000,(t4:0.1000000,(t3:0.1000000,t2:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[6]][1],petersen[[8]][1]),c(petersen[[6]][2],petersen[[8]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t4:0.1000000,(t5:0.1000000,(t2:0.1000000,t3:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[6]][1],petersen[[9]][1]),c(petersen[[6]][2],petersen[[9]][2]),col=2)
    }
    else if(dist.topo(read.tree(text=toString(e$V1[input$Iteration])),
                      read.tree(text=toString(
                        "(t5:0.1000000,(t2:0.1000000,(t3:0.1000000,t4:0.1000000):0.1):0.1,t1:0.1000000);")))==0){
      lines(c(petersen[[8]][1],petersen[[10]][1]),c(petersen[[8]][2],petersen[[10]][2]),col=2)
    }
    else{print("Not 5 taxa tree")}
  })
  
  output$acceptancerate <- renderPlot({
    if(input$animation){
    plot(e$V7[1:input$Iteration],type='l',ylab="acceptance probability",main="Acceptance Probability",cex.main=0.75)
    }
    else{
      plot(e$V7,type='l',ylab="acceptance probability",main="Acceptance Probability",cex.main=0.75)
    }
      })
  output$likelihood <- renderPlot({
    if(input$animation){
    plot(e$V6[1:input$Iteration],type='l',ylim=c(-4600,-4400),ylab="loglikelihood",main="loglikelihood",cex.main=0.75)
    }
    else{
      plot(e$V6,type='l',ylim=c(-4600,-4400),ylab="loglikelihood",main="loglikelihood",cex.main=0.75)
    } 
     })
  
  output$densityplot <- renderPlot({
    if(input$animation){
    dat <- as.data.frame(array(c(e$V2[1:input$Iteration],e$V3[1:input$Iteration]),dim=c(length(e$V2[1:input$Iteration]),2)))
    m <- ggplot(dat, aes(x = V1, y = V2)) +
      geom_point() + xlim(0, 1) + ylim(0, 1)
    m + geom_density2d()
    }
    else{
      dat <- as.data.frame(array(c(e$V2,e$V3),dim=c(length(e$V2),2)))
      m <- ggplot(dat, aes(x = V1, y = V2)) +
        geom_point() + xlim(0, 1) + ylim(0, 1)
      m + geom_density2d()
    }
  })
  
}

shinyApp(ui = ui, server = server)