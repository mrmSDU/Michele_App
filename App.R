library(shiny)
library(dplyr)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(reshape2)
library(shiny)
library(feather)
library(grid)

d.counts <- as.data.frame(read_feather("./Data/d_counts.feather"))
d.average <- as.data.frame(read_feather("./Data/d_average.feather"))

ui <- fluidPage(
  
  titlePanel("Look up your favorite gene"),
  
  sidebarLayout(
    
    sidebarPanel(
      selectInput("gene", "Chose a gene:", 
                  choices=d.average$Gene %>% 
                    as.character() %>% 
                    unique() %>% sort() ),
      hr(),
      helpText("Chose one of the 18932 genes detectable in my RNA-seq dataset")
      
      
    ),
    
    mainPanel(
      plotOutput("genePlot") %>% withSpinner(color = "#0dc5c1")
    )
  )
  
)


# Define a server for the Shiny app
server <- function(input, output,session) {

  
  # Fill in the spot we created for a plot
  output$genePlot <- renderPlot({
    
    
    tmp.average <- d.average[ d.average$Gene == input$gene,]
    tmp.average$ID <- factor(tmp.average$ID, levels=tmp.average$ID)
    
    tmp.counts <- d.counts[ d.counts$Gene == input$gene,]
    tmp.counts$ID <- factor(tmp.counts$ID, levels=tmp.average$ID)
    
    
    p1 <- ggplot(tmp.average, aes(ID, Count, fill=Lineage))+
      geom_bar(stat="identity", alpha=0.65,
               col="grey4")+
      scale_fill_manual(values = c("Msc"="grey",
                                   "Ob"="blue",
                                   "Ad" = "red"))+
      geom_point(data=tmp.counts, aes(ID, Count), size=3, alpha=0.5)+
      theme_minimal()+
      theme(axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(1,1,3,2, "cm"))+
      ggtitle(paste0(tmp.average$Gene," mRNA levels"))+
      labs(y="Normalized Counts")
    
    
    # Code to override clipping
    gt <- ggplotGrob(p1)
    gt$layout$clip[gt$layout$name=="panel"] <- "off"
    grid.draw(gt)
    
    #Getting x position information within the viewport drawing bars
    grid.force()
    
    names <- grid.ls()$name
    HMmatch = grep("geom_rect.rect", names, value = TRUE)
    hm = grid.get(HMmatch)
    
    x_coordinates <- as.character(hm$x) %>% strsplit(split ="native") %>% unlist() %>% as.numeric()
    x_widths <- as.character(hm$width) %>% strsplit(split ="native") %>% unlist() %>% as.numeric()
    
    #Now drawing the text within the bar viewport
    seekViewport("panel.7-5-7-5")
    grid.draw(textGrob(label=tmp.average$Lineage, 
                       x = unit(x_coordinates+(x_widths/2),"native"),
                       y = unit(0, "native")))
    
    grid.draw(textGrob(label=tmp.average$TimePoint, 
                       x = unit(x_coordinates+(x_widths/2),"native"),
                       y = unit(-0.05, "native")))
    
    grid.draw(textGrob(label=tmp.average$Dex_Concentration, 
                       x = unit(x_coordinates+(x_widths/2),"native"),
                       y = unit(-0.1, "native")))
    
    grid.draw(textGrob(label=tmp.average$VitaminD_Concentration, 
                       x = unit(x_coordinates+(x_widths/2),"native"),
                       y = unit(-0.15, "native")))
    
    grid.draw(textGrob(label=c("Lineage", "Timepoint","[Dex]nM", "[VitD]nM"),
                       y=unit(c(0,-0.05,-0.1,-0.15), "native"),
                       x= -0.015, just = "right", gp=gpar(fontface="bold")))
    
    popViewport()
    
  })
}

shinyApp(ui = ui, server = server)
