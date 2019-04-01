
ggmaplot <- function(tdf, 
                 xvar,
                 yvar,
                 adjpvalvar,
                 Gene.names,
                 fdr.range.cut = c(0,0.001,0.01,0.05,1)) {
  require(ggplot2)
  require(RColorBrewer)
  require(plotly)
  source('I:/Ball_Kerri/Staurosporine_pDrive/txt_all_temp_staurosporine/DataAnalysis/fdrfunctions.R')
#  browser()
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  adjpvalvar <- enquo(adjpvalvar) 
  adjpvalvarcut <- fdrcut(tdf %>% select(!! adjpvalvar) %>% pull(), fdr.range.cut = fdr.range.cut)
  rcols <- rev(brewer.pal(n = 5, name = "BuPu"))[1:(length(fdr.range.cut) - 1)]
  gcol <- brewer.pal(n = 8, name = "Greys")[3]
  
  ggy <- tdf %>% ggplot(mapping = aes_string( quo_name(xvar), quo_name(yvar))) +
    geom_point(aes(text = Gene.names,
                   color = adjpvalvarcut)) +
    scale_colour_manual(values = c(gcol, rev(rcols))) +
    geom_vline(xintercept = 0,color = "grey") +
    theme_classic()
  
  
  tdf.thresh10pct <- tdf %>% filter(!! adjpvalvar <= 0.1 & !! adjpvalvar > 0.001)
    
 
  
  tdf.thresh1pct <- tdf %>% filter(!! adjpvalvar <= 0.0001 | (!! adjpvalvar <= 0.01 & (!! xvar > percentile90[1,1] | !! xvar < -percentile90[1,1])))
    

  ggply <- ggplotly(ggy, tooltip = "text")
  
  ggply <- ggply %>% add_annotations(x = tdf.thresh10pct %>% select(!! xvar) %>% pull(),
                                     y  = tdf.thresh10pct %>% select(!! yvar) %>% pull(),
                                     text = tdf.thresh10pct$Gene.name,
                                     font = list(family = 'Arial',
                                                 size = 10,
                                                 color = 'black'),
                                     showarrow = T,
                                     visible = F, 
                                     arrowhead = 4,
                                     arrowwidth = 0.5,
                                     arrowcolor = 'grey',
                                     arrowsize = .5,
                                     ax = 20,
                                     ay = -20,
                                     clicktoshow = "onoff") %>% 
    add_annotations(x = tdf.thresh1pct %>% select(!! xvar) %>% pull(),
                    y  = tdf.thresh1pct %>% select(!! yvar) %>% pull(),
                    text = tdf.thresh1pct$Gene.name,
                    font = list(family = 'Arial',
                                size = 10,
                                color = 'black'),
                    showarrow = T,
                    visible = T, 
                    arrowhead = 4,
                    arrowwidth = 0.5,
                    arrowcolor = 'grey',
                    arrowsize = .5,
                    ax = 20,
                    ay = -20,
                    clicktoshow = "onoff")
  # 
  # ggply <- ggply %>% add_annotations(x = tdf.thresh10pct$abcam,
  #                                     y  = tdf.thresh10pct$avgabcam,
  #                                     text = tdf.thresh10pct$Gene.name,
  #                                     font = list(family = 'Arial',
  #                                                 size = 10,
  #                                                 color = 'black'),
  #                                     showarrow = T,
  #                                     visible = F, 
  #                                     arrowhead = 4,
  #                                     arrowwidth = 0.5,
  #                                     arrowcolor = 'grey',
  #                                     arrowsize = .5,
  #                                     ax = 20,
  #                                     ay = -20,
  #                                     clicktoshow = "onoff") %>% add_annotations(x = tdf.thresh1pct$abcam,
  #                                                                                y  = tdf.thresh1pct$avgabcam,
  #                                                                                text = tdf.thresh1pct$Gene.name,
  #                                                                                font = list(family = 'Arial',
  #                                                                                            size = 10,
  #                                                                                            color = 'black'),
  #                                                                                showarrow = T,
  #                                                                                visible = T, 
  #                                                                                arrowhead = 4,
  #                                                                                arrowwidth = 0.5,
  #                                                                                arrowcolor = 'grey',
  #                                                                                arrowsize = .5,
  #                                                                                ax = 20,
  #                                                                                ay = -20,
  #                                                                                clicktoshow = "onoff")
  # 
  
  return(list(maplot.ggplot=ggy,
              maplot.plotly = ggply))
  # ggplot(mpg, aes(displ, hwy)) +
  #   geom_point(aes(color = class)) +
  #   geom_smooth(se = FALSE) +
  #   theme_bw()
  
}

# Try:
# myggma <- tophits %>% ggmaplot(xvar = abcam, yvar = avgabcam, adjpvalvar = abcam.adj.P.Val)
# ggplotly(myggma, tooltip = "text")
