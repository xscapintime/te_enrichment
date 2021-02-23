colors <- c("#FF3200", "#E9A17C", "#E9E4A6",
                "#1BB6AF", "#0076BB", "#172869")
tmp <- (grDevices::colorRampPalette(colors))(16)

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
            "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
            "#E34234", "#9A2A2A", "#630330", "#342D7E",
            "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

morecols <- (grDevices::colorRampPalette(mycols))(56) ### expand the platte!


normal <- c("#ebac23","#b80058", "#008cf9", "#006e00", "#00bbad", "#d163e6",
            "#b24502", "#ff9287", "#5954d6", "#00c6f8", "#00c6f8", "#00a76c",
            "#bdbdbd")

brew <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
          "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")

library(ggsci)
jco <-  pal_jco("default")(10)
d3_20 <- pal_d3("category20")(20)
simpsons <- pal_simpsons("springfield")(16)


## show palette
show_col(pal_d3("category20")(20))
show_col(pal_simpsons("springfield")(16))



## plot palette 
# https://bookdown.org/xiangyun/masr/sec-colors.html#subsec:color-palettes
pal <- function(n = 20, colors = colors, border = "light gray", ...) {
  colorname <- (grDevices::colorRampPalette(colors))(n)
  plot(0, 0,
       type = "n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, ...
  )
  rect(0:(n - 1) / n, 0, 1:n / n, 1, col = colorname, border = border)
}
par(mar = rep(0, 4))
pal(n = 16, colors = pal_simpsons("springfield")(16), xlab = "Colors from Peach to Pear", ylab = "")
    
