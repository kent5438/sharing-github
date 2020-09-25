# Rscript Scatter.r output/ebseq/AM1-D7,AM3-D7-AM1-2D-1.genMat.result.1 ./
library(scatterD3)
library(htmlwidgets)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
output<-args[2]

data <- read.table(fn, header = TRUE, sep = "\t", check.names=F)
C1Mean <-data$`sample1`
C2Mean <-data$`sample2`
FDR <-data$FDR
log2_Ctrl <- log2(C1Mean+1)
log2_Treat <- log2(C2Mean+1)

# R-squared calculation
correlation <- lm(log2_Treat~log2_Ctrl)
slope <- coef(correlation)["log2_Ctrl"]
intercept <- coef(correlation)["Intercept"]

# Plot three coefficient line
three_line <- data.frame(slope1 = c(slope = slope, intercept = intercept), slope2 = c(slope = 1, intercept = -1), slope3 = c(slope = 1, intercept = 1))
three_line_tr <- data.frame(t(three_line))

# Caption (custom tooltips)
tooltips <- paste("ID: <strong>", rownames(data),"</strong><br />FDR", data$FDR)

#output = paste0(dir, "/Scatter.html")
# Scatter Plot
scatter <- scatterD3(data = NULL, x = log2_Ctrl, y = log2_Treat, point_size = 10, xlim=c(0,20), ylim=c(0,20), hover_size = 2, col_var = FDR, lines = data.frame(slope = three_line_tr$slope, intercept = three_line_tr$intercept, stroke = "blue", stroke_width = 2, stroke_dasharray = 5), tooltip_text = tooltips)

# HTML output
saveWidget(scatter, file=output, selfcontained=FALSE)


