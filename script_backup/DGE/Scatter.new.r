# Rscript Scatter.r output/ebseq/AM1-D7,AM3-D7-AM1-2D-1.genMat.result.1 ./
library(scatterD3)
library(htmlwidgets)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
dir<-args[2]

data <- read.table(fn, header = TRUE, sep = "\t")
C1Mean <-data$C1Mean
C2Mean <-data$C2Mean
PValue <-data$PPEE
log2_C1Mean <- log2(C1Mean+1)
log2_C2Mean <- log2(C2Mean+1)

# R-squared calculation
correlation <- lm(log2(C1Mean+1)~log2(C2Mean+1))
slope <- coef(correlation)["log2(C2Mean + 1)"]
intercept <- coef(correlation)["(Intercept)"]

# Plot three coefficient line
three_line <- data.frame(slope1 = c(slope = slope, intercept = intercept), slope2 = c(slope = 1, intercept = -1), slope3 = c(slope = 1, intercept = 1))
three_line_tr <- data.frame(t(three_line))

# Caption (custom tooltips)
tooltips <- paste("Gene: <strong>", data$Gene,"</strong><br />p-value", PValue)

output = paste0(dir, "/Scatter.html")
# Scatter Plot
scatter <- scatterD3(data = NULL, x = log2_C2Mean, y = log2_C1Mean, point_size = 10, xlim=c(0,20), ylim=c(0,20), hover_size = 2, col_var = PValue, lines = data.frame(slope = three_line_tr$slope.log2.C2Mean...1., intercept = three_line_tr$intercept..Intercept., stroke = "blue", stroke_width = 2), tooltip_text = tooltips)

# HTML output
saveWidget(scatter, file=output, selfcontained=FALSE)


