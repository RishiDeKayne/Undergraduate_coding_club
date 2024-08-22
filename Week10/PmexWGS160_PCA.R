#160 p.mex PCAs - De-Kayne 2024

#load the necessary packages
#if you have not installed these already you will need to do so first!
library(ggplot2)
library("ggrepel")
library(tidyverse)

#set the working directory and load in background information 
#these filepaths are specific to my computer so you will have to change your own
#to match where you put the files
setwd("Re-sequencing_project_2024/PCA/")
background <- read.csv("Samples_selected_background_shape_col.csv", header = T, sep = ",")
str(background)

#make objects for file names (I typically do this to help with re-using this code)
eigenvec <- as.character("160_unfilt_out.eigenvec")
eigenval <- as.character("160_unfilt_out.eigenval")

#start by making pca df file from the input eigenvec file
pca <- read_table(eigenvec, col_names = FALSE)
#remove the first column of the pca DF leaving only one ID column
pca <- pca[,-1]

#and load in the eigenval values
eigenval <- scan(eigenval)

#we want to add a names column as well as any info we will use to determine shape/colour of points
#this is especially useful for differentiating populations as we will do here
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#we start by making empty vectors which we will later fill in a loop
#we want one per characteristic we wish to differentiate individuals by later
colour_plot <- rep(NA, length(pca$ind))
shape_plot <- rep(NA, length(pca$ind))

#here we go line by line in the PCA df and extract metadata from our background file for that indiv
#we then populate our empty vectors with that information
#in this case we are adding the corresponding colour and symbol to each individual
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$Sample.ID) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$colour)
  shape_plot[i] <- as.character(col_name$shape)
}

#now we make a new df called pca, effectively updating our df
#by adding the corresponding  colour and shape values vectors as columns
pca <- as.data.frame(data.frame(pca, colour_plot, shape_plot))

#we also get the pve values i.e. variance explained by each PC axis
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

#now we make plot - we start by looking at the loadings of each PC axis
#this script will plot the loadings as a nice barchart
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

#we can also calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#let's start with a simple plot
plot(pca$PC1, 
     pca$PC2)

#that looks good but we cant tell who is who - so now we can update our plot

#now we build on this - we specify the colour/type of point
#we can also write our own axis labels
plot(pca$PC1, 
     pca$PC2, 
     #now we have plotted the values we can add design details
     #this includes the colour of each point which is from our vector
     col = "black", 
     #the shape of the points
     pch = as.integer(pca$shape_plot),
     #the background colour which will be from our vector
     bg = as.character(pca$colour_plot),
     #the line width
     lwd = 2,
     xlab = "PC1 - 23%", ylab = "PC2 - 9%"
)

#even better we can automate the axis labels to include the correct pve fore each PC we plot
plot(pca$PC1, 
     pca$PC2, 
     #now we have plotted the values we can add design details
     #this includes the colour of each point which is from our vector
     col = "black", 
     #the shape of the points
     pch = as.integer(pca$shape_plot),
     #the background colour which will be from our vector
     bg = as.character(pca$colour_plot),
     #the line width
     lwd = 2,
     #and remembering to label our axes with the variance explained for each
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     cex = 1.75,
     #and a plot title
     main = "PC1 vs. PC2")

#next we plot PC1 vs PC3
plot(pca$PC1, 
     pca$PC3, 
     col = "black", 
     pch = as.integer(pca$shape_plot),
     bg = as.character(pca$colour_plot),
     lwd = 2,
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC1 vs. PC3")

#and finally PC2 vs. PC3
plot(pca$PC2, 
     pca$PC3, 
     col = "black", 
     pch = as.integer(pca$shape_plot),
     bg = as.character(pca$colour_plot), 
     lwd = 2,
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "PC2 vs. PC3")


#we can also add a custom legend to our plot
#we want to show each population and the colours and shapes that correspond to each
pop_list <- c("Pich", "Ixta", "Puya", "Taco", "non-sulfidic", "sulfidic")
pop_list_shape <- c(21, 23, 24, 22, "16", "16") 
pop_list_col <- c("black", "black", "black", "black", "#5DA5DA", "#DECF3F")

legend('topright', 
 legend = pop_list, 
 col = pop_list_col,
 pch = as.integer(pop_list_shape),
 pt.cex = 1,
 cex = 0.5) 

#if we want to label thd individuals we can add it with the following
#dont forget to change the 'PC' bit if you add this to another combination of PC axes
#you'll know if its wrong because the text won't be located over the points
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)

#a good habit to get into is exporting plots rather than saving from the plot pane
#we can do this a number of ways in R but the following is my favourite

#first we initiate the plot - you'll notice the plot no longer shows up in the pane
#it is a good idea to have a plot how you like it first and then save it as a .jpg or .tiff with the following
#first we initiate the plot - from here on every plotting command will be added to this plot
tiff("PmexWGS160_PC1_PC2.tiff", 
     height=8, 
     width=8, 
     units="in", 
     res=300, 
     compression="lzw")

#next we plot our plot
plot(pca$PC1, 
     pca$PC2, 
     #now we have plotted the values we can add design details
     #this includes the colour of each point which is from our vector
     col = "black", 
     #the shape of the points
     pch = as.integer(pca$shape_plot),
     #the background colour which will be from our vector
     bg = as.character(pca$colour_plot),
     #the line width
     lwd = 2,
     #and remembering to label our axes with the variance explained for each
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     #and a plot title
     main = "PC1 vs. PC2")

#and then we tell R we are done and it should essentially close the plot and 'save' the file
dev.off()

#### general data explanation
#plot a histogram of standard lengths to see distribution
hist(background$SL)
#now the same with mass
hist(background$Mass)

#now we can see if there is a relationship between the two
plot(background$SL, background$Mass)

#lets make the plot a bit nicer
plot(background$SL, 
     background$Mass,
     pch = 21,
     col = "black",
     bg = background$colour,
     xlab = "Standard Length",
     ylab = "Mass",
     cex = 1)

#let's make a linear model to see the relationship 
#this will be a very basic one and somethign like a logistic curve may fit better
#look into how to add that
our_model <- lm(background$Mass~background$SL)
#this will tell us the model parameters including the R2 and p-value
summary(our_model)

#we can also extract our p-value and summary stats like this:
summary(our_model)$r.squared
anova(our_model)$'Pr(>F)'[1]

#Now we will add a legend containing these
legend('topleft', 
       legend = c("R2 = 0.66", "P-value = 3.28x10-35"), 
       cex = 1)

#we can add this line to our plot
abline(lm(background$Mass~background$SL), col="red", lwd = 2)

#what if we want to see if the slope/line differs between sulfidic and non-sulfidic
#we can subset our data into two groups
NS <- subset(background, background$Sulfur.Fresh == "Fresh")
S <- subset(background, background$Sulfur.Fresh == "Sulfur")

#first for NS
our_model_NS <- lm(NS$Mass~NS$SL)
summary(our_model_NS)

summary(our_model_NS)$r.squared
anova(our_model_NS)$'Pr(>F)'[1]

#now we can plot the line 
abline(lm(NS$Mass~NS$SL), col="darkblue", lwd = 2, lty = 2)

#and now for S
our_model_S <- lm(S$Mass~S$SL)
summary(our_model_S)

summary(our_model_S)$r.squared
anova(our_model_S)$'Pr(>F)'[1]

#now we can plot the line 
abline(lm(S$Mass~S$SL), col="goldenrod", lwd = 2, lty = 2)

legend('topleft', 
       legend = c("R2 = 0.66", "P-value = 3.28x10-35", "R2 = 0.50", "P-value = 1.60x10-13", "R2 = 0.76", "P-value = 1.02x10-20"), 
       col = c("red", "white", "darkblue", "white", "goldenrod", "white"),
       lty = c(1,1, 2, 2, 2, 2),
       cex = 1)

#what does this plot tell us about the size distributions of sulfidic and non-sulfidic fish?

