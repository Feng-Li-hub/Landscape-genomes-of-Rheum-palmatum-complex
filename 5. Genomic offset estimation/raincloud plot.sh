setwd("F:\\王聪颖\\景观_feng\\云雨图")
# Read the CSV file and convert 'Order2' column to an ordered factor
data = read.csv("core_GFGO126.csv")
data %>%
  mutate(Order2 = factor(Order2, levels = c("GO_5026_GFwest", "GO_5026_GFeast", "GO_9026_GFwest", "GO_9026_GFeast"))) -> df

# Manually define colors for each category
ordercolors <- c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")

# Create a ggplot visualization
ggplot(data = df, aes(x = Order2, y = GO * 100, fill = Order2)) +
  # Create a half violin plot (right side), with transparency and no border color
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  # Create a half boxplot (right side), remove error bars, set width and line thickness
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  # Create a half point plot (left side), set point shape, size, and color
  geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
  # Set fill color mapping based on the predefined colors
  scale_fill_manual(values = ordercolors) +
  # Set y-axis limits from 0 to 10 and remove expansion space
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
  # Set y-axis label and remove x-axis label
  labs(y = "GO (%)", x = NULL) +
  # Use a classic theme without background grid lines
  theme_classic() +
  # Customize the plot theme
  theme(
    # Set legend position at the bottom
    legend.position = "bottom",
    # Set axis title font size and color
    axis.title = element_text(size = 16, color = "black"),
    # Set y-axis text style (font size and color)
    axis.text.y = element_text(size = 13, color = "black"),
    # Set x-axis text style (font size and color)
    axis.text.x = element_text(size = 11, color = "black")
  )




data = read.csv("core_GFGO585.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5085_GFwest", "GO_5085_GFeast", "GO_9085_GFwest", "GO_9085_GFeast"))) -> df  ###重复上面操作
ordercolors <- c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"), 
axis.text.x = element_text(size = 11, color = "black")  
)


data = read.csv("pan_GFGO126.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5026_GFwest", "GO_5026_GFeast", "GO_9026_GFwest", "GO_9026_GFeast"))) -> df
ordercolors <- c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"),  
axis.text.x = element_text(size = 11, color = "black")  
)


data = read.csv("pan_GFGO585.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5085_GFwest", "GO_5085_GFeast", "GO_9085_GFwest", "GO_9085_GFeast"))) -> df  ###重复上面操作
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"), 
axis.text.x = element_text(size = 11, color = "black")  
)


data = read.csv("core_RDAGO126.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5026_RDAwest", "GO_5026_RDAeast", "GO_9026_RDAwest", "GO_9026_RDAeast"))) -> df
ordercolors <- c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"), 
axis.text.x = element_text(size = 11, color = "black") 
)


data = read.csv("core_RDAGO585.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5085_RDAwest", "GO_5085_RDAeast", "GO_9085_RDAwest", "GO_9085_RDAeast"))) -> df
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"), 
axis.text.x = element_text(size = 11, color = "black")  
)


data = read.csv("pan_RDAGO126.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5026_RDAwest", "GO_5026_RDAeast", "GO_9026_RDAwest", "GO_9026_RDAeast"))) -> df
ordercolors <- c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"),  
axis.text.x = element_text(size = 11, color = "black")  
)


data = read.csv("pan_RDAGO585.csv")
data %>%
mutate(Order2 = factor(Order2, levels = c("GO_5085_RDAwest", "GO_5085_RDAeast", "GO_9085_RDAwest", "GO_9085_RDAeast"))) -> df
ggplot(data = df,
aes(x = Order2, y = GO * 100, fill = Order2)) +
geom_half_violin(side = "r", color = NA, alpha = 0.35) +
geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
scale_fill_manual(values = ordercolors) +
scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
labs(y = "GO (%)", x = NULL) +
theme_classic() +
theme(
legend.position = "bottom",
axis.title = element_text(size = 16, color = "black"),
axis.text.y = element_text(size = 13, color = "black"), 
axis.text.x = element_text(size = 11, color = "black")  
)

