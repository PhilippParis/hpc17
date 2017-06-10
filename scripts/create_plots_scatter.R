library(ggplot2)
library(dplyr)

read_file_data <- function(fname, expid) {
  df <- read.table(fname, header=TRUE)
  df$expid <- expid
  df
}

df1 <- read_file_data("./data_scatter_min_alloc/scatter_out1.dat", 1)
df2 <- read_file_data("./data_scatter_min_alloc/scatter_out2.dat", 2)
df3 <- read_file_data("./data_scatter_min_alloc/scatter_out3.dat", 3)
df4 <- read_file_data("./data_scatter_min_alloc/scatter_out4.dat", 4)
df5 <- read_file_data("./data_scatter_min_alloc/scatter_out5.dat", 5)

df <- rbind(df1, df2, df3, df4, df5)

# compute median for each case
df1 <- df %>% group_by(test, count, expid) %>% summarize(median=median(runtime_sec)*1E6)
df1.1 <- df1 %>% summarize(mmedian=median(median))

p1 <- ggplot(df1.1, aes(x=count, y=mmedian, group=test, color=test)) + 
  geom_point() +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        legend.position="top") +
  ylab("latency [microseconds]") +
  xlab("message size [Bytes]")
plot(p1)


normalize_medians <- function(df) {
  df$rel_lat    <- df$median  / min(df[df$test=="MPI_Scatter",]$median)
  df
}

# compute relative latency, use min(median) of MPI_Gather as reference
df2 <- df1 %>% ungroup() %>% group_by(count) %>% do( normalize_medians(.) )
p2 <- ggplot(df2, aes(x=factor(count), y=rel_lat, fill=test)) + 
  geom_boxplot(position = position_dodge(width = .9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        legend.position="top") +
  ylab("relative latency") +
  xlab("message size [Bytes]")
plot(p2)

pdf("scatter32x16_runtime.pdf", width=6, height=4)
plot(p1)
dev.off()

pdf("scatter32x16_rel_runtime.pdf", width=6, height=4)
plot(p2)
dev.off()
