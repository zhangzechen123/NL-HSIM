dat <- read.table("C:/research/paper/review/memory and time/a10000memory_A.csv",
                  sep = "|", header = TRUE, stringsAsFactors = FALSE)

dat$MaxRSS_MB <- as.numeric(sub("M$", "", dat$MaxRSS))
dat$AveRSS_MB <- as.numeric(sub("M$", "", dat$AveRSS))

summary_mem <- data.frame(
  n_rep            = sum(!is.na(dat$MaxRSS_MB)),
  mean_MaxRSS_MB   = mean(dat$MaxRSS_MB, na.rm = TRUE),
  sd_MaxRSS_MB     = sd(dat$MaxRSS_MB, na.rm = TRUE),
  median_MaxRSS_MB = median(dat$MaxRSS_MB, na.rm = TRUE),
  max_MaxRSS_MB    = max(dat$MaxRSS_MB, na.rm = TRUE),
  mean_AveRSS_MB   = mean(dat$AveRSS_MB, na.rm = TRUE),
  sd_AveRSS_MB     = sd(dat$AveRSS_MB, na.rm = TRUE)
)

print(summary_mem)



TIME <- read.table(
  "C:/research/paper/review/memory and time/a5000time_A.csv",
  sep = ",",
  header = FALSE,
  stringsAsFactors = FALSE
)


x <- as.numeric(TIME$V4)


summary_time <- data.frame(
  n        = length(x),
  mean     = mean(x),
  sd       = sd(x),
  median   = median(x),
  IQR      = IQR(x),
  min      = min(x),
  max      = max(x)
)

print(summary_time)
