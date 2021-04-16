# Create test dataset with three clusters of 100 cells each
test.df <- data.frame(cluster = factor(rep(c(1, 2, 3), each = 100)))

# Create 6 donors that are cases or controls
donors.df <- data.frame(donor = rep(paste("Donor", LETTERS[1:6], sep = "_"), each = 50),
                        sex = rep(c("M", "F", "M", "F", "F", "M"), each = 50),
                        status = factor(rep(c("Case", "Case", "Control", "Control", "Case", "Control"), each = 50)))

# Now make cluster 1 mostly case, cluster 2 mostly controls, etc
cases <- donors.df[donors.df$status == "Case",]
cases <- cases[sample(nrow(cases)),]
controls <- donors.df[donors.df$status == "Control",]
controls <- controls[sample(nrow(controls)),]

test.df <- cbind(rbind(cases[1:75,], controls[1:25,], cases[76:115,],
                       controls[26:85,], cases[116:150,], controls[86:150,]),
                 test.df)

# Test set call
library(lme4)
MASC(dataset = test.df,
     cluster = test.df$cluster,
     contrast = "status",
     random_effects = "donor",
     fixed_effects = "sex",
     verbose = TRUE, save_models = FALSE)
