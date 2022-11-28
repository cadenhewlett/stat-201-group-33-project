#install.packages("pracma")

### IMPORT 201 LIBRARIES ###
library(tidyverse)
library(broom)
library(repr)
library(digest)
library(infer)
library(gridExtra)
library(dplyr)
### Import Additional Libraries ###
library(MASS)
library(pracma)

#getwd()
df = read.csv("https://raw.githubusercontent.com/cadenhewlett/stat-201-group-33-project/main/monthly_20221102T0350.csv")

# Rename the columns according to original order (which was lost)
colnames(df) = c(
  "ID","PARAM",	"TYPE",	"Year",	"Jan",	"SYM",
  "Feb",	"SYM",	"Mar",	"SYM",	"Apr",	"SYM",
  "May",	"SYM",	"Jun",	"SYM",	"Jul",	"SYM",	"Aug",	
  "SYM",	"Sep",	"SYM",	"Oct",	"SYM",	"Nov",	"SYM",
  "Dec",	"SYM",	"Mean"
)


rowMax = function(i, threeDF){ 
  v = c( threeDF[i, 1], # get each month
         threeDF[i, 2],
         threeDF[i, 3])
  max(as.numeric(v))} # return the max of the three

# Preview the data frame's relevant columns.
head(df[colnames(df)!='SYM' & colnames(df)!='PARAM' & colnames(df)!='TYPE'])

# Create another cross-column filter, like we did before.
datafilter = (   df$Jul != "" & df$Mean != "" & 
                   df$Jun != "" & df$Aug != "" & 
                   df$Jun != "Jun" &  df$Aug != "Aug" & 
                   df$Jul != "Jul" &  df$Mean != "Mean"& 
                   df$TYPE != "8" )

# Filter the original dataframe according to this list
jun_filtered  =  df$Jun[datafilter]
jul_filtered  =  df$Jul[datafilter]
aug_filtered  =  df$Aug[datafilter]
mean_filtered =  df$Mean[datafilter] # with the addition of the mean col.

# Create new dataframe with relevant columns
computeDF = data.frame(jun_filtered, jul_filtered, aug_filtered, mean_filtered)


# Define a new column: result of the max function applied to compute df
computeDF[, "max"] <- sapply(1:nrow(computeDF), rowMax, threeDF = computeDF)

# Reast relevant columns to numeric
computeDF$max = as.numeric(computeDF$max)
computeDF$mean_filtered = (as.numeric(computeDF$mean_filtered))

run = 2002
###################### START LOOP HERE #####################
for(i in 1:39){
  run = 2002 + i
  set.seed(seed = run)
# take a random sample of size 20 from the population
trainDF = computeDF %>% rep_sample_n(size = 30) %>% ungroup() %>% dplyr::select(-replicate)

# multiply x by a factor of 0.1 to reduce matrix magnitude (next section)
train_x = 0.1*trainDF$max
train_y = 0.5*trainDF$mean_filtered

################ LOG REGRESSION ################
x = train_x
y = train_y

A_log = matrix(c(train_x^0, log(train_x)), ncol = 2)

# Solve the matrix equation c = (A^-1)y to find coefficients a,b,c and d
c_log = tryCatch(
  # try to solve with QR Decomposition
  { qr.solve(A_log, y) },
  #i f the matrix is singular, compute pseudoinverse
  error = function(e) {
    message('Warn: Matrix is Singular, using Pseudoinverse Solution to LSE')
    # compute the pseudoinverse c = A^+ y
    pinv(A_log) %*% y
  })


logtoCurve <- function(x_val){
  # this is equivalent to alog(x) + b
  (c_log[2]*(log(x_val)) + c_log[1])
}

################ CUBE REGRESSION ################
# Build Matrix from system of Minimization Equations
A_cube = matrix( c( sum(x^6), sum(x^5), sum(x^4), sum(x^3), 
                    sum(x^5), sum(x^4), sum(x^3), sum(x^2),
                    sum(x^4), sum(x^3), sum(x^2), sum(x^1),
                    sum(x^3), sum(x^2), sum(x^1), sum(x^0)), 
                 nrow = 4, ncol = 4)

# Build "Results" Vector, y, which is 4x1.
y_cube = matrix( c( sum(x^3 * y), 
                    sum(x^2 * y), 
                    sum(x * y), 
                    sum(y)), # equivalent to x^0*y
                 nrow = 4, ncol = 1)


# Solve the matrix equation c = (A^-1)y to find coefficients a and b
c_cube = tryCatch(
  # try to solve with QR Decomposition, R1 x = Q1^T y
  { qr.solve(A_cube, y_cube) },
  #i f the matrix is singular, compute pseudoinverse
  error = function(e) {
    message('Warn: Matrix is Singular, using Pseudoinverse Solution to LSE')
    # compute the pseudoinverse c = A^+ y
    pinv(A_cube) %*% y_cube
  })

# This function determines the values of the regression line at input x
cubetoCurve <- function(x_val){
  # this is equivalent to ax^3 + bx^2 + cx + d
  (c_cube[1]*(x_val^3) + c_cube[2]*(x_val^2) + c_cube[3]*(x_val) + c_cube[4])
}

################ LINEAR SIN REGRESSION ################
A_sin = matrix(c (train_x^0, train_x, sin(train_x)), ncol = 3) 

c_sin = qr.solve(A_sin, y)

sintoCurve <- function(x_val){
  # this is equivalent to alog(x) + b
  (c_sin[3]*sin(x_val) + c_sin[2]*(x_val) + c_sin[1])
}


################ LINEAR REGRESSION ################
# determine slope estimate from sample data
beta_hat = sum((x - mean(x))*(y - mean(y))) / sum((x - mean(x))^2)
# determine intercept estimate from sample data
alpha_hat = mean(y) - beta_hat*mean(x)

# Determines the values of the regression line at input x
lineartoCurve <- function(x_val){
  beta_hat*(x_val) + alpha_hat}

#### Regression Summary ####
# Reshape dataframe columns according to computational simplifications
modified_fn_input =  0.10*as.numeric(computeDF$max)
modified_df_actual = 0.5*(as.numeric(computeDF$mean_filtered))

resid_sin = (sintoCurve(modified_fn_input) - modified_df_actual)
resid_cub = (cubetoCurve(modified_fn_input) - modified_df_actual)
resid_log = (logtoCurve(modified_fn_input) - modified_df_actual)
resid_lin = (lineartoCurve(modified_fn_input) - modified_df_actual)

# Summarize our findings in a Data Frame, and make it look like a table
pt2_Results = data.frame(
  # by definition RSS is the sum of the squared residuals
  Mean_Residual = c(
    mean(resid_sin),  
    mean(resid_cub), 
    mean(resid_log), 
    mean(resid_lin),
          run))
  # also find the avg. residual for further comparisons


# Set Row Names to Improve Clarity
rownames(pt2_Results) = c("Sinusoidal Linear","Cubic",
                          "Logarithmic","Linear", "Seed")
#?write.table
#t(pt2_Results[,1])
# # pt2_Results$RSS
write.table(t(pt2_Results), file = "average_residual_regression_results.csv",
            append = T, sep = ",", col.names = F)}
#write.table(t(pt2_Results), file = "average_residual_regression_results.csv", sep = ",")
#}
# write.table(row, file = csv_fname, sep = ",",
#             append = TRUE, quote = FALSE,
#             col.names = FALSE, row.names = FALSE)