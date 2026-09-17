# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05. Inferential statistics
# Compare betas to the residuals
# Has the GLM cpatured the variance in the signal well?
# Is there any remaining variance in the residuals?
# RQ: Is there a significant difference between residuals and betas?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a.) Scatter Plot: Plot Mean Beta vs. Mean Residual for each ROI.
#	  If there is a correlation, your model might be misspecified (e.g., missing a regressor).
# b.) Histogram of Residuals: Check if residuals are normally distributed (Gaussian).
#     GLM assumes residuals are normally distributed.
# c.) Spatial Map of Residual Variance: Plot the standard deviation of residuals across the brain.
#     High variance in specific regions might indicate artifacts not captured by your model.