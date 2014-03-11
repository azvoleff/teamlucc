context("cloud_remove")

# Test that areas that are NOT cloud masked have identical reflectance to those 
# in the input image

# Test that areas that ARE cloud masked have reflectance that do not vary 
# largely from the clear reflectance

plotRGB(cloudy, 3, 2, 1)
plotRGB(clear, 3, 2, 1)
plotRGB(clout_filled, 3, 2, 1)
