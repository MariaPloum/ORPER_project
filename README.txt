Read me before running individualRegressions_ORPER.R script.

####### Project Description ###########

In this project, human subjects were presented with visual stimuli inside the MRI scanner. 

The stimuli are 4 colors that form a gradient (4 shades of blue, from dark blue to light blue).

Based on the experiment manipulations (not of interest here), we expect brain responses in specific brain regions 
to also follow a gradient.
This means that we expect more brain activation for the first stimulus in the gradient and gradually less activation 
for the next stimuli along the gradient.

Participants are split in two experimental groups. We expect that the brain activation gradient is more pronounced
in one group and more flat in the other group. 

####### Code Purpose #########

The purpose of this script is to fit a linear regression line on each participant's brain activation data, and then compare 
the fitted lines between the two groups of participants to verify our hypothesis.

####### Important concepts ########

- Brain activation is measured using a General Linear Model that calculates a coefficient for each condition 
(each stimulus in our case). This coefficient, otherwise called beta value, is what we use as a measure of 
brain activation.
So our dependent variable in this script is called "betas", which you can think of as brain activation magnitude.

- The brain regions that we expect to follow the activation gradient are all together our "Regions Of Interest", or ROIs.
ROIs are traditionally coded with an acronym that corresponds to a specific brain region, e.g., aINS = anterior insula, etc.

