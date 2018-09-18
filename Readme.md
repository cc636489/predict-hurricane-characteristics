# Predict Hurricane Characteristics using a Bayesian approach

A Fortran solution to predict future hurricane characteristics using joint probability methods with optimal sampling (JPMOS).

## Introduction

The JPMOS approach combines the following three inputs:
* The annual occurance rate λ of hurricanes on a site of interest.
* The joint probability distribution of all hurricane characteristics.
* The hurricane-generated surge η(x) at the site of interest, given the hurricane characteristics.

And formulated an annual exceedence probability function, based on the previous three inputs.
<p align="center">
<img src="/doc/Annual_occurrence_rate_probability.png">
</p>

This probability function basically can be formulated by asking the following three sub-questions:
- What is the probability for 1 hurricane to happen at reference point? (Integral f(x) in x space)
- What is the probability for 1 hurricane which exceed the designed surge η, to happen at reference point? (Integral f(x)*P(η) in x space)
- What is the probability of hurricanes which exceed the designed surge η, to happen at reference point, within a year? (Multiple by λ)

And it gives the following outputs:
* Annual exceedence probability of hurricanes.
* A set of synthetic hurricanes (optimizers), represented by predicted hurricane characteristics.


## Dependencies

The following libraries can be found at /library folder.

- nlopt
- quadpack

## Execution

* compile & run normal_optimization.f90
* compile & run normal_weights.f90
* compile & run physical_rosenblatt.f90

## Results
### Training data sets.
1. landfall location.
<p align="center">
<img src="/doc/Plot_landfall_location.png">
</p>

2. central pressure.
<p align="center">
<img src="/doc/Plot_central_pressure.png">
</p>

3. Radius of maximum wind.
<p align="center">
<img src="/doc/Plot_radius_max_wind.png">
</p>

4. forward speed
<p align="center">
<img src="/doc/Plot_forward_speed.png">
</p>

5. heading angle.
<p align="center">
<img src="/doc/Plot_heading_angle.png">
</p>

Where landfall location, central pressure, radius of max wind are given by historical data directly, whilst forward speed and heading angle has been calculated based on geophysical relations correspondingly.

### Statistical Modeling for 5 hurricane characteristics.
1. Annual occurrence rate λ. (Kernel estimation)

Using Gaussian shape kernels, and defined annual occurrence rate as:
<p align="center">
<img src="/doc/Annual_occurrence_rate_equation.png">
</p>
<p align="center">
<img src="/doc/Kernel_formula.png">
</p>
which gives, the number of hurricanes which makes landfall within 1 km range from the reference point in each year.

2. Heading angle. (Beta distribution)

This has been obtained through regression analysis.
<p align="center">
<img width="440" height="300" src="/doc/Heading_1.png"> <img width="440" height="300" src="/doc/Heading_2.png">
</p>


3. Radius of maximum wind. (negtively correlated to pressure deficit)

The probability density function of RMW, conditioned on pressure deficit, has been obtained through regression analysis.
It takes the form of log(rmw) = a+b*log(dp) 
<p align="center">
<img width="440" height="300" src="/doc/Radius_1.png"> <img width="440" height="300" src="/doc/Radius_2.png">
</p>

4. Forward speed. (Log-normal distribution)

This has been obtained through regression analysis. 
<p align="center">
<img width="440" height="300" src="/doc/Forward_1.png"> <img width="440" height="300" src="/doc/Forward_2.png">
</p>

5. Central pressure. (Weibull distribution )

This characteristics is of great importance. We formulated a log-likelihood cost function, and maximized it to obtain the parameters. Bootstrapping strategy has also been used to quantify the statistical uncertainty on the optimized parameters.
<p align="center">
<img src="/doc/objective_function_optimization3.png">
</p>

### JPMOS methods with comparing to FEMA annual report.

After solving the optimization problem, we obtained 5 optimized hurricane characteristics, which is comparable to FEMA annual report.
<p align="center">
<img width="600" height="280"   src="/doc/Optimization_process.png">
</p>

<p align="center">
<img width="440" height="300" src="/doc/Comparable_result_1.png"> <img width="440" height="180" src="/doc/Comparable_result_2.png">
</p>

