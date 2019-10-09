# CCS-Calculation

Hello All!

This algorithm intends to simplify CCS calculations based on raw data measured on Waters Synapt G2 files.

You can directly upload raw data files into the user interface and note the temperature, pressure and type of drift gas you used
to measure your data. You also need to give the mass and charge of the peaks you want to analyse as input.

Afterwards the algorithm will follow these simple steps:

- translate the raw data files into xyz text files
- search for mass and charge of the peaks you want
- if it finds a peak, it will copy the data and extract the arrival time distribution for all measured voltages
- it will try to fit the data with a gaussian fit:
  either one, two or three fits depending on an error limit which is an empirical value
- each gaussian fit results in a ceontroid drift time value
- a typical CCS measurements contains up to eight drift voltages, therefore we get up to eight different centroids
- the algorithm will plot the reversed voltage against the ceontroid drift times to extract the mobility
  for each peak out of the linear slope
- using the Mason-Schamp equation, we can then calculate the CCS out of all given parameters
- this loop repeats for each given mass 
- at the end, all results are written to an csv file which can be watched by excel

For further information, please have a look at the manual.pdf or write me
