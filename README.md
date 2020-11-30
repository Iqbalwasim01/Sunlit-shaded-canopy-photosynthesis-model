# Sunlit-shaded-canopy-photosynthesis-model
For the R scripts/code to work you must have the following R packages installed:
required:'Bigleaf','pracma','foreach','dplyr'
optional:'ggplot2'- if running example simulation

# 1 Loading model functions
Once all necessary R packages have been installed. Run the following file to load all model functions:
'Canopy model functions.R'
This can be done with 'Ctrl+a' and pressing run:

# 2 Running model example
Download and open the script: 
'Example simulation-sugar beet.R'

Download the site environmental data from the file:
'Final.daily.data.csv'
Import the file as a 'text (base)' and specify the first row as headings. If the file is NOT imported as a 'text(base)' the functions may have problems recognising field names.

Using the 'Example simulation-sugar beet.R' run each line of code to run the simulation. 

NOTE: Make sure all files are within the same file directory and your R session is linked to the same directory. 

# 3 Brief explanation of canopy model function:
A lot of work has been put into condensing the canopy model into a simple excutable function: 

'Two.big.leaf.concept(File.length = nrow('filename'),No.kinetic.sets=,LAI=,crop.type="C3",leaf.angle="Spherical",pathway="C3",Tair =,wind =,humidity = x/100,Precip =
                             ,SWC_1 =,Rn=,H=,LE=,LMA=,pressure =,SZA =,SW_down =,Ca=,Sensor.height =,Kcat =,Kc =,Sc.o =,KcHA=,VmaxHA=,GamHA=)'

The function above takes an environnmental observation as a row within the dataframe whether this is a daily,hourly,yearly or monthly observation. No.kinetic.sets must be specified if you have imported another file with more than one set of Rubisco kinetics (i.e. Kcat, Kc Sc.o etc). If running only one set of kinetics then simply set it to '1' and enter the one set of kinetics directly into the function. 

File.length and No.kinetic.sets are used to give users a progress/counter from 0-100% on the R console. Simulations can sometimes take a few minutes to complete so its convient to have an indication of progress. 

A full list input meanings can be found:
'Canopy model function inputs.txt'






