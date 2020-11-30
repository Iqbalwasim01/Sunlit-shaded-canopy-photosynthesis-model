# Sunlit/shaded canopy photosynthesis-model
For the R scripts/code to work you must have the following R packages installed:

required:'foreach','ggplot2'

Important notes:

Ignore this warning during model execution: 'Warning message: executing %dopar% sequentially: no parallel backend registered'. Seems to be an unknown problem with the foreach package. However, has not caused any problems with my modelling. 

Make sure all files are within the same file directory and your R session is linked to the same directory. 

The 'Canopy model functions.R' file has many redundant functions because of trial and error. For the final functions see the functions within the 'Two.big.leaf.concept()' function.  

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


# 3 Brief explanation of canopy model function:
The canopy model is a single executable function: 

'Two.big.leaf.concept(File.length = nrow('filename'),No.kinetic.sets=,LAI=,crop.type="C3",leaf.angle="Spherical",pathway="C3",Tair =,wind =,humidity = x/100,Precip =
                             ,SWC_1 =,Rn=,H=,LE=,LMA=,pressure =,SZA =,SW_down =,Ca=,Sensor.height =,Kcat =,Kc =,Sc.o =,KcHA=,VmaxHA=,GamHA=)'

The function above takes an environnmental observation as a row within the dataframe whether this is a daily,hourly,yearly or monthly observation. No.kinetic.sets must be specified if you have imported another file with more than one set of Rubisco kinetics (i.e. Kcat, Kc Sc.o etc). If running only one set of kinetics then simply set it to '1' and enter the one set of kinetics directly into the function. 

File.length and No.kinetic.sets are used to give users a progress/counter from 0-100% on the R console. Simulations can sometimes take a few minutes to complete so its convient to have an indication of progress. 

A full list input meanings can be found:
'Canopy model function inputs.txt'






