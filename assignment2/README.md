# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

#### Cat:<br/>
![alt text](Images/Q_1/Cat.GIF "Title")

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).

#### Cat:<br/>
![alt text](Images/Q_2/Cat.JPG "Title")

#### Luigi:<br/>
![alt text](Images/Q_2/Luigi.JPG "Title")

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).

#### Sphere:<br/>
##### Parameters: Res_15_WR_0.1_Epsilon_0.01_PD_1_SF_Unset <br/>
![alt text](Images/Q_3/Sphere_Res_15_WR_0.1_Epsilon_0.01_PD_1_SF_Unset.JPG "Title")
##### Parameters: Res_15_WR_0.1_Epsilon_0.01_PD_1_SF_1 <br/>
![alt text](Images/Q_3/Sphere_Res_15_WR_0.1_Epsilon_0.01_PD_1_SF_1.JPG "Title")

#### Cat:<br/>
##### Parameters: Res_20_WR_0.1_Epsilon_0.01_PD_2_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_20_WR_0.1_Epsilon_0.01_PD_2_SF_150.JPG "Title")
##### Parameters: Res_20_WR_0.1_Epsilon_0.01_PD_1_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_20_WR_0.1_Epsilon_0.01_PD_1_SF_150.JPG "Title")
##### Parameters: Res_20_WR_0.1_Epsilon_0.01_PD_0_SF_300 <br/>
![alt text](Images/Q_3/Cat_Res_20_WR_0.1_Epsilon_0.01_PD_0_SF_300.JPG "Title")
##### Parameters: Res_20_WR_0.1_Epsilon_0.01_PD_0_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_20_WR_0.1_Epsilon_0.01_PD_0_SF_150.JPG "Title")
##### Parameters: Res_30_WR_0.1_Epsilon_0.01_PD_0_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_30_WR_0.1_Epsilon_0.01_PD_0_SF_150.JPG "Title")
##### Parameters: Res_40_WR_0.1_Epsilon_0.01_PD_0_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_40_WR_0.1_Epsilon_0.01_PD_0_SF_150.JPG "Title")
##### Parameters: Res_50_WR_0.1_Epsilon_0.01_PD_0_SF_150 <br/>
![alt text](Images/Q_3/Cat_Res_50_WR_0.1_Epsilon_0.01_PD_0_SF_150.JPG "Title")

#### Luigi:<br/>
##### Parameters: Res_30_WR_0.05_Epsilon_0.01_PD_0_SF_10 <br/>
![alt text](Images/Q_3/Luigi_Res_30_WR_0.05_Epsilon_0.01_PD_0_SF_10.JPG "Title")
##### Parameters: Res_40_WR_0.05_Epsilon_0.01_PD_0_SF_10 <br/>
![alt text](Images/Q_3/Luigi_Res_40_WR_0.05_Epsilon_0.01_PD_0_SF_10.JPG "Title")

#### Conclusions:<br/>
##### 1. As it can be seen, from the Sphere, increasing the Smoothing Factor from 0.35 to 1 has improved significantly the quality but as be seen, from the Cat, to high Smoothing Factor is also bad <br/>
##### 2. As expected, from the Cat, increasing the Resolution from 20 to 50 has improved significantly the quality, using high resolution is possible only because of the Spatial Index <br/>
##### 3. As surprising it is, from the Cat, Polinom Degree of 0 is bather than 1 and 2 since it leaves a lot of artifacts and cats the ears of the cat <br/>

4) Theory question: Save your notes to assignment2/results and add a link to this page.
##### Answer can be found here
![a relative link](Results/Q_4)

### Optional Tasks

1) Save your notes and add a link to this page.

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with this method, and try to understand the differences. Report your findings.
