reconnection_project
====================
Matlab files for the project in the reconnection course at IRF Uppsala. 

Authors: Sergio Toledo-Redondo and Andreas Johlander

Installation
------------
Go to your Matlab working directory and write in the terminal:
> git clone git://github.com/ajohlander/reconnection_project.git

[irfu-matlab](https://github.com/irfu/irfu-matlab "IRFU's github") is required to run this function.

Test
-----
Run "example1.m", which downloads magnetic field data and position data from Cluster Science Archive (or CAA?) and runs the function "c_4_v_timing_mva.m".

Usage
-----
**c\_4\_v\_timing\_mva** Performs timing and minimum variance analysis on four s/c
field data with a graphical user interface.
   **c\_4\_v\_timing\_mva**(b1,b2,b3,b4,R,column) interactive discontinuity
   analysis analyzer on magnetic field b1,...b4 with position R using
   column number 'column'. R has the form R.R1,...R.R4.
   
   **c\_4\_v\_timing\_mva**('B?',R) uses B1, B2, B3 and B4 from workspace.
   
   **c\_4\_v\_timing\_mva**('B?',R, column) also uses column. 2=x, 3=y, 4=z.

   The user is asked to define an interval in which the discontinuity
   analysis should be made. Two new figures open when this is done.
   In the first window, the normal vectors obtained from MVA and timing is
   shown in 3D space. Normal vectors that do not fulfill l2/l3>5 is shown
   as dotted arrows. The dotted arrows are not used in determining the
   maximum angle.
   In the second figure, the magnetic field is plotted in the
   LMN-system, where L=maximum, M=intermiediate and N = minimum.

   A prompt asks the user to save the variables: velocity of discontinuity V,
   uncertainty in velocity dV, all 5 normal vectors in n, where
   n.nTiming is the normal vector from timing and n.n1,...n4 is normal
   vectors obtained from MVA, all vectors from MVA v, and eigenvalues from
   MVA. The variables do not appear until the main figure is closed.
