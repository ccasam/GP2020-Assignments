# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

The image below reports the constrained points for the cat. Epsilon has been adjusted iteratively in such a way that normals do not intersect, meaning that the nearest point to every off surface constraint is the corresponding on surface constraint.

![Cat constraints](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/constraints.png)

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).

Red points represent positive values of the implicit function, green points negative values. The grid has been slightly expanded around the mesh in all dimensions, to avoid artifacts in the mesh reconstruction.

Cat:
![Cat grid](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/catgrid.png)

Luigi:

![Luigi grid](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigigrid.png)

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).

Cat:
![Cat](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/cat.png)


Luigi:
![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigi.png)

### Comments

* Wendland radius
* Resolution
* Degree
* Axis alignment


### Theory question: Save your notes to assignment2/results and add a link to this page.

1) Save your notes and add a link to this page.

Proof that normal to surface is proportional to the gradient:
![normal](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/.png)

Moving Least Squares explicit gradient computation:
![MLS gradient](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/gradient.jpg)

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.
