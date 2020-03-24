# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

The image below reports the constrained points for the cat. Epsilon has been adjusted iteratively in such a way that normals do not intersect, meaning that the nearest point to every off surface constraint is the corresponding on surface constraint.

![Cat constraints](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/constraints.png)

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).

Red points represent positive values of the implicit function, green points negative values. The grid has been slightly expanded around the mesh in all dimensions, to avoid artifacts in the mesh reconstruction.

* Cat

![Cat grid](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/catgrid.png)

* Luigi

![Luigi grid](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigigrid.png)

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).

* Cat

![Cat](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/cat.png)


* Luigi

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigi.png)

### Comments

* Wendland radius: increasing the radius gives a smoother surface and generally cause loss of details. On the contrary, a too small radius does not allow to reconstruct the mesh, since otherwise no points on the grid would have neighboring constrained points and we could not fit a polynomial locally. We use a radius deoendent on the dimension of the bounding box diagonal (in the two meshes shown here we used r= 0.07\*diagonal).

* Resolution: increasing the resolution improves the quality of the reconstruction, but makes the computations more expensive. The cat was reconstructed with a resolution of 20; for Luigi we used a resolution of 30, wich was necessary because the legs are near to each other (lower resolution caused the legs to be united).

* Degree: with degree 0 we obtained better reconstructions; in fact, when increasing the degree, some extra matter around the true mesh appears.

* Axis alignment: to save computation time, we load the mesh of Luigi in such a way that it is axis aligned. We implemented our simple axis alignment algorithm with PCA, namely taking the direction with highest variance in the point cloud and rotating the pointcloud (and corresponding normals!) in such a way that such direction is aligned with one of the axes (axis y for example).


### Theory question: Save your notes to assignment2/results and add a link to this page.

1) Save your notes and add a link to this page.

* Proof that normal to surface is proportional to the gradient:
![normal](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/gradient_orthogonal.jpeg)

* Moving Least Squares explicit gradient computation:
![MLS gradient](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/mlsgrad.jpeg)

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

* Standard reconstruction

![normal](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/hound_standard.png)

* Paper reconstruction

![normal](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/hound_paper.png)

Note that the back of the neck of the hound, which is a sharp feature, is better preserved in the paper reconstruction. Their method based on normals is in fact more suitable in presence of such features. The other important advantage of the papaer implementation is that it is sensibly faster than the standard implementation. This is due to the fact that in the papaer implementation we only have n constraints, versus the 3n constraints of the standard implementation.

3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.

![Poisson](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/poissonhound.png )

![RIMLS](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/marchingcubesRIMLS.png )
