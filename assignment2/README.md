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

* Cat (radius = 0.07 \*diagonal, resolution = 20, polydegree = 0)

![Cat](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/cat.png)

* Cat (radius = 0.07 \*diagonal, resolution = 60, polydegree = 0)

![Cat](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/catresolution.png)


* Luigi (radius = 0.07\*diagonal, resolution = 30, polydegree = 0)

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigi.png)


* Luigi (radius = 0.12\*diagonal, resolution = 30, polydegree = 1)

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigiradius.png)


* Sphere (radius = 0.045\*diagonal, resolution = 20, polydegree = 0)

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/sphere00.png)


* Sphere (radius = 0.045\*diagonal, resolution = 20, polydegree = 1)

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/sphere1.png)



### Comments

* Wendland radius: increasing the radius gives a smoother surface and generally cause loss of details. On the contrary, a too small radius does not allow to reconstruct the mesh, since otherwise no points on the grid would have neighboring constrained points and we could not fit a polynomial locally. We use a radius dependent on the dimension of the bounding box diagonal. The two previous pictures of Luigi show the impact of the radius on the reconstruction.

* Resolution: increasing the resolution improves the quality of the reconstruction, but makes the computations more expensive. The cat was reconstructed with a resolution of 20; for Luigi we used a resolution of 30, wich was necessary because the legs are near to each other (lower resolution caused the legs to be united). The two previous pictures of the cat show the impact of changing the resolution.

* Degree: with degree 0 we obtained better reconstructions; in fact, when increasing the degree, some extra matter around the true mesh appears. The best example to see the impact of the degree is the sphere (which in our case does not have the problem of extra matter around the true mesh). Keeping fixed the other parameters we can see that increasing polydegree from 0 to 1 has a smoothing effect. Note that such an effect can be obtained in general also increasing the wendland radius, which in some cases can be a better choice because of the extra matter artifacts caused by the increased degree (we give an example of some artifacts below, on luigi). Here on the sphere we attentively chose a wendland radius to demonstrate the effect of polydegree.

* Axis alignment: to save computation time, we load the mesh of Luigi in such a way that it is axis aligned. We implemented our simple axis alignment algorithm with PCA, namely taking the direction with highest variance in the point cloud and rotating the pointcloud (and corresponding normals!) in such a way that such direction is aligned with one of the axes (axis y for example).


After these comments we report another couple of nice reconstructions, and an example of the extra matter problem caused bu increasing polydegree:

* Horse

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/horse.png)

* Bunny

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/bunny.png)

* Luigi (polydeg1)

![Luigi](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/luigideg1.png)



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

Note that the back of the neck of the hound, which is a sharp feature (it is not really a sharp feature as we will point out in the next section, but it could be interpreted as such), is better preserved in the paper reconstruction. Their method based on normals is in fact more suitable in presence of such features. The other important advantage of the papaer implementation is that it is sensibly faster than the standard implementation. This is due to the fact that in the papaer implementation we only have n constraints, versus the 3n constraints of the standard implementation.

3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.

* Poisson

![Poisson](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/poissonhound.png )

* RIMLS

![RIMLS](https://github.com/ccasam/GP2020-Assignments/blob/master/assignment2/results/RIMLS.png)

Tha main difference we notice is that the "sharp feature" of the neck that was reconstructed by our MLS method is now not present anymore. Visually, is is clear that such a sharp feature is in fact a mistake and both poisson and RIMLS are able to detect the anomaly and correct for it, returning a smoother result on the neck. Such a sharp feature is in fact caused by a lack of data points on the neck but RIMLS is a method whose improvement is indeed robustness to outliers and situations analogous to this one. The screen we reported for RIMLS has resolution = 200 (we left all the Meshlab parameters unchanged), and even if this is much higher than the resolution we used for our implementation (usually 20 to 30), the computing time is much lower. Similarly, Poisson surface reconstruction is extremely fast, but gives less detailed reconstrunctions when compared to RIMLS. To summarize, paper implementation versus standard marching cubes behaves better on sharp features; poisson is super fat and RIMLS is super fast but also much more robust.
