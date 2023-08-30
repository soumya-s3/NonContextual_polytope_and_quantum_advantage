# Construction of NC polytope and their violation in Quantum Theory
(i) For a prepare and measure scenario with both preparation and measurement equivalence obtain the extremal points of both the convex polytopes using [lcon2vert](https://in.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra) function.

(ii) Obtain the extremal points of Noncontextual polytope by multiplying the vertices of the two polytopes obtained in previous step. `Extremal_points.ipynb` performs this task. 

(iii) Convert the vertex representation to halfspace representation of the polytope. [Polymake](https://polymake.org/doku.php/user_guide/tutorials/apps_polytope) is a useful candidate for this task. The halfspace inequalities are the NC behaviours. Few such inequalities corresponding to different scenarios are listed in the folder `Facet inequalities`.

(iv) Find similar facet inequalities those can be converted into each other by exchanging labels, respecting the symmetry of the given scenario. Represent all similar facets by a single facet from that set. The folder `Similar facets identification` contains few such examples of different scenarios. 

(v) Check quantum behaviours using `SeeSaw` and `Q1_bound`.  
