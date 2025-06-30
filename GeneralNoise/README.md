# Perlin Noise

There already exist great tutorials on implementing the original (simple) Perlin noise algorithm, and its extension to improved Perlin noise and octaves [[1],[2],[3]]. Furthermore, digitalfreepen formally proves in [[4]] that the range of 2-dimensional Perlin noise is bounded in $[-\frac{\sqrt{2}}{2}, \frac{\sqrt{2}}{2}]$, which generalizes to a range bounded in $[-\frac{\sqrt{n}}{2}, \frac{\sqrt{n}}{2}]$ for n-dimensional (simple) Perlin noise.

Most tutorials rightly take the path of simple visuals, intuitive explanations and pseudocode to present the Perlin noise algorithm. Because there already exists plenty of quality reading material on the topic, all that follows is me having fun figuring out my own, unintelligible notation for expressing Perlin noise semi-formally. This is supplemented by visualizations and animations whenever fitting.

Prefer to depend on the sources cited in the [references](#references) to anything I write below.


## Uses of multi-dimensional noise

Perlin noise elegantly generalizes to $n$ noise dimensions. That is, it is possible to sample noise values from an $n$-dimensional cartesian space. We mention some possible uses of high-dimensional noise sampling.

1. 2-dimensional Perlin noise can represent a height map for a 2-dimensional terrain mesh.
1. 3-dimensional Perlin noise can represent a density map for a 3-dimensional terrain mesh. For example, for a 3-dimensional voxel cube, any noise value in the lower half of the noise range, i.e. noise in $[-1, 0] \subset [-1, 1]$, could be interpreted as an air voxel while all remaining noise values are interpreted as soil. This approach generates pockets of air which represent caves inside of the voxel cube.
1. 3-dimensional Perlin noise can represent the evolution of a 2-dimensional height map over time. In other words, if we sample a 3D cube of noise values, then each slice of the cube parallel to the xy-plane corresponds to a distinct height map. The z-axis then corresponds to the specific point in time for which the slice is a height map.
1. 4-dimensional Perlin noise can represent the evolution of a 3-dimensional density map over time. Similarly to the voxel cave example, the fourth noise dimension can model the changing of the caves over time.



## Variable step-sizes

Initially I was confused why Perlin noise is specifically defined with respect to an n-dimensional lattice (nD grid) of points where all sides of the hypercubes that compose it are length $1$. In other words, in 1D this means a number line with segments of length $1$, in 2D this means a grid of squares with sides of length $1$, in 3D this means cubes with sides of length $1$, etc. If my intuition is correct, then it is very much possible to define a unique step size $0 \lt s_k \ne 1$ for each $k$-th dimension of the lattice and then generate Perlin noise based on the step sizes $s_1, \dots, s_n$. However, this can be reduced to the case where all step sizes are one, $s_1 = \dots = s_n = 1$, which corresponds to the usual definition of the lattice for Perlin noise.


### Defining the lattice

Perlin noise elegantly generalizes to arbitrary step sizes in each dimension. The gradients of the grid must not necessarily be located at integer coordinates $(x, y) \in \mathbb{N}^2$. Given $s_1, \dots, s_n \in \mathbb{R}$ where $0 \lt s_k$ is the step size between grid cells in the $k$-th dimension, and $m_1, \dots, m_n \in \mathbb{N}$ where $0 \lt m_k$ is the number of grid *cells* in the $k$-th dimension (meaning there are $m_k+1$ grid *points* in the $k$-th dimension). We define the grid points $G_p$ as the set of n-dimensional cartesian coordinates that each correspond to one of the gradient vectors.
$$
    G_p = \{ (c_1 \cdot s_1, \dots, c_n \cdot s_n) \mid k = 1, \dots, n : c_k \in \{ 0, 1, \dots, m_k \} \}
$$
In other words, each gridpoint is an n-dimensional vector where the $k$-th vector component $c_k \cdot s_k$ is a multiple of the $k$-th step size $s_k$. But, it is fairly trivial to map a cartesian space that uses step sizes $s_k \ne 1$ for $k = 1, \dots, n$ to a space that does use all step sizes $s_k = 1$ for $k = 1, \dots, n$. We map $G_p$ to $G'_p$.
$$\begin{align*}
    G'_p &= \{ \vec{g} \cdot \left(\frac{1}{s_1}, \dots, \frac{1}{s_n}\right) \mid \vec{g} \in G_p \}\\
    &= \{ (c_1, \dots, c_n) \mid k = 1, \dots, n : c_k \in \{ 0, 1, \dots, m_k \} \}
\end{align*}$$
Then, any point $\vec{x} = (x_1, \dots, x_n)$ for which we wish to compute a noise value must first undergo the same transformation, $\vec{y} = \left( \dfrac{x_1}{s_1}, \dots, \dfrac{x_n}{s_n} \right)$. Essentially, you independently stretch out (if $s_k < 1$) or shrink (if $s_k > 1$) each dimension of the old cartesian space, to obtain the modified grid $G'_p$ and modified sample points $\vec{y}$. Thus, we should prefer a Perlin noise implementation that exclusively uses grid step sizes $s_k = 1$ because this simplifies the implementation and the math. We can simulate step sizes different from $1$ by simply scaling the sample points $\vec{x}$ appropriately.


### Example: 2D grid with variable step size

The layout of a 2-dimensional grid is as follows, where $0 \leq i \leq m_1$ and $0 \leq j \leq m_2$.

$G_p : \begin{Bmatrix}
    (0, m_2 \cdot s_2) & \cdots & (i \cdot s_1, m_2 \cdot s_2) & \cdots & (m_1 \cdot s_1, m_2 \cdot s_2)\\
    \vdots & \ddots & \vdots & & \vdots \\
    (0, j \cdot s_2) & \cdots & (i \cdot s_1, j \cdot s_2) & \cdots & (m_1 \cdot s_1, j \cdot s_2) \\
    \vdots & & \vdots & \ddots & \vdots\\
    (0,0) & \cdots & (i \cdot s_1, 0) & \cdots & (m_1 \cdot s_1, 0)\\
\end{Bmatrix}$

However, it is easier implementation-wise to assume that for all $k = 1, \dots, n$ it holds that $s_k = 1$. This implies the grid points only have integer coordinates. For example, a 2-dimensional grid then consists of grid points $(x, y) \in \mathbb{N}^2$.

$G_p : \begin{Bmatrix}
    (0, m_2) & \cdots & (i, m_2) & \cdots & (m_1, m_2)\\
    \vdots & \ddots & \vdots & & \vdots \\
    (0, j) & \cdots & (i, j) & \cdots & (m_1, j) \\
    \vdots & & \vdots & \ddots & \vdots\\
    (0,0) & \cdots & (i, 0) & \cdots & (m_1, 0)\\
\end{Bmatrix}$

We can uniquely identify any gridpoint $(i \cdot s_1, j \cdot s_2)$ by its 2-dimensional tuple of indices $(i, j)$ that specify its location within the grid.

To make the transformation more explicit, consider a 2-dimensional $m_1 \times m_2 = 3 \times 2$ grid with step sizes $\vec{s} = (s_1, s_2) = (2, 2)$ and grid $G_p$. Both plots below denote the four corner grid points.
$$
    A = (0, 0), B = (0, m_2 \cdot s_2), C = (m_1 \cdot s_1, 0), D = (m_1 \cdot s_1, m_2 \cdot s_2)
$$
For grid $G_p$ these corner points are as follows.
$$
    A = (0, 0), B = (0, 4), C = (6, 0), D = (4, 6)
$$
To obtain $G'_p$ from $G_p$ we simply scale all elements in $G_p$ by $\vec{r} = (\frac{1}{s_1}, \frac{1}{s_2}) = (0.5, 0.5)$ so that $\vec{s}' = \vec{r} \cdot \vec{s} = (1, 1)$. Note that any sample point $\vec{x}$ must undergo the same scaling, $\vec{y} = \vec{r} \cdot \vec{x} = (0.5, 0.5) \cdot (3, 3) = (1.5, 1.5)$. The scaled corner grid points $G'_p$ are as follows,
$$
    A = (0, 0), B = (0, 2), C = (3, 0), D = (2, 3)
$$
which indeed results in a grid where the sides of the hypercubes are always length one.

![](/artefacts/perlin/grid-a.png)


### Ordering grid cell corners

We define the Perlin noise grid *cells* $G_c$ as the set of $n$-dimensional cartesian coordinates where each of its elements uniquely identifies one of the hypercubes in the lattice. As before, we assume step sizes $s_k = 1$ for $k = 1, \dots, n$.

$$\begin{align*}
    G_c &= \{ (c_1 \cdot s_1, \dots, c_n \cdot s_n) \mid k = 1, \dots, n : c_k \in \{ 0, 1, \dots, m_k - 1 \} \}\\
    \Rightarrow G_c &= \{ (c_1, \dots, c_n) \mid k = 1, \dots, n : c_k \in \{ 0, 1, \dots, m_k - 1 \} \}
\end{align*}$$
The sample point $\vec{x} = (x_1, \dots, x_n)$ is contained in an n-dimensional hypercube, a gridcell. An n-dimensional hypercube by definition has $2^n$ corners, and we associate one gradient vector with each corner. Thus, we must interpolate $2^n$ noise samples to obtain the final, singular Perlin noise value for the sample point $\vec{x}$.

Implementation wise, suppose that each gridcell in $G_c$ (a hypercube, not a *grid point*) is uniquely identified by an n-tuple of indices $(i_1, \dots, i_n)$ into the n-dimensional grid. That is, a 2-dimensional gridcell with $2^2 = 4$ corners is identified by its "bottom left" corner's indices $(i_1, i_2)$ in the grid. The same holds for a 3-dimensional gridcell with $2^3 = 8$ corners, which is identified by its "bottom left" corner's indices $(i_1, i_2, i_3)$. Then the sequence that specifies the *grid point* indices of all corners of the hypercube that contains $\vec{x}$ can be ordered in ascending order by those indices, where the importance of the indices to the ordering is $i_n > \dots > i_1$. This ordering will simplify evaluating the interpolation of all gradients that correspond to the chosen corner grid points. Given an n-dimensional gridcell (hypercube) $S_{i_1, \dots, i_n}$ uniquely identified by its "bottom left" corner indices $(i_1, \dots, i_n)$, we describe the set containing all of its corners' indices as follows.

$$
    S_{i_1, \dots, i_n} = \{ i_1, i_1+1 \} \times \dots \times \{ i_n, i_n+1 \}
$$


For example, assume $n=3$ dimensions and $\vec{x} = (x_1, x_2, x_3)$ is contained in the gridcell (hypercube) uniquely described by indices $(i_1, i_2, i_3) = (2, 1, 3)$ corresponding to the hypercube's "bottom left" corner. Then the sequence $S$ of $2^n = 2^3 = 8$ grid points is ordered by indices according to the ordering $i_3 > i_2 > i_1$.


$$
    S_{i_1, i_2, i_3} = S_{2, 1, 3} = \{ 2, 3 \} \times \{ 1, 2 \} \times \{ 3, 4 \}
$$

This gives rise to the following, ordered sequence of gridpoint indices. Note that the sequence is sorted top-to-bottom, where we first sort by index $i_3$, then by index $i_2$ and finally by index $i_1$.

```
(2, 1, 3)   # The start of the sorted sequence.
(3, 1, 3)
(2, 2, 3)
(3, 2, 3)
(2, 1, 4)
(3, 1, 4)
(2, 2, 4)
(3, 2, 4)   # The end of the sorted sequence.
```

Note, when testing with Julia version `1.11.4`, this is exactly the ordering of indices produced by the following Julia expression. The caveat is that the elements of the lists passed to `Iterators.product` are sorted already.

```julia-repl
julia> indices = Iterators.product([2,3],[1,2],[3,4])
julia> reduce(vcat, indices)
8-element Vector{Tuple{Int64, Int64, Int64}}:
 (2, 1, 3)
 (3, 1, 3)
 (2, 2, 3)
 (3, 2, 3)
 (2, 1, 4)
 (3, 1, 4)
 (2, 2, 4)
 (3, 2, 4)
```

### Ordering noise interpolations
Suppose $p_{i_1, i_2, i_3}$ is the noise value that corresponds to the hypercube corner (gridpoint) uniquely identified by the indices $(i_1, i_2, i_3)$.
The usefulness of the sequence ordering discussed previously becomes apparent when we describe how to compose the  linear interpolations in the figure below, to reduce $2^n = 2^3 = 8$ noise values to a single noise value $p$.
- Computing $p_{i_2, i_3}$ requires linearly interpolating noise two values with different $i_1$ index but the same $i_2$ and $i_3$ indices, implying this interpolation happens with respect to the **x-axis** because $i_1$ corresponds to the x-axis.
- Next, computing $p_{i_3}$ requires linearly interpolating noise values with different $i_2$ index but the same $i_3$ index, implying this interpolation happens with respect to the **y-axis**  because $i_2$ corresponds to the y-axis.
- Finally, computing computing $p$ requires linearly interpolating noise values with different $i_3$ index, implying this interpolation happens with respect to the **z-axis** because $i_3$ corresponds to the z-axis.

$$
    \left.\begin{aligned}
        \left.\begin{aligned}
            \left.\begin{aligned}
                p_{2, 1, 3}\\
                p_{3, 1, 3}\\
            \end{aligned}\right\} = lerp(t_1, p_{2, 1, 3}, p_{3, 1, 3}) = p_{1,3}\\
            \left.\begin{aligned}
                p_{2, 2, 3}\\
                p_{3, 2, 3}\\
            \end{aligned}\right\} = lerp(t_1, p_{2, 2, 3}, p_{3, 2, 3}) = p_{2,3}\\
        \end{aligned}\right\} = lerp(t_2, p_{1,3}, p_{2,3}) = p_3\\
        \left.\begin{aligned}
            \left.\begin{aligned}
                p_{2, 1, 4}\\
                p_{3, 1, 4}\\
            \end{aligned}\right\} = lerp(t_1, p_{2, 1, 4}, p_{3, 1, 4}) = p_{1,4}\\
            \left.\begin{aligned}
                p_{2, 2, 4}\\
                p_{3, 2, 4}\\
            \end{aligned}\right\} = lerp(t_1, p_{2, 2, 4}, p_{3, 2, 4}) = p_{2,4}
        \end{aligned}\right\} = lerp(t_2, p_{1,4}, p_{2,4}) = p_4
    \end{aligned}\right\} = lerp(t_3, p_3, p_4) = p
$$

Note that the linear interpolation also uses the three fractions $t_1, t_2, t_3 \in [0, 1]$. These fractions correspond to the variable $t$ in the linear interpolation equation 

$$
    lerp(t, a, b) = a + t \cdot (b - a)
$$

where $t \in [0, 1]$ linearly interpolates between $a$ and $b$. Instead of linear interpolation, you could use a smoothing function to make the interpolation more gradual. An example of a family of smoothing functions are sigmoid-like functions. 

$$
    slerp(t, a, b) = a + smooth(t) \cdot (b - a)
$$

We can determine $t_k$ by affine transformation given sample point $\vec{x} = (x_1, \dots, x_n)$ and knowing that the $k$-th dimension of the gridcell (hypercube) uniquely identified by the indices $(i_1, \dots, i_n)$ ranges $x_k \in D_k =  [i_k \cdot s_k, (i_k + 1) \cdot s_k]$ where $s_k$ is the $k$-th dimension's step size and $0 \leq i_k \leq m_k - 1$ holds.

$$
    t_k = \dfrac{x_k - i_k \cdot s_k}{((i_k + 1) \cdot s_k) - (i_k \cdot s_k)} = \dfrac{x_k - i_k \cdot s_k}{s_k}
$$

In other words, this affine transformation applies a translation and scaling onto $x_k$ based on its domain, to map $x_k \in D_k$ onto $t_k \in [0, 1]$. This can be verified via interval arithmetic.

$$
    t_k = \dfrac{x_k - i_k \cdot s_k}{s_k} \quad \Rightarrow \quad t_k \in \dfrac{[i_k \cdot s_k, (i_k + 1) \cdot s_k] - i_k \cdot s_k}{s_k} = \dfrac{[0, s_k]}{s_k} = [0, 1]
$$

This choice of interpolation fractions $(t_1, t_2, t_3)$ ensures that the closer the sample point is in the k-th dimension to a specific corner, the more the lerped noise value resembles the noise value of that corner.

We illustrate the linear interpolation with the gif animations below. Note that the x-axis is red, the y-axis is green and the z-axis is blue. The hypercube (gridcell) of interest here is identified by indices $(i_1, i_2, i_3) = (0, 0, 0)$. **The main takeaway is that linearly interpolating between two corners (n-dimensional points in space) is analogous to linearly interpolating between the noise values that those corners correspond to.**

- **Fig 1.** and **Fig 2.** show that the sample point $\vec{x} = P$ travels along a spiral path through the chosen gridcell, so that we showcase the linear interpolation values for a variety of coordinates.
- **Fig 3.** shows that the x-coordinate $x_1$ of the sample point $\vec{x} = P = (x_1, x_2, x_3)$ directly impacts the value of x-axis interpolation fraction $t_1$. Namely, if the hypercube ranges $x_1 \in D_1 = [0, 4]$ on the x-axis and $x_1 = 3.0$, then $t_1 = \dfrac{x_k - i_k \cdot s_k}{s_k} = \dfrac{3.0 - 0 \cdot 4.0}{4.0} = 0.75$. Said differently, if the sample point $\vec{x}$ is $75\%$ of the way through the hypercube on the x-axis, then x-axis lerping will use a factor $t_1=0.75$. This is reflected by the green points in the animation. If the green point $p_{2,2} = C_{1,2,2} + t_{1} (C_{2,2,2} - C_{1,2,2})$ is $75\%$ of the way to corner $C_{2,2,2}$ from corner $C_{1,2,2}$, meaning it is lerped $75\%$ of the way between corners, then we compute the lerped noise value exactly as $p_{2,2} = p_{1,2,2} + 0.75 (p_{2,2,2} - p_{1,2,2})$ where $p_{1,2,2}$ is the noise value of corner $C_{1,2,2}$ and $p_{2,2,2}$ is the noise value of corner $C_{1,2,2}$. The green points correspond to the first round of linear interpolation.
- **Fig 4.** shows the same, but for the y-axis. The yellow points are linear interpolations between two green points each, but their y-coordinate depends on that of the sample point $\vec{x}$. This reflects that the sample y-coordinate $x_2$ directly determines $t_2$.
- **Fig 5.** shows the same, but for the z-axis. Note that the red point $p$, which corresponds to the single, final Perlin noise value $p$, overlaps exactly with the sample point $\vec{x}$. This is significant, because it validates that the interpolation fractions $(t_1, t_2, t_3)$ were chosen correctly, and that the corner sequence ordering discussed previously indeed evaluates to the lerped noise value of the sample point $\vec{x}$.

Fig 1. The path followed by point $P$ to showcase $lerp$.|Fig 2. Project point $P$ on the sides of the hypercube.
-|-
![Lerp Cube showing the path point P takes.](/artefacts/perlin/cube-0-cut.gif) | ![Lerp Cube showing the projection of point P onto the sides of the cube.](/artefacts/perlin/cube-1-cut.gif)


Fig 3. The x-axis $lerp$ produces green points by lerping corners.|Fig 4. The y-axis $lerp$ produces yellow points by lerping green points.|Fig 5. The z-axis $lerp$ produces a red point by lerping yellow points.
-|-|-
![Lerp Cube showing the x-axis lerping.](/artefacts/perlin/cube-2-cut.gif) | ![Lerp Cube showing the y-axis lerping.](/artefacts/perlin/cube-3-cut.gif) | ![Lerp Cube showing the z-axis lerping.](/artefacts/perlin/cube-4-cut.gif)




## Dot product

The dot product is a central operation with regards to Perlin noise. Thus, we discus some of its properties and its relation to the noise algorithm.

The dot product $\vec{u} . \vec{v}$ between two vectors has two equivalent definitions, given vectors $\vec{u} = (u_1, \dots, u_n)$ and $\vec{v} = (v_1, \dots, v_n)$.
$$\begin{align*}
    \vec{u} . \vec{v} &= u_1 \cdot v_1 + \dots + u_n \cdot v_n\\
    &= |\vec{u}|_2 \cdot |\vec{v}|_2 \cdot cos(\theta)\\
\end{align*}$$
where $\theta$ is the angle between the vectors $\vec{u}$ and $\vec{v}$. The dot product $\vec{u} . \vec{v} = \vec{v} . \vec{u}$ is commutative:
$$\begin{align*}
    \vec{u} . \vec{v} &= u_1 \cdot v_1 + \dots + u_n \cdot v_n\\
    &= v_1 \cdot u_1 + \dots + v_n \cdot u_n\\
    &= \vec{v} . \vec{u}
\end{align*}$$


### Gridcell dot products

With respect to Perlin noise, we always compute the dot product between an $n$-dimensional gradient vector and an $n$-dimensional sample point. Let us demonstrate the range that such dot products, and consequently each of the $2^n$ unlerped noise values, can take.
We give a 2-dimensional example, where we call $\vec{g}$ the gradient vector, we regard all four squares that share the gridcell corner at which $\vec{g}$ is located, and $\vec{x} = (x, y)$ is any sample point within those four squares (gridcells).

Given a square that ranges across $x,y \in [-1, 1]$ on the x-axis and y-axis, and a vector $\vec{g}$ starting at the origin, with $L_2$ norm $|\vec{g}|_2 = 1$. Then we can generate a heatmap of the dot product values $\vec{g} . \vec{x}$ for the unit vector $\vec{g}$ with every vector $\vec{x} = (x, y) \in [-1, 1] \times [-1, 1]$ in the square. The square region $[-1, 1] \times [-1, 1]$ corresponds to four distinct quadrants, each of which reprents a distinct grid cell.
- $x,y \in [-1, 0] \times [-1, 0]$ is the bottom left quadrant, and has $\vec{g}$ located in its top right corner.
- $x,y \in [0, 1] \times [-1, 0]$ is the bottom right quadrant, and has $\vec{g}$ located in its top left corner.
- $x,y \in [-1, 0] \times [0, 1]$ is the top left quadrant, and has $\vec{g}$ located in its bottom right corner.
- $x,y \in [0, 1] \times [0, 1]$ is the top right quadrant, and has $\vec{g}$ located in its bottom left corner.

We now determine the global minimum and maximum for the dot product restricted to the $x,y$ values in the entire square.

The longest vector $\vec{x}$ is found at the corners of the space, namely $\vec{x} = (\pm 1, \pm 1)$ so that $|\vec{x}|_2 = \sqrt{(\pm 1)^2 + (\pm 1)^2} = \sqrt{2}$. Additionally, $\cos(0) = 1$ and $\cos(\theta) \in [1, -1]$ is strictly, monotonically decreasing for $\theta \in [0, \pi]$ so that $\cos(\pi) = -1$. Mirroring this behavior, $\cos(\pi) = -1$ and $\cos(\theta) \in [-1, 1]$ is strictly, monotonically increasing for $\theta \in [\pi, 2\pi]$ so that $\cos(2 \pi) = 1$. The period of cosine is $2\pi$, so this behavior repeats for multiples of the period.

We have that the dot product given the angle $\theta$ between $\vec{g}$ and $\vec{x}$,
* $\vec{g} . \vec{x} = |\vec{g}|_2 \cdot |\vec{x}|_2 \cdot \cos(\theta) = 1 \cdot \sqrt{2} \cdot 1 = \sqrt{2}$, is maximized given $|\vec{g}|_2 = 1$, $|\vec{x}|_2 = \sqrt{2}$ and $\theta = 0 \degree$.
* $\vec{g} . \vec{x} = |\vec{g}|_2 \cdot |\vec{x}|_2 \cdot \cos(\theta) = 1 \cdot \sqrt{2} \cdot -1 = -\sqrt{2}$, is minimized given $|\vec{g}|_2 = 1$, $|\vec{x}|_2 = \sqrt{2}$ and $\theta = 180 \degree$.

![](/artefacts/perlin/dot-a.gif)

The white lines denote the contour plot overlayed on top of the heatmap. A contour plot consists of lines where each point on the same line shares the same value. Namely, in our plot each white line corresponds to a set of points $\vec{x} = (x, y)$ that share the exact same dot product value $c$, namely $\{ \vec{x} \mid \vec{g} . \vec{x} = c \}$. Note that all contour lines are perpendicular to the vector $\vec{g}$. This matches the geometric interpretation of the dot product. Namely, the dot product is the length of the projection of $\vec{x}$ onto $\vec{g}$. The projection of all points on such a perpendicular line w.r.t. $\vec{g}$ is the same vector $\vec{x}_p$ on $\vec{g}$.

![](/artefacts/perlin/dot-b.gif)
![](/artefacts/perlin/dot-c.gif)




## Perlin noise as a weighted sum

Interpreting Perlin noise as a chain of linear interpolations of a sequence of noise values is somewhat difficult to picture. Thus, it may be easier to write out the sequence of interpolations, and see that it resolves to a weighted sum of the noise values.

Suppose we consider the 1D gridcell represented by index $(i_1) = (1)$ and the 2D gridcell represented by indexes $(i_1, i_2) = (1, 1)$. Recall that for $n$ dimensions there are $2^n$ unlerped noise values. For 1 dimension, we have
$$\begin{align*}
    p = lerp(t_1, p_1, p_2) = p_1 (1 - t_1) + p_2 t_1
\end{align*}$$
and for 2 dimensions we have
$$\begin{align*}
    p_1 &= lerp(t_1, p_{1,1}, p_{2,1}) = p_{1,1} (1 - t_1) + p_{2,1} t_1\\
    p_2 &= lerp(t_1, p_{1,2}, p_{2,2}) = p_{1,2} (1 - t_1) + p_{2,2} t_1\\
    p &= lerp(t_2, p_1, p_2)\\
    &= p_1 (1 - t_2) + p_2 t_2\\
    &= (p_{1,1} (1 - t_1) + p_{2,1} t_1) (1 - t_2) + (p_{1,2} (1 - t_1) + p_{2,2} t_1) t_2\\
    &= p_{1,1} (1 - t_1) (1 - t_2) + p_{2,1} t_1 (1 - t_2) + p_{1,2} (1 - t_1) t_2 + p_{2,2} t_1 t_2\\
    &= p_{1,1} w_1 + p_{2,1} w_2 + p_{1,2} w_3 + p_{2,2} w_4
\end{align*}$$
where the weights $w_1, w_2, w_3, w_4$ sum to one.
$$\begin{align*}
    &w_1 + w_2 + w_3 + w_4\\
    &= (1 - t_1) (1 - t_2) + t_1 (1 - t_2) + (1 - t_1) t_2 + t_1 t_2\\
    &= 1 - t_2 - t_1 + t_1 t_2 + t_1 - t_1 t_2 + t_2 - t_1 t_2 + t_1 t_2\\
    &= 1
\end{align*}$$
We can generally describe the $n$-dimensional, simple Perlin noise function $P(\vec{x})$ for a sample point $\vec{x}$ as a weighted sum, where $p_i$ is some unlerped noise value and $w_i$ its corresponding weight.
$$
    P(\vec{x}) = \sum_{i=1}^{2^n} w_i p_i \quad\text{where}\quad w_i = \prod_{j=1}^n T_j \quad\text{with}\quad T_j \in \{t_j, (1-t_j)\}
$$

Note that the weighted sum $P(\vec{x})$ is non-linear! Each of its terms is a product of a weight and a noise value, both of which are dependent on the variable sample point $\vec{x} \in \mathbb{R}^n$. As shown in [[4]], the weighted sum is a dome-shaped function, where the extreme values occur at exactly $\vec{x} = (\frac{1}{2}, \dots, \frac{1}{2})$ the middle of the gridcell. The maximum occurs when all gradient vectors point exactly towards that center, and the minimum occurs when the gradient vectors point exactly away from that center.



# References


- [[1]] Eevee's Perlin noise tutorial
- [[2]] Raouf's Perlin noise tutorial
- [[3]] Adrian's Perlin noise tutorial
- [[4]] digitalfreepen's Perlin noise range proof
- [[5]] forum: gamedev.net Perlin noise range discussion
- [[6]] forum: SO Perlin noise range discussion
- [[7]] wikipedia: Perlin noise
- [[8]] wikipedia: smoothstep
- [[9]] youtube: Perlin noise visualization
- [[10]] Perlin noise & optimizations pseudocode

[1]: https://eev.ee/blog/2016/05/29/perlin-noise/
[2]: https://rtouti.github.io/graphics/perlin-noise-algorithm
[3]: https://adrianb.io/2014/08/09/perlinnoise.html
[4]: https://digitalfreepen.com/2017/06/20/range-perlin-noise.html
[5]: https://www.gamedev.net/forums/topic/285533-2d-perlin-noise-gradient-noise-range--/
[6]: https://stackoverflow.com/questions/18261982/output-range-of-perlin-noise
[7]: https://en.wikipedia.org/wiki/Perlin_noise
[8]: https://en.wikipedia.org/wiki/Smoothstep
[9]: https://www.youtube.com/watch?v=9B89kwHvTN4
[10]: https://www.cs.montana.edu/courses/spring2005/525/students/Thiesen1.pdf
