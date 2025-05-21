# Perlin Noise


## Generalized Perlin Noise

### Noise Grid Dimensions

Perlin noise elegantly generalizes to $n$ noise dimensions.

1. 2-dimensional perlin noise can represent a height map for a 2-dimensional terrain mesh.
1. 3-dimensional perlin noise can represent a density map for a 3-dimensional terrain mesh. For example, for a 3-dimensional voxel cube, any noise value in the lower half of the noise range, i.e. noise in $[0, 1/2] \subset [0, 1]$, could be interpreted as an air voxel. This generates pockets of air which represent caves inside of the voxel cube.
1. 3-dimensional perlin noise can **represent the evolution of a 2-dimensional height map over time**. In other words, we can interpret the third dimension of the noise map as the time dimension, so that each step along the third axis corresponds to a different height map composed of the first two noise map dimensions. This produces smooth results because that third dimension is still generated via perlin noise, which makes the time steps smooth via the interpolation of the noise along that third, time dimension.

### Noise Grid Dimension Step Sizes

Perlin noise elegantly generalizes to arbitrary step sizes in each dimension. The gradients of the grid must not necessarily be located at integer points $(x, y) \in \mathbb{N}^2$. Given $s_1, \dots, s_n \in \mathbb{R}^+$ where $0 \lt s_k$ is the step size between gridcells in the $k$-th dimension, and $m_1, \dots, m_n \in \mathbb{N}$ where $0 \lt m_k$ is the number of grid cells in the $k$-th dimension. We define the perlin noise gridpoints $G$ as the set of n-dimensional cartesian coordinates that each correspond to one of the gradient vectors.

$
    G = \{ (c_1 \cdot s_1, \dots, c_n \cdot s_n) \in \mathbb{R}^2 \mid k = 1, \dots, n : c_k \in \{ 0, 1, \dots, m_k \} \}
$

<p style="color:red">Read the bold note below; is it preferable to always use all step sizes 1? This simplifies the computations AND the math? And you just need to scale the input space to match the all 1 step size space?</p>

**BUT, it is fairly trivial to map an input space that uses step sizes differing from 1 to that space that does use all step sizes 1? I.e. given a coordinate $\vec{x} = (x_1, x_2, x_3)$ with step sizes $\vec{s} = (s_1, s_2, s_3)$ you transform the sample point $\vec{x}$ to the new point $\vec{y}$ in the space that uses step sizes (1, 1, 1) as follows: $\vec{y} = \left( \dfrac{x_1}{s_1}, \dfrac{x_2}{s_2}, \dfrac{x_3}{s_3} \right)$. Essentially, you independently stretch out (if $s_i < 1$) or shrink (if $s_i > 1$) each dimension of the old space so that one unit on dimension $i$ in the new space equals $1$ instead of $s_i$ as it did in the old space. Thus, we should prefer a perlin noise implementation that exclusively uses grid step sizes $s_i = 1$ because this really simplifies the implementation AND the math. We can simulate step sizes different from $1$ by simply scaling the sample points.**

In other words, each gridpoint is an n-dimensional vector where the $k$-th vector component is a multiple of the $k$-th step size.

#### Example: 2D grid with variable step size

The layout of a 2-dimensional grid is as follows, where $0 \leq i \leq m_1$ and $0 \leq j \leq m_2$.

$G = \begin{Bmatrix}
    (0,0) & \cdots & (i \cdot s_1, 0) & \cdots & (m_1 \cdot s_1, 0)\\
    \vdots & \ddots & \vdots & & \vdots\\
    (0, j \cdot s_2) & \cdots & (i \cdot s_1, j \cdot s_2) & \cdots & (m_1 \cdot s_1, j \cdot s_2) \\
    \vdots & & \vdots & \ddots & \vdots \\
    (0, m_2 \cdot s_2) & \cdots & (i \cdot s_1, m_2 \cdot s_2) & \cdots & (m_1 \cdot s_1, m_2 \cdot s_2)\\
\end{Bmatrix}$

Some perlin noise implementations make the assumption that all step sizes are exactly $1$, meaning for all $k = 1, \dots, n$ it holds that $s_k = 1$. This implies the gridpoints only have integer coordinates. For example, a 2-dimensional grid then consists of gridpoints $(x, y) \in \mathbb{N}^2$.

$G = \begin{Bmatrix}
    (0,0) & \cdots & (i, 0) & \cdots & (m_1, 0)\\
    \vdots & \ddots & \vdots & & \vdots\\
    (0, j) & \cdots & (i, j) & \cdots & (m_1, j) \\
    \vdots & & \vdots & \ddots & \vdots \\
    (0, m_2) & \cdots & (i, m_2) & \cdots & (m_1, m_2)\\
\end{Bmatrix}$

Note, we can uniquely identify any gridpoint $(i \cdot s_1, j \cdot s_2)$ by its 2-dimensional tuple of indices $(i, j)$ that specify its location within the grid.


### Perlin Noise Sequence

The sample point $\vec{x} = (x_1, \dots, x_n)$ is contained in an n-dimensional hypercube, a gridcell part of the noise map. An n-dimensional hypercube by definition has $2^n$ corners, and we associate one gradient vector with each corner. Thus, we must interpolate $2^n$ noise samples to obtain the final, singular perlin noise value for the sample point $\vec{x}$.

Implementation wise, suppose that each gridcell (a hypercube, not a *gridpoint*) is uniquely identified by an n-tuple of indices $(i_1, \dots, i_n)$ into the n-dimensional grid. That is, a 2-dimensional gridcell with $2^2 = 4$ corners is identified by its "bottom left" corner's indices $(i_1, i_2)$ in the grid. The same holds for a 3-dimensional gridcell with $2^3 = 8$ corners, which is identified by its "bottom left" corner's indices $(i_1, i_2, i_3)$. Then the sequence that specifies the *gridpoint* indices of all corners of the hypercube that contains $\vec{x}$ can be ordered in ascending order by those indices, where the importance of the indices to the ordering is $i_n > \dots > i_1$. This ordering will simplify evaluating the interpolation of all gradients that correspond to the chosen corner gridpoints. Given an n-dimensional gridcell (hypercube) $S_{i_1, \dots, i_n}$ uniquely identified by its "bottom left" corner indices $(i_1, \dots, i_n)$, we describe the set containing all of its corners' indices as follows.

$$
    S_{i_1, \dots, i_n} = \{ i_1, i_1+1 \} \times \dots \times \{ i_n, i_n+1 \}
$$


For example, assume $n=3$ dimensions and $\vec{x} = (x_1, x_2, x_3)$ is contained in the gridcell (hypercube) uniquely described by indices $(i_1, i_2, i_3) = (2, 1, 3)$ corresponding to the hypercube's "bottom left" corner. Then the sequence $S$ of $2^n = 2^3 = 8$ gridpoints is ordered by indices according to the ordering $i_3 > i_2 > i_1$.


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

### Perlin Noise Sequence Interpolation

Suppose $p_{j_1, j_2, j_3}$ is the noise value that corresponds to the hypercube corner (gridpoint) uniquely identified by the indices $(j_1, j_2, j_3)$.
The usefulness of the sequence ordering discussed previously becomes apparent when we describe how to compose the  linear interpolations in the figure below, to reduce $2^n = 2^3 = 8$ noise values to a single noise value $p$.
- Computing $p_{j_2, j_3}$ requires linearly interpolating noise two values with different $j_1$ index but the same $j_2$ and $j_3$ indices, implying this interpolation happens with respect to the **x-axis**.
- Next, computing $p_{j_3}$ requires linearly interpolating noise values with different $j_2$ index but the same $j_3$ index, implying this interpolation happens with respect to the **y-axis**.
- Finally, computing computing $p$ requires linearly interpolating noise values with different $j_3$ index, implying this interpolation happens with respect to the **z-axis**.

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

We can determine $t_k$ by affine transformation given sample point $\vec{x} = (x_1, \dots, x_,)^T$ and knowing that the $k$-th dimension of the gridcell (hypercube) uniquely identified by the indices $(i_1, \dots, i_n)$ ranges $x_k \in D_k =  [i_k \cdot s_k, (i_k + 1) \cdot s_k]$ where $s_k$ is the $k$-th dimension's step size and $0 \leq i_k \leq m_k - 1$ holds.

$$
    t_k = \dfrac{x_k - i_k \cdot s_k}{((i_k + 1) \cdot s_k) - (i_k \cdot s_k)} = \dfrac{x_k - i_k \cdot s_k}{s_k}
$$

In other words, this affine transformation applies a translation and scaling onto $x_k$ based on its domain, to map $x_k \in D_k$ onto $t_k \in [0, 1]$. This can be verified via interval arithmetic.

$$
    t_k = \dfrac{x_k - i_k \cdot s_k}{s_k} \quad \Rightarrow \quad t_k \in \dfrac{[i_k \cdot s_k, (i_k + 1) \cdot s_k] - i_k \cdot s_k}{s_k} = \dfrac{[0, s_k]}{s_k} = [0, 1]
$$

This choice of interpolation fractions $(t_1, t_2, t_3)$ ensures that the closer the sample point is in the k-th dimension to a specific corner, the more the lerped noise value resembles the noise value of that corner.

We illustrate the linear interpolation with the gif animations below. Note that the x-axis is red, the y-axis is green and the z-axis is blue. Note that the hypercube (gridcell) of interest here is identified by indices $(i_1, i_2, i_3) = (0, 0, 0)$. **The main takeaway is that linearly interpolating between two corners (n-dimensional points in space) is analogous to linearly interpolating between the noise values that those corners correspond to.**

- **Fig 1.** and **Fig 2.** show that the sample point $\vec{x} = P$ travels along a spiral path through the chosen gridcell, so that we showcase the linear interpolation values for a variety of coordinates.
- **Fig 3.** shows that the x-coordinate $x_1$ of the sample point $\vec{x} = P = (x_1, x_2, x_3)^T$ directly impacts the value of x-axis interpolation fraction $t_1$. Namely, if the hypercube ranges $x_1 \in D_1 = [0, 4]$ on the x-axis and $x_1 = 3.0$, then $t_1 = \dfrac{x_k - i_k \cdot s_k}{s_k} = \dfrac{3.0 - 0 \cdot 4.0}{4.0} = 0.75$. Said differently, if the sample point $\vec{x}$ is $75\%$ of the way through the hypercube on the x-axis, then x-axis lerping will use a factor $t_1=0.75$. This is reflected by the green points in the animation. If the green point $p_{2,2} = C_{1,2,2} + t_{1} (C_{2,2,2} - C_{1,2,2})$ is $75\%$ of the way to corner $C_{2,2,2}$ from corner $C_{1,2,2}$, meaning it is lerped $75\%$ of the way between corners, then we compute the lerped noise value exactly as $p_{2,2} = p_{1,2,2} + 0.75 (p_{2,2,2} - p_{1,2,2})$ where $p_{1,2,2}$ is the noise value of corner $C_{1,2,2}$ and $p_{2,2,2}$ is the noise value of corner $C_{1,2,2}$. The green points correspond to the first round of linear interpolation.
- **Fig 4.** shows the same, but for the y-axis. The yellow points are linear interpolations between two green points each, but their y-coordinate depends on that of the sample point $\vec{x}$. This reflects that the sample y-ccordinate $x_2$ directly determines $t_2$.
- **Fig 5.** shows the same, but for the z-axis. Note that the red point $p$, which corresponds to the single, final perlin noise value $p$, overlaps exactly with the sample point $\vec{x}$. This is significant, because it validates that the interpolation fractions $(t_1, t_2, t_3)$ were chosen correctly, and that the corner sequence ordering discussed previously indeed evaluates to the lerped noise value of the sample point $\vec{x}$.

Fig 1. The path followed by point $P$ to showcase $lerp$.|Fig 2. Project point $P$ on the sides of the hypercube.
-|-
![Lerp Cube showing the path point P takes.](/artefacts/cube-0-cut.gif) | ![Lerp Cube showing the projection of point P onto the sides of the cube.](/artefacts/cube-1-cut.gif)


Fig 3. The x-axis $lerp$ produces green points by lerping corners.|Fig 4. The y-axis $lerp$ produces yellow points by lerping green points.|Fig 5. The z-axis $lerp$ produces a red point by lerping yellow points.
-|-|-
![Lerp Cube showing the x-axis lerping.](/artefacts/cube-2-cut.gif) | ![Lerp Cube showing the y-axis lerping.](/artefacts/cube-3-cut.gif) | ![Lerp Cube showing the z-axis lerping.](/artefacts/cube-4-cut.gif)


## Unit Circle

There are multiple ways to express the unit circle.

### Cartesian Definition

First, its classical, cartesian definition is $x^2 + y^2 = 1$ with respect to cartesian coordinates $x$ and $y$. That is, the unit circle is a set of cartesian points $(x, y)$ at $L_2$ (cartesian) distance $1$ from the origin.

$
    C = \{ (x, y) \in \mathbb{R}^2 \mid x^2 + y^2 = 1 \}
$

Note that the expression $x^2 + y^2 = 1$ really just specifies the cartesian distance requirement of the unit circle for a vector $\vec{v} = (x, y)^T$.

$
    x^2 + y^2 = 1\\
    \Rightarrow x^2 + y^2 = 1^2\\
    \Rightarrow \sqrt{x^2 + y^2} = 1\\
    \Rightarrow |\vec{v}|_2 = 1
$

Here $|\vec{v}|_2$ expresses the $L_2$ (cartesian) norm of the vector $\vec{v}$. Thus, we can rephrase the unit circle set definition.

$
    C = \{ \vec{v} \in \mathbb{R}^2 \mid |\vec{v}|_2 = 1 \}
$


### Trigonometric Definition

Second, we can equivalently express the unit circle using trigonometric functions $sin$ and $cos$. It is based on the goniometric identity $cos^2(\alpha) + sin^2(\alpha) = 1$. This identity is reminiscent of the cartesian unit circle definition $x^2 + y^2 = 1$. From this, we conclude $x = cos(\alpha)$ and $y = sin(\alpha)$ specifically when defining the unit circle with regards to the $L_2$ norm of $\vec{v} = (x, y)^T$.

$
    C = \{ (cos(\alpha), sin(\alpha)) \in \mathbb{R}^2 \mid \alpha \in [0, 2\pi] \}
$

Note that basing the unit circle definition on the goniometric identity makes it implicit that $|\vec{v}|_2 = 1$; the set definition does not explicitly state any requirements on the $L_2$ norm of $\vec{v} = (cos(\alpha), sin(\alpha))$ since the goniometric identity implicitly satisfies the norm requirement of the unit circle already.

### Vector Math

Multiplying a vector $\vec{v} = (v_1, \dots, v_n)^T$ by a scalar $c$ also multiplies its $L_2$ (cartesian) length $|\vec{v}|_2$ by that scalar.

$$
    c \cdot |\vec{v}|_2 = |c \cdot \vec{v}|_2
$$

Note that $c \cdot \vec{v} = c \cdot (v_1, \dots, v_n)^T = (c \cdot v_1, \dots, c \cdot v_n)^T$. Then we can algebraically verify the claim.

$
    c \cdot |\vec{v}|_2\\
    = c \cdot \sqrt{v_1^2 + \dots + v_n^2}\\
    = \sqrt{c^2 \cdot (v_1^2 + \dots + v_n^2)}\\
    = \sqrt{c^2 \cdot v_1^2 + \dots + c^2 \cdot v_n^2}\\
    = \sqrt{(c \cdot v_1)^2 + \dots + (c \cdot v_n)^2}\\
    = |(c \cdot v_1 + \dots + c \cdot v_n)|_2\\
    = |c \cdot \vec{v}|_2\\
$

### Dot Product

The dot product between two vectors $dot(\vec{u}, \vec{v}) = \vec{u} . \vec{v}$ has two equivalent definitions, given vectors $\vec{u} = (u_1, \dots, u_n)^T$ and $\vec{v} = (v_1, \dots, v_n)^T$.

$
    \vec{u} . \vec{v}\\
    = u_1 \cdot v_1 + \dots + u_n \cdot v_n\\
    = |\vec{u}|_2 \cdot |\vec{v}|_2 \cdot cos(\theta)
$

where $\theta$ is the angle between the vectors $\vec{u}$ and $\vec{v}$. Note that $cos(\theta) = cos(\theta + 2\pi) = cos(2\pi - \theta)$. Consequently the dot product $\vec{u} . \vec{v} = \vec{v} . \vec{u}$ is commutative because if $\vec{u}.\vec{v}$ measures the inner angle $\theta$ then $\vec{v}.\vec{u}$ measures the outer angle $2\pi - \theta$.




## Perlin Noise Range

### Rewrite Lerp using cos

Let us first determine the value range of the initially sampled noise values, computed by the dot product of gradient vector $\vec{g}$ and delta vector $\vec{d}$ where $\vec{g}, \vec{d} \in \mathbb{R}^n$, before any linear interpolations are applied.

$
    \vec{g}.\vec{d} = |\vec{g}|_2 \cdot |\vec{d}|_2 \cdot \cos(\vec{g},\vec{d})\\
    = |\vec{g}|_2 \cdot |\vec{d}|_2 \cdot \cos(\theta)\\
    = [1, 1] \cdot [0, \sqrt{n}] \cdot [-1, 1]\\
    = [-\sqrt{n}, \sqrt{n}]
$

where $|\vec{g}|_2$ denotes the length range of the gradient vector and $|\vec{d}|_2$ denotes the length range of the delta vector.
1. The gradient vector is always normalized, so its length range is $[1, 1]$ degenerate.
1. The delta / offset vector has all zero components at its shortest, and all one / negative one components at its longest. Namely, $|\vec{x}|_2 = \sqrt{x_1^2 +  \dots + x_n^2}$ so that $|\vec{x}|_2 = \sqrt{0^2 + \dots + 0^2} = 0$ is the shortest and $|\vec{x}|_2 = \sqrt{(\pm 1)^2 + \dots + (\pm 1)^2} = \sqrt{n}$ is the longest.
1. **BUT**, the previous point only holds for the classic perlin noise implementations where the noise hypercubes always have cube edge length $1$. So, if you want to use a "step size" different from $1$, then you instead get that $|\vec{x}|_2 = \sqrt{s_1^2 + \dots + s_n^2}$ at its longest.
1. $cos(A,B)$ is equivalent to $cos(\theta)$. Both mean the cosine of the angle between vectors $A$ and $B$. By definition, the range of the cosine is always $[-1, 1]$. Note, you can compute the angle $\theta$ between two $n$-dimensional vectors, the vectors' dimensionality can be $> 2$.


Given two adjacent corners $\vec{a}, \vec{b} \in \mathbb{R}^n$. We call $\vec{g}_a, \vec{g}_b \in \mathbb{R}^n$ the gradient vectors, $\vec{x} \in \mathbb{R}^n$ the sample point, $\vec{d}_a = \vec{x} - \vec{a}$ and $\vec{d}_b = \vec{x} - \vec{b}$ the delta vectors between the corners and the sample point, $\theta_a$ the angle between $\vec{g}_a$ and $\vec{d}_a$, and $\theta_b$ the angle between $\vec{g}_b$ and $\vec{d}_b$. Furthermore, require that $x_i \in [a_i, b_i]$ where $a_i \le b_i$ for all $i = 1, \dots, n$.

The unlerped noise values $p_a, p_b \in \mathbb{R}$ are computed as follows. They are **not** vectors, but real numbers!

$
    p_a = \vec{g}_a . \vec{d}_a = |\vec{g}_a|_2 \cdot |\vec{d}_a|_2 \cdot \cos(\theta_a)   \quad \text{and} \quad   p_b = \vec{g}_b . \vec{d}_b = |\vec{g}_b|_2 \cdot |\vec{d}_b|_2 \cdot \cos(\theta_b)
$

The interpolation fractions are computed as follows.

$
    t_1 = \dfrac{x_1 - a_1}{b_1 - a_1} \quad \dots \quad t_i = \dfrac{x_i - a_i}{b_i - a_i} \quad \dots \quad t_n = \dfrac{x_n - a_n}{b_n - a_n}\\
$

Assume the gradient vectors $g_a$ and $g_b$ are unit vectors, so that $|\vec{g}_a|_2 = |\vec{g}_b|_2 = 1$.

$$\begin{align*}
    lerp(t_1, p_a, p_b) &= p_a + t_1 \cdot (p_b - p_a)\\
    &= p_a \cdot (1 - t_1) + p_b \cdot t_1\\
    &= |\vec{g}_a|_2 \cdot |\vec{d}_a|_2 \cdot cos(\theta_a) \cdot (1 - t_1) + |\vec{g}_b|_2 \cdot |\vec{d}_b|_2 \cdot cos(\theta_b) \cdot t_1\\
    &= |\vec{d}_a|_2 \cdot cos(\theta_a) \cdot (1 - t_1) + |\vec{d}_b|_2 \cdot cos(\theta_b) \cdot t_1\\
\end{align*}$$


<p style="color:red">TODO: Can the fraction $t/d_a$ OR $t/d_b$ give any indication of the relation between t and d_a and d_b? I.e. can I plot this so that it is visually intuitive?</p>



### Notes

We can determine the range of values produced by the perlin noise algorithm by depending on the trigonometric defintion of the unit circle.


**Notes from here on, rephrase them!**

1.  $[-sqrt(n)/2, sqrt(n)/2]$ (https://stackoverflow.com/questions/18261982/output-range-of-perlin-noise)
1. Try to reformulate the dot product results in function of $sin$ and $cos$ and stuff? I.e. w.r.t. the trigonometric unit circle definition? How does the unit sphere (3D) affect the dot product results, since that is one additional term to the dot product's sum compared to the unit circle (2D)?

Given vectors $\vec{u}, \vec{v} \in \mathbb{R}^n$ where $n = 1, \dots$ . Denote the dot product of $\vec{u} = (u_1, \dots, u_n)^T$ and $\vec{v} = (v_1, \dots, v_n)^T$ as $\vec{u}.\vec{v} = u_1 \cdot v_1 + \dots + u_n \cdot v_n$

1. $\vec{u}.\vec{v}$ is largest when the angle $\alpha$ between $\vec{u}$ and $\vec{v}$ is $0\degree$ or $360\degree$.
1. $\vec{u}.\vec{v}$ is $0$ when the angle $\alpha$ between $\vec{u}$ and $\vec{v}$ is $90\degree$ or $-90\degree$.
1. $\vec{u}.\vec{v}$ is smallest when the angle $\alpha$ between $\vec{u}$ and $\vec{v}$ is $180\degree$ or $-180\degree$.

Note that the offset vector between a given pair of grid point (to which a gradient vector corresponds) and sample point is fixed! Only the gradient can realistically vary. Thus, the range of noise values is influenced by two factors:
1. When is the dot product smallest or largest? (Because all gradient vectors are of the same length (unit length?), so only the orientation of the gradients w.r.t. the sample point is important)
1. What is the largest possible offset vector? (Because the offset vectors can still vary in length!)


https://www.gamedev.net/forums/topic/285533-2d-perlin-noise-gradient-noise-range--/

Step 1:
Nothing interesting.

Step 2:
The smallest possible delta vector for each corner is of length 0. The largest possible delta vector is when all the components are +/- 1, giving a vector of length Sqrt(n), with n being the dimentionality of the noise. The gradient vectors are always of length 1, since they are normalized during initialization.

Step 3:
Since A.B = |A|*|B|*cos(A,B), this gives us values in the range -Sqrt(n)..Sqrt(n).

Step 4:
This one is actually trickier than it looks. To simplify the analysis, let's consider the 1D example first:
A---P------B

P is the sample point, A and B are the corners of the containing unit segment, t is the distance from A to P, and s(t) is the interpolation factor. Let's call the dot product results for A and B, a and b, respectively. The value at P is then:
p = b*s + a*(1-s)

This gives us the follwing behavior. What's important to remember is that a and b are also functions of t:

- At t=0.5, p is exactly the average of a and b, which is 0.5 for the 1D case.
- At t=0 and t=1, p is zero, because then you either have "p = b\*0 + 0 \* 1" or "p = 0\*1 + a\*0". <p style="color: red">p = 0 at t=0 and t=1 because this implies the sample point is overlaps exactly with the corner point! This means the delta / offset vector would be the zero vector, resulting in a dot product $|\vec{x}|_2 = 0$ which is zero! Thus, the interpolation would also result in a zero noise value?</p>
    - <p style="color: red">Given formula "p = b*s + a*(1-s)". The values a, b are the noise values that result from the dot products, for corners A, B respectively.</p>
    - <p style="color: red">t=0 means s=0 so that "p = b*0 + a*(1-0) = b*0 + a*1" is self evident. But, since we interpolate from a to b, then t=0 means sample x overlaps with the corner that produces a. Thus, the offset vector for corner A is the zero vector, meaning a=0. Thus, "p = b*0 + a*(1-0) = b*0 + a*1 = b*0 + 0*1 = 0".</p>
    - <p style="color: red">Similarly, t=1 means s=1 and the offset vector for corner B is the zero vector, meaning b=0 so that "p = b*1 + a*(1-1) = b*1 + a*0 = 0*1 + a*0 = 0".</p>
- As we approach from t=0.5 to t=0, a decreases linearly, b increases linearly, but s and 1-s do not change linearly. Depending on how s and s-1 interact with their counterparts, this may result in one of two cases: a) p drops monotonically towards zero; or b) p slightly rises to some maximum, and then drops back to zero monotonically.
- As we approach from t=0.5 to t=1, the same thing occurs as when approaching t=0.


### Perlin Noise Range

<p style="color:red">My own writeup.</p>

Suppose the sample point $\vec{x}$ is contained in the hypercube uniquely identified by indices $(i_1, \dots, i_n)$. Then the set of corners is

$$
    S_{i_1, \dots, i_n} = \{ i_1, i_1+1 \} \times \dots \times \{ i_n, i_n+1 \}
$$

Note that the perlin noise equation works recursively. <p style="color:red">This is not entirely true! nD noise uses nD coordinates, meaning sample vector x, corner coordinates c and the offset vectors (x-c) are all nD.</p>

1. For 1-dimensional perlin noise, you compute $2^1=2$ noise values via dot product and lerp them to produce $p$.
1. For 2-dimensional perlin noise, you sample $2^2=4$ noise values and lerp them pairwise to produce $2^1=2$ intermediate noise values. Then lerp those again to produce $p$. The intermediate noise values are computed exactly as in 1-dimensional perlin noise.
1. For 3-dimensional perlin noise, you sample $2^3=8$ noise values and lerp them pairwise to produce $2^2=4$ intermediate noise values. Then lerp those again to produce $2^1=2$ intermdiate noise values. Then lerp those to produce $p$. The first set of intermediate noise values are computed exactly as in 1-dimensional perlin noise. The second set of intermediate noise values are computed exactly as in 2-dimensional perlin noise.
1. For n-dimensional perlin noise, this works analogously to 3-dimensional perlin noise.


Now, we try to write out the equations for each successive noise value.

The initial noise values always result from the linear interpolation of dot products.
$$\begin{align*}
    p_{i_1, j_2 \dots, j_n} &= \vec{g}_{i_1, j_2 \dots, j_n} . (\vec{x} - \vec{c}_{i_1, j_2 \dots, j_n})\\
    p_{i_1+1, j_2 \dots, j_n} &= \vec{g}_{i_1+1, j_2 \dots, j_n} . (\vec{x} - \vec{c}_{i_1+1, j_2 \dots, j_n})\\
    p_{j_2 \dots, j_n} &= lerp(t_1, p_{i_1, j_2 \dots, j_n}, p_{i_1+1, j_2 \dots, j_n})
\end{align*}$$

Here we computed the first linear interpolation, which for 1-dimensional perlin noise would be the final, output value. Note that we interpolate from corner $(i_1, i_2, \dots, i_n)$ to corner $(i_1+1, i_2, \dots, i_n)$ which **only** differ in the very first (x-axis) component. In other words, we will lerp across the x-axis only. Here $\vec{x}$ is the sample point, $\vec{g}_{j_1, j_2 \dots, j_n}$ is a gradient vector and $\vec{c}_{j_1, j_2 \dots, j_n} = (c_{j_1, j_2 \dots, j_n}^1, \dots, c_{j_1, j_2 \dots, j_n}^n)^T$ is the absolute coordinate vector of the corner that corresponds to that gradient vector. Thus $(\vec{x} - \vec{c}_{j_1, j_2 \dots, j_n})$ is the delta vector from the corner to the sample point.

Recall $x_k \in D_k =  [i_k \cdot s_k, (i_k + 1) \cdot s_k] = [c_{i_1, j_2 \dots, j_n}^k, c_{i_1+1, j_2 \dots, j_n}^k]$.

$$\begin{align*}
    lerp(t, a, b) &= a + t \cdot (b - a)\\
    &= a + b \cdot t - a \cdot t\\
    &= a \cdot (1 - t) + b \cdot t\\

    t_1 &= \dfrac{x_1 - i_1 \cdot s_1}{s_1} = \dfrac{x_1 - c_{i_1, j_2 \dots, j_n}^1}{c_{i_1+1, j_2 \dots, j_n}^1 - c_{i_1, j_2 \dots, j_n}^1}\\

    p_{j_2 \dots, j_n} &= lerp(t_1, p_{i_1, j_2 \dots, j_n}, p_{i_1+1, j_2 \dots, j_n})\\

    &= p_{i_1, j_2 \dots, j_n} \cdot (1 - t_1) + p_{i_1+1, j_2 \dots, j_n} \cdot t_1\\
    &= \vec{g}_{i_1, j_2 \dots, j_n} . (\vec{x} - \vec{c}_{i_1, j_2 \dots, j_n}) \cdot (1 - t_1) + \vec{g}_{i_1+1, j_2 \dots, j_n} . (\vec{x} - \vec{c}_{i_1+1, j_2 \dots, j_n}) \cdot t_1\\
\end{align*}$$


#### 1D Lerp

Here $n=1$ so each gridcell has $2^n=2^1=2$ corners $\vec{c}_{i_1}$ and $\vec{c}_{i_1+1}$.

$$\begin{align*}
    \vec{c}_{i_1} &= (c_{i_1}^1)^T\\
    \vec{c}_{i_1+1} &= (c_{i_1+1}^1)^T\\
    t_1 &= \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1}\\
    p &= lerp(t_1, p_{i_1}, p_{i_1+1})\\
    &= p_{i_1} + t_1 \cdot (p_{i_1+1} - p_{i_1})\\

    \\

    &= (\vec{g}_{i_1} . (\vec{x} - \vec{c}_{i_1})) + \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1} \cdot ((\vec{g}_{i_1+1} . (\vec{x} - \vec{c}_{i_1+1})) -(\vec{g}_{i_1} . (\vec{x} - \vec{c}_{i_1})))\\
    &= g_{i_1}^1 \cdot (x_1 - c_{i_1}^1) + \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1} \cdot ((\vec{g}_{i_1+1} . (\vec{x} - \vec{c}_{i_1+1})))\\

    \\

    &= p_{i_1} \cdot (1 - t_1) + p_{i_1+1} \cdot t_1\\
    &= p_{i_1} \cdot (1 - \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1}) + p_{i_1+1} \cdot \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1}\\
    &= p_{i_1} \cdot \dfrac{(c_{i_1}^1 - c_{i_1+1}^1) - (x_1 - c_{i_1}^1)}{c_{i_1}^1 - c_{i_1+1}^1} + p_{i_1+1} \cdot \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1}\\
    &= \dfrac{p_{i_1} \cdot ((c_{i_1}^1 - c_{i_1+1}^1) - (x_1 - c_{i_1}^1)) + p_{i_1+1} \cdot (x_1 - c_{i_1}^1)}{c_{i_1}^1 - c_{i_1+1}^1}\\
    &= g_{i_1}^1 \cdot (x_1 - c_{i_1}^1) \cdot \dfrac{(c_{i_1}^1 - c_{i_1+1}^1) - (x_1 - c_{i_1}^1)}{c_{i_1}^1 - c_{i_1+1}^1} + g_{i_1+1}^1 \cdot (x_1 - c_{i_1+1}^1) \cdot \dfrac{x_1 - c_{i_1}^1}{c_{i_1}^1 - c_{i_1+1}^1}\\
\end{align*}$$







# References

* https://en.wikipedia.org/wiki/Perlin_noise
* https://en.wikipedia.org/wiki/Smoothstep
* https://stackoverflow.com/questions/18261982/output-range-of-perlin-noise
* https://www.cs.montana.edu/courses/spring2005/525/students/Thiesen1.pdf
* https://www.youtube.com/watch?v=9B89kwHvTN4