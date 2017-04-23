
https://en.wikipedia.org/wiki/Flexibility_method
https://en.wikipedia.org/wiki/Direct_stiffness_method

We have set of $N$ nodes with positions $r_i=(x_i,y_i,z_i)$ connectedby $M$ sticks of length $l_k = |r_k - r_k|$ and stiffness $s_k$. The structure is under stres due to external forces $g_i$ applied in nodes. We want to find deflections of node positions $\delta r_i = r_i - r_i^0$. 

### First order

Assuming that deflections are very small, angles and directions between stics does not change due to due to strain. Than we solve system of linear equations in matrix form.
$ A q_i = g_i $ 
with stiffness matrix $A$. The quest is to express $A$ in terms of stick stiffness $s_k$.  Force acting along the stick is $f_k = s_k \delta l_k$ where $\delta l_k = l_k - l_k^0$ is a small change of stick lenght from . Equlibrium is achieved when sum of all stick forces acting on node $i$ equals to external force $g_i$. 



Or by variation of stick energy with position of vectors $ E_k = s_k \delta l_k^2/2 = s_{ij} (|r_i-r_j| - l^0_{ij})^2/2$. 

Derivative of total energy $E = \sum E_k = \sum s_{ij} (|r_i-r_j| - l^0_{ij})^2/2 $

$\partial E_k / \partial x_i = s_{ij} (x_i-x_j) (|r_i-r_j| - l^0_{ij}) /|r_i-r_j| 
\approx s_{ij} \delta l_{ij} (x_i-x_j)/l_{ij}^0 
= s_{ij} \delta l_{ij} d^x_{ij} $

for total energy $E$ we obtain 
$  \partial E / \partial x_i = \sum_j s_{ij} \delta l_{ij} d^x_{ij}  $

but we have to find relation exressed explicitly in terms of $\delta x_i$ rather than in terms of $\delta l_{ij}$. We utilize relation beween $\delta l_k$ and $\delta r_i$ assuming direction encoded in unitary vecor $d = (d^x,d^y,d^z)$  is constant.

$\delta l_k 
= \delta (x_i-x_j) d^x_k 
= \delta (y_i-y_j) d^y_k 
= \delta (z_i-z_j) d^z_k$

Therefore 
$ g_i =   \partial E / \partial x_i = \sum_j s_{ij} \delta (x_i-x_j) (d^x_{ij})^2  $









