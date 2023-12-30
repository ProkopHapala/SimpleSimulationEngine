
## Truss Simulation Performance

### Elastic constrains: Potential vs Velocity reflection

There are two basic mechanisms how to solve constrains (e.g. rope, stick, or particle-particle collision):
 1. Calculate distance how much the constrain is violated (e.g. how much stick lengh $l$ differs from relaxed $l_0$, how much is the distance between the particle $r_{ij}$ closer than sum of their radius $R_{ij} = R_i +R_j$ ). Then calculate the force as $\vec f_{ij}$ as derivative of potential energy $V(r_{i,j})$. Typically we assume that close to the optimim it is qudratic (harmonic potential)  $V(r_{i,j}) \approx k(r_{ij} - R_{ij})^2$.
 2. Momentum approach where we assume elastic or inelastic of the two partilces, or reflection of the partilce from a hard (i.e. inpenetrable) constrain. In this approach we do not consider how much violated is the constrain (the potential), but rather we reflect the component of velocity $\vec v$ which is oriented along the constrained direction (e.g. normal to the inpenetrable surface, or along the stick or rope bonding the two points). The change of the velocity can be calculated from momentum preservation: $ m_i {\vec v}_i + m_j {\vec v}_j = m_i {\vec v'}_i + m_j {\vec v'}_j $ and form convervation of a part of kinetic energy $E_k = (m_i v_i^2 + m_j v_j^2)/2 $. We do this by first extracting the velocity component along the constrain, and evaluating velocity of center of mass $v_{cog} = (m_i v_i + m_j v_j)/(m_i +m_j)$. Then we consider that valocity with respect to center of mass $u_i = v_i - v_{cog}$ was inverted and rescaled by coeficient of restitution $\sqrt{\epsilon}$ in order to change corresponding kinetic energy by $\epsilon$. Therefore $ u'_i = \sqrt{\epsilon}(-v_i + v_{cog}) $ and after return to original reference fram we end up with $ v'_i = v_{cog} - \sqrt{\epsilon} (v_i - v_{cog}) = (1+\varepsilon) v_{cog} - \varepsilon v_i $.


 $ (1+e) v_{cog} - e v_i $
 $( (1+e) (m_i v_i + m_j v_j)- e (m_i + m_j) v_i    ) /(m_i +m_j)  $
 $( (1+e) (m_i v_i + m_j v_j)- e (m_i + m_j) v_i    ) /(m_i +m_j)  $

 $ m_i v_i + m_j v_j  + e m_i v_i + e m_j v_j - e m_i  v_i  - e m_j v_i  $
 $ m_i v_i + m_j v_j              + e m_j v_j               - e m_j v_i  $
 $ v_i' =  v_{cog} +                  e m_j( v_j - v_i )/m_{cog} $


$ \Delta v_i = v_i'-v_i =   (1+e)m_j( v_j - v_i )/m_{cog} $
$ m_i v_i + m_j v_j              + e m_j v_j               - e m_j v_i  -  m_i  v_i - m_j  v_i $
$         + m_j v_j              + e m_j v_j               - e m_j v_i              - m_j  v_i $
$     (1+e) m_j v_j              +                       -(1+e) m_j v_i                          $

 Both approaches produce a change of velocity:
  1. For **potential derivative** (resp. force) based update it is $ \Delta {\vec v}_i = ({\vec f}_{ij}/m_i) \Delta t $.
  2. For **velocity reflection** (resp. elastic restitution) it is $ \Delta {\vec v}_i =  (1+e)m_j( v_j - v_i )/m_{cog} $.

We should consider when these two velocity changes equals each other. If the momentum approach is limiting case of the potential approach we can switch betwee the two regimes when the two matches. 

$ \vec f_{ij}/m_i) \Delta t  =  (1+e)m_j( v_j - v_i )/m_{cog} $

$ \vec f_{ij} (m_{cog}/m_i) \Delta t  =  (1+e)m_j( v_j - v_i ) $