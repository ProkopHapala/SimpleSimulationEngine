
# Linear Motion

## Leap Frog 
 - **Pseudo-code** ```
    v(t+0.5*dt) = v(t-0.5*dt) + a(t)*dt
    x(t+dt)     = x(t) + v(t+0.5*dt)*dt
    ```
 - simplest integration scheme 
 - stable (conserve energy)
 - access to velocity, but in differet time   v(t+0.5dt)
 - we know that dx = v(t+0.5*dt)*dt, there is no need to store old_x 
 
### Verlet (Stormer)
 - **Pseudo-code** ``` x(t+dt) = 2*x(t) - x(t-dt) + a(t)*(dt*dt*0.5)  ```
 - does not require velocity
 - store previous position instead, this is very **usefull for collision detection**
 - problematic initialization ( set x(t-dt) ) when time-step is changed
 - velocity can be easily obtained as v(t-0.5*dt)=(x(t)-x(t-dt))/dt
 
### Verlet Difference
 - **Pseudo-code** ```
    dx(t+0.5*dt) = dx(t-dt*0.5) + a(t)*(dt*dt*0.5)
    x (t+dt)     = x(t) + dx(dt)
    ```
 - dx is more usefull than x(t-dt) in many situation (e.g. for sweeps in collision detection)
 - store previous position instead, this is very **usefull for collision detection**
 - problematic initialization ( set x(t-dt) ) when time-step is changed
 - velocity can be easily obtained as v(t-0.5*dt)=(x(t)-x(t-dt))/dt
 
### Using v-dir
 - We may use normalized direction vector vdir instead ov velocity or dx
 - This has advantage to easily compute many vecotr equations, without normalization to find hdir (i.e. saving sqrt() call)
    - collisions ( sphere-sweep, terrain-grid marching ) 
    - aerodynamics ( Fdrag =  )
