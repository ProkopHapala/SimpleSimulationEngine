
# Fast Moving Targets
 - **Problem** - Shots are represented as line segments in collision-detection in order to work properly for fast moving targers. This is however generally not done for other objects (e.g. vehicles) since they has more complex collision shapes. However in some situations, e.g. Air combat, Space Combat speed of target is non-negligible. How to reliably and accurately determined hit in such sutuation which reasonable performance cost and algorithmic complexity?
 - **velocity gauge transform** - Physics is inveriant with repect to any intertial transformation. We may calculate the collision in coordinate system ov moving body (the target) as well. To do so we have to transform line segments of projectiels by velocity.
    **Pseudo-code** ```
        target.move(dt)->{pos,old_pos}
        shot  .move(dt)->{pos,old_pos}
        collide( shot.sweep_shape( shot.old_pos, shot.pos-(target.pos-target.old_pos) ), target.shape(target.old_pos) )
        ```
    or if we use leap-frog ```
        collide( shot.sweep_shape( shot.pos, (shot.vel-target.vel)*dt ), target.shape(target.pos) )
        ```
 - **Sphere sweep Broad phase** - In 
