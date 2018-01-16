
## ToDo 
- Formation BBox needs willpower-feedback which prohibits it to leave soldiers - it should be not possible that formation bbox moves while bulk of soldiers got stuck behind obstacle.
- soldiers perhaps should not rotate weapon with infinite speed


## What is done

- Loading different types of soldiers
- write soldier type parameters with captions to the string

## Alternative Solutions

#### poly-line formations
- Currently stick which is in sketches_SDL/3D/test_Stick is implemented using MMFF (cpp/common/dynamics/Molecular) it would be more efficient do 2D version
- two polylines can interact by circles in control points, but this is not very efficient, since than controlpoints must be very close. 
   - More efficient is define colision/interaction as distance between point and line-segment. Big problem if two lines accidently cross
   - also it would be good efficiently implement mid-angle-line between two segments which is normalize(normalize(perp(a))+normalize(perp(a))) (i.e. to avoid normalize, it is also cosinus of half angle so (ahat+bhat)/sqrt(1-dot(ahat,bhat)))
- Soldiers are splited into smaller groups (~8-16?) assigned to smaller bounding boxes attached to each line segment. Larger boundig box (with margin) is used for neighbor-search acceleration. If soldier leaves bounding box it is assigned to big bounding box of the whole formation. From that it is iteratively re-assigned to the nearest not-saturated segment-bbox. There is also mutual exchange of soldiers between neighboring segments.

