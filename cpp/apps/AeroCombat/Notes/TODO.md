
# Simple

* make visible projectiles, so player can see some effects of cannons (curve on which they move)
  * ad some targets ( like moving on some curve, or aircrafts with stupid AI )
* Aircraft control mouse-aim like in warthunder
* test set of manuevers which Long-Shot is usually doing

# more work

* loading lua scripts to perform test manuevers


# Tactical Fighting

#### Can interceptor catch bomber?

1. Evaluate flight envelope of both aircrafts independently
	* steady state speed vs angle of climb
	* x = horizontal distance
	* y = vertical distance (climb, dive)
2. overlap the two envelopes
	* envelopes are shifted by initial advantage
		* distance between the aircrafts
		* relative height of the aircrafts

$$x_5 = $$