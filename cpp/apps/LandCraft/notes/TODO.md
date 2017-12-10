
# Simple

- Floading terrain
- critical angle of repose
  - some materials does not have angle of repose
- Edit terraain height - save, load


- build power plant, watermil
- Slopes with respect to sun control plants which can grow there
- Wind controled by terrain (coarse grained)



#### Which facilities depend on terrain?
- Harbor - cover ships from waves
- dam/powerplant - max height and gradient with minimal barrier
- irrigation
- channel/watterway with channel lock
- road/train

- Basins - map of all pixels associated to particular sink
    - each sink can set water level independently
    - each sink knows height of saddlepoints to other sinks (basin connectivity network)
    - How to Do it:
        - random droplet tracking which writes to know[] resp. contour
        - when we have pixel map of basis, we can just trace boundaries of these basins

- River Building
	- downward rainfall gathering
    		1. make sort pixels by ground height
    		2. trace rainfall on each pixel downwards according to height sorted order. After this operation each pixel knows total amout of rainwater which flow throught it per unit of time.
    		3. 
	- upward river tracking 
		-  identify sink-less pixelss (those which does not have lower laying neighbor
		- from each sinkless pixel build 1D river upwards
		- for each pixel find neighbor with highest flow, this is the main branch, if there are other neighbors with flow higher than treshhold, they are other side branches (feeders)
