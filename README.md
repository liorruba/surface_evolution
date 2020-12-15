# REGOLIT: REworking and Gardening of Lunar Impacted Terrains
## A surface evolution model for impacted surfaces to study mechanical gardening of ice and regolith
The model is based on [Richardson et al. 2009 CTEM model](https://www.sciencedirect.com/science/article/pii/S0019103509003194?casa_token=j6uDz1cdmAkAAAAA:UsgleZ2OBuARNfT8Gj0a2jaye59Fh9o4tzBj2rApSYEn_61GKxn3XCTfej-JPxHY2O2Un595JA)
but adds an efficient 3-D description of the subsurface using layers.

### Installation
1. Install [Boost](https://www.boost.org/).
2. To compile, simply run `make` in the `REGOLIT` directory
3. To run, execute `./build/apps/regolit_main.run`
4. The output will saved in the `output` directory: an ascii matrix of surface elevation (see animation below). The layers are saved in a binary file that can be read using the attached Matlab script.

### An example topography evolution simulation:
![Surface evolution](https://github.com/liorruba/surface_evolution/blob/master/craters.gif)
