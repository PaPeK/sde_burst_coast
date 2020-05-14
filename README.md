# SDE Burst-Coast-Model – simulation code for agent-based burst-coast behavior (similar to Zebra Fish "Danio Rerio")
Version 0.1 described and applied in https://www.biorxiv.org/content/10.1101/809442v1.abstract

(Python-wrapper package calling c++ running)

[Github](https://github.com/PaPeK/sde_burst_coast)

## General Notes

SDE_burst_coast is a C++ code which runs fast numerical simulations of the agent based stochastic differential equations.
The motivation of this code is the swimming behavior of Zebra Fish (danio rerio) which and the base-parameters are also set such that the models resembles experimetal data.
The agents decide at each burst-initiation if they consider social- or environmental information, parametrized by the parameter "prob_social".
If agents respond to social cues they behave according to the three-zone-model (repulsion-alignment-attraction).
Otherwise agents respond to a predator if it is detected and present otherwise they chose a random burst-direction.
The random burst-direction reflects a false environmental cue received by the agents (reflection on water surface).
In fact the stochasticity is introduced in the randomness of the burst-decision if it is based on environmental cues.
Different scenarios are implemented:

- burst-coast: as in the experimental setup swim the agents in a circular tank with inelastic boundary conditions and a wall-avoidance mechanism
- natPred: a natural predator which follows the closest agent hunts the agents in periodic boundary conditions
- sinFisher: a single fisher (as the natural predator but randomly moving) hunts the agents
- mulFisher: multiple aligned fisher agents representing a trawling-net hunt the agents


## Required C++ libraries:

next to standard C libraries the following libraries need to be installed via ``apt-get install``'

- libhdf5-serial-dev
- libgsl-dev
- libcgal-dev

### Docker-alternative

The ```Dockerfile``` contains the build instruction for an image
which can be used to run the c++ code. In this case the above libraries do not need to be installed.
For building the docker:

- install docker following https://docs.docker.com/engine/install/.
- build the image by running 
    ```docker build -t gcc_docker .```
    in the directory of this git-repository
    - note that ```gcc\_docker``` is an arbitrary name for the image and can be chose freely
- check if the image exists by running ```docker images``` ("gcc_docker" should be listed)
- change in "RunSingle.py" the line ```dockerName = None``` to ```dockerName = 'gcc_docker'```

## Required python packages

- numpy
- h5py
- configparser
- time
- os

## Compilation(MacOS/Linux):

On MacOS/Linux system:
The included Makefile is for compilation of the source code (see tutorial: http://mrbook.org/blog/tutorials/make/).

```make``` 

in the terminal in the source directory should compile the source code.

### Docker-alternative

run in the directory of this repository
```docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp gcc_docker make```

## Running

After successful compilation the executable file swarmdyn will be generated.
The executable is meant to be called via running:

- RunSingle.py - starting a single simulation of swarmdyn and defining parameters.

RunSingle.py will also call a visualization of the simulations by calling:

- AnimateRun.py - animating the results using pythons matplotlib

## User Agreement

By downloading SDE_burst_coast you agree with the following points: SDE_burst_coast is provided without any warranty or conditions of any kind. We assume no responsibility for errors or omissions in the results and interpretations following from application of SDE_burst_coast.

## License

Copyright (C) 2016-2020 Pascal Klamser, Pawel Romanczuk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
