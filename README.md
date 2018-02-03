# Unscented Kalman Filter

[image1]: ./sim_image.png "Simulator view"

This project is part of the Udacity Self-Driving Car Engineer Nanodegre (SDCND) Program (. The project starter code can be found [here](https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project). To run the project, the Term 2 simulator is required from [here](https://github.com/udacity/self-driving-car-sim/releases).

### Goals
In this project, a Kalman filter is implemented to combine sensor measurements of lidar and radar for localization of a vehicle. Often, each measurement comes with a small, random error. By utilizing a Kalman filter, one can combine both sensors and reduce each errors to achieve a better result. 

![alt text][image1]

Here you can see the view of the simulator with two error-prone input sources (lidar and radar) in blue and red used for localization of the vehicle. The green triangles is the result of applying the Kalman filter to combine both inputs.


### Background
For an introduction to Kalman filters, please see the [wikipedia page](https://en.wikipedia.org/wiki/Kalman_filter). Often, the state transition function is non-linear, as is the case with any reasonable motion model in vehicles. For this purpose, the extended Kalman filter was developed (see [here](https://en.wikipedia.org/wiki/Extended_Kalman_filter)). It uses a Jacobian Matrix to linearize the state transition matrix, which can get very complex as more detailed motion models are used. An alternative approach is the unscented Kalman filter, which translates a few sample points precisely with the non-linear function to calculate a new mean and covariance. Therefore, errors stemming from the linearization are reduced.

### Dependencies

You need a websocket interface to communicate with the simulator, [uWebSockets](https://github.com/uWebSockets/uWebSockets). Installation of this tool can be done with the two scripts install-ubuntu.sh for Ubuntu Linux and install-mac.sh for MacOS. Under Windows, using a VM, Docker or Ubuntu Bash for Windows is recommended.

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

### Installation

1. Clone this repository
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt` or `./UnscentedKF` to use the standard `./data/obj_pose_laser_radar_synthetic-input.txt` data file