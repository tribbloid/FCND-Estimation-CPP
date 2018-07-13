Project outline:

 - [Step 1: Sensor Noise](#step-1-sensor-noise)
 - [Step 2: Attitude Estimation](#step-2-attitude-estimation)
 - [Step 3: Prediction Step](#step-3-prediction-step)
 - [Step 4: Magnetometer Update](#step-4-magnetometer-update)
 - [Step 5: Closed Loop + GPS Update](#step-5-closed-loop--gps-update)
 - [Step 6: Adding Your Controller](#step-6-adding-your-controller)



### Step 1: Sensor Noise ###

I simply use Excel to read the generated csv file and calculate the mean and variance.

### Step 2: Attitude Estimation ###

The previous attitude of the vehicle is converted into quaternion, integrated with body rates measured by gyroscope, then converted back to Euler angles.

Yaw angle is normalized to be between -pi and pi.

### Step 3: Prediction Step ###

PredictState() and the main body of Predict() implements prediction of mean and covariance of the state respectively.

PredictState() use quaternion rotation to convert acceleration in body frame into inertial frame, then integrated into velocity in inertial frame, the previous speed is integrated into position.

Predict() first attempts to solve the Jacobian of the transition function first: 3 elements (the partial derivatives of 3 velocities in world frame w.r.t. yaw angle) of all the 49 in the matrix are solved numerically. Since dt is small enough, it can be used as a perturbation to numerically approximate partial derivatives without incurring too much round-up error. After this, the Jacobian is used directly to update the covariance matrix of the state.

### Step 4: Magnetometer Update ###

For Magnetometer, the observation matrix H is simply set to extract yaw from state vector, if the predicted yaw is too far away from the observed yaw, it is rounded up by +- 2pi to get closer.

### Step 5: Closed Loop + GPS Update ###

The implementation of GPS update is quite similar to Magnetometer update, except that the 6x7 observation matrix now extracts the first 6 elements of the state vector, now round up was performed.

### Step 6: Adding Your Controller ###

The QuadControl is replaced with my implementation for previous project, with is full quaternion based and has a much more sensitive set of parameters.

After migration the parameters turns out to be too sensitive to follow the estimated path and frequently caused recoverable flipping, which may indicates that Kalman gain is still too high. My solution is to increase the attitude control gain while decrease linear position/velocity gain, such that the vehicle can respond faster to jagged path but won't accelerate to high speed/overcommit to a waypoint even it is far from the vehicle.