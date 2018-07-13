#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"
#include <matrix/math.hpp>
#include <iostream>
#include <cfloat>

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

using matrix::SquareMatrix;

SquareMatrix<float, 4> QuadControl::getInvMixing()
{
  auto coll = matrix::Vector<double, 4>(std::vector<double>({1,1,1,1}).data());
  auto roll = matrix::Vector<double, 4>(std::vector<double>({1,-1,1,-1}).data()) * this->L;
  auto pitch = matrix::Vector<double, 4>(std::vector<double>({1,1,-1,-1}).data()) * this->L;
  auto yaw = matrix::Vector<double, 4>(std::vector<double>({-1,1,1,-1}).data()) * this->kappa;

  auto mixing = SquareMatrix<double, 4>();
  mixing.setRow(0, coll);
  mixing.setRow(1, roll);
  mixing.setRow(2, pitch);
  mixing.setRow(3, yaw);

  auto invMixing = mixing.I();

  auto final = SquareMatrix<float, 4>();

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      final._data[i][j] = (float) invMixing._data[i][j];
    }
  }

  return final;
}

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;

#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();

  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);

  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);

#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif

  demixing = getInvMixing();

  maxAcc = maxMotorThrust*4.f / mass;
//  // CAUTION: must be smaller than maxAcc: when all thrusters are close to limit you lose maneuverability
//  maxSafeAcc = maxMotorThrust*2.f / mass;
//  maxSafeAccSq = maxSafeAcc * maxSafeAcc;
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to
  //   individual motor thrust commands
  // INPUTS:
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS:
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // just left-multiply by invMixing
  float proto[] = {collThrustCmd, momentCmd.x, momentCmd.y, momentCmd.z};
  auto demixed = demixing * matrix::Vector<float, 4>(proto);
  auto data = (demixed.data());

  cmd.desiredThrustsN[0] = CONSTRAIN(data[0], minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[1] = CONSTRAIN(data[1], minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[2] = CONSTRAIN(data[2], minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[3] = CONSTRAIN(data[3], minMotorThrust, maxMotorThrust);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS:
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS:
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

//  float Iarray[3][3] = {{Ixx,0,0,}, {0,Iyy,0}, {0,0,Izz}};
//  auto IM = Matrix3f(*Iarray);

  auto errors = pqrCmd - pqr;
  auto proportions = errors * kpPQR;

  auto IMat = Mat3x3F::Zeros();
  IMat[0] = Ixx;
  IMat[4] = Iyy;
  IMat[8] = Izz;

//  auto Ivec = V3F(Ixx, Iyy, Izz);

  auto momentBase = IMat * proportions;
  auto momentEulerCorrection = pqr.cross(IMat * pqr);
  momentCmd = momentBase + momentEulerCorrection;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// source:
// https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
Quaternion<float> getRotationQuaternion(V3F v1, V3F v2)
{
  auto cross = v1.cross(v2);
  auto residual = std::sqrt(v1.magSq() * v2.magSq()) + v1.dot(v2);
  auto quat = Quaternion<float>(residual, cross.x, cross.y, cross.z).Normalise();
  return quat;
}

// returns a desired roll and pitch rate
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS:
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  auto R33 = R[8];

  auto accThrust = - collThrustCmd / mass;
  auto accThrust_z = accThrust * R33;

  auto a_cmd = V3F(accelCmd[0], accelCmd[1], accThrust_z);
  auto a_cmd_optimized = optimizeAccCmd(a_cmd);

  auto a_cmd_body = attitude.Rotate_ItoB(a_cmd_optimized);
  auto a_current_body = V3F(0, 0, accThrust);

  auto d_q = getRotationQuaternion(a_current_body, a_cmd_body);

  auto d_angle = d_q.ToEulerRPY();
  pqrCmd = V3F((float) d_angle.x * kpBank, (float) d_angle.y * kpBank, 0);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

//yield a more feasible banking angle
V3F QuadControl::optimizeAccCmd(V3F cmd)
{
  auto result = cmd;

  auto maxTiltSin = sin(maxTiltAngle);
  auto length = cmd.mag();

  auto maxZ = - std::max(length * maxTiltSin, -0.00001f); //z is down
  if (cmd.z >= maxZ) {
    result.z = maxZ;
  }

//  if (result.z >= -FLT_MIN){
//    //avoid flipping
////    cout << "flipping!\n";
//    result.z = -FLT_MIN;
//  }
//
//  auto a_cmd_L2Sq = result.magSq();
//  if (a_cmd_L2Sq > maxSafeAccSq) {
//
//    auto max_xyL2Sq = maxSafeAccSq - result.z * result.z;
//    if (max_xyL2Sq < 0) {
//      result = {0, 0, - maxSafeAcc};
//    }
//    else {
//      auto xyL2Sq = result.x * result.x + result.y * result.y;
//      auto ratio = std::sqrt(max_xyL2Sq / xyL2Sq);
//      result.x = result.x * ratio;
//      result.y = result.y * ratio;
//    }
//  }

  return result;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical
  //   acceleration feed-forward command
  // INPUTS:
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  auto error_z = posZCmd - posZ;
  auto error_dot_z = velZCmd - velZ;
  integratedAltitudeError += error_z *dt;

  auto p = kpPosZ * kpVelZ;
  auto d = kpVelZ;
  auto i = KiPosZ;

  auto acc_z = p*error_z + i*integratedAltitudeError + d*error_dot_z + accelZCmd;
  auto acc_z_withG = acc_z - (float) CONST_GRAVITY;

  auto R33 = R[8];
  float acc_thrust;
  if (R33 <=0) {
    cout << "flipped!\n";
    acc_thrust = maxAcc;
  }
  else {
    // counter-projection
    acc_thrust = CONSTRAIN(- acc_z_withG / R33, 0, maxAcc);
  }
  thrust = acc_thrust * mass;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS:
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations.
  //     the Z component should be 0
  // HINTS:
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  auto error_xy = posCmd - pos;
  auto error_dxy = velCmd - vel;

  auto p = kpPosXY * kpVelXY;
  auto d = kpVelXY;

  accelCmd += p*error_xy + d*error_dxy;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS:
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS:
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b].
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  auto error = yawCmd - yaw;
  yawRateCmd = kpYaw * error;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);

  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);

  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}

