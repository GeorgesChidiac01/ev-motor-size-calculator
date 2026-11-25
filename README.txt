EV Motor and Battery Sizing Calculator
======================================

A first-order sizing tool for the traction motor and battery of a small
electric car. Enter the vehicle, tires, driving conditions and performance
targets, and it estimates:

  - Road load (air drag, grade, rolling resistance) and wheel power
  - Motor speed, torque and power through a single-speed gear
  - Peak torque and power from acceleration, top speed and hill start
  - Approximate motor dimensions and mass
  - Inverter figures (phase voltage, current, electrical frequency)
  - Battery size, mass and peak discharge rate for a target range

Results are early-design estimates (roughly +/- 30% on motor size), not a
replacement for detailed electromagnetic or thermal design.


FILES
-----
  index.html             Interactive version with sliders and charts.
                         Open it in any web browser, or use it live at
                         https://georgeschidiac01.github.io/ev-motor-size-calculator/
  ev_motor_sizing.py     Command-line version. Same model, same results.


USING THE PYTHON SCRIPT
-----------------------
Requires Python 3. matplotlib is optional (only needed for --plot).

  python ev_motor_sizing.py                       run with default values
  python ev_motor_sizing.py --tire 185/60R15      try a different tire
  python ev_motor_sizing.py --gear 8.5 --vdc 800  change gear ratio and voltage
  python ev_motor_sizing.py --range-km 500        size the battery for 500 km
  python ev_motor_sizing.py --plot                show the torque-speed plot
  python ev_motor_sizing.py --help                list every parameter


DEFAULT CASE
------------
Kia Rio-sized car, 1800 kg, Cd 0.31, frontal area 2.2 m^2, P225/45R17 tires,
9:1 gear. Design cruise case: 135 km/h into a 25 m/s headwind on a 10% grade.
Range target: 400 km at 100 km/h.


MODEL
-----
  Road load        F = 0.5 rho Cd A (v + v_wind)^2 + m g sin(theta) + Crr m g cos(theta)
  Tire radius      r = (25.4 D + 2 W AR/100) / 2000   [m]
  Motor speed      n = 60 G v / (2 pi r)              [rpm]
  Motor torque     T = F r / (G eta_gear)             [N m]
  Accel. power     P = delta m (vf^2 + vb^2) / (2t) + 2/3 m g Crr vf + 1/5 rho Cd A vf^3
                   (Ehsani et al., Modern Electric, Hybrid Electric and Fuel
                   Cell Vehicles)
  Motor envelope   constant torque to base speed, then constant power to max speed
  Motor size       rotor volume = peak torque / torque per rotor volume (TRV)
  Battery          energy = range x consumption / usable fraction / (1 - margin)


LICENSE
-------
MIT License. Free to use, modify and share with attribution.
