"""
EV motor and battery sizing calculator
======================================

Estimates, for a small electric car:
  1. Road load (aero drag, grade, rolling resistance) at a design cruise case
  2. Motor operating point through a single-speed gear (speed, torque, power)
  3. Peak motor requirements from acceleration, top speed and hill start
  4. Approximate motor dimensions and mass (torque-per-rotor-volume sizing)
  5. Inverter / electrical figures (phase voltage, current, electrical frequency)
  6. Battery size and mass for a target range

All numbers are first-order engineering estimates for early design work.

Usage
-----
    python ev_motor_sizing.py                          # defaults
    python ev_motor_sizing.py --tire 185/60R15 --gear 8.5
    python ev_motor_sizing.py --range-km 500 --range-speed 110
    python ev_motor_sizing.py --plot                   # torque-speed plot (needs matplotlib)
    python ev_motor_sizing.py --help                   # every parameter

Model notes
-----------
Road load:        F = 1/2 rho Cd A (v + v_wind)^2 + m g sin(theta) + Crr m g cos(theta)
Tire radius:      r = (25.4 D + 2 W AR/100) / 2000   [m]
Motor speed:      n_m = 60 G v / (2 pi r)            [rpm]
Motor torque:     T_m = F r / (G eta_gear)           [N m]
Acceleration power (Ehsani, Modern EV/HEV/FCV, 2nd ed.):
                  P = delta m / (2 t) (v_f^2 + v_b^2) + 2/3 m g Crr v_f + 1/5 rho Cd A v_f^3
Motor envelope:   constant torque T_pk up to base speed, then constant power P_pk
                  up to max speed (base speed = max speed / CPSR)
Rotor volume:     V_r = T_pk / TRV,  D_r = (4 V_r / (pi k))^(1/3),  L = k D_r
Stator OD:        D_s = D_r / split ratio
Battery:          E = range * (Wh/km at range speed, flat, no wind) / usable fraction / (1 - margin)
"""

import argparse
import math
import re

G_ACCEL = 9.81  # m/s^2


# ----------------------------------------------------------------- defaults
DEFAULTS = dict(
    # tire
    tire_width=225, tire_aspect=45, rim_in=17,
    # vehicle
    mass=1800.0, cd=0.31, area=2.2, crr=0.008, rho=1.225,
    # drivetrain
    gear=9.0, eta_gear=0.97, eta_motor=0.92, poles=8,
    # design cruise case (sustained, sets continuous power)
    cruise_kmh=135.0, wind=25.0, grade=10.0,
    # performance (sets peak torque / power)
    accel_time=9.0, top_kmh=150.0, cpsr=3.5, hill_start=30.0, delta=1.05,
    # motor construction
    trv=140.0, aspect_ld=1.0, split=0.65,
    # electrical
    vdc=400.0, pf=0.85,
    # battery
    range_km=400.0, range_kmh=100.0, aux_kw=1.0, usable=0.90, margin=0.15, pack_whkg=160.0,
)


# ----------------------------------------------------------------- tire
def parse_tire(code):
    """'P225/45R17' -> (225, 45, 17)."""
    m = re.fullmatch(r"\s*[A-Za-z]*\s*(\d{3})\s*/\s*(\d{2})\s*[A-Za-z]*\s*(\d{2})\s*", code)
    if not m:
        raise ValueError(f"Can't read tire code '{code}'. Use a format like P225/45R17.")
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def tire_radius(width_mm, aspect, rim_in):
    h = width_mm * aspect / 100.0      # sidewall height, mm
    d = 25.4 * rim_in + 2.0 * h        # overall diameter, mm
    return h, d, d / 2000.0            # radius, m


# ----------------------------------------------------------------- physics
def road_forces(p, v, v_wind, grade_pct):
    """Road forces (N) at steady vehicle speed v (m/s)."""
    theta = math.atan(grade_pct / 100.0)
    v_air = v + v_wind
    F_aero = 0.5 * p["rho"] * p["cd"] * p["area"] * v_air * abs(v_air)
    F_grade = p["mass"] * G_ACCEL * math.sin(theta)
    F_roll = p["crr"] * p["mass"] * G_ACCEL * math.cos(theta)
    return F_aero, F_grade, F_roll


def compute(p):
    """Run the whole model. p is a dict with the keys in DEFAULTS."""
    out = {}
    rpm = 60.0 / (2.0 * math.pi)

    # --- tire
    h, d, r = tire_radius(p["tire_width"], p["tire_aspect"], p["rim_in"])
    out.update(sidewall_mm=h, tire_d_mm=d, r=r)

    # --- design cruise case (continuous requirement)
    v = p["cruise_kmh"] / 3.6
    Fa, Fg, Fr = road_forces(p, v, p["wind"], p["grade"])
    F = Fa + Fg + Fr
    out.update(v=v, F_aero=Fa, F_grade=Fg, F_roll=Fr, F_total=F,
               P_aero=Fa * v, P_grade=Fg * v, P_roll=Fr * v, P_wheel=F * v)
    w_wheel = v / r
    out.update(n_wheel=w_wheel * rpm, T_wheel=F * r,
               n_cruise=p["gear"] * w_wheel * rpm,
               T_cruise=F * r / (p["gear"] * p["eta_gear"]),
               P_cont=F * v / p["eta_gear"])

    # --- acceleration 0-100 km/h (peak requirement)
    v_f = 100.0 / 3.6
    v_top = p["top_kmh"] / 3.6
    v_b = v_top / p["cpsr"]                        # vehicle speed at motor base speed
    m_eff = p["delta"] * p["mass"]
    if v_b <= v_f:
        P_acc_inertia = m_eff / (2.0 * p["accel_time"]) * (v_f ** 2 + v_b ** 2)
    else:                                          # constant torque all the way to 100 km/h
        P_acc_inertia = m_eff * v_f * v_b / p["accel_time"]
    P_acc_resist = (2.0 / 3.0) * p["mass"] * G_ACCEL * p["crr"] * v_f \
        + 0.2 * p["rho"] * p["cd"] * p["area"] * v_f ** 3
    P_accel_wheel = P_acc_inertia + P_acc_resist

    # --- top speed on flat, no wind
    Fa_t, Fg_t, Fr_t = road_forces(p, v_top, 0.0, 0.0)
    P_top_wheel = (Fa_t + Fg_t + Fr_t) * v_top

    # --- peak power covers every case
    P_peak_wheel = max(P_accel_wheel, P_top_wheel, F * v)
    P_peak = P_peak_wheel / p["eta_gear"]
    drivers = {"acceleration": P_accel_wheel, "top speed": P_top_wheel, "cruise case": F * v}
    out["power_driver"] = max(drivers, key=drivers.get)

    # --- motor speeds
    w_base = p["gear"] * v_b / r
    w_max = p["gear"] * v_top / r
    out.update(v_base=v_b, n_base=w_base * rpm, n_max=w_max * rpm)

    # --- peak torque: constant-power corner, or hill start, whichever is larger
    th = math.atan(p["hill_start"] / 100.0)
    F_hill = p["mass"] * G_ACCEL * (math.sin(th) + p["crr"] * math.cos(th))
    T_hill = F_hill * r / (p["gear"] * p["eta_gear"])
    T_corner = P_peak / w_base
    T_peak = max(T_corner, T_hill)
    P_peak = max(P_peak, T_peak * w_base) if T_hill > T_corner else P_peak
    out.update(P_accel=P_accel_wheel / p["eta_gear"], P_top=P_top_wheel / p["eta_gear"],
               T_hill=T_hill, T_corner=T_corner, T_peak=T_peak, P_peak=P_peak,
               torque_driver="hill start" if T_hill > T_corner else "power at base speed")
    out["cruise_ok"] = v <= v_top * 1.0001

    # --- motor dimensions (TRV sizing)
    V_r = T_peak / (p["trv"] * 1000.0)             # m^3
    D_r = (4.0 * V_r / (math.pi * p["aspect_ld"])) ** (1.0 / 3.0)
    L = p["aspect_ld"] * D_r
    D_s = D_r / p["split"]
    end_winding = 0.18 * D_s                       # axial overhang per side (rough)
    housing_od = D_s + 2 * 0.025                   # 25 mm wall + cooling jacket
    housing_len = L + 2 * end_winding + 2 * 0.02   # 20 mm end covers
    active_mass = math.pi / 4 * D_s ** 2 * L * 6000.0
    out.update(V_rotor=V_r, D_rotor=D_r, L_stack=L, D_stator=D_s,
               housing_od=housing_od, housing_len=housing_len,
               mass_active=active_mass, mass_motor=2.0 * active_mass,
               power_density=P_peak / (2.0 * active_mass))

    # --- electrical
    V_ph_pk = p["vdc"] / math.sqrt(3.0)            # SVPWM, linear region
    V_ll_rms = p["vdc"] / math.sqrt(2.0)
    P_elec = P_peak / p["eta_motor"]
    I_rms = P_elec / (math.sqrt(3.0) * V_ll_rms * p["pf"])
    out.update(V_ph_pk=V_ph_pk, V_ll_rms=V_ll_rms, I_ph_rms=I_rms,
               f_base=p["poles"] * out["n_base"] / 120.0,
               f_max=p["poles"] * out["n_max"] / 120.0)

    # --- battery
    v_rng = p["range_kmh"] / 3.6
    Fa_r, Fg_r, Fr_r = road_forces(p, v_rng, 0.0, 0.0)
    wh_km_wheel = (Fa_r + Fr_r) / 3.6              # 1 N = 1 kJ/km = 0.2778 Wh/km
    wh_km = wh_km_wheel / (p["eta_gear"] * p["eta_motor"]) + p["aux_kw"] * 1000.0 / p["range_kmh"]
    E_needed = p["range_km"] * wh_km / 1000.0      # kWh
    E_pack = E_needed / p["usable"] / (1.0 - p["margin"])
    P_batt_peak = P_peak / p["eta_motor"] + p["aux_kw"] * 1000.0
    wh_km_cruise = F / 3.6 / (p["eta_gear"] * p["eta_motor"]) + p["aux_kw"] * 1000.0 / p["cruise_kmh"]
    out.update(wh_km=wh_km, E_needed=E_needed, E_pack=E_pack,
               batt_mass=E_pack * 1000.0 / p["pack_whkg"],
               P_batt_peak=P_batt_peak, c_rate=P_batt_peak / 1000.0 / E_pack,
               I_batt_peak=P_batt_peak / p["vdc"],
               range_cruise=E_pack * p["usable"] * 1000.0 / wh_km_cruise)
    return out


# ----------------------------------------------------------------- report
def report(p, o):
    k = lambda w: f"{w / 1000:.1f} kW"
    print("=" * 70)
    print(f"TIRE  {p['tire_width']}/{p['tire_aspect']}R{p['rim_in']}: "
          f"sidewall {o['sidewall_mm']:.1f} mm, diameter {o['tire_d_mm']:.1f} mm, radius {o['r']:.4f} m")

    print(f"\nDESIGN CRUISE CASE  {p['cruise_kmh']:.0f} km/h, {p['wind']:.0f} m/s headwind, {p['grade']:.0f}% grade")
    print(f"  Air drag            {o['F_aero']:8.1f} N   {k(o['P_aero'])}")
    print(f"  Grade               {o['F_grade']:8.1f} N   {k(o['P_grade'])}")
    print(f"  Rolling resistance  {o['F_roll']:8.1f} N   {k(o['P_roll'])}")
    print(f"  Total at wheels     {o['F_total']:8.1f} N   {k(o['P_wheel'])}")
    print(f"  Motor at this point {o['n_cruise']:8.0f} rpm, {o['T_cruise']:.1f} N m, "
          f"{k(o['P_cont'])} (continuous rating)")

    print("\nPEAK MOTOR REQUIREMENT")
    print(f"  Power for 0-100 km/h in {p['accel_time']:.1f} s     {k(o['P_accel'])}")
    print(f"  Power for {p['top_kmh']:.0f} km/h top speed       {k(o['P_top'])}")
    print(f"  Peak power (set by {o['power_driver']})   {k(o['P_peak'])}")
    print(f"  Torque for {p['hill_start']:.0f}% hill start         {o['T_hill']:.0f} N m")
    print(f"  Peak torque (set by {o['torque_driver']})  {o['T_peak']:.0f} N m")
    print(f"  Base speed {o['n_base']:.0f} rpm (car at {o['v_base'] * 3.6:.0f} km/h), "
          f"max speed {o['n_max']:.0f} rpm")
    if not o["cruise_ok"]:
        print("  WARNING: cruise speed is above top speed; raise top speed.")

    print("\nAPPROXIMATE MOTOR DIMENSIONS (+/- 30%)")
    print(f"  Rotor diameter      {o['D_rotor'] * 1000:6.0f} mm")
    print(f"  Stack length        {o['L_stack'] * 1000:6.0f} mm")
    print(f"  Stator outer dia.   {o['D_stator'] * 1000:6.0f} mm")
    print(f"  Housing (OD x len)  {o['housing_od'] * 1000:6.0f} mm x {o['housing_len'] * 1000:.0f} mm")
    print(f"  Motor mass          {o['mass_motor']:6.0f} kg  ({o['power_density'] / 1000:.1f} kW/kg peak)")

    print(f"\nELECTRICAL  ({p['vdc']:.0f} V DC bus, {p['poles']} poles)")
    print(f"  Max phase voltage   {o['V_ph_pk']:.0f} V peak, line-line {o['V_ll_rms']:.0f} V rms")
    print(f"  Phase current at peak power  {o['I_ph_rms']:.0f} A rms (PF {p['pf']})")
    print(f"  Electrical frequency  {o['f_base']:.0f} Hz at base speed, {o['f_max']:.0f} Hz at max speed")

    print(f"\nBATTERY  ({p['range_km']:.0f} km at {p['range_kmh']:.0f} km/h, flat, no wind)")
    print(f"  Consumption         {o['wh_km']:.0f} Wh/km (incl. {p['aux_kw']} kW auxiliaries)")
    print(f"  Energy needed       {o['E_needed']:.1f} kWh")
    print(f"  Pack size           {o['E_pack']:.1f} kWh  ({p['usable'] * 100:.0f}% usable, "
          f"{p['margin'] * 100:.0f}% margin)")
    print(f"  Pack mass           {o['batt_mass']:.0f} kg at {p['pack_whkg']:.0f} Wh/kg")
    print(f"  Peak battery power  {k(o['P_batt_peak'])}  ->  {o['c_rate']:.1f} C, {o['I_batt_peak']:.0f} A")
    print(f"  Range if driven at the design cruise case: {o['range_cruise']:.0f} km")
    print("=" * 70)


def plot(p, o):
    import matplotlib.pyplot as plt
    n = [o["n_max"] * i / 300 for i in range(1, 301)]
    w = [x / (60 / (2 * math.pi)) for x in n]
    env = [min(o["T_peak"], o["P_peak"] / wi) for wi in w]
    def demand(wind, grade):
        T = []
        for ni in n:
            v = ni / (60 / (2 * math.pi)) / p["gear"] * o["r"]
            F = sum(road_forces(p, v, wind, grade))
            T.append(F * o["r"] / (p["gear"] * p["eta_gear"]))
        return T
    plt.plot(n, env, lw=2, label="Motor peak envelope")
    plt.plot(n, demand(0, 0), lw=2, label="Steady speed, flat, no wind")
    plt.plot(n, demand(p["wind"], p["grade"]), lw=2, label="Steady speed, design grade + headwind")
    plt.scatter([o["n_cruise"]], [o["T_cruise"]], s=60, zorder=5, color="k", label="Design cruise point")
    plt.xlabel("Motor speed (rpm)"); plt.ylabel("Motor torque (N m)")
    plt.ylim(0, o["T_peak"] * 1.15); plt.grid(alpha=0.3); plt.legend()
    plt.title("Motor torque-speed requirement")
    plt.show()


def main():
    ap = argparse.ArgumentParser(description="EV motor and battery sizing")
    add = lambda name, key, help_: ap.add_argument(name, dest=key, type=type(DEFAULTS[key]),
                                                   default=DEFAULTS[key], help=f"{help_} (default {DEFAULTS[key]})".replace("%", "%%"))
    ap.add_argument("--tire", help="tire code like P225/45R17 (overrides width/aspect/rim)")
    add("--tire-width", "tire_width", "tire width, mm")
    add("--tire-aspect", "tire_aspect", "tire aspect ratio, %")
    add("--rim", "rim_in", "rim diameter, in")
    add("--mass", "mass", "vehicle mass incl. battery, kg")
    add("--cd", "cd", "drag coefficient")
    add("--area", "area", "frontal area, m^2")
    add("--crr", "crr", "rolling resistance coefficient")
    add("--rho", "rho", "air density, kg/m^3")
    add("--gear", "gear", "reduction gear ratio")
    add("--eta-gear", "eta_gear", "gearbox efficiency")
    add("--eta-motor", "eta_motor", "motor + inverter efficiency")
    add("--poles", "poles", "number of motor poles")
    add("--cruise-kmh", "cruise_kmh", "design cruise speed, km/h")
    add("--wind", "wind", "headwind at design cruise, m/s")
    add("--grade", "grade", "grade at design cruise, %")
    add("--accel-time", "accel_time", "0-100 km/h time, s")
    add("--top-kmh", "top_kmh", "top speed, km/h")
    add("--cpsr", "cpsr", "constant power speed ratio (max speed / base speed)")
    add("--hill-start", "hill_start", "steepest grade to start on, %")
    add("--delta", "delta", "rotating mass factor")
    add("--trv", "trv", "torque per rotor volume, kN m/m^3")
    add("--aspect-ld", "aspect_ld", "stack length / rotor diameter")
    add("--split", "split", "rotor dia. / stator outer dia.")
    add("--vdc", "vdc", "battery / DC bus voltage, V")
    add("--pf", "pf", "motor power factor at peak")
    add("--range-km", "range_km", "target range, km")
    add("--range-speed", "range_kmh", "speed for range estimate, km/h")
    add("--aux-kw", "aux_kw", "auxiliary load (HVAC etc.), kW")
    add("--usable", "usable", "usable fraction of pack")
    add("--margin", "margin", "extra margin for cold weather / ageing (0-1)")
    add("--pack-whkg", "pack_whkg", "pack-level energy density, Wh/kg")
    ap.add_argument("--plot", action="store_true", help="show torque-speed plot (matplotlib)")
    a = ap.parse_args()

    p = {k: getattr(a, k) for k in DEFAULTS}
    if a.tire:
        p["tire_width"], p["tire_aspect"], p["rim_in"] = parse_tire(a.tire)
    o = compute(p)
    report(p, o)
    if a.plot:
        plot(p, o)


if __name__ == "__main__":
    main()
