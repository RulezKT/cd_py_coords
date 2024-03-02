import math
from cd_py_math import atan2_in_degrees


"""	
# переводим в эклиптические прямоугольные (декартовы координаты) из экваториальных
#If eps_angle represents an angle in the range [0, π/2], 
then math.cos(eps_angle) will always be non-negative, 
and using abs() would be redundant. 
However, if eps_angle can be outside that range, then abs() ensures that 
the sign of c is preserved, which might be necessary depending on the context 
of your calculations.
"""


def equ_to_ecl(equ, eps):
    ecl = {}
    # print(equ)
    ecl["x"] = equ["x"]

    cos_eps = math.cos(eps)
    sin_eps = math.sin(eps)

    ecl["y"] = equ["y"] * abs(cos_eps) + equ["z"] * sin_eps

    ecl["z"] = -equ["y"] * sin_eps + equ["z"] * cos_eps

    # то же самое делаем для скорости
    # ecl.velocity_x = equ.velocity_x
    # ecl.velocity_y = equ.velocity_y * abs(cos_eps) + equ.velocity_z * sin_eps
    # ecl.velocity_z = -equ.velocity_y * sin_eps + equ.velocity_z * cos_eps

    # print(f"eps = {eps},\nequ ={equ},\necl = {ecl}\n")
    return ecl


def from_cartesian_to_polar(coords):
    """
    conversion of cartesian coordinates x,y,z
    into polar r, theta, phi
    theta in [-90 deg, +90 deg]
    phi in [0 deg, +360 deg]
    """

    x = coords["x"]
    y = coords["y"]
    z = coords["z"]
    # print(x, y, z)

    r = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))
    phi = atan2_in_degrees(y, x)
    if phi < 0:
        phi = phi + 360
    theta = atan2_in_degrees(z, (math.sqrt(math.pow(x, 2) + math.pow(y, 2))))

    print(f"longitude= {phi/57.29577951308232} theta = {theta} r = {r}")
    return {
        "radius": r,
        "latutude": theta,
        "longitude": phi,
    }


def cart_to_polarRAD(coords):
    """
    Conversion of cartesian coordinates (x, y, z) to polar coordinates (r, theta, phi).
    Theta is in the range [-pi/2, pi/2].
    Phi is in the range [0, 2*pi].
    """

    x = coords["x"]
    y = coords["y"]
    z = coords["z"]

    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.atan2(z, math.sqrt(x**2 + y**2))
    phi = math.atan2(y, x)
    if phi < 0:
        phi += 2 * math.pi

    # print(f"longitude= {phi} theta = {theta} r = {r}")
    return {
        "radius": r,
        "latitude": theta,
        "longitude": phi,
    }


def cart_to_polar_with_speed(x, y, z, vx, vy, vz):
    """
    Conversion of cartesian coordinates (x, y, z) and velocities (vx, vy, vz)
    to polar coordinates (r, theta, phi) and velocities (vr, vtheta, vphi).
    Theta is in the range [-90 deg, +90 deg].
    Phi is in the range [0 deg, +360 deg].
    """
    r = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))
    phi = atan2_in_degrees(y, x)
    if phi < 0:
        phi = phi + 360
    theta = atan2_in_degrees(z, (math.sqrt(math.pow(x, 2) + math.pow(y, 2))))

    vr = (x * vx + y * vy + z * vz) / r
    vtheta = (z * (x * vx + y * vy) - r * vz * math.sqrt(x**2 + y**2)) / (r**2 + z**2)
    vphi = (x * vy - y * vx) / (r**2 * math.sin(theta) ** 2)

    return r, theta, phi, vr, vtheta, vphi


if __name__ == "__main__":
    eps = 0.708242238014157
    equ = {"x": -27.182109028554287, "y": -11.058890899570304, "z": 4.934455939952573}
    ecl = equ_to_ecl(equ, eps)

    eps = 0.43994978906527565
    equ = {"x": -6.17769276839023, "y": -27.79971990559918, "z": -11.215292203937771}
    ecl = equ_to_ecl(equ, eps)

    eps = 0.4225764500783662
    equ = {"x": -12.607489704795924, "y": -12.329185215246486, "z": -5.207204495443294}
    ecl = equ_to_ecl(equ, eps)

    ecl = {"x": -27.182109028554287, "y": -5.189436885606041, "z": 10.941556935018141}
    cart_to_polarRAD(ecl)
