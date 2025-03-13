import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os  # For creating the output folder

# Constants
G = 6.674 * 10**-11  # Gravitational constant (m^3 kg^-1 s^-2)
M_sun = 1.989 * 10**30  # Mass of the Sun (kg)
AU = 1.496 * 10**11  # Astronomical unit (m) - average Earth-Sun distance

# Planet data (initial distance from Sun in AU, initial velocity in m/s)
planets = {
    "Mercury": {"r": 0.39 * AU, "v": 47.87 * 1000},
    "Venus": {"r": 0.72 * AU, "v": 35.02 * 1000},
    "Earth": {"r": 1.00 * AU, "v": 29.78 * 1000},
    "Mars": {"r": 1.52 * AU, "v": 24.07 * 1000},
    "Jupiter": {"r": 5.20 * AU, "v": 13.07 * 1000},
    "Saturn": {"r": 9.58 * AU, "v": 9.69 * 1000},
    "Uranus": {"r": 19.22 * AU, "v": 6.81 * 1000},
    "Neptune": {"r": 30.05 * AU, "v": 5.43 * 1000},
}

# Add the Moon (orbiting Earth)
moon = {
    "r": 0.00257 * AU,  # Distance from Earth (m)
    "v": 1.022 * 1000,  # Orbital velocity around Earth (m/s)
}

# Add an asteroid belt (simple representation)
asteroids = [{"r": (2.2 + i * 0.1) * AU, "v": 20 * 1000} for i in range(20)]  # 20 asteroids

# Time step and simulation duration
dt = 60 * 60 * 24  # 1 day in seconds
t_max = 60 * 60 * 24 * 365 * 165  # 165 years (enough for Neptune's orbit)

# Initialize positions and velocities
for planet in planets:
    planets[planet]["x"] = planets[planet]["r"]
    planets[planet]["y"] = 0
    planets[planet]["z"] = 0
    planets[planet]["vx"] = 0
    planets[planet]["vy"] = planets[planet]["v"]
    planets[planet]["vz"] = 0

# Initialize Moon (orbiting Earth)
moon["x"] = planets["Earth"]["x"] + moon["r"]
moon["y"] = planets["Earth"]["y"]
moon["z"] = 0
moon["vx"] = 0
moon["vy"] = planets["Earth"]["vy"] + moon["v"]
moon["vz"] = 0

# Initialize asteroids
for asteroid in asteroids:
    asteroid["x"] = asteroid["r"]
    asteroid["y"] = 0
    asteroid["z"] = 0
    asteroid["vx"] = 0
    asteroid["vy"] = asteroid["v"]
    asteroid["vz"] = 0

# Lists to store orbits for plotting
orbits = {planet: {"x": [], "y": [], "z": []} for planet in planets}
orbits["Moon"] = {"x": [], "y": [], "z": []}
orbits["Asteroids"] = [{"x": [], "y": [], "z": []} for _ in asteroids]

# Simulation loop
for t in np.arange(0, t_max, dt):
    # Update planets
    for planet in planets:
        # Distance from Sun
        r = np.sqrt(planets[planet]["x"]**2 + planets[planet]["y"]**2 + planets[planet]["z"]**2)
        
        # Gravitational force (components)
        Fx = -G * M_sun * planets[planet]["x"] / r**3
        Fy = -G * M_sun * planets[planet]["y"] / r**3
        Fz = -G * M_sun * planets[planet]["z"] / r**3
        
        # Update velocity
        planets[planet]["vx"] += Fx * dt
        planets[planet]["vy"] += Fy * dt
        planets[planet]["vz"] += Fz * dt
        
        # Update position
        planets[planet]["x"] += planets[planet]["vx"] * dt
        planets[planet]["y"] += planets[planet]["vy"] * dt
        planets[planet]["z"] += planets[planet]["vz"] * dt
        
        # Store orbit for plotting
        orbits[planet]["x"].append(planets[planet]["x"])
        orbits[planet]["y"].append(planets[planet]["y"])
        orbits[planet]["z"].append(planets[planet]["z"])

    # Update Moon (orbiting Earth)
    r_moon = np.sqrt((moon["x"] - planets["Earth"]["x"])**2 + 
                     (moon["y"] - planets["Earth"]["y"])**2 + 
                     (moon["z"] - planets["Earth"]["z"])**2)
    Fx_moon = -G * planets["Earth"]["x"] / r_moon**3
    Fy_moon = -G * planets["Earth"]["y"] / r_moon**3
    Fz_moon = -G * planets["Earth"]["z"] / r_moon**3
    
    moon["vx"] += Fx_moon * dt
    moon["vy"] += Fy_moon * dt
    moon["vz"] += Fz_moon * dt
    
    moon["x"] += moon["vx"] * dt
    moon["y"] += moon["vy"] * dt
    moon["z"] += moon["vz"] * dt
    
    orbits["Moon"]["x"].append(moon["x"])
    orbits["Moon"]["y"].append(moon["y"])
    orbits["Moon"]["z"].append(moon["z"])

    # Update asteroids
    for i, asteroid in enumerate(asteroids):
        r = np.sqrt(asteroid["x"]**2 + asteroid["y"]**2 + asteroid["z"]**2)
        Fx = -G * M_sun * asteroid["x"] / r**3
        Fy = -G * M_sun * asteroid["y"] / r**3
        Fz = -G * M_sun * asteroid["z"] / r**3
        
        asteroid["vx"] += Fx * dt
        asteroid["vy"] += Fy * dt
        asteroid["vz"] += Fz * dt
        
        asteroid["x"] += asteroid["vx"] * dt
        asteroid["y"] += asteroid["vy"] * dt
        asteroid["z"] += asteroid["vz"] * dt
        
        orbits["Asteroids"][i]["x"].append(asteroid["x"])
        orbits["Asteroids"][i]["y"].append(asteroid["y"])
        orbits["Asteroids"][i]["z"].append(asteroid["z"])

# Plot the orbits in 3D
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111, projection="3d")

# Plot planets
for planet in planets:
    ax.plot(orbits[planet]["x"], orbits[planet]["y"], orbits[planet]["z"], label=planet, linewidth=2)

# Plot Moon
ax.plot(orbits["Moon"]["x"], orbits["Moon"]["y"], orbits["Moon"]["z"], label="Moon", linestyle="--", linewidth=1.5, color="gray")

# Plot asteroids
for asteroid in orbits["Asteroids"]:
    ax.plot(asteroid["x"], asteroid["y"], asteroid["z"], color="red", linewidth=0.5, alpha=0.7)

# Plot the Sun at the center
ax.scatter([0], [0], [0], color="yellow", s=200, label="Sun")

# Set axis limits to zoom in on the inner solar system
ax.set_xlim([-2 * AU, 2 * AU])
ax.set_ylim([-2 * AU, 2 * AU])
ax.set_zlim([-2 * AU, 2 * AU])

# Add labels and title
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_title("3D Solar System Model with Moon and Asteroids (Inner Solar System)")
ax.legend(loc="upper right")

# Adjust the viewing angle
ax.view_init(elev=20, azim=30)  # Tilt the view for better perspective

# Add a grid
ax.grid(True)

# Save the plot as a PNG file
output_folder = "output_images"  # Folder name where the image will be saved
output_filename = "solar_system_plot.png"  # Name of the output file
output_path = f"{output_folder}/{output_filename}"  # Full path to the output file

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Save the figure
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Plot saved to {output_path}")

# Show the plot
plt.show()