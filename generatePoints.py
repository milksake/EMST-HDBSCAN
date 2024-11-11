import random

num_points = 1000
points = [(random.uniform(-50, 50), random.uniform(-50, 50)) for _ in range(num_points)]

with open("points.txt", "w") as f:
    for x, y in points:
        f.write(f"{x} {y}\n")
