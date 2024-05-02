import numpy as np
from scipy.stats import norm


class Position:
    """
    Represents a point in 3D space and supports basic vector arithmetic.

    Attributes:
        x (float): x-coordinate in 3D space.
        y (float): y-coordinate in 3D space.
        z (float): z-coordinate in 3D space.
    """

    def __init__(self, x=0, y=0, z=0):
        """
        Initializes a Position instance with given coordinates.

        Parameters:
            x (float): x-coordinate in 3D space.
            y (float): y-coordinate in 3D space.
            z (float): z-coordinate in 3D space.
        """
        self.x, self.y, self.z = x, y, z

    def __add__(self, other):
        """ Enables vector addition with another Position object. """
        return Position(self.x + other.x, self.y + other.y, self.z + other.z)

    def __truediv__(self, scalar):
        """ Enables division of the position's coordinates by a scalar. """
        return Position(self.x / scalar, self.y / scalar, self.z / scalar)

    def distance_to(self, other):
        """ Calculates Euclidean distance between this position and another. """
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2) ** 0.5


class Macromolecule:
    """
    Simulates a macromolecule using the freely jointed chain model.

    Attributes:
        degree_of_polymerization (int): Number of monomer units in the polymer.
        segment_length (float): Length of each segment in the chain model.
        positions (list of Position): List of positional objects representing the monomer units.
    """

    def __init__(self, degree_of_polymerization, segment_length=0.154):
        """
        Initializes a Macromolecule with a specified degree of polymerization and segment length.

        Parameters:
            degree_of_polymerization (int): Number of monomer units.
            segment_length (float): Length of each segment in the chain model.
        """
        self.degree_of_polymerization = degree_of_polymerization
        self.segment_length = segment_length
        self.positions = [Position()]
        self.simulate()

    def simulate(self):
        """ Generates the positions of each segment using a random direction for freely jointed model. """
        current_pos = self.positions[0]
        for _ in range(1, self.degree_of_polymerization):
            direction = np.random.normal(0, 1, 3)
            direction /= np.linalg.norm(direction)
            direction *= self.segment_length
            new_position = Position(current_pos.x + direction[0],
                                    current_pos.y + direction[1],
                                    current_pos.z + direction[2])
            self.positions.append(new_position)
            current_pos = new_position

    def center_of_mass(self):
        """ Calculates the center of mass of the polymer. """
        total_pos = Position()
        for pos in self.positions:
            total_pos += pos
        return total_pos / self.degree_of_polymerization

    def end_to_end_distance(self):
        """ Calculates the end-to-end distance of the polymer chain. """
        return self.positions[0].distance_to(self.positions[-1])

    def radius_of_gyration(self):
        """ Calculates the radius of gyration of the polymer. """
        cm = self.center_of_mass()
        return np.sqrt(sum(pos.distance_to(cm) ** 2 for pos in self.positions) / self.degree_of_polymerization)


def simulate_polymers(N, count):
    """
    Simulates a set of macromolecules and computes statistical properties.

    Parameters:
        N (int): Target degree of polymerization.
        count (int): Number of macromolecules to simulate.
    """
    degrees_of_polymerization = norm.rvs(loc=N, scale=0.1 * N, size=count)
    polymers = [Macromolecule(int(degree)) for degree in degrees_of_polymerization]

    centers_of_mass = [poly.center_of_mass() for poly in polymers]
    end_to_end_distances = [poly.end_to_end_distance() for poly in polymers]
    radii_of_gyration = [poly.radius_of_gyration() for poly in polymers]

    avg_center_of_mass = np.mean([[pos.x, pos.y, pos.z] for pos in centers_of_mass], axis=0)
    std_end_to_end = np.std(end_to_end_distances)
    avg_end_to_end = np.mean(end_to_end_distances)
    std_radius_of_gyration = np.std(radii_of_gyration)
    avg_radius_of_gyration = np.mean(radii_of_gyration)

    pdi = np.var(degrees_of_polymerization) / np.mean(degrees_of_polymerization)

    print(f"Metrics for {count} molecules of degree of polymerization = {N}")
    print(
        f"Avg. Center of Mass (nm) = {avg_center_of_mass[0]:.3f}, {avg_center_of_mass[1]:.3f}, {avg_center_of_mass[2]:.3f}")
    print(f"End-to-end distance (μm):\n\tAverage = {avg_end_to_end:.3f}\n\tStd. Dev. = {std_end_to_end:.3f}")
    print(
        f"Radius of gyration (μm):\n\tAverage = {avg_radius_of_gyration:.3f}\n\tStd. Dev. = {std_radius_of_gyration:.3f}")
    print(f"PDI = {pdi:.3f}")


def main():
    N = int(input("Degree of polymerization (default 1000)?: ") or 1000)
    count = int(input("How many molecules (default 50)?: ") or 50)
    simulate_polymers(N, count)


if __name__ == "__main__":
    main()
