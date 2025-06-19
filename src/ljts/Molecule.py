class Molecule:
    """
    Just a molecule with a 3D position.
    """

    def __init__(self, x, y, z):
        # store position
        self._x = x
        self._y = y
        self._z = z

    def get_position(self):
        """
        Returns current position as (x, y, z).
        """
        return (self._x, self._y, self._z)
