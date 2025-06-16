class Molecule:
    def __init__(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z

    def get_position(self):
        return (self._x, self._y, self._z)