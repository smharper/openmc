from numbers import Integral
from xml.etree import ElementTree as ET

import openmc.checkvalue as cv
from .mixin import EqualityMixin, IDManagerMixin


class TallyDerivative(EqualityMixin, IDManagerMixin):
    """A material perturbation derivative to apply to a tally.

    Parameters
    ----------
    derivative_id : int
        Unique identifier for the tally derivative. If none is specified, an
        identifier will automatically be assigned
    variable : str
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    materials : Iterable of int
        The perturbed material IDs
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    Attributes
    ----------
    id : int
        Unique identifier for the tally derivative
    variable : str
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    materials : int
        The perturubed material IDs
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    """

    next_id = 1
    used_ids = set()

    def __init__(self, derivative_id=None, variable=None, materials=None,
                 nuclide=None):
        # Initialize Tally class attributes
        self.id = derivative_id
        self.variable = variable
        self.materials = materials
        self.nuclide = nuclide

    def __repr__(self):
        string = 'Tally Derivative\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tVariable', self.variable)
        string += '{: <16}=\t{}\n'.format('\tMaterials', self.materials)
        if self.variable == 'nuclide_density':
            string += '{: <16}=\t{}\n'.format('\tNuclide', self.nuclide)

        return string

    @property
    def variable(self):
        return self._variable

    @property
    def materials(self):
        return self._materials

    @property
    def nuclide(self):
        return self._nuclide

    @variable.setter
    def variable(self, var):
        if var is not None:
            cv.check_type('derivative variable', var, str)
            cv.check_value('derivative variable', var,
                           ('density', 'nuclide_density', 'temperature'))
        self._variable = var

    @materials.setter
    def materials(self, mats):
        if mats is not None:
            cv.check_iterable_type('derivative materials', mats, Integral)
        self._materials = mats

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('derivative nuclide', nuc, str)
        self._nuclide = nuc

    def to_xml_element(self):
        """Return XML representation of the tally derivative

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing derivative data

        """

        element = ET.Element("derivative")
        element.set("id", str(self.id))
        element.set("variable", self.variable)
        element.set("materials", ' '.join(str(m) for m in self.materials))
        if self.variable == 'nuclide_density':
            element.set("nuclide", self.nuclide)
        return element
