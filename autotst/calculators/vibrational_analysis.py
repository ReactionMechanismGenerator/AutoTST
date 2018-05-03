from autotst.reaction import *
from cclib.io import ccread


def percent_change(original,new):
    percent_change = (abs(new - original) / original) * 100
    return percent_change


class Vibrational_Analysis():
    """
    A class designed to perform vibrational analysis
    """

    def __init__(self, reaction):
        """
        reaction = AutoTST_Reaction
        """

        self.reaction = reaction

    def get_log_file(self, reaction):

        self.log_file = reaction.label + "_overall.log"

    def parse_vibrations(self):

        log_file_info = ccread(self.log_file)

        self.vibrations = zip(log_file_info.vibfreqs, log_file_info.vibdisps)

        return self.vibrations

    def obtain_geometries(self, reaction=None):

        assert reaction

        self.before_geometry = reaction.ts.ase_ts.copy()
        self.post_geometry = reaction.ts.ase_ts.copy()

        for vib, displacements in reaction.vibrations:
            if vib < 0: # Finding the imaginary frequency
                got_imaginary_frequency = True
                self.post_geometry.arrays["positions"] -= displacements

        return self.before_geometry, self.post_geometry

    def obtain_percent_changes(self, reaction):

        results = []
        for torsion in reaction.ts.torsions:
            i,j,k,l = torsion.indices
            before = self.before_geometry.get_dihedral(i,j,k,l)
            after = self.post_geometry.get_dihedral(i,j,k,l)
            results.append(["Tor", torsion.indices, torsion.reaction_center, percent_change(before, after)])

        for angle in reaction.ts.angles:
            i,j,k = angle.indices
            before = self.before_geometry.get_angle(i,j,k)
            after = self.post_geometry.get_angle(i,j,k)
            results.append(["Ang", angle.indices, angle.reaction_center, percent_change(before, after)])

        for bond in reaction.ts.bonds:
            i,j = bond.indices
            before = self.before_geometry.get_distance(i,j)
            after = self.post_geometry.get_distance(i,j)
            results.append(["Bond", bond.indices, bond.reaction_center, percent_change(before, after)])

        results = pd.DataFrame(results)
        results.columns = ["type", "index", "center", "percent_change"]

        self.percent_changes = results

    def validate_ts(self):

        self.get_log_file(self.reaction)

        self.parse_vibrations()

        self.obtain_geometries(self.reaction)

        self.obtain_percent_changes(self.reaction)

        if (np.log10(((self.percent_changes[self.percent_changes.center == "Yes"].mean()))) > np.log10(((self.percent_changes[self.percent_changes.center != "Yes"].mean()))) + 1).all():
            logging.info("Vibrational analysis was successful")
            return True

        else:
            logging.info("Cannot reasonably say that we have arrived at a TS through vibrational analysis.\nRunning an IRC calc.")
            return False
