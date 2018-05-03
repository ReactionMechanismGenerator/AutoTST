import os
from autotst.reaction import *
from autotst.molecule import *
from autotst.geometry import *
from ase.io.gaussian import read_gaussian, read_gaussian_out

from autotst.calculators.vibrational_analysis import *

from ase.calculators.gaussian import *

from ase.optimize import BFGS

class AutoTST_Gaussian:

    def __init__(self, autotst_reaction, scratch="."):

        self.reaction = autotst_reaction

        self.scratch = scratch

        self.get_reactant_and_product_calcs(self.scratch)

        self.shell_calc = self.get_shell_calc(self.scratch)
        self.center_calc = self.get_center_calc(self.scratch)
        self.overall_calc = self.get_overall_calc(self.scratch)
        self.irc_calc = self.get_irc_calc(self.scratch)

        self.completed_irc = False

    def reactants_or_products_calc(self, autotst_mol, scratch="."):

        autotst_mol.rmg_molecule.updateMultiplicity()

        label = autotst_mol.rmg_molecule.toSMILES().replace("(", "left").replace(")", "right")

        calc = Gaussian(mem="5GB",
                        nprocshared="20",
                        label=label,
                        scratch=scratch,
                        extra="opt=(verytight,gdiis) freq IOP(2/16=3)",
                        multiplicity = autotst_mol.rmg_molecule.multiplicity
                        )
        del calc.parameters['force']

        return calc

    def get_reactant_and_product_calcs(self, scratch="."):

        self.reactant_calcs = {}
        self.product_calcs = {}

        for reactant in self.reaction.reactant_mols:
            calc = self.reactants_or_products_calc(reactant, scratch)
            self.reactant_calcs[reactant] = calc

        for product in self.reaction.product_mols:
            calc = self.reactants_or_products_calc(product, scratch)
            self.product_calcs[product] = calc


    def get_shell_calc(self, scratch="."):

        indicies = []
        for i, atom in enumerate(self.reaction.ts.rmg_ts.atoms):
            if not (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a,b = combo
            combos += "{0} {1} F\n".format(a+1,b+1)

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace("(", "left").replace(")", "right") + "_shell"

        calc = Gaussian(mem="5GB",
                        nprocshared="20",
                        label=label,
                        scratch=scratch,
                        extra="Opt=(ModRedun,Loose) Int(Grid=SG1)",
                        multiplicity = self.reaction.ts.rmg_ts.multiplicity,
                        addsec = [combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_center_calc(self, scratch="."):

        indicies = []
        for i, atom in enumerate(self.reaction.ts.rmg_ts.atoms):
            if (atom.label == ""):
                indicies.append(i)

        combos = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a,b = combo
            combos += "{0} {1} F\n".format(a+1,b+1)

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace("(", "left").replace(")", "right") + "_center"

        calc = Gaussian(mem="5GB",
                        nprocshared="20",
                        label=label,
                        scratch=scratch,
                        extra="opt=(ts,calcfc,noeigentest,ModRedun)",
                        multiplicity = self.reaction.ts.rmg_ts.multiplicity,
                        addsec = [combos[:-1]])

        del calc.parameters['force']
        return calc

    def get_overall_calc(self, scratch="."):

        self.reaction.ts.rmg_ts.updateMultiplicity()

        label = self.reaction.label.replace("(", "left").replace(")", "right") + "_overall"

        calc = Gaussian(mem="5GB",
                        nprocshared="20",
                        label=label,
                        scratch=scratch,
                        extra="opt=(ts,calcfc,noeigentest) freq",
                        multiplicity = self.reaction.ts.rmg_ts.multiplicity)

        del calc.parameters['force']
        return calc

    def get_irc_calc(self, scratch="."):

        self.reaction.ts.rmg_ts.updateMultiplicity()
        label = self.reaction.label.replace("(", "left").replace(")", "right") + "_irc"

        calc = Gaussian(mem="5GB",
                        nprocshared="20",
                        label=label,
                        scratch=scratch,
                        extra="irc=(calcall)",
                        multiplicity = self.reaction.ts.rmg_ts.multiplicity)

        del calc.parameters['force']
        return calc


    def calculate(self, autotst_object, calc):

        if isinstance(autotst_object, AutoTST_Molecule):
            if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
                calc.calculate(autotst_object.ase_molecule)
            autotst_object.ase_molecule = read_gaussian_out(os.path.join(calc.scratch, calc.label + ".log"))
            autotst_object.update_from_ase_mol()

        elif isinstance(autotst_object, AutoTST_Reaction):
            if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
                calc.calculate(autotst_object.ts.ase_ts)
            autotst_object.ts.ase_ts = read_gaussian_out(os.path.join(calc.scratch, calc.label + ".log"))
            autotst_object.ts.update_from_ase_ts()

        elif isinstance(autotst_object, AutoTST_TS):
            if not os.path.exists(os.path.join(calc.scratch, calc.label + ".log")):
                calc.calculate(autotst_object.ase_ts)
            autotst_object.ase_ts = read_gaussian_out(os.path.join(calc.scratch, calc.label + ".log"))
            autotst_object.update_from_ase_ts()

        return autotst_object

    def run_reactants_and_products(self):
        for mol, calc in self.reactant_calcs.iteritems():
            mol = self.calculate(mol, calc)

        for mol, calc in self.product_calcs.iteritems():
            mol = self.calculate(mol, calc)

    def run_shell(self):
        logging.info("Running shell optimization with center frozen...")
        self.reaction = self.calculate(self.reaction, self.shell_calc)
        logging.info("Shell optimization complete!")

    def run_center(self):
        logging.info("Running center optimization with shell frozen...")
        self.reaction = self.calculate(self.reaction, self.center_calc)
        logging.info("Center optimization complete!")

    def run_overall(self):
        logging.info("Running overall optimization...")
        self.reaction = self.calculate(self.reaction, self.overall_calc)
        logging.info("Overall optimization complete!")

    def run_irc(self):
        logging.info("Running IRC calculation")
        try:
            self.irc_calc.calculate(self.reaction.ts.ase_ts)
        except:
            # This normally fails because of an issue with ase's `read_results` method.
            pass
        logging.info("IRC calc complete!")

    def validate_irc(self):
        logging.info("Validating IRC file...")
        irc_path = os.path.join(self.irc_calc.scratch, self.irc_calc.label + ".log")
        if not os.path.exists(irc_path):
            logging.info("It seems that the file was `fixed`, reading in the `fixed` version.")
            irc_path = irc_path.replace("left", "(").replace("right", ")")

            if not os.path.exists(irc_path):
                logging.info("It seems that the IRC claculation has not been run.")
                irc_path = False

        if irc_path:
            f = open(irc_path, "r")
            file_lines = f.readlines()[-5:]

            for file_line in file_lines:
                if " Normal termination" in file_line:
                    logging.info("IRC successfully validated")
                    self.validated_irc = True
                else:
                    logging.info("IRC could not be validated")
                    self.validated_irc = False
        else:
            self.validated_irc = False

        return self.validated_irc


    def run_all(self, vibrational_analysis=True):
        self.run_shell()
        self.run_center()
        self.run_overall()

        vib = Vibrational_Analysis(self.reaction)
        logging.info("Performing Vibrational Analysis...")
        if vibrational_analysis and vib.validate_ts():
            logging.info("Vibrational analysis successful! Successfully arrived at a TS.")
            result = True
        elif vibrational_analysis and not vib.validate_ts():
            logging.info("Could not validate via vibrational analysis... \nRunning IRC instead...")
            self.run_irc()
            result = self.validate_irc()
        else:
            logging.info("Running without vibrational analysis... \nRunning IRC instead...")
            self.run_irc()
            result = self.validate_irc()

        logging.info("Fixing file names...")
        self.fix_io_files()

        if result:
            logging.info("Arrived at a TS!")
            return result

        else:
            logging.info("Could not arrive at a TS!")
            return result

    def fix_io_file(self, calc):
        old_log_file = calc.label + ".log"
        old_log_path = os.path.join(calc.scratch, old_log_file)
        if os.path.exists(old_log_path):
            new_log_path = old_log_path.replace("left", "(").replace("right", ")")
            os.rename(old_log_path, new_log_path)

        old_ase_file = calc.label + ".ase"
        old_ase_path = os.path.join(calc.scratch, old_ase_file)
        if os.path.exists(old_ase_path):
            new_ase_path = old_ase_path.replace("left", "(").replace("right", ")")
            os.rename(old_ase_path, new_ase_path)

        old_com_file = calc.label + ".com"
        old_com_path = os.path.join(calc.scratch, old_com_file)
        if os.path.exists(old_com_path):
            new_com_path = old_com_path.replace("left", "(").replace("right", ")")
            os.rename(old_com_path, new_com_path)


    def fix_io_files(self):

        for calc in self.reactant_calcs.values:
            self.fix_io_file(calc)

        for calc in self.product_calcs.values:
            self.fix_io_file(calc)

        self.fix_io_file(self.shell_calc)
        self.fix_io_file(self.center_calc)
        self.fix_io_file(self.overall_calc)
        self.fix_io_file(self.irc_calc)
