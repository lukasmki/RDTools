import cmd

from RDTools import SMILES

known_molecules = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "benzene": "C1=CC=CC=C1",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}


class App(cmd.Cmd):
    intro = "RDTools"
    prompt = "(RDTools) "

    molecule = None

    def do_query(self, line):
        pass

    def do_load(self, line):
        "Load a molecule to context"
        self.molecule = SMILES(line)
        self.stdout.write(f"Loaded {line}\n")

    def complete_load(self, text, line, start_index, end_index):
        if text:
            return [sm for mol, sm in known_molecules.items() if mol.startswith(text)]
        else:
            return known_molecules.keys()

    def do_add(self, line):
        "Add a new molecule to context"
        self.molecule.add(line.strip().split())
        self.stdout.write(f"Added {line} -> {self.molecule}\n")

    def do_remove(self, line):
        "Remove a SMARTS pattern from all molecules in context"
        self.molecule.remove(line)
        self.stdout.write(f"Removed {line} -> {self.molecule}\n")

    def do_replace(self, line):
        "Replace a SMARTS pattern with a SMARTS fragment in all molecules in context"
        target, fragment = line.split()
        self.molecule.replace(target_smarts=target, fragment_smarts=fragment)
        self.stdout.write(f"Replaced {line} -> {self.molecule}\n")

    def do_react(self, line):
        "Perform a reaction defined by reaction SMARTS to all molecules in context.\nUsage: `react {REACTION_SMARTS}`\n"
        self.molecule.react(line)
        self.stdout.write(f"Reacted {line} -> {self.molecule}\n")

    def do_list(self, line):
        "List molecules currently loaded in context"
        self.columnize(str(self.molecule).split("."))

    def do_bye(self, line):
        return True


if __name__ == "__main__":
    App().cmdloop()
