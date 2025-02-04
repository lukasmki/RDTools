import cmd

known_molecules = ["aspirin", "benzene", "caffeine"]


class App(cmd.Cmd):
    intro = "RDTools"
    prompt = "(RDTools) "

    def do_query(self, line):
        pass

    def do_load(self, line):
        self.stdout.write(f"{line}\n")
        pass

    def complete_load(self, text, line, start_index, end_index):
        if text:
            return [mol for mol in known_molecules if mol.startswith(text)]
        else:
            return known_molecules

    def do_add(self, line):
        pass

    def do_remove(self, line):
        pass

    def do_replace(self, line):
        pass

    def do_exit(self, line):
        return True


if __name__ == "__main__":
    App().cmdloop()
