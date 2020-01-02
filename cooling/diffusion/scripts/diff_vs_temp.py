    # Parse the number of atoms
    atoms = re.split('(\d+)', comp)
    atoms = sum([int(i) for i in atoms if i.isdigit()])

