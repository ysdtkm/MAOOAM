with open("IC_atm66_ocn99.txt") as f:
    a = f.read().split()
with open("IC_atm66_ocn99mod.txt", "w") as f:
    for c in a:
        f.write(c)
        f.write("\n")
