def print_args(args):
    print("")
    print("+++++++++++++++++++++++++++++++")
    print("ARGUMENTS")
    print("---------")
    for arg, value in vars(args).items():
        print(f"{arg}: {value}")
    print("+++++++++++++++++++++++++++++++")
    print("")
