import logging

from mykrobe.species_data import DataDir

def describe(parser, args):
    args = parser.parse_args()
    ddir = DataDir(args.panels_dir)
    print(f"Gathering data from {ddir.root_dir}")
    ddir.print_panels_summary()

def update_metadata(parser, args):
    args = parser.parse_args()
    ddir = DataDir(args.panels_dir)
    ddir.update_manifest(filename=args.filename)

def update_species(parser, args):
    args = parser.parse_args()
    ddir = DataDir(args.panels_dir)
    logging.info(f"Loaded panels metdata from {ddir.root_dir}")

    if args.remove:
        if args.species == "all":
            raise NotImplementedError("Can only delete individual species")
        ddir.remove_species(args.species)
        logging.info(f"Removed species {args.species}")
    else:
        if args.species == "all":
            ddir.update_all_species()
        else:
            ddir.update_species(args.species)

