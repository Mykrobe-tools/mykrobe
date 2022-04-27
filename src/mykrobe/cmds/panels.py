import logging
logger = logging.getLogger(__name__)

from mykrobe.species_data import DataDir

def describe(parser, args):
    ddir = DataDir(args.panels_dir)
    print(f"Gathering data from {ddir.root_dir}")
    ddir.print_panels_summary()

def update_metadata(parser, args):
    ddir = DataDir(args.panels_dir)
    ddir.update_manifest(filename=args.filename)

def update_species(parser, args):
    ddir = DataDir(args.panels_dir)
    logger.info(f"Loaded panels metdata from {ddir.root_dir}")

    if args.remove:
        if args.species == "all":
            raise NotImplementedError("Can only delete individual species")
        ddir.remove_species(args.species)
        logger.info(f"Removed species {args.species}")
    else:
        if args.species == "all":
            ddir.update_all_species()
        else:
            ddir.update_species(args.species)

